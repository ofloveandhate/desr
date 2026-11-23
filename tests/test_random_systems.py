'''
Statements every reduction must satisfy, swept over systems nobody chose.

The other test modules use a dozen systems picked by hand, and a hand-picked system is
one somebody already understood.  Here the systems come from a seeded generator
(`random_systems.py`), each is reduced every way that applies, and each result is held
to statements that do not depend on how it was computed:

* every row of the scaling matrix is a scaling symmetry, by the definition of one;
* every reduced equation is true, by the chain rule;
* where two schemes apply, they give the same system up to names;
* what the translator reports about a reduction is what the reduction contains;
* reducing a reduced system gives a true reduction with as many symbols.

A failure names a shape and a seed, and `random_system(seed, **SHAPES[shape])` is the
system.  Nothing here may be refused either: an exception from any of these is a failure.
'''
import functools
import itertools
from unittest import TestCase

import sympy

from desr.ode_translation import ODETranslator, ReductionScheme

from oracles import false_equations, forget
from random_systems import SHAPES, random_system
from test_scheme_equivalence import disagreements

#: How many seeds of each shape, counted from zero.
SAMPLE = {'usual': 24, 'small': 40, 'non-autonomous': 15, 'large': 6,
          'reaction network': 40}

#: Seeds that once showed a fault, kept whatever the sample size becomes.  Each of these
#: is a dependent-variable system with an auxiliary that is a monomial in the parameters,
#: for which `expected_variable_count` did not count the auxiliary as normalised away.
ONCE_FAILED = [('small', 66), ('small', 214), ('small', 254), ('small', 312),
               ('large', 24), ('non-autonomous', 76), ('usual', 458)]

CASES = sorted(set(ONCE_FAILED) | {(shape, seed) for shape, count in SAMPLE.items()
                                   for seed in range(count)})

#: Every setting of the general scheme's two further steps.
STEPS = list(itertools.product((True, False), repeat=2))

#: A second reduction needs names the first did not use.
SECOND_NAMES = {'auxiliaries': 'w', 'dependents': 'v', 'parameters': 'q'}


@functools.lru_cache(maxsize=None)
def case(shape, seed):
    ''' The system for a case, and its translator. '''
    system = random_system(seed, **SHAPES[shape])
    return system, ODETranslator.from_ode_system(system)


@functools.lru_cache(maxsize=None)
def general(shape, seed, normalise, identify):
    system, translator = case(shape, seed)
    return translator.translate_general(system, normalise_parameters=normalise,
                                        identify_indep_var=identify)


@functools.lru_cache(maxsize=None)
def dep_var(shape, seed, normalise):
    system, translator = case(shape, seed)
    return translator.translate_dep_var(system, normalise_parameters=normalise)


@functools.lru_cache(maxsize=None)
def parameter(shape, seed):
    system, translator = case(shape, seed)
    return translator.translate_parameter(system)


def cases_under(scheme):
    return [(shape, seed) for shape, seed in CASES
            if case(shape, seed)[1].scheme() == scheme]


def rows_that_are_not_symmetries(system, scaling_matrix):
    '''
    The rows of a scaling matrix that are not scaling symmetries of the system.

    A row `a` scales each symbol `z_i` to `lambda^(a_i) z_i`, the clock among them.  It is
    a symmetry when the scaled symbols satisfy the same equations, which for
    `dz_i/dt = f_i(z)` says `f_i(lambda^a z) = lambda^(a_i - a_t) f_i(z)`.
    '''
    scale = sympy.Symbol('lambda_', positive=True)
    clock = system.independent_variable_index
    found = []
    for row in range(scaling_matrix.rows):
        weights = scaling_matrix[row, :]
        scaled = {symbol: scale ** weights[column] * symbol
                  for column, symbol in enumerate(system.symbols)}
        for column, rate in enumerate(system.derivatives):
            if rate is None or column == clock:
                continue
            difference = (rate.subs(scaled, simultaneous=True)
                          - scale ** (weights[column] - weights[clock]) * rate)
            if sympy.simplify(difference) != 0:
                found.append(row)
                break
    return found


class TestTheSampleIsWorthSweeping(TestCase):
    ''' So that none of the sweeps below passes for want of anything to check. '''

    def test_it_has_every_scheme(self):
        for scheme in ReductionScheme:
            with self.subTest(scheme=scheme.name):
                self.assertGreaterEqual(len(cases_under(scheme)), 10)

    def test_it_has_parameter_scheme_systems_with_symmetries_to_remove(self):
        # Without the reaction networks there are next to none: the parameter scheme
        # applies to the other shapes almost only where the scaling matrix has no rows.
        with_symmetry = [each for each in cases_under(ReductionScheme.PARAMETER)
                         if case(*each)[1].r > 0]
        self.assertGreaterEqual(len(with_symmetry), 15)

    def test_it_has_systems_with_no_symmetry_and_with_several(self):
        ranks = {case(*each)[1].r for each in CASES}
        self.assertLessEqual({0, 1, 2, 3, 4}, ranks)

    def test_a_seed_is_one_system(self):
        self.assertEqual(random_system(7), random_system(7))
        self.assertNotEqual(random_system(7), random_system(8))


class TestEveryRowOfTheScalingMatrixIsASymmetry(TestCase):

    def test_by_the_definition_of_a_scaling_symmetry(self):
        for shape, seed in CASES:
            with self.subTest(shape=shape, seed=seed):
                system, translator = case(shape, seed)
                self.assertEqual([], rows_that_are_not_symmetries(
                    system, translator.scaling_matrix))

    def test_a_row_that_is_not_a_symmetry_is_noticed(self):
        # du/dt = p u has the symmetries of A = [[1, 0, -1]] on (t, u, p) and no other
        # direction: scaling u alone is one, scaling p alone is not.
        t, u, p = sympy.symbols('t u p')
        system = random_system(0)  # for its class; rebuilt below
        system = type(system)([t, u, p], [sympy.Integer(1), p * u, None], indep_var=t)
        self.assertEqual([], rows_that_are_not_symmetries(
            system, sympy.Matrix([[1, 0, -1], [0, 1, 0]])))
        self.assertEqual([0], rows_that_are_not_symmetries(
            system, sympy.Matrix([[0, 0, 1], [0, 1, 0]])))


class TestEveryReducedEquationIsTrue(TestCase):

    def test_under_the_scheme_translate_chooses(self):
        for shape, seed in CASES:
            with self.subTest(shape=shape, seed=seed):
                system, translator = case(shape, seed)
                self.assertEqual([], false_equations(system, translator.translate(system)))

    def test_under_the_general_scheme_with_every_setting_of_its_steps(self):
        for shape, seed in CASES:
            system, translator = case(shape, seed)
            # With no symmetry there is no auxiliary to normalise and no invariant but
            # the symbols themselves, so the steps have nothing to act on.
            for normalise, identify in (STEPS if translator.r else STEPS[:1]):
                with self.subTest(shape=shape, seed=seed, normalise=normalise,
                                  identify=identify):
                    self.assertEqual([], false_equations(
                        system, general(shape, seed, normalise, identify)))

    def test_under_the_dependent_variable_scheme(self):
        for shape, seed in cases_under(ReductionScheme.DEP_VAR):
            for normalise in (True, False):
                with self.subTest(shape=shape, seed=seed, normalise=normalise):
                    system, _ = case(shape, seed)
                    self.assertEqual([], false_equations(
                        system, dep_var(shape, seed, normalise)))

    def test_under_the_parameter_scheme(self):
        for shape, seed in cases_under(ReductionScheme.PARAMETER):
            with self.subTest(shape=shape, seed=seed):
                system, _ = case(shape, seed)
                self.assertEqual([], false_equations(system, parameter(shape, seed)))

    def test_a_false_equation_is_noticed(self):
        # The same reduction with one right-hand side doubled.
        shape, seed = next(each for each in cases_under(ReductionScheme.PARAMETER)
                           if case(*each)[1].r > 0)
        system, _ = case(shape, seed)
        reduced = parameter(shape, seed)
        position = next(index for index, rate in enumerate(reduced.derivatives)
                        if rate is not None and rate != 0
                        and reduced.symbols[index] != reduced.independent_variable)
        rates = list(reduced.derivatives)
        rates[position] = 2 * rates[position]
        corrupted = reduced.copy()
        corrupted._derivatives = tuple(rates)
        self.assertEqual([reduced.symbols[position]], false_equations(system, corrupted))


class TestTheSchemesAgree(TestCase):
    ''' Up to names, compared through the monomial each reduced symbol denotes. '''

    def split_by_clock(self):
        '''
        The parameter-scheme cases, by whether the general scheme ends with one clock
        or keeps `t` beside a rescaled one.
        '''
        one, two = [], []
        for shape, seed in cases_under(ReductionScheme.PARAMETER):
            system, _ = case(shape, seed)
            monomials = set(general(shape, seed, True, True).named_symbols.values())
            rescaled = [monomial for monomial in monomials
                        if monomial != system.independent_variable
                        and system.independent_variable in monomial.free_symbols]
            (two if rescaled and system.independent_variable in monomials
             else one).append((shape, seed))
        return one, two

    def test_the_general_scheme_with_one_clock_is_the_parameter_scheme(self):
        one, _ = self.split_by_clock()
        self.assertGreaterEqual(len(one), 10)
        for shape, seed in one:
            with self.subTest(shape=shape, seed=seed):
                self.assertEqual([], disagreements(general(shape, seed, True, True),
                                                   parameter(shape, seed)))

    def test_the_general_scheme_with_two_clocks_is_the_parameter_scheme_and_t(self):
        # The parameter scheme takes the rescaled clock k*t as its clock and has no t.
        # The general scheme identifies k*t with t only when normalising the auxiliaries
        # makes k equal to 1.  When k is left as a combination of the invariant
        # parameters it keeps both: the parameter scheme's symbols, and t.  Both are
        # true reductions, which `TestEveryReducedEquationIsTrue` checks; this checks
        # that t is the whole of the difference.
        _, two = self.split_by_clock()
        self.assertGreaterEqual(len(two), 3)
        for shape, seed in two:
            with self.subTest(shape=shape, seed=seed):
                system, _ = case(shape, seed)
                by_general = set(general(shape, seed, True, True).named_symbols.values())
                by_parameter = set(parameter(shape, seed).named_symbols.values())
                self.assertEqual(by_parameter | {system.independent_variable}, by_general)
                self.assertNotIn(system.independent_variable, by_parameter)

    def test_the_general_scheme_is_the_dependent_variable_scheme(self):
        # Whether or not it is asked to identify a clock invariant: the scaling leaves
        # the clock alone, so the clock is already its own invariant.
        for shape, seed in cases_under(ReductionScheme.DEP_VAR):
            for normalise, identify in STEPS:
                with self.subTest(shape=shape, seed=seed, normalise=normalise,
                                  identify=identify):
                    self.assertEqual([], disagreements(
                        general(shape, seed, normalise, identify),
                        dep_var(shape, seed, normalise)))


class TestTheReportsDescribeTheReduction(TestCase):
    '''
    The translator answers questions about a reduction without running it.  Each answer
    is checked against the reduction it is about, under every scheme that applies.
    '''

    def reductions(self):
        ''' Every (case, scheme, reduction) with both steps on, as the reports assume. '''
        for shape, seed in CASES:
            scheme = case(shape, seed)[1].scheme()
            yield shape, seed, ReductionScheme.GENERAL, general(shape, seed, True, True)
            if scheme == ReductionScheme.PARAMETER:
                yield shape, seed, scheme, parameter(shape, seed)
            if scheme == ReductionScheme.DEP_VAR:
                yield shape, seed, scheme, dep_var(shape, seed, True)

    def test_the_named_symbols(self):
        for shape, seed, scheme, reduced in self.reductions():
            with self.subTest(shape=shape, seed=seed, scheme=scheme.name):
                _, translator = case(shape, seed)
                reported = translator.named_symbols(scheme=scheme)
                self.assertEqual(tuple(reduced.symbols), tuple(reported))
                self.assertEqual(dict(reduced.named_symbols), dict(reported))

    def test_the_names_by_role_are_the_reduced_systems_symbols_in_order(self):
        # The dependent-variable scheme's names are held to its own layout in
        # `test_vocabulary.py`.
        for shape, seed, scheme, reduced in self.reductions():
            if scheme == ReductionScheme.DEP_VAR:
                continue
            with self.subTest(shape=shape, seed=seed, scheme=scheme.name):
                _, translator = case(shape, seed)
                if scheme == ReductionScheme.PARAMETER:
                    expected = tuple(translator.invariant_names(scheme))
                else:
                    expected = ((reduced.independent_variable,)
                                + tuple(translator.auxiliary_names(scheme))
                                + tuple(translator.invariant_names(scheme)))
                self.assertEqual(expected, tuple(reduced.symbols))

    def test_exactly_one_symbol_is_the_one_to_integrate_against(self):
        for shape, seed, scheme, reduced in self.reductions():
            with self.subTest(shape=shape, seed=seed, scheme=scheme.name):
                _, translator = case(shape, seed)
                kinds = translator.symbol_kinds(scheme=scheme)
                clocks = [name for name, kind in kinds.items() if kind == 'independent']
                self.assertEqual([reduced.independent_variable], clocks)

    def test_the_count_under_the_general_scheme(self):
        for shape, seed in CASES:
            for normalise, identify in STEPS:
                with self.subTest(shape=shape, seed=seed, normalise=normalise,
                                  identify=identify):
                    system, translator = case(shape, seed)
                    self.assertEqual(
                        len(general(shape, seed, normalise, identify).symbols),
                        translator.expected_variable_count(
                            system, scheme=ReductionScheme.GENERAL,
                            normalise_parameters=normalise, identify_indep_var=identify))

    def test_the_count_under_the_dependent_variable_scheme(self):
        for shape, seed in cases_under(ReductionScheme.DEP_VAR):
            for normalise in (True, False):
                with self.subTest(shape=shape, seed=seed, normalise=normalise):
                    system, translator = case(shape, seed)
                    self.assertEqual(
                        len(dep_var(shape, seed, normalise).symbols),
                        translator.expected_variable_count(
                            system, scheme=ReductionScheme.DEP_VAR,
                            normalise_parameters=normalise))

    def test_the_count_under_the_parameter_scheme(self):
        for shape, seed in cases_under(ReductionScheme.PARAMETER):
            with self.subTest(shape=shape, seed=seed):
                system, translator = case(shape, seed)
                self.assertEqual(len(parameter(shape, seed).symbols),
                                 translator.expected_variable_count(system))
                self.assertEqual(len(system.symbols) - translator.r,
                                 len(parameter(shape, seed).symbols))

    def test_the_dependent_variable_scheme_can_normalise_an_auxiliary(self):
        # The premise of the count above: some of these systems lose a symbol to it.
        lost = [each for each in cases_under(ReductionScheme.DEP_VAR)
                if len(dep_var(*each, True).symbols) < len(dep_var(*each, False).symbols)]
        self.assertGreaterEqual(len(lost), 5)


class TestReducingAReducedSystem(TestCase):

    def again(self, shape, seed):
        system, translator = case(shape, seed)
        once = forget(translator.translate(system))
        second = ODETranslator.from_ode_system(once, naming_scheme=SECOND_NAMES)
        return once, second, second.translate(once)

    def test_gives_a_true_reduction_with_as_many_symbols(self):
        for shape, seed in CASES:
            with self.subTest(shape=shape, seed=seed):
                once, _, twice = self.again(shape, seed)
                self.assertEqual([], false_equations(once, twice))
                self.assertEqual(len(once.symbols), len(twice.symbols))

    def test_a_parameter_reduction_has_no_symmetry_left(self):
        for shape, seed in cases_under(ReductionScheme.PARAMETER):
            with self.subTest(shape=shape, seed=seed):
                once, second, twice = self.again(shape, seed)
                self.assertEqual(0, second.r)
                self.assertEqual(once, forget(twice))
