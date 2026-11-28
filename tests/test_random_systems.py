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
#: for which `expected_symbol_count` did not count the auxiliary as normalised away.
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
    return system, ODETranslator(system)


@functools.lru_cache(maxsize=None)
def general(shape, seed, normalise, identify):
    system, translator = case(shape, seed)
    return translator.translate_general(normalise_parameters=normalise,
                                        rescale_clock=identify)


@functools.lru_cache(maxsize=None)
def dep_var(shape, seed, normalise):
    system, translator = case(shape, seed)
    return translator.translate_dep_var(normalise_parameters=normalise)


@functools.lru_cache(maxsize=None)
def parameter(shape, seed):
    system, translator = case(shape, seed)
    return translator.translate_parameter()


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
        system = type(system)._from_variables([t, u, p], [sympy.Integer(1), p * u, None],
                                              indep_var=t)
        self.assertEqual([], rows_that_are_not_symmetries(
            system, sympy.Matrix([[1, 0, -1], [0, 1, 0]])))
        self.assertEqual([0], rows_that_are_not_symmetries(
            system, sympy.Matrix([[0, 0, 1], [0, 1, 0]])))


class TestEveryReducedEquationIsTrue(TestCase):

    def test_under_the_scheme_translate_chooses(self):
        for shape, seed in CASES:
            with self.subTest(shape=shape, seed=seed):
                system, translator = case(shape, seed)
                self.assertEqual([], false_equations(system, translator.translate()))

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

    def test_the_general_scheme_with_both_steps_is_the_parameter_scheme(self):
        with_symmetry = 0
        for shape, seed in cases_under(ReductionScheme.PARAMETER):
            with self.subTest(shape=shape, seed=seed):
                self.assertEqual([], disagreements(general(shape, seed, True, True),
                                                   parameter(shape, seed)))
                with_symmetry += bool(case(shape, seed)[1].r)
        self.assertGreaterEqual(with_symmetry, 15)

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


class TestTheClockAndTheParametersOfAGeneralReduction(TestCase):
    '''
    One clock.  If some invariant is the independent variable times parameters and
    nothing else, the reduced system is written against it and the original clock is
    gone; otherwise the original clock stays.  Whether the auxiliaries were normalised
    does not come into it.  And an invariant built from parameters alone is a parameter
    of the reduced system, named as one.
    '''

    def test_a_rescaled_time_is_the_clock_and_t_is_not(self):
        changed = kept = autonomous = 0
        for shape, seed in CASES:
            system, translator = case(shape, seed)
            t = system.independent_variable
            parameters = set(system.parameters)
            # k*t with k a monomial in the parameters: t to the first power, by
            # differentiating, and nothing else that moves.
            rescaled = [monomial for monomial in translator.invariant_monomials()
                        if monomial != t and sympy.simplify(
                            t * sympy.diff(monomial, t) - monomial) == 0
                        and monomial.free_symbols <= parameters | {t}]
            for normalise in (True, False):
                with self.subTest(shape=shape, seed=seed, normalise=normalise):
                    reduced = general(shape, seed, normalise, True)
                    denotes = reduced.named_symbols
                    clock = reduced.independent_variable
                    if not rescaled:
                        # The original clock, under its own name.  One exception, which
                        # needs the normalisation: an auxiliary t/p whose p is set to 1
                        # is then the clock itself, and the clock is named for it.
                        if denotes[clock] != t:
                            self.assertTrue(normalise)
                            self.assertIn(denotes[clock],
                                          list(translator.auxiliary_monomials()))
                            self.assertEqual('tau', str(clock))
                        else:
                            self.assertEqual(t, clock)
                        kept += 1
                        continue
                    changed += 1
                    self.assertEqual(rescaled[0], denotes[clock])
                    self.assertEqual('tau', str(clock))
                    # The original clock is no invariant of the reduced system.  It can
                    # still be an auxiliary the system carries.
                    carried = set(reduced.auxiliary_names)
                    self.assertNotIn(t, {monomial for symbol, monomial in denotes.items()
                                         if symbol not in carried})
                    # Where the original is autonomous and the clock is the only reduced
                    # symbol built with t, the reduced system is autonomous too.
                    others = [monomial for symbol, monomial in denotes.items()
                              if symbol != clock]
                    if (not any(rate.has(t) for rate in system.derivatives
                                if rate is not None)
                            and not any(monomial.has(t) for monomial in others)):
                        autonomous += 1
                        for rate in reduced.derivatives:
                            if rate is not None:
                                self.assertFalse(rate.has(clock))
        self.assertGreaterEqual(changed, 40)
        self.assertGreaterEqual(kept, 40)
        self.assertGreaterEqual(autonomous, 20)

    def test_declining_the_change_keeps_the_original_clock(self):
        for shape, seed in CASES:
            with self.subTest(shape=shape, seed=seed):
                system, _ = case(shape, seed)
                reduced = general(shape, seed, True, False)
                self.assertEqual(system.independent_variable,
                                 reduced.named_symbols[reduced.independent_variable])

    def test_an_invariant_in_the_parameters_alone_is_a_kappa(self):
        seen = 0
        for shape, seed in CASES:
            system, translator = case(shape, seed)
            parameters = set(system.parameters)
            reductions = [general(shape, seed, True, True)]
            if translator.scheme() == ReductionScheme.DEP_VAR:
                reductions.append(dep_var(shape, seed, True))
            for reduced in reductions:
                with self.subTest(shape=shape, seed=seed, scheme=reduced.scheme.name):
                    for symbol, monomial in reduced.named_symbols.items():
                        constant = bool(monomial.free_symbols) and \
                            monomial.free_symbols <= parameters
                        named_as_one = str(symbol).startswith('kappa')
                        if symbol == monomial:
                            continue  # an original symbol left unchanged keeps its name
                        if symbol in reduced.auxiliary_names:
                            continue  # an auxiliary is named as an auxiliary
                        self.assertEqual(constant, named_as_one,
                                         '{} denotes {}'.format(symbol, monomial))
                        seen += constant
                    self.assertEqual(
                        {symbol for symbol in reduced.symbols
                         if str(symbol).startswith('kappa')} - set(reduced.parameters),
                        set())
        self.assertGreaterEqual(seen, 40)


class TestIsolatingTheClock(TestCase):
    '''
    `isolate_clock` puts a translator in another basis of the same invariants, in which
    the independent variable appears in at most one.  It is the same scaling, so
    everything true of a reduction stays true; and where the invariant holding the clock
    is `k*t`, an autonomous system reduces to an autonomous one.
    '''

    def isolated(self, shape, seed):
        ''' The system, its translator as built, and a second one with the clock isolated. '''
        system, translator = case(shape, seed)
        isolated = ODETranslator(system)
        isolated.isolate_clock()
        return system, translator, isolated

    def test_the_clock_is_in_at_most_one_invariant_with_a_positive_power(self):
        rearranged = 0
        for shape, seed in CASES:
            with self.subTest(shape=shape, seed=seed):
                system, translator, isolated = self.isolated(shape, seed)
                t = system.independent_variable
                holding = [monomial for monomial in isolated.invariant_monomials()
                           if monomial.has(t)]
                self.assertLessEqual(len(holding), 1)
                for monomial in holding:
                    self.assertGreater(sympy.degree(monomial, t), 0)
                before = [m for m in translator.invariant_monomials() if m.has(t)]
                self.assertEqual(bool(before), bool(holding))
                rearranged += len(before) > 1
        self.assertGreaterEqual(rearranged, 15)

    def test_it_is_the_same_scaling_and_the_same_invariants(self):
        # Each new invariant is invariant, and each old one is a product of integer
        # powers of the new: the multiplier of the one in terms of the other is an
        # integer matrix of determinant +-1.
        for shape, seed in CASES:
            with self.subTest(shape=shape, seed=seed):
                system, translator, isolated = self.isolated(shape, seed)
                self.assertEqual(translator.scaling_matrix, isolated.scaling_matrix)
                self.assertEqual(translator.herm_form, isolated.herm_form)
                self.assertEqual(1, abs(isolated.herm_mult.det()))
                for monomial in isolated.invariant_monomials():
                    self.assertTrue(translator.is_invariant_expr(monomial))
                # V' = V C.  The new invariants are products of the old invariants alone,
                # by a unimodular block; each new auxiliary is the old one times
                # invariants, so it scales as it did.
                change = translator.inv_herm_mult * isolated.herm_mult
                r = translator.r
                self.assertTrue(change[:r, r:].is_zero_matrix)
                self.assertEqual(1, abs(change[r:, r:].det()))
                self.assertEqual(sympy.eye(r), change[:r, :r])

    def test_it_changes_the_translator_it_is_called_on_and_nothing_is_left_stale(self):
        for shape, seed in CASES[::3]:
            with self.subTest(shape=shape, seed=seed):
                system, _ = case(shape, seed)
                translator = ODETranslator(system)
                translator.inv_herm_mult        # computed, and kept, before the change
                translator.isolate_clock()
                n = translator.n
                self.assertEqual(sympy.eye(n),
                                 translator.inv_herm_mult * translator.herm_mult)
                # What it reports is what it reduces with.
                reduced = translator.translate()
                self.assertEqual(dict(translator.named_symbols()),
                                 dict(reduced.named_symbols))

    def test_doing_it_twice_is_doing_it_once(self):
        for shape, seed in CASES[::3]:
            with self.subTest(shape=shape, seed=seed):
                system, _, isolated = self.isolated(shape, seed)
                once = isolated.herm_mult.copy()
                isolated.isolate_clock()
                self.assertEqual(once, isolated.herm_mult)

    def test_another_translator_for_the_same_system_is_not_touched(self):
        shape, seed = CASES[0]
        system, translator = case(shape, seed)
        before = translator.herm_mult.copy()
        self.isolated(shape, seed)
        self.assertEqual(before, translator.herm_mult)

    def test_every_reduced_equation_is_still_true(self):
        for shape, seed in CASES:
            with self.subTest(shape=shape, seed=seed):
                system, translator, isolated = self.isolated(shape, seed)
                self.assertEqual([], false_equations(system, isolated.translate()))
                self.assertEqual([], false_equations(
                    system, isolated.translate_general(normalise_parameters=False)))

    def test_an_autonomous_system_with_a_rescaled_time_reduces_to_an_autonomous_one(self):
        # The statement the operation is for.  The premise is read off the isolated
        # multiplier: the invariant holding t is k*t, with k in the parameters.
        checked = were_not = 0
        for shape, seed in CASES:
            system, translator, isolated = self.isolated(shape, seed)
            t = system.independent_variable
            if any(rate.has(t) for rate in system.derivatives if rate is not None):
                continue
            holding = [monomial for monomial in isolated.invariant_monomials()
                       if monomial.has(t)]
            if not holding or holding[0].free_symbols - set(system.parameters) - {t} \
                    or sympy.degree(holding[0], t) != 1:
                continue
            with self.subTest(shape=shape, seed=seed):
                reduced = isolated.translate_general()
                clock = reduced.independent_variable
                self.assertEqual(holding[0], reduced.named_symbols[clock])
                for rate in reduced.derivatives:
                    if rate is not None:
                        self.assertFalse(rate.has(clock))
                checked += 1
                plain = general(shape, seed, True, True)
                were_not += any(rate is not None and rate.has(plain.independent_variable)
                                and symbol != plain.independent_variable
                                for symbol, rate in zip(plain.symbols, plain.derivatives))
        self.assertGreaterEqual(checked, 40)
        # And it changes something: some of these were not autonomous without it.
        self.assertGreaterEqual(were_not, 5)


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
                        translator.expected_symbol_count(scheme=ReductionScheme.GENERAL,
                            normalise_parameters=normalise, rescale_clock=identify))

    def test_the_count_under_the_dependent_variable_scheme(self):
        for shape, seed in cases_under(ReductionScheme.DEP_VAR):
            for normalise in (True, False):
                with self.subTest(shape=shape, seed=seed, normalise=normalise):
                    system, translator = case(shape, seed)
                    self.assertEqual(
                        len(dep_var(shape, seed, normalise).symbols),
                        translator.expected_symbol_count(scheme=ReductionScheme.DEP_VAR,
                            normalise_parameters=normalise))

    def test_the_count_under_the_parameter_scheme(self):
        for shape, seed in cases_under(ReductionScheme.PARAMETER):
            with self.subTest(shape=shape, seed=seed):
                system, translator = case(shape, seed)
                self.assertEqual(len(parameter(shape, seed).symbols),
                                 translator.expected_symbol_count())
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
        once = forget(translator.translate())
        second = ODETranslator(once, naming_scheme=SECOND_NAMES)
        return once, second, second.translate()

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
