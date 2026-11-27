'''
The scaling matrix has full row rank, including when it has no rows.

Hubert & Labahn assume a scaling matrix of full row rank from their section 2.1 on, and
their Proposition 5.1 shows the matrix computed from a system always has it: with
:math:`U K` in row Hermite normal form "having exactly r zero rows", :math:`A` is the
last :math:`r` rows of :math:`U`.  When there are no zero rows there are no rows of
:math:`A`, the system has no scaling symmetry, and its reduction is the identity.

The statements tested here are the proposition's, the definition of a symmetry, and the
chain rule.  None is a comparison of the library with its own earlier output.
'''
from unittest import TestCase

import sympy

from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.matrix_normal_forms import normal_hnf_col
from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator, ReductionScheme

from oracles import false_equations, forget


WITH_SYMMETRY = {
    'Michaelis-Menten': lambda: ODESystem.from_equations(
        ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
         'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']),
    'predator-prey': lambda: ODESystem.from_equations(
        ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)']),
    'Hubert & Labahn 6.4': lambda: ODESystem.from_equations(
        ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)']),
    'Hubert & Labahn 6.6': lambda: ODESystem.from_equations(
        ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
         'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)']),
    'enzyme chain': lambda: ChemicalReactionNetwork.from_diagram(
        'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system(),
    'decay': lambda: ODESystem.from_equations(['dn/dt = -n/tc']),
}

WITHOUT_SYMMETRY = {
    'one variable': lambda: ODESystem.from_equations(['dz/dt = z + 1 + t']),
    'two variables': lambda: ODESystem.from_equations(
        ['dx/dt = x + y + 1', 'dy/dt = x*y + t + 2']),
    'with a parameter': lambda: ODESystem.from_equations(['dz/dt = a*z + a + z + t']),
}


class TestProposition51(TestCase):
    '''
    What Hubert & Labahn's Proposition 5.1 says of the matrix computed from a system.
    '''

    def every_system(self):
        for label, build in dict(WITH_SYMMETRY, **WITHOUT_SYMMETRY).items():
            yield label, build()

    def test_it_annihilates_the_exponent_matrix(self):
        # The defining condition: F is invariant under the scaling of A iff A K = 0.
        for label, system in self.every_system():
            with self.subTest(system=label):
                A, K = system.maximal_scaling_matrix(), system.exponent_matrix()
                self.assertEqual(sympy.zeros(A.rows, K.cols), A * K)

    def test_it_has_full_row_rank(self):
        # "The first property implies in particular that A is of full row rank."
        for label, system in self.every_system():
            with self.subTest(system=label):
                A = system.maximal_scaling_matrix()
                self.assertEqual(A.rows, A.rank())

    def test_its_column_hermite_normal_form_is_the_identity_and_zeros(self):
        # (i) The column Hermite normal form of A is [I_r, 0].
        for label, system in self.every_system():
            with self.subTest(system=label):
                A = system.maximal_scaling_matrix()
                form, multiplier = normal_hnf_col(A)
                expected = sympy.eye(A.rows).row_join(sympy.zeros(A.rows, A.cols - A.rows))
                self.assertEqual(expected, form)
                self.assertEqual(1, abs(multiplier.det()))

    def test_it_is_maximal(self):
        # (ii) Every integer B with B K = 0 is M A, so the rows of A span the whole left
        # kernel of K: its dimension is n minus the rank of K.
        for label, system in self.every_system():
            with self.subTest(system=label):
                A, K = system.maximal_scaling_matrix(), system.exponent_matrix()
                self.assertEqual(K.rows - K.rank(), A.rows)

    def test_it_has_as_many_columns_as_symbols_even_with_no_rows(self):
        for label, build in WITHOUT_SYMMETRY.items():
            with self.subTest(system=label):
                system = build()
                self.assertEqual((0, len(system.symbols)),
                                 system.maximal_scaling_matrix().shape)


class TestASystemWithNoScalingSymmetry(TestCase):
    '''
    Its exponent matrix has full row rank, so nothing but the zero vector annihilates it,
    the scaling matrix has no rows, and there is nothing to reduce by.
    '''

    def test_the_premise(self):
        # No scaling leaves these invariant: K has rank n, so its left kernel is trivial.
        for label, build in WITHOUT_SYMMETRY.items():
            with self.subTest(system=label):
                K = build().exponent_matrix()
                self.assertEqual(K.rows, K.rank())

    def test_it_reports_no_symmetries(self):
        for label, build in WITHOUT_SYMMETRY.items():
            with self.subTest(system=label):
                system = build()
                translator = ODETranslator(system)
                self.assertEqual(0, translator.num_scaling_symmetries)
                self.assertEqual([], system.scaling_symmetries())
                self.assertEqual([], list(translator.auxiliary_monomials()))
                self.assertEqual({}, translator.normalisations())

    def test_every_symbol_is_its_own_invariant(self):
        # With nothing scaled, each symbol is invariant, and the n of them generate.
        for label, build in WITHOUT_SYMMETRY.items():
            with self.subTest(system=label):
                system = build()
                translator = ODETranslator(system)
                self.assertEqual(set(system.symbols),
                                 set(translator.invariant_monomials()))
                for symbol in system.symbols:
                    self.assertTrue(translator.is_invariant_expr(symbol))

    def test_every_symbol_is_dimensionless(self):
        # A structural dimension is a weight under the scaling, and there is no scaling.
        for label, build in WITHOUT_SYMMETRY.items():
            with self.subTest(system=label):
                dimensions = build().dimensions_of()
                self.assertEqual({sympy.Integer(1)}, set(dimensions.values()))

    def test_the_reduction_is_the_identity(self):
        for label, build in WITHOUT_SYMMETRY.items():
            with self.subTest(system=label):
                system = build()
                reduced = ODETranslator(system).translate()
                # The equations are true, there are as many symbols as before, and each
                # reduced symbol denotes one original symbol: a relabelling.
                self.assertEqual([], false_equations(system, reduced))
                self.assertEqual(len(system.symbols), len(reduced.symbols))
                self.assertEqual(set(system.symbols), set(reduced.named_symbols.values()))
                self.assertEqual(0, reduced.num_scaling_symmetries)
                # And it is the system itself, names included: nothing was changed, so
                # nothing is renamed.
                self.assertEqual({symbol: symbol for symbol in system.symbols},
                                 dict(reduced.named_symbols))
                self.assertEqual(system, forget(reduced))

    def test_the_general_scheme_agrees(self):
        # Hubert & Labahn's general reduction adds the retained clock, so n + 1, and with
        # no auxiliaries to carry.
        for label, build in WITHOUT_SYMMETRY.items():
            with self.subTest(system=label):
                system = build()
                reduced = ODETranslator(system).translate_general(normalise_parameters=False, rescale_clock=False)
                self.assertEqual(ReductionScheme.GENERAL, reduced.scheme)
                self.assertEqual([], false_equations(system, reduced))
                self.assertEqual(len(system.symbols) + 1, len(reduced.symbols))
                self.assertEqual((), tuple(reduced.auxiliary_names))

    def test_the_reports_all_answer(self):
        system = WITHOUT_SYMMETRY['two variables']()
        translator = ODETranslator(system)
        reduced = translator.translate()
        self.assertEqual(len(system.symbols), len(translator.unit_table()))
        self.assertEqual({'value'}, {row['role'] for row in translator.unit_table()})
        self.assertEqual(0, translator.num_units_and_dimensions()['units'])
        self.assertEqual(set(system.symbols),
                         {image.free_symbols.pop() and symbol
                          for symbol, image in translator.substitutions().items()})
        self.assertEqual(set(reduced.symbols),
                         set().union(*[image.free_symbols
                                       for image in translator.substitutions().values()]))


class TestReductionIsIdempotent(TestCase):
    '''
    Reducing a reduced system changes nothing.

    A reduction quotients by every scaling symmetry the system has.  What is left is
    either a system with none, whose scaling matrix has no rows, or a system whose only
    symmetries rescale an auxiliary it still carries -- the free constant of a
    quadrature.  In both cases a second reduction is a relabelling.
    '''

    def test_a_parameter_reduction_leaves_no_symmetry(self):
        for label in ('Michaelis-Menten', 'predator-prey'):
            with self.subTest(system=label):
                system = WITH_SYMMETRY[label]()
                translator = ODETranslator(system)
                self.assertEqual(ReductionScheme.PARAMETER, translator.scheme())
                once = forget(translator.translate())
                self.assertEqual(0, once.maximal_scaling_matrix().rows)

    def test_reducing_it_again_returns_it_unchanged(self):
        # Same names, same equations: the second reduction is the first one's output.
        for label in ('Michaelis-Menten', 'predator-prey'):
            with self.subTest(system=label):
                system = WITH_SYMMETRY[label]()
                once = forget(ODETranslator(system).translate())
                twice = ODETranslator(once).translate()
                self.assertEqual([], false_equations(once, twice))
                self.assertEqual({symbol: symbol for symbol in once.symbols},
                                 dict(twice.named_symbols))
                self.assertEqual(once, forget(twice))

    def test_and_a_third_time(self):
        system = WITH_SYMMETRY['Michaelis-Menten']()
        once = forget(ODETranslator(system).translate())
        twice = forget(ODETranslator(once).translate())
        thrice = forget(ODETranslator(twice).translate())
        self.assertEqual(once, twice)
        self.assertEqual(twice, thrice)

    def test_what_is_left_when_an_auxiliary_is_carried(self):
        # The only scalings remaining move a carried auxiliary, or the clock of a system
        # written as t dz/dt = z * F with F free of t (eq 20): the remaining rows of A are unit
        # vectors on those symbols.
        for label in ('Hubert & Labahn 6.4', 'Hubert & Labahn 6.6', 'enzyme chain',
                      'decay'):
            with self.subTest(system=label):
                system = WITH_SYMMETRY[label]()
                reduced = ODETranslator(system).translate()
                once = forget(reduced)
                A = once.maximal_scaling_matrix()
                moved = {once.symbols[column] for column in range(A.cols)
                         if not A[:, column].is_zero_matrix}
                allowed = set(reduced.auxiliary_names) | {reduced.independent_variable}
                self.assertTrue(moved)
                self.assertLessEqual(moved, allowed)
                for row in range(A.rows):
                    self.assertEqual(1, sum(abs(entry) for entry in A.row(row)))

    def test_reducing_that_again_is_a_relabelling(self):
        fresh = {'parameters': 'q', 'dependents': 'w', 'auxiliaries': 'u'}
        for label in ('Hubert & Labahn 6.4', 'Hubert & Labahn 6.6', 'enzyme chain',
                      'decay'):
            with self.subTest(system=label):
                system = WITH_SYMMETRY[label]()
                once = forget(ODETranslator(system).translate())
                twice = ODETranslator(once, naming_scheme=fresh).translate()
                self.assertEqual([], false_equations(once, twice))
                self.assertEqual(len(once.symbols), len(twice.symbols))
                # Each symbol of the second reduction denotes one symbol of the first.
                self.assertEqual(set(once.symbols), set(twice.named_symbols.values()))


class TestAScalingMatrixMustHaveFullRowRank(TestCase):
    '''
    A matrix handed to the translator directly is held to the paper's hypothesis.  A row
    that adds no scaling would be counted as a symmetry, the reduction would quotient by
    one that is not there, and the equations it returned would be false.
    '''

    def setUp(self):
        self.system = WITH_SYMMETRY['Michaelis-Menten']()
        self.real = self.system.maximal_scaling_matrix()
        self.n = self.real.cols

    def deficient(self):
        zero = sympy.zeros(1, self.n)
        return {
            'one zero row': zero,
            'two zero rows': sympy.zeros(2, self.n),
            'a real row and a zero row': self.real.row(0).col_join(zero),
            'the same row twice': self.real.row(0).col_join(self.real.row(0)),
            'two rows and their sum': self.real.col_join(
                self.real.row(0) + self.real.row(1)),
            'the real rows and a zero row': self.real.col_join(zero),
        }

    def test_each_is_rank_deficient(self):
        for label, matrix in self.deficient().items():
            with self.subTest(matrix=label):
                self.assertLess(matrix.rank(), matrix.rows)

    def test_each_is_refused(self):
        for label, matrix in self.deficient().items():
            with self.subTest(matrix=label):
                with self.assertRaises(ValueError) as caught:
                    ODETranslator(self.system, scaling_matrix=matrix)
                message = str(caught.exception)
                self.assertIn('full row rank', message)
                self.assertIn('{} row'.format(matrix.rows), message)
                self.assertIn('rank {}'.format(matrix.rank()), message)

    def test_a_matrix_of_full_row_rank_is_accepted(self):
        translator = ODETranslator(self.system, scaling_matrix=self.real)
        self.assertEqual(self.real.rows, translator.num_scaling_symmetries)

    def test_a_single_real_row_is_accepted(self):
        # A sub-group of the maximal scaling is still a scaling symmetry, and reducing
        # by it alone is a true reduction that leaves the other symmetry in place.
        translator = ODETranslator(self.system, scaling_matrix=self.real.row(0))
        reduced = translator.translate_general(normalise_parameters=False, rescale_clock=False)
        self.assertEqual(1, translator.num_scaling_symmetries)
        self.assertEqual([], false_equations(self.system, reduced))

    def test_a_matrix_with_no_rows_is_accepted(self):
        # Its rank is its row count.  It is the trivial scaling, and reducing by it is
        # the identity whatever symmetries the system has.
        translator = ODETranslator(self.system, scaling_matrix=sympy.zeros(0, self.n))
        self.assertEqual(0, translator.num_scaling_symmetries)
        reduced = translator.translate_general(normalise_parameters=False, rescale_clock=False)
        self.assertEqual([], false_equations(self.system, reduced))
        self.assertEqual(len(self.system.symbols) + 1, len(reduced.symbols))


class TestASharedNameThatMeansTheSameThing(TestCase):
    '''
    The name-clash refusal is about one spelling for two quantities.  A reduced symbol
    that denotes the original symbol of the same name is one spelling for one quantity.
    '''

    def test_a_name_that_would_denote_something_else_is_still_refused(self):
        # Predator-prey with its parameter `d` spelled `kappa0`.  The reduction names an
        # invariant `kappa0` too, and that invariant is `kappa0/K`: one spelling, two
        # quantities.
        system = ODESystem.from_equations(
            ['dn/dt = n*( r*(1 - n/K) - k*p/(n+kappa0) )', 'dp/dt = s*p*(1 - h*p / n)'])
        with self.assertRaises(ValueError) as caught:
            ODETranslator(system).translate()
        self.assertIn('kappa0', str(caught.exception))
        self.assertIn('two different quantities', str(caught.exception))

    def test_a_name_that_denotes_itself_is_not(self):
        # A system with no symmetry, already written in the names a reduction generates.
        nu0, tau = sympy.symbols('nu0 tau')
        system = ODESystem({nu0: nu0 + 1 + tau}, indep_var=tau)
        self.assertEqual(0, system.maximal_scaling_matrix().rows)
        reduced = ODETranslator(system).translate()
        self.assertEqual({symbol: symbol for symbol in system.symbols},
                         dict(reduced.named_symbols))
