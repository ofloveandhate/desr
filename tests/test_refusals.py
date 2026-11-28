'''
Refusals that nothing had ever triggered.

A guard clause that is never reached is untested twice over: that it fires, and that the
message it builds can be built at all.  `is_invariant_expr` once ended in a line that
called `.format()` with no arguments, which no test or doctest had run.  Each test here
provokes one refusal and reads its message.
'''
from unittest import TestCase

import numpy as np
import sympy

from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.matrix_normal_forms import is_hnf_col
from desr.numerics import InsufficientKnownValues, NumericTranslator
from desr.ode_system import ODESystem, parse_de
from desr.ode_translation import ODETranslator, ReductionScheme

MICHAELIS_MENTEN = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                    'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
HL_6_4 = ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)']
HL_6_6 = ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
          'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)']


class TestTranslatorRefusals(TestCase):

    def setUp(self):
        self.system = ODESystem.from_equations(MICHAELIS_MENTEN)
        self.translator = ODETranslator(self.system)

    def test_a_naming_scheme_that_is_not_a_dict(self):
        with self.assertRaises(ValueError) as caught:
            ODETranslator(self.system, naming_scheme='kappa')
        self.assertIn('dict over the keys', str(caught.exception))
        self.assertIn('str', str(caught.exception))

    def test_a_scaling_matrix_that_does_not_match_the_symbols(self):
        with self.assertRaises(ValueError) as caught:
            ODETranslator(self.system, scaling_matrix=self.translator.scaling_matrix[:, :3])
        self.assertIn('has 3 columns and the system has 7 symbols', str(caught.exception))

    def test_the_dependent_variable_scheme_on_a_system_whose_clock_is_scaled(self):
        # Michaelis-Menten's scaling moves t, so that scheme does not apply.
        column = self.system.independent_variable_index
        self.assertFalse(self.translator.scaling_matrix[:, column].is_zero_matrix)
        with self.assertRaises(ValueError) as caught:
            self.translator.translate_dep_var()
        self.assertIn('is acted on by', str(caught.exception))

    def test_swapping_an_auxiliary_column(self):
        # Only invariant columns may be exchanged; an auxiliary's column of A*V is not
        # zero, and moving it would change the Hermite form.
        auxiliary = list(self.translator.auxiliary_monomials())[0]
        invariant = list(self.translator.invariant_monomials())[0]
        for first, second in ((auxiliary, invariant), (invariant, auxiliary)):
            with self.subTest(first=str(first)):
                with self.assertRaises(ValueError) as caught:
                    self.translator.multiplier_swap_columns(first, second)
                self.assertIn('Cannot swap non-zero column', str(caught.exception))

    def test_swapping_a_column_with_itself_changes_nothing(self):
        invariant = list(self.translator.invariant_monomials())[0]
        before = self.translator.herm_mult.copy()
        self.translator.multiplier_swap_columns(invariant, invariant)
        self.assertEqual(before, self.translator.herm_mult)


class TestSystemRefusals(TestCase):

    def test_a_line_that_is_not_a_differential_equation(self):
        with self.assertRaises(ValueError) as caught:
            parse_de('x = 3')
        self.assertIn('Invalid differential equation', str(caught.exception))

    def test_an_equation_in_another_independent_variable(self):
        with self.assertRaises(ValueError) as caught:
            parse_de('dx/ds = x', indep_var='t')
        self.assertIn('ordinary DEs in t', str(caught.exception))

    def test_reordering_with_the_wrong_symbols(self):
        system = ODESystem.from_equations(HL_6_4)
        with self.assertRaises(ValueError) as caught:
            system.reorder_variables(['t', 'z1', 'w'])
        self.assertIn('Mismatching variables', str(caught.exception))

    def test_reordering_accepts_a_string(self):
        system = ODESystem.from_equations(HL_6_4)
        system.reorder_variables('t z2 z1')
        self.assertEqual(('t', 'z2', 'z1'), tuple(map(str, system.symbols)))


class TestSystemsThatAreNotEqual(TestCase):
    ''' `==` had only ever been asked of systems that were equal. '''

    def test_different_equations(self):
        first = ODESystem.from_equations(['dx/dt = x'])
        second = ODESystem.from_equations(['dx/dt = 2*x'])
        self.assertNotEqual(first, second)

    def test_different_symbols(self):
        first = ODESystem.from_equations(['dx/dt = x'])
        second = ODESystem.from_equations(['dy/dt = y'])
        self.assertNotEqual(first, second)

    def test_a_parameter_against_a_conserved_variable(self):
        # dw/dt = 0 says w is conserved; no equation says w is a parameter.  Not the same.
        x, w, t = sympy.symbols('x w t')
        conserved = ODESystem({w: sympy.Integer(0), x: w * x})
        parameter = ODESystem({x: w * x})
        self.assertNotEqual(conserved, parameter)

    def test_something_that_is_not_a_system(self):
        self.assertNotEqual(ODESystem.from_equations(['dx/dt = x']), 'dx/dt = x')


class TestNumericRefusals(TestCase):

    def test_nothing_to_recover_when_every_auxiliary_is_carried(self):
        system = ODESystem.from_equations(HL_6_4)
        numeric = NumericTranslator(ODETranslator(system).translate())
        self.assertTrue(numeric.carries_auxiliaries)
        with self.assertRaises(ValueError) as caught:
            numeric.recover_auxiliaries({})
        self.assertIn('nothing to recover', str(caught.exception))

    def test_known_values_and_auxiliaries_together(self):
        system = ChemicalReactionNetwork.from_diagram(
            'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
        numeric = NumericTranslator(ODETranslator(system).translate_general(include_aux_vars=False))
        with self.assertRaises(ValueError) as caught:
            numeric.recover_auxiliaries({}, known_values={'k_0_1': 1.0},
                                        auxiliaries=[1.0, 1.0, 1.0, 1.0])
        self.assertIn('not both', str(caught.exception))


class TestMoreTranslatorRefusals(TestCase):

    def test_a_multiplier_that_does_not_give_the_hermite_form(self):
        # A V = [1, -1] for V = I, which is not [h, 0].
        with self.assertRaises(ValueError) as caught:
            ODETranslator(ODESystem.from_equations(HL_6_4),
                          scaling_matrix=sympy.Matrix([[1, -1]]),
                          hermite_multiplier=sympy.eye(2))
        self.assertIn('is not in HNF', str(caught.exception))
        self.assertFalse(is_hnf_col(sympy.Matrix([[1, -1]])))

    def test_the_parameter_scheme_on_a_system_it_does_not_fit(self):
        # The parameter scheme refuses; the substitutions, which hold under every scheme,
        # still answer.
        system = ODESystem.from_equations(HL_6_6)
        translator = ODETranslator(system)
        self.assertNotEqual(ReductionScheme.PARAMETER, translator.scheme())
        with self.assertRaises(ValueError) as caught:
            translator.translate_parameter()
        self.assertIn('not compatible for parameter translation', str(caught.exception))
        self.assertEqual(set(system.symbols), set(translator.substitutions()))

    def test_reversing_with_the_wrong_number_of_solutions(self):
        system = ODESystem.from_equations(HL_6_6)
        translator = ODETranslator(system)
        one = sympy.symbols('a:1')
        for reverse in (lambda: translator.reverse_translate(one),
                        lambda: translator.reverse_translate_general(one),
                        lambda: translator.reverse_translate_dep_var(one, 0)):
            with self.assertRaises(ValueError) as caught:
                reverse()
            self.assertIn('Incorrect number of values', str(caught.exception))

class TestMoreNumericRefusals(TestCase):

    def carrying(self):
        system = ODESystem.from_equations(HL_6_4)
        return NumericTranslator(ODETranslator(system).translate())

    def dropping(self):
        system = ODESystem.from_equations(MICHAELIS_MENTEN)
        return NumericTranslator(ODETranslator(system).translate_general(include_aux_vars=False))

    def test_a_carried_auxiliary_with_no_value_given(self):
        numeric = self.carrying()
        self.assertTrue(numeric.carries_auxiliaries)
        values = {symbol: 1.0 for symbol in numeric.reduced.symbols
                  if symbol not in numeric.auxiliary_names}
        with self.assertRaises(InsufficientKnownValues):
            numeric.reverse(values)

    def test_auxiliaries_given_as_series(self):
        numeric = self.dropping()
        values = {symbol: np.ones(3) for symbol in numeric.reduced.symbols}
        with self.assertRaises(ValueError) as caught:
            numeric.reverse(values, auxiliaries=[np.ones(3)] * numeric.r)
        self.assertIn('single numbers', str(caught.exception))

    def test_recovering_with_no_times_or_with_an_invariant_missing(self):
        numeric = self.dropping()
        clock = numeric.reduced.independent_variable
        values = {symbol: np.ones(3) for symbol in numeric.reduced.symbols}
        without_clock = {symbol: series for symbol, series in values.items()
                         if symbol != clock}
        with self.assertRaises(ValueError) as caught:
            numeric.recover_auxiliaries(without_clock, auxiliaries=[1.0] * numeric.r)
        self.assertIn('No values given for {}'.format(clock), str(caught.exception))
        with self.assertRaises(ValueError) as caught:
            numeric.recover_auxiliaries({clock: values[clock]},
                                        auxiliaries=[1.0] * numeric.r)
        self.assertIn('No values given for the invariant', str(caught.exception))


class TestMultiplierAndSystemEdits(TestCase):

    def setUp(self):
        self.system = ODESystem.from_equations(MICHAELIS_MENTEN)
        self.translator = ODETranslator(self.system)

    def test_adding_a_column_to_itself(self):
        invariant = list(self.translator.invariant_monomials())[0]
        with self.assertRaises(ValueError) as caught:
            self.translator.multiplier_add_columns(invariant, invariant, 1)
        self.assertIn('to itself', str(caught.exception))

    def test_negating_an_auxiliary_column(self):
        auxiliary = list(self.translator.auxiliary_monomials())[0]
        with self.assertRaises(ValueError) as caught:
            self.translator.multiplier_negate_column(auxiliary)
        self.assertIn('Cannot negate non-zero column', str(caught.exception))

    def test_negating_an_invariant_column_inverts_that_invariant(self):
        before = list(self.translator.invariant_monomials())
        self.translator.multiplier_negate_column(before[0])
        after = list(self.translator.invariant_monomials())
        self.assertEqual(1, sympy.simplify(before[0] * after[0]))
        self.assertEqual(before[1:], after[1:])
        self.assertTrue(self.translator.is_invariant_expr(after[0]))

    def test_a_constraint_against_zero(self):
        s, c = sympy.symbols('s c')
        with self.assertRaises(ValueError) as caught:
            self.system.add_constraint(s - c, 0)
        self.assertIn('equality constraint with 0', str(caught.exception))

    def test_a_substitution_that_is_neither_a_value_nor_a_renaming(self):
        s, c = sympy.symbols('s c')
        with self.assertRaises(ValueError) as caught:
            self.system.subs({s: c ** 2 + 1})
        self.assertIn('unable to make substitution', str(caught.exception))

    def test_translating_no_values_forward(self):
        numeric = NumericTranslator(self.translator.translate())
        with self.assertRaises(ValueError):
            numeric.forward({})
