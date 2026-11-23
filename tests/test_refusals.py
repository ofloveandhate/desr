'''
Refusals that nothing had ever triggered.

A guard clause that is never reached is untested twice over: that it fires, and that the
message it builds can be built at all.  `is_invariant_expr` once ended in a line that
called `.format()` with no arguments, which no test or doctest had run.  Each test here
provokes one refusal and reads its message.
'''
from unittest import TestCase

import sympy

from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.numerics import NumericTranslator
from desr.ode_system import ODESystem, parse_de
from desr.ode_translation import ODETranslator

MICHAELIS_MENTEN = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                    'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
HL_6_4 = ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)']


class TestTranslatorRefusals(TestCase):

    def setUp(self):
        self.system = ODESystem.from_equations(MICHAELIS_MENTEN)
        self.translator = ODETranslator.from_ode_system(self.system)

    def test_a_naming_scheme_that_is_not_a_dict(self):
        with self.assertRaises(ValueError) as caught:
            ODETranslator.from_ode_system(self.system, naming_scheme='kappa')
        self.assertIn('dict over the keys', str(caught.exception))
        self.assertIn('str', str(caught.exception))

    def test_symbols_that_do_not_match_the_columns(self):
        with self.assertRaises(ValueError) as caught:
            ODETranslator(self.translator.scaling_matrix,
                          symbols_domain=self.system.symbols[:3])
        self.assertIn('3 variables given', str(caught.exception))

    def test_the_dependent_variable_scheme_on_a_system_whose_clock_is_scaled(self):
        # Michaelis-Menten's scaling moves t, so that scheme does not apply.
        column = self.system.independent_variable_index
        self.assertFalse(self.translator.scaling_matrix[:, column].is_zero_matrix)
        with self.assertRaises(ValueError) as caught:
            self.translator.translate_dep_var(self.system)
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
        conserved = ODESystem([t, w, x], [sympy.Integer(1), sympy.Integer(0), w * x],
                              indep_var=t)
        parameter = ODESystem([w, t, x], [None, sympy.Integer(1), w * x], indep_var=t)
        self.assertNotEqual(conserved, parameter)

    def test_something_that_is_not_a_system(self):
        self.assertNotEqual(ODESystem.from_equations(['dx/dt = x']), 'dx/dt = x')


class TestNumericRefusals(TestCase):

    def test_nothing_to_recover_when_every_auxiliary_is_carried(self):
        system = ODESystem.from_equations(HL_6_4)
        numeric = NumericTranslator(ODETranslator.from_ode_system(system).translate(system))
        self.assertTrue(numeric.carries_auxiliaries)
        with self.assertRaises(ValueError) as caught:
            numeric.recover_auxiliaries({})
        self.assertIn('nothing to recover', str(caught.exception))

    def test_known_values_and_auxiliaries_together(self):
        system = ChemicalReactionNetwork.from_diagram(
            'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
        numeric = NumericTranslator(ODETranslator.from_ode_system(
            system).translate_general(system, include_aux_vars=False))
        with self.assertRaises(ValueError) as caught:
            numeric.recover_auxiliaries({}, known_values={'k_0_1': 1.0},
                                        auxiliaries=[1.0, 1.0, 1.0, 1.0])
        self.assertIn('not both', str(caught.exception))
