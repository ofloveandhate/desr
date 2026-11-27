'''
The inverse multiplier `W = V^-1`, its blocks, and completing a matrix to a unimodular one.
'''
from unittest import TestCase

import sympy

from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator, _int_inv, extend_rectangular_matrix

MICHAELIS_MENTEN = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                    'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
HL_6_4 = ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)']
HL_6_6 = ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
          'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)']


class TestTheBlocksOfTheInverseMultiplier(TestCase):
    '''
    `W = V^-1`, with `W_u` its first `r` rows and `W_d` the rest, so `W V = I` read
    block by block says `W_u V_i = I`, `W_u V_n = 0`, `W_d V_i = 0`, `W_d V_n = I`.
    '''

    def test_for_systems_under_each_scheme(self):
        for equations in (MICHAELIS_MENTEN, HL_6_4, HL_6_6):
            with self.subTest(equations=equations[0]):
                system = ODESystem.from_equations(equations)
                translator = ODETranslator(system)
                r, n = translator.r, translator.n
                self.assertGreater(r, 0)
                W_u, W_d = translator.inv_herm_mult_u, translator.inv_herm_mult_d
                V_i, V_n = translator.herm_mult_i, translator.herm_mult_n
                self.assertEqual(sympy.eye(r), W_u * V_i)
                self.assertTrue((W_u * V_n).is_zero_matrix)
                self.assertTrue((W_d * V_i).is_zero_matrix)
                self.assertEqual(sympy.eye(n - r), W_d * V_n)
                self.assertEqual(translator.inv_herm_mult, sympy.Matrix.vstack(W_u, W_d))

    def test_the_dependent_variable_inverse_is_the_inverse_each_time_it_is_asked(self):
        # The second call is answered from a cache.
        system = ODESystem.from_equations(HL_6_4)
        translator = ODETranslator(system)
        index = system.independent_variable_index
        for _ in range(2):
            self.assertEqual(sympy.eye(2), translator.dep_var_inv_herm_mult(index)
                             * translator.dep_var_herm_mult(index))

    def test_an_integer_inverse_needs_determinant_one(self):
        with self.assertRaises(ValueError) as caught:
            _int_inv(sympy.Matrix([[2, 0], [0, 1]]))
        self.assertIn('det 2', str(caught.exception))


class TestExtendingToAUnimodularMatrix(TestCase):
    '''
    `extend_rectangular_matrix` completes a matrix to a square one of determinant +-1,
    keeping the given matrix as its first columns, or rows when it is wider than tall.
    '''

    def test_a_tall_matrix_is_kept_as_the_first_columns(self):
        tall = sympy.Matrix([[1, 2], [0, 1], [3, 4]])
        extended = extend_rectangular_matrix(tall, check_unimodular=True)
        self.assertEqual((3, 3), extended.shape)
        self.assertEqual(tall, extended[:, :2])
        self.assertEqual(1, abs(extended.det()))

    def test_a_wide_matrix_is_kept_as_the_first_rows(self):
        wide = sympy.Matrix([[1, 2, 3], [0, 1, 4]])
        extended = extend_rectangular_matrix(wide, check_unimodular=True)
        self.assertEqual((3, 3), extended.shape)
        self.assertEqual(wide, extended[:2, :])
        self.assertEqual(1, abs(extended.det()))

    def test_a_square_matrix_is_returned_as_it_is(self):
        square = sympy.Matrix([[1, 2], [3, 5]])
        self.assertEqual(square, extend_rectangular_matrix(square))

    def test_a_matrix_with_no_such_extension_is_refused_when_asked_to_check(self):
        # Every 2x2 minor of this is even, so no integer column completes it to +-1.
        with self.assertRaises(ValueError) as caught:
            extend_rectangular_matrix(sympy.Matrix([[2, 0], [0, 1], [0, 0]]),
                                      check_unimodular=True)
        self.assertIn('Unable to extend', str(caught.exception))


class TestRewriteRulesAndTheMovingFrame(TestCase):
    '''
    Both are monomial maps read off the multiplier.  The rewrite rules send `z` to
    `z^(V_n W_d)`, and the moving frame is `z^(-V_i)`.
    '''

    def translators(self):
        for equations in (MICHAELIS_MENTEN, HL_6_4, HL_6_6):
            system = ODESystem.from_equations(equations)
            yield system, ODETranslator(system)

    def test_the_rewrite_is_a_projection_onto_the_invariants(self):
        # P = V_n W_d satisfies P P = P, because W_d V_n = I; and A P = 0, because the
        # columns of V_n are invariants.  So rewriting twice is rewriting once, and
        # everything it produces is invariant.
        for system, translator in self.translators():
            with self.subTest(system=str(system.symbols)):
                P = translator.herm_mult_n * translator.inv_herm_mult_d
                self.assertEqual(P, P * P)
                self.assertTrue((translator.scaling_matrix * P).is_zero_matrix)

    def test_the_rewrite_rules_leave_every_invariant_as_it_is(self):
        for system, translator in self.translators():
            with self.subTest(system=str(system.symbols)):
                rules = translator.rewrite_rules(system.symbols)
                self.assertEqual(set(system.symbols), set(rules))
                for invariant in translator.invariant_monomials():
                    rewritten = invariant.subs(rules, simultaneous=True)
                    self.assertEqual(0, sympy.simplify(rewritten - invariant))
                for image in rules.values():
                    self.assertTrue(translator.is_invariant_expr(image))

    def test_the_moving_frame_undoes_the_scaling(self):
        # Scaling z by lambda^A scales z^(-V_i) by lambda^(-A V_i) = lambda^(-H), with H
        # the nonzero block of the Hermite form.
        for system, translator in self.translators():
            with self.subTest(system=str(system.symbols)):
                r = translator.r
                frame = translator.moving_frame(system.symbols)
                self.assertEqual(r, len(frame))
                self.assertEqual(-translator.herm_form[:, :r],
                                 translator.scaling_matrix * (-translator.herm_mult_i))
                for component, auxiliary in zip(frame, translator.auxiliary_monomials()):
                    self.assertEqual(1, sympy.simplify(component * auxiliary))

    def test_the_tex_of_a_translator_names_its_three_matrices(self):
        _, translator = next(self.translators())
        tex = translator.to_tex()
        self.assertLess(tex.index('A='), tex.index('V='))
        self.assertLess(tex.index('V='), tex.index('W='))
