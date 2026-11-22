'''
`ODETranslator.is_invariant_expr`, checked against the definition of invariance.

An expression is invariant under a scaling when its value at a point equals its value at
every rescaled image of that point.  That can be asked numerically, with no reference to
how the method decides, and is the oracle here.
'''
from unittest import TestCase

import sympy

from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator
from desr.sympy_helper import sympify_expression


def invariant_by_definition(expression, symbols, scaling_matrix):
    '''
    Evaluate at a point and at its image under each of several rescalings.  Exact
    rational arithmetic throughout, so equality means equality.
    '''
    point = {symbol: sympy.Rational(prime)
             for symbol, prime in zip(symbols, (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31))}
    here = expression.subs(point)
    for scales in ((2, 3, 5), (7, 1, 3), (1, 5, 2)):
        image = {}
        for column, symbol in enumerate(symbols):
            factor = sympy.Integer(1)
            for row in range(scaling_matrix.rows):
                factor *= sympy.Rational(scales[row]) ** scaling_matrix[row, column]
            image[symbol] = point[symbol] * factor
        if sympy.simplify(expression.subs(image) - here) != 0:
            return False
    return True


class TestIsInvariantExpr(TestCase):

    MICHAELIS_MENTEN = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                        'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
    PREDATOR_PREY = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )',
                     'dp/dt = s*p*(1 - h*p / n)']

    #: For Michaelis-Menten.  Chosen to include sums, quotients of sums, and expressions
    #: that are invariant under one of the two scalings and not the other.
    EXPRESSIONS = [
        's', 'c', 't', 'e_0', 'k_1',
        'k_2/k_m1', 'k_2*t', 'e_0*k_1*t', 's/e_0', 'c/s',
        'c + s', '(c + s)/e_0', 'k_2 + k_m1', '(k_2 + k_m1)*t', '(k_2 + k_m1)/(e_0*k_1)',
        'c/e_0 + k_2*t', 's/e_0 + t', 'k_1*s*t', 'k_1*s', '1', 'c*s/e_0**2 - k_m1/k_2',
        '(c - s)/(c + s)', 'e_0*k_1/(k_2 + k_m1) + s/c',
    ]

    def setUp(self):
        self.system = ODESystem.from_equations(self.MICHAELIS_MENTEN)
        self.translator = ODETranslator.from_ode_system(self.system)

    def test_it_agrees_with_the_definition(self):
        answers = set()
        for text in self.EXPRESSIONS:
            with self.subTest(expression=text):
                expression = sympify_expression(text)
                expected = invariant_by_definition(
                    expression, self.system.symbols, self.translator.scaling_matrix)
                self.assertEqual(expected, self.translator.is_invariant_expr(expression))
                answers.add(expected)
        # Both answers occur, so agreeing is not agreeing on a constant.
        self.assertEqual({True, False}, answers)

    def test_the_invariants_the_reduction_finds_are_invariant(self):
        for equations in (self.MICHAELIS_MENTEN, self.PREDATOR_PREY):
            system = ODESystem.from_equations(equations)
            translator = ODETranslator.from_ode_system(system)
            for monomial in translator.invariant_monomials():
                with self.subTest(monomial=str(monomial)):
                    self.assertTrue(translator.is_invariant_expr(monomial))

    def test_an_auxiliary_is_never_invariant(self):
        # The auxiliaries are the directions the scaling moves; that is what they are.
        for equations in (self.MICHAELIS_MENTEN, self.PREDATOR_PREY):
            system = ODESystem.from_equations(equations)
            translator = ODETranslator.from_ode_system(system)
            for monomial in translator.auxiliary_monomials():
                with self.subTest(monomial=str(monomial)):
                    self.assertFalse(translator.is_invariant_expr(monomial))

    def test_any_rational_function_of_invariants_is_invariant(self):
        # Hubert & Labahn's Theorem 4.2 runs the other way -- every rational invariant is
        # a rational function of the generating ones -- and this is the easy direction.
        a, b, c, d = list(self.translator.invariant_monomials())[:4]
        for expression in (a + b, a * b - c, (a + 2 * b) / (1 + c * d), a ** 3 / (b + d)):
            with self.subTest(expression=str(expression)):
                self.assertTrue(self.translator.is_invariant_expr(expression))

    def test_a_symbol_the_scaling_leaves_alone_is_invariant(self):
        # `q` enters only as `q - 1`, so no scaling can move it: its column is zero.
        system = ODESystem.from_equations(['dp/dt = p*(q - 1)'])
        translator = ODETranslator.from_ode_system(system)
        column = list(system.symbols).index(sympy.Symbol('q'))
        self.assertTrue(translator.scaling_matrix[:, column].is_zero_matrix)
        self.assertTrue(translator.is_invariant_expr('q'))
        self.assertTrue(translator.is_invariant_expr('q**2 + 1'))

    def test_a_string_and_an_expression_are_the_same_question(self):
        self.assertEqual(self.translator.is_invariant_expr('(c + s)/e_0'),
                         self.translator.is_invariant_expr(sympify_expression('(c + s)/e_0')))

    def test_a_symbol_the_scaling_does_not_know_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            self.translator.is_invariant_expr('s/q')
        self.assertIn('q', str(caught.exception))
        self.assertIn('does not act on', str(caught.exception))

    def test_a_translator_built_from_a_bare_matrix(self):
        bare = ODETranslator(self.translator.scaling_matrix)
        with self.assertRaises(ValueError) as caught:
            bare.is_invariant_expr('s')
        self.assertIn('variable_order', str(caught.exception))

        # Given the order, it answers as the system-built one does.
        order = self.system.symbols
        for text in ('s', 'k_2/k_m1', '(c + s)/e_0'):
            with self.subTest(expression=text):
                self.assertEqual(self.translator.is_invariant_expr(text),
                                 bare.is_invariant_expr(text, variable_order=order))

    def test_an_order_of_the_wrong_length_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            self.translator.is_invariant_expr('s', variable_order=self.system.symbols[:3])
        self.assertIn('columns', str(caught.exception))
