'''
The `multiplier_*` operations, told which columns by their monomials.

A column of the Hermite multiplier is an exponent vector, which is a monomial in the
system's symbols.  Adding column j to column i, `power` times, multiplies monomial i by
monomial j to that power; negating a column inverts its monomial; swapping two exchanges
them.  Those are the definitions the tests hold the operations to, along with the two
facts that make an operation admissible: the multiplier stays unimodular, and its
product with the scaling matrix does not change.
'''
from unittest import TestCase

import sympy

from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator

MICHAELIS_MENTEN = ['dE/dt = - k_1*E*S + k_m1*C + k_2*C',
                    'dS/dt = - k_1*E*S + k_m1*C',
                    'dC/dt = k_1*E*S - k_m1*C - k_2*C',
                    'dP/dt = k_2*C']


class TestColumnsByMonomial(TestCase):

    def setUp(self):
        self.system = ODESystem.from_equations(MICHAELIS_MENTEN)
        self.symbols = {str(symbol): symbol for symbol in self.system.symbols}

    def translator(self):
        return ODETranslator(self.system)

    def monomials(self, translator):
        return list(translator.forward_monomials())

    def assert_admissible(self, before, after):
        self.assertEqual(1, abs(after.herm_mult.det()))
        self.assertEqual(before.scaling_matrix * before.herm_mult,
                         after.scaling_matrix * after.herm_mult)

    def two_invariants(self, translator):
        ''' Two invariant columns, as (index, monomial) pairs. '''
        r = translator.r
        monomials = self.monomials(translator)
        return (r, monomials[r]), (r + 2, monomials[r + 2])

    def test_adding_multiplies_the_monomials(self):
        for power in (1, -1, 2, -3):
            with self.subTest(power=power):
                reference, translator = self.translator(), self.translator()
                (i, first), (j, second) = self.two_invariants(translator)
                translator.multiplier_add_columns(first, second, power=power)

                after = self.monomials(translator)
                self.assertEqual(sympy.simplify(first * second ** power),
                                 sympy.simplify(after[i]))
                # Every other column is as it was.
                for column, monomial in enumerate(self.monomials(reference)):
                    if column != i:
                        self.assertEqual(monomial, after[column])
                self.assert_admissible(reference, translator)

    def test_negating_inverts_the_monomial(self):
        reference, translator = self.translator(), self.translator()
        (i, first), _ = self.two_invariants(translator)
        translator.multiplier_negate_column(first)
        self.assertEqual(sympy.simplify(1 / first),
                         sympy.simplify(self.monomials(translator)[i]))
        self.assert_admissible(reference, translator)

    def test_swapping_exchanges_them(self):
        reference, translator = self.translator(), self.translator()
        (i, first), (j, second) = self.two_invariants(translator)
        translator.multiplier_swap_columns(first, second)
        after = self.monomials(translator)
        self.assertEqual((second, first), (after[i], after[j]))
        self.assert_admissible(reference, translator)

    def test_an_expression_a_string_and_an_index_are_one_column(self):
        # Three ways to say it, one multiplier afterwards.
        by_expression, by_string, by_index = (self.translator() for _ in range(3))
        (i, first), (j, second) = self.two_invariants(by_expression)

        by_expression.multiplier_add_columns(first, second, power=-1)
        by_string.multiplier_add_columns(str(first), str(second), power=-1)
        by_index.multiplier_add_columns(i, j, power=-1)

        self.assertEqual(by_index.herm_mult, by_expression.herm_mult)
        self.assertEqual(by_index.herm_mult, by_string.herm_mult)

    def test_the_monomials_can_be_written_in_the_systems_own_symbols(self):
        # Not taken from the translator: typed, as a user with the model in front of
        # them would.
        translator = self.translator()
        k_1, k_2, P, C = (self.symbols[name] for name in ('k_1', 'k_2', 'P', 'C'))
        monomials = self.monomials(translator)
        self.assertIn(P * k_1 / k_2, monomials)
        self.assertIn(C * k_1 / k_2, monomials)

        translator.multiplier_add_columns(P * k_1 / k_2, C * k_1 / k_2, power=-1)
        self.assertIn(P / C, [sympy.simplify(m) for m in self.monomials(translator)])

    def test_a_column_is_named_by_what_it_holds_now(self):
        # After an operation the column's monomial has changed, and the new one is its
        # name.  The old one no longer names anything.
        translator = self.translator()
        (i, first), (j, second) = self.two_invariants(translator)
        translator.multiplier_add_columns(first, second)
        product = sympy.simplify(first * second)

        translator.multiplier_negate_column(product)
        self.assertEqual(sympy.simplify(1 / product),
                         sympy.simplify(self.monomials(translator)[i]))
        with self.assertRaises(ValueError):
            translator.multiplier_negate_column(first)

    def test_a_monomial_that_is_no_column_is_refused_with_the_ones_there_are(self):
        translator = self.translator()
        E, S = self.symbols['E'], self.symbols['S']
        with self.assertRaises(ValueError) as caught:
            translator.multiplier_negate_column(E * S)
        message = str(caught.exception)
        self.assertIn('E*S', message)
        for monomial in self.monomials(translator):
            self.assertIn(str(monomial), message)

    def test_an_index_out_of_range_is_refused(self):
        translator = self.translator()
        columns = translator.herm_mult.cols
        for index in (columns, -columns - 1):
            with self.subTest(index=index):
                with self.assertRaises(ValueError) as caught:
                    translator.multiplier_negate_column(index)
                self.assertIn('out of range', str(caught.exception))

    def test_a_negative_index_counts_from_the_end(self):
        by_negative, by_positive = self.translator(), self.translator()
        columns = by_negative.herm_mult.cols
        by_negative.multiplier_negate_column(-1)
        by_positive.multiplier_negate_column(columns - 1)
        self.assertEqual(by_positive.herm_mult, by_negative.herm_mult)

    def test_true_is_not_a_column(self):
        # A bool is an int in Python, and `True` would otherwise be column 1.
        with self.assertRaises(ValueError):
            self.translator().multiplier_negate_column(True)

    def test_an_auxiliary_column_still_cannot_be_added(self):
        # Naming it by its monomial does not make it admissible.
        translator = self.translator()
        auxiliary = list(translator.auxiliary_monomials())[0]
        (_, invariant), _ = self.two_invariants(translator)
        with self.assertRaises(ValueError) as caught:
            translator.multiplier_add_columns(invariant, auxiliary)
        self.assertIn('auxiliary column', str(caught.exception))
