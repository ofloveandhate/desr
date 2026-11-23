'''
The normal forms held to their definitions, over integer matrices drawn from a seed.

`smith_normal_form` and `hnf_col_lll` were reached only by docstring examples, which
compare printed output.  Here what each returns is checked against what it is defined to
be, for matrices nobody chose.
'''
import random
from unittest import TestCase

import sympy

from desr.matrix_normal_forms import (hnf_col_lll, is_hnf_col, is_hnf_row,
                                      is_non_negative, is_smith_normal_form,
                                      smith_normal_form)


def random_matrix(seed):
    ''' An integer matrix of between one and four rows and columns, some of low rank. '''
    draw = random.Random(seed)
    rows, columns = draw.randint(1, 4), draw.randint(1, 4)
    matrix = sympy.Matrix(rows, columns, lambda i, j: draw.randint(-6, 6))
    if rows > 1 and draw.random() < 0.3:
        # A repeated row, so that some of these are rank deficient.
        matrix[rows - 1, :] = 2 * matrix[0, :]
    return matrix


SEEDS = range(300)


class TestSmithNormalForm(TestCase):
    '''
    `S = U M V` with `U` and `V` unimodular, and `S` diagonal with non-negative entries
    each dividing the next.
    '''

    def test_by_definition(self):
        ranks = set()
        for seed in SEEDS:
            with self.subTest(seed=seed):
                matrix = random_matrix(seed)
                ranks.add(matrix.rank())
                S, U, V = smith_normal_form(matrix)
                self.assertEqual(S, U * matrix * V)
                self.assertEqual(1, abs(U.det()))
                self.assertEqual(1, abs(V.det()))
                self.assertEqual(matrix.shape, S.shape)
                diagonal = [S[i, i] for i in range(min(S.shape))]
                for i in range(S.rows):
                    for j in range(S.cols):
                        if i != j:
                            self.assertEqual(0, S[i, j])
                self.assertTrue(all(entry >= 0 for entry in diagonal))
                for first, second in zip(diagonal, diagonal[1:]):
                    if first == 0:
                        self.assertEqual(0, second)
                    else:
                        self.assertEqual(0, second % first)
                self.assertTrue(is_smith_normal_form(S))
        self.assertLessEqual({1, 2, 3}, ranks)

    def test_the_rank_is_the_number_of_nonzero_diagonal_entries(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                matrix = random_matrix(seed)
                S = smith_normal_form(matrix)[0]
                self.assertEqual(matrix.rank(),
                                 sum(1 for i in range(min(S.shape)) if S[i, i] != 0))

    def test_what_is_not_a_smith_normal_form(self):
        self.assertFalse(is_smith_normal_form(sympy.diag(2, 1, 1)))    # 2 does not divide 1
        self.assertFalse(is_smith_normal_form(sympy.diag(-1, 1, 2)))   # negative entry
        self.assertFalse(is_smith_normal_form(sympy.diag(0, 1)))       # zero before nonzero
        self.assertFalse(is_smith_normal_form(sympy.Matrix([[1, 1], [0, 1]])))  # off diagonal


class TestColumnHermiteNormalFormByLLL(TestCase):
    ''' `H = M V` with `V` unimodular and `H` in column Hermite normal form. '''

    def test_by_definition(self):
        checked = 0
        for seed in SEEDS:
            matrix = random_matrix(seed)
            if matrix.rows < 2 or matrix.is_zero_matrix:
                continue  # its stated domain: nonzero, at least two rows
            with self.subTest(seed=seed):
                H, V = hnf_col_lll(matrix)
                self.assertEqual(H, matrix * V)
                self.assertEqual(1, abs(V.det()))
                self.assertTrue(is_hnf_col(H))
                checked += 1
        self.assertGreaterEqual(checked, 60)

    def test_a_row_form_is_the_transpose_of_a_column_form(self):
        for seed in SEEDS:
            matrix = random_matrix(seed)
            if matrix.rows < 2 or matrix.is_zero_matrix:
                continue
            with self.subTest(seed=seed):
                H, _ = hnf_col_lll(matrix)
                self.assertEqual(is_hnf_col(H), is_hnf_row(H.T))


class TestSmallMatrixHelpers(TestCase):

    def test_non_negative(self):
        self.assertTrue(is_non_negative(sympy.Matrix([[0, 1], [2, 3]])))
        self.assertFalse(is_non_negative(sympy.Matrix([[0, 1], [-2, 3]])))
