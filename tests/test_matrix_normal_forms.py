from unittest import TestCase

import sympy

from desr.matrix_normal_forms import (is_hnf_row, hnf_row_lll, is_hnf_col,
                                      is_normal_hermite_multiplier, normal_hnf_col)


class TestHermiteMethods(TestCase):

    def test_example_sage(self):
        ''' Example from the Sage website
            See: http://doc.sagemath.org/html/en/reference/matrices/sage/matrix/matrix_integer_dense.html
        '''
        A = sympy.Matrix([[1,2],[3,4]])
        H, V = hnf_row_lll(A)
        self.assertEqual(H, sympy.Matrix([[1,0],[0,2]]))
        self.assertEqual(H, V * A)
        self.assertTrue(is_hnf_row(H))

        A = sympy.Matrix(range(25)).reshape(5,5)
        H, V = hnf_row_lll(A)
        H_A = sympy.Matrix.vstack(sympy.Matrix([[5,   0,  -5, -10, -15],
                                                [0,   1,   2,   3,   4]]),
                                  sympy.zeros(3, 5))
        self.assertEqual(H, H_A)
        self.assertEqual(H, V * A)
        self.assertTrue(is_hnf_row(H))

        ##TODO This test case is broken, the answer should be [[0, 1]]
        # A = sympy.Matrix([[0, -1]])
        # H, V = hnf_row_lll(A)
        # self.assertEqual(H, sympy.Matrix([[0, -1]]))
        # self.assertEqual(H, V * A)
        # self.assertTrue(is_hnf_row(H))

    def test_canonicalise_hnf_row(self):
        ''' _canonicalise_hnf_row turns a valid but non-canonical row HNF into the unique
            canonical one (positive pivots; above-pivot entries reduced into [0, pivot)),
            while preserving U * A == H.  This guards the post-pass that makes hnf_row_lll
            backend-independent: FLINT output is already canonical, but other HNF algorithms
            (e.g. an LLL-based Hermite method) can return a non-canonical triangularisation.
        '''
        from desr.matrix_normal_forms import _canonicalise_hnf_row
        A = sympy.Matrix([[3, 1, 4], [1, 5, 9], [2, 6, 5]])
        H, U = hnf_row_lll(A)
        self.assertTrue(is_hnf_row(H))

        # Perturb into a valid but non-canonical pair by adding multiples of lower pivot
        # rows to the rows above them (valid unimodular row ops), pushing above-pivot
        # entries out of [0, pivot).
        H_bad, U_bad = H.copy(), U.copy()
        H_bad[0, :] += 3 * H[2, :]; U_bad[0, :] += 3 * U[2, :]
        H_bad[1, :] -= 2 * H[2, :]; U_bad[1, :] -= 2 * U[2, :]
        H_bad[0, :] -= 5 * H[1, :]; U_bad[0, :] -= 5 * U[1, :]
        self.assertEqual(U_bad * A, H_bad)       # still a valid factorisation
        self.assertFalse(is_hnf_row(H_bad))      # but no longer canonical

        _canonicalise_hnf_row(H_bad, U_bad)
        self.assertTrue(is_hnf_row(H_bad))       # now canonical
        self.assertEqual(H_bad, H)               # equals the unique canonical HNF
        self.assertEqual(U_bad, U)               # row ops mirrored -> multiplier restored
        self.assertEqual(U_bad * A, H_bad)       # invariant U * A == H preserved

    def test_example1(self):
        ''' Example from Extended gcd and Hermite nomal form via lattice basis reduction '''

        A = sympy.Matrix(10, 10, lambda i, j: (i + 1) ** 3 * (j + 1) ** 2 + i + j + 2)

        H_answer_paper = sympy.Matrix([[1, 0, 7, 22, 45, 76, 115, 162, 217, 280],
                                       [0, 1, 4, 9, 16, 25, 36, 49, 64, 81],
                                       [0, 0, 12, 36, 72, 120, 180, 252, 336, 432]])
        H_answer_webcalc = sympy.Matrix([[1, 0, 7, 22, 45, 76, 115, 162, 217, 280],
                                          [0, 1, 4, 9, 16, 25, 36, 49, 64, 81],
                                          [0, 0, 12, 36, 72, 120, 180, 252, 336, 432]])

        self.assertTrue(is_hnf_row(H_answer_webcalc))
        self.assertTrue(is_hnf_row(sympy.Matrix.vstack(H_answer_webcalc, sympy.zeros(7, 10))))

        V_answer_paper = sympy.Matrix([[-10, -8, -5, 1, 2, 3, 5, 3, 0, -4],
                                       [-2, -1, 0, 1, -1, 0, 1, 0, 1, -1],
                                       [-15, -11, -4, 0, 4, 5, 4, 3, 1, -5],
                                       [1, -1, -1, 0, 2, -1, 0, 0, 0, 0],
                                       [0, 1, -1, -1, 1, -1, 2, -1, 0, 0],
                                       [1, 0, -1, -1, -1, 2, 0, 1, -1, 0],
                                       [1, 0, -1, 1, -1, 1, -1, 1, 1, -1],
                                       [-1, 0, 1, 0, 1, 1, -1, -2, 0, 1],
                                       [1, -1, 0, -1, 1, 0, 0, -1, 2, -1],
                                       [1, -2, 1, 1, -2, 0, 2, -1, 0, 0]])

        V_answer_webcalc = sympy.Matrix([[-10, -8, -5, 1, 2, 3, 5, 3, 0, -4],
                                         [-2, -1, 0, 1, -1, 0, 1, 0, 1, -1],
                                         [-15, -11, -4, 0, 4, 5, 4, 3, 1, -5],
                                         [-1, 2, -1, -1, 2, 0, -2, 1, 0, 0],
                                         [-1, 1, 0, 1, -1, 0, 0, 1, -2, 1],
                                         [1, 0, -1, 0, -1, -1, 1, 2, 0, -1],
                                         [-1, 0, 2, -1, 1, -1, 1, -1, -1, 1],
                                         [-1, 0, 1, 1, 1, -2, 0, -1, 1, 0],
                                         [0, -1, 1, 1, -1, 1, -2, 1, 0, 0],
                                         [-1, 1, 1, 0, -2, 1, 0, 0, 0, 0]])

        H, V = hnf_row_lll(A)
        H_nz, H_z = H[:3, :], H[3:, :]
        self.assertTrue(H_z.is_zero_matrix)
        self.assertEqual(H_nz, H_answer_webcalc)
        self.assertEqual(H, V * A)

    def test_wiki_example(self):
        ''' Test the examples from wikipedia
            https://en.wikipedia.org/wiki/Hermite_normal_form
        '''

        A1 = sympy.Matrix([[3,3,1,4],
                           [0,1,0,0],
                           [0,0,19,16],
                           [0,0,0,3]])
        H1 = sympy.Matrix([[3,0,1,1],
                           [0,1,0,0],
                           [0,0,19,1],
                           [0,0,0,3]])

        A2 = sympy.Matrix([[5, 0, 1, 4],
                           [0, -1, -4, 99],
                           [0, 20, 19, 16],
                           [0, 0, 2, 1],
                           [0, 0, 0, 3]])
        A2 = sympy.Matrix.hstack(sympy.zeros(5, 2), A2)
        A2 = sympy.Matrix.vstack(A2, sympy.zeros(1, 6))
        H2 = sympy.Matrix([[5, 0, 0, 2],
                           [0, 1, 0, 1],
                           [0, 0, 1, 2],
                           [0, 0, 0, 3]])
        H2 = sympy.Matrix.hstack(sympy.zeros(4, 2), H2)
        H2 = sympy.Matrix.vstack(H2, sympy.zeros(2, 6))

        A3 = sympy.Matrix([[2, 3, 6, 2],
                          [5, 6, 1, 6],
                          [8, 3, 1, 1]])
        H3 = sympy.Matrix([[1, 0, 50, -11],
                           [0, 3, 28,-2],
                           [0, 0, 61, -13]])

        for A, H in ((A1, H1), (A2, H2), (A3, H3)):
            H_calc, V = hnf_row_lll(A)
            self.assertTrue(is_hnf_row(H))
            self.assertEqual(H, H_calc)
            self.assertEqual(H, V * A)

    def test_normal_hermite_multiplier_example(self):
        ''' Example from Hubert Labahn '''
        # Broken due to different definitions of HNF
        A = sympy.Matrix([[8, 2, 15, 9, 11],
                          [6, 0, 6, 2, 3]])

        H_answer = sympy.Matrix.hstack(sympy.eye(2), sympy.zeros(2, 3))

        self.assertTrue(is_hnf_col(H_answer))

        H, V = normal_hnf_col(A)
        self.assertEqual(H, H_answer)
        self.assertTrue(is_hnf_col(H))

        # Check V_n is in column hnf.
        self.assertTrue(is_hnf_col(V[:, 2:]))
        self.assertTrue(is_normal_hermite_multiplier(V, A))

        # We can't compare Hermite mmultipliers since we are using different definitions
        V_answer = sympy.Matrix([[-1, -2, -2, -2, -1],
                                 [-3, -14, -7, -13, -7],
                                 [1, 1, 2, 1, 0],
                                 [0, 2, 0, 3, 0],
                                 [0, 1, 0, 0, 2]])
