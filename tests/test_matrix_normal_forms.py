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

        # This case was commented out with a note that it was broken and the answer should
        # be [[0, 1]].  It is right now: a canonical row HNF has positive pivots, so a
        # negative one is negated and the factorisation absorbs the sign.  The
        # canonicalising post-pass added with the FLINT backend is what fixed it.
        A = sympy.Matrix([[0, -1]])
        H, V = hnf_row_lll(A)
        self.assertEqual(H, sympy.Matrix([[0, 1]]))
        self.assertEqual(H, V * A)
        self.assertTrue(is_hnf_row(H))

    def test_a_negative_pivot_is_made_positive(self):
        ''' The case above, generalised: a canonical row HNF has positive pivots, whatever
            the sign of the input.  This was broken before the canonicalising post-pass, and
            silently so -- a non-canonical answer is still a valid factorisation, so only
            the pivot sign gives it away.
        '''
        for entries in ([[0, -1]], [[-1, 0]], [[0, 0, -3]], [[-2, -4]],
                        [[0, -1], [0, 0]]):
            with self.subTest(matrix=entries):
                A = sympy.Matrix(entries)
                H, V = hnf_row_lll(A)
                self.assertTrue(is_hnf_row(H))
                self.assertEqual(H, V * A)
                self.assertIn(V.det(), (1, -1))  # unimodular, so the lattice is unchanged

    def test_canonicalise_hnf_row(self):
        ''' _canonicalise_hnf_row_hmm turns a valid but non-canonical row HNF into the
            unique canonical one (positive pivots; above-pivot entries reduced into
            [0, pivot)), while preserving U * A == H.  This guards the post-pass that makes
            the Hermite machinery backend-independent: FLINT output is already canonical,
            but other HNF algorithms (e.g. an LLL-based Hermite method) can return a
            non-canonical triangularisation.

            It works in the Havas-Majewski-Matthews convention, which is what the FLINT
            path produces internally -- not the convention desr reports.  See
            desr.matrix_normal_forms._CONVENTION.
        '''
        from desr.matrix_normal_forms import (_canonicalise_hnf_row_hmm, _hnf_row_hmm,
                                              _is_hnf_row_hmm)
        A = sympy.Matrix([[3, 1, 4], [1, 5, 9], [2, 6, 5]])
        H, U = _hnf_row_hmm(A)
        self.assertTrue(_is_hnf_row_hmm(H))

        # Perturb into a valid but non-canonical pair by adding multiples of lower pivot
        # rows to the rows above them (valid unimodular row ops), pushing above-pivot
        # entries out of [0, pivot).
        H_bad, U_bad = H.copy(), U.copy()
        H_bad[0, :] += 3 * H[2, :]; U_bad[0, :] += 3 * U[2, :]
        H_bad[1, :] -= 2 * H[2, :]; U_bad[1, :] -= 2 * U[2, :]
        H_bad[0, :] -= 5 * H[1, :]; U_bad[0, :] -= 5 * U[1, :]
        self.assertEqual(U_bad * A, H_bad)       # still a valid factorisation
        self.assertFalse(_is_hnf_row_hmm(H_bad)) # but no longer canonical

        _canonicalise_hnf_row_hmm(H_bad, U_bad)
        self.assertTrue(_is_hnf_row_hmm(H_bad))  # now canonical
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
        ''' Hubert & Labahn 2013, example 2.4, p. 486.

            The paper's own normal Hermite multiplier is asserted in
            tests/test_hubert_labahn.py, where the difference of convention between desr
            and the paper is recorded.  What is checked here is that whatever desr computes
            satisfies the conditions of Prop 2.3(ii) under desr's own convention.
        '''
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

    def test_a_valid_but_unnormalised_multiplier_is_refused(self):
        ''' Hubert & Labahn Prop 2.3(i): a Hermite multiplier is unique up to right
            multiplication by [[I, 0], [M, U]], which is to say V_n -> V_n*U and
            V_i -> V_i + V_n*M.  Perturbing by such an M leaves a perfectly valid Hermite
            multiplier that is no longer the *normal* one, so it is what condition (c)
            exists to reject.

            Worth asserting explicitly: condition (c) was previously checked over the first
            r rows rather than V_n's pivot rows, and a predicate that is never shown to say
            "no" is not known to work.  On this example the old test accepted 6 of these 18
            perturbations -- every one along the third invariant column, whose pivot row is
            2 and so was never looked at.
        '''
        A = sympy.Matrix([[8, 2, 15, 9, 11],
                          [6, 0, 6, 2, 3]])
        r = A.rows
        H, V = normal_hnf_col(A)
        self.assertTrue(is_normal_hermite_multiplier(V, A))

        # Add a multiple of an invariant column onto an auxiliary column: the M block.
        for column in range(r):
            for kernel_column in range(r, V.cols):
                for multiple in (1, -1, 3):
                    with self.subTest(column=column, kernel_column=kernel_column,
                                      multiple=multiple):
                        perturbed = V.copy()
                        perturbed[:, column] += multiple * V[:, kernel_column]

                        # Still a valid Hermite multiplier: unimodular, same A*V.
                        self.assertEqual(abs(perturbed.det()), 1)
                        self.assertEqual(A * perturbed, H)
                        # The invariant columns are untouched, so (b) still holds ...
                        self.assertTrue(is_hnf_col(perturbed[:, r:]))
                        # ... and only (c) can tell it is not the normal multiplier.
                        self.assertFalse(is_normal_hermite_multiplier(perturbed, A))
