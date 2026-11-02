'''
The worked examples of the paper this library implements, asserted against its printed
values.

Reference: E. Hubert and G. Labahn, "Scaling Invariants and Symmetry Reduction of
Dynamical Systems", Foundations of Computational Mathematics 13(4):479-516, 2013.

These are **characterisation tests**: they say what desr must produce for a reader to be
able to compare its output with the paper.  Until the Hermite normal form convention was
corrected, desr reported the mirror image of the paper's normal form.  The two
conventions, in full:

===============  ================================  ============  ===================
                 pivot of a column                 H             reduce on pivot row
===============  ================================  ============  ===================
HL Def 2.1       **last** nonzero, zeros below it   upper tri.    to the **right**
Havas-Majewski-  **first** nonzero, zeros above it  lower tri.    to the **left**
Matthews
===============  ================================  ============  ===================

desr used to implement the second.  It arrived not as a choice about output but as a side effect:
the original Hermite implementation was Havas-Majewski-Matthews' lattice-reduction method,
whose own paper states that convention, and the library adopted it along with the
algorithm.  That backend has since been replaced by FLINT, so the convention is now a
vestige of a removed dependency, while the theory desr implements -- and the output a
reader compares against -- is Hubert & Labahn's.

Nothing in the library could detect the difference while it lasted: condition (c) of
Prop 2.3(ii) was checked over the first `r` rows rather than over `V_n`'s pivot rows, so
`is_normal_hermite_multiplier` rejected the paper's own published multiplier while
accepting multipliers that were not normal.  Fixing that was the companion to this change.
'''
from unittest import TestCase

import sympy

from desr.matrix_normal_forms import is_hnf_col, is_normal_hermite_multiplier, normal_hnf_col
from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator

t, z1, z2 = sympy.symbols('t z1 z2')


class TestExample24(TestCase):
    '''
    Example 2.4, p. 486, with the normal Hermite multiplier of Prop 2.3, pp. 484-485.
    '''

    A = sympy.Matrix([[8, 2, 15, 9, 11],
                      [6, 0, 6, 2, 3]])

    #: The paper's normal Hermite multiplier.  Prop 2.3(ii) makes the normal multiplier
    #: *unique*, so satisfying (a), (b) and (c) is itself proof that this is the printed
    #: one -- there is nothing else it could be.  Verified clause by clause: V_n's columns
    #: have pivots at rows 2, 3, 4 with values 2, 3, 2, zeros below each, and V_i reduced
    #: into [0, pivot) on exactly those rows.
    V = sympy.Matrix([[-1, -2, -2, -2, -1],
                      [-3, -14, -7, -13, -7],
                      [1, 1, 2, 1, 0],
                      [0, 2, 0, 3, 0],
                      [0, 1, 0, 0, 2]])

    H = sympy.Matrix.hstack(sympy.eye(2), sympy.zeros(2, 3))

    def test_the_paper_multiplier_is_a_valid_hermite_multiplier(self):
        # True under either convention, and it is what makes the disagreement a matter of
        # normal form rather than of correctness.
        self.assertEqual(abs(self.V.det()), 1)
        self.assertEqual(self.A * self.V, self.H)

    def test_desr_computes_the_paper_multiplier(self):
        H, V = normal_hnf_col(self.A)
        self.assertEqual(H, self.H)
        self.assertEqual(V, self.V)

    def test_the_paper_multiplier_is_recognised_as_normal(self):
        # This failed twice over before: the convention, and the broken condition (c).
        self.assertTrue(is_normal_hermite_multiplier(self.V, self.A))

    def test_the_invariant_columns_are_in_hermite_normal_form(self):
        # Condition (b) of Prop 2.3(ii).
        self.assertTrue(is_hnf_col(self.V[:, 2:]))


class TestExample66(TestCase):
    '''
    Example 6.6, pp. 503-504 -- the one the library's own examples and documentation
    reproduce, which makes it the test that matters most.

    The paper prints the scaling matrix [3 -1 5] over (t, z1, z2), the invariants
    y1 = t*z1**3 and y2 = z1**2*z2/t, and the auxiliary variable x = t*z1**2.
    '''

    EQUATIONS = ['dz1/dt = z1/t * (-2/3 + 1/3 *z1**5 *z2 )',
                 'dz2/dt = z2/t * ( 10/3 - 2/3 *z1**5 *z2 + z1**2*z2/t )']

    SCALING_MATRIX = sympy.Matrix([[3, -1, 5]])

    INVARIANTS = (t * z1 ** 3, z1 ** 2 * z2 / t)

    AUXILIARIES = (t * z1 ** 2,)

    def setUp(self):
        self.system = ODESystem.from_equations(self.EQUATIONS)
        self.assertEqual(self.system.symbols, (t, z1, z2))

    def test_the_scaling_matrix_matches(self):
        # This one already agrees.  A single-row matrix has no cross-row conditions to
        # violate, so it happens to satisfy both conventions at once -- which is part of
        # why the disagreement stayed hidden for so long.
        self.assertEqual(self.system.maximal_scaling_matrix(), self.SCALING_MATRIX)

    def test_the_invariants_match(self):
        # desr used to give (t*z1**3, z1**5*z2).  The first agreed; the second was the
        # paper's y2 multiplied by y1, which is why the documentation used to reach the
        # paper's answer only after a hand-applied multiplier_add_columns(2, 1, -1).
        translator = ODETranslator.from_ode_system(self.system)
        self.assertEqual(tuple(translator.invariant_monomials()), self.INVARIANTS)

    def test_the_auxiliary_variable_matches(self):
        translator = ODETranslator.from_ode_system(self.system)
        self.assertEqual(tuple(translator.auxiliary_monomials()), self.AUXILIARIES)

    def test_no_hand_applied_column_operations_are_needed(self):
        # The point of the change: the paper's answer should be what desr says by default.
        translator = ODETranslator.from_ode_system(self.system)
        self.assertEqual(tuple(translator.auxiliary_monomials()) + tuple(translator.invariant_monomials()),
                         self.AUXILIARIES + self.INVARIANTS)


class TestTheConventionItself(TestCase):
    '''
    The convention stated directly, so that a future reader can see which one desr follows
    without reconstructing it from an example.
    '''

    #: Hubert & Labahn Def 2.1: in each column the pivot is the last nonzero entry, so the
    #: entries below it are zero and the leading block is upper triangular.
    HL_FORM = sympy.Matrix([[1, -1],
                            [3, 2],
                            [0, 1]])

    #: The mirror, which is what desr used to report: the pivot is the first nonzero entry,
    #: the entries above it are zero, and the leading block is lower triangular.
    HMM_FORM = sympy.Matrix([[1, 0],
                             [3, 5],
                             [0, 1]])

    def test_the_two_forms_span_the_same_lattice(self):
        # Both are bases of the kernel of [3 -1 5], so the choice between them changes
        # which invariants are reported, never which invariants exist.
        A = sympy.Matrix([[3, -1, 5]])
        self.assertTrue((A * self.HL_FORM).is_zero_matrix)
        self.assertTrue((A * self.HMM_FORM).is_zero_matrix)
        self.assertEqual(sympy.Matrix.hstack(self.HL_FORM, self.HMM_FORM).rank(), 2)

    def test_desr_reports_the_paper_form(self):
        self.assertTrue(is_hnf_col(self.HL_FORM))

    def test_desr_no_longer_accepts_the_mirror(self):
        # Both are valid bases of the same lattice, so this is not a judgement about
        # correctness -- only about which of the two canonical forms desr reports.
        self.assertFalse(is_hnf_col(self.HMM_FORM))
