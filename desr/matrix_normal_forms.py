
##### Richard Tanburn
import sympy
import itertools

from flint import fmpz_mat as _fmpz_mat


def _to_fmpz(m):
    """sympy.Matrix (integer entries) → fmpz_mat"""
    return _fmpz_mat([[int(v) for v in row] for row in m.tolist()])


def _from_fmpz(m):
    """fmpz_mat → sympy.Matrix"""
    return sympy.Matrix(m.tolist())

#: Which Hermite normal form desr reports.
#:
#: Hubert & Labahn 2013, Definition 2.1: an integer matrix is in **column** Hermite normal
#: form when, writing :math:`i_1 < \\cdots < i_r` for the pivot rows,
#:
#: 1. the first :math:`r` columns are nonzero;
#: 2. :math:`h_{k,j} = 0` for :math:`k > i_j` -- so the pivot of column :math:`j` is its
#:    **last** nonzero entry, and the entries below it are zero;
#: 3. :math:`0 \\le h_{i_j,k} < h_{i_j,j}` for :math:`j < k` -- the pivot row is reduced to
#:    the **right** of the pivot, and the pivot is positive.
#:
#: Equivalently: the submatrix on the pivot rows and the first :math:`r` columns is upper
#: triangular, with each pivot row's maximum on the diagonal.  The row form is the same
#: definition with the row and column indices exchanged.
#:
#: desr used to report the mirror of this -- Havas, Majewski and Matthews' convention,
#: where the pivot is the *first* nonzero entry of a column, the entries above it are zero,
#: and the leading block is lower triangular.  That was not a decision about output: the
#: original Hermite routine was their lattice-reduction method, and the library took the
#: convention along with the algorithm.  FLINT has since replaced that backend, so the
#: convention was a vestige of a dependency no longer present, while the theory desr
#: implements -- and the output a reader compares against -- is Hubert & Labahn's.
#:
#: The two forms are related by reversing the order of the rows and of the nonzero columns,
#: which is how the functions below are built: the Havas-Majewski-Matthews machinery over
#: FLINT is untouched and reached through that reversal.  Zero rows stay at the bottom and
#: zero columns at the right under *both* conventions, so only the nonzero block reverses.
_CONVENTION = 'Hubert & Labahn 2013, Definition 2.1'


def _reversal(size):
    '''The permutation matrix that reverses the order of `size` rows or columns.'''
    return sympy.Matrix.hstack(*[sympy.eye(size)[:, size - 1 - i] for i in range(size)])


def _num_nonzero_rows(matrix_):
    '''How many rows of `matrix_` are nonzero, requiring that they come first.'''
    nonzero = [i for i in range(matrix_.rows) if any(matrix_[i, :])]
    if nonzero != list(range(len(nonzero))):
        raise ValueError('Zero row above a nonzero row in {}'.format(matrix_))
    return len(nonzero)


def _first_nonzero_indices(matrix_):
    '''
    For each row, the column index of its **first** nonzero entry.

    The pivot rule of the Havas-Majewski-Matthews convention, used by the machinery that
    this module's public functions wrap.  See :data:`_CONVENTION`: what desr *reports* is
    Hubert & Labahn's mirror of it, for which see :func:`get_pivot_row_indices`.
    '''
    indices = []
    for row_ind in range(matrix_.shape[0]):
        row = matrix_[row_ind, :]
        for ind, col in enumerate(row):
            if col != 0:
                indices.append(ind)
                break
        # If we haven't found an index. throw an error, as we're assuming the rows are non-zero
        if len(indices) != row_ind + 1:
            raise ValueError('Zero row in {}'.format(matrix_))
    return indices


def get_pivot_row_indices(matrix_):
    ''' Return the pivot indices of the matrix

        The pivot of a row is its **last** nonzero entry; see :data:`_CONVENTION`.

        Args:
            matrix_ (sympy.Matrix): The matrix in question.

        Returns:
            indices (list): A list of indices, where the pivot entry of row i is [i, indices[i]]

        >>> matrix_ = sympy.eye(3)
        >>> get_pivot_row_indices(matrix_)
        [0, 1, 2]

        >>> matrix_ = sympy.Matrix([[1, 0, 0],[0, 0, 1]])
        >>> get_pivot_row_indices(matrix_)
        [0, 2]

        The pivot is the last nonzero entry, not the first, so a row with several nonzero
        entries has its pivot on the right.

        >>> get_pivot_row_indices(sympy.Matrix([[3, -1, 5]]))
        [2]
    '''
    indices = []
    for row_ind in range(matrix_.shape[0]):
        row = matrix_[row_ind, :]
        found = [ind for ind, col in enumerate(row) if col != 0]
        if not found:
            raise ValueError('Zero row in {}'.format(matrix_))
        indices.append(found[-1])
    return indices

def get_pivot_col_indices(matrix_):
    '''
    Return the pivot row of each column of a matrix in column Hermite normal form.

    The counterpart of :func:`get_pivot_row_indices`, and the one place that knows where a
    column's pivot sits.  Which entry of a column is the pivot is a property of the
    convention -- Havas-Majewski-Matthews take the first nonzero entry, Hubert & Labahn the
    last -- so anything that needs a pivot's position asks here rather than deciding for
    itself.

    Args:
        matrix_ (sympy.Matrix): A matrix in column Hermite normal form.

    Returns:
        indices (list): A list of indices, where the pivot entry of column j is
        [indices[j], j].

    >>> get_pivot_col_indices(sympy.eye(3))
    [0, 1, 2]

    >>> get_pivot_col_indices(sympy.Matrix([[1, 0], [0, 0], [0, 1]]))
    [0, 2]

    The pivot is the last nonzero entry of the column, so it sits at the bottom.

    >>> get_pivot_col_indices(sympy.Matrix([[1, -1], [3, 2], [0, 1]]))
    [1, 2]
    '''
    return get_pivot_row_indices(matrix_.T)

def is_hnf_col(matrix_):
    '''
    Decide whether a matrix is in row Hermite normal form, when acting on rows.

    Args:
        matrix_ (sympy.Matrix): The matrix in question.

    Returns:
        bool
    '''
    return is_hnf_row(matrix_.T)

def is_hnf_row(matrix_):
    ''' Decide whether a matrix is in Hermite normal form, when acting on rows.

        This is Hubert & Labahn 2013 Definition 2.1 with the row and column indices
        exchanged, as that definition's closing remark directs: the pivot of a row is its
        **last** nonzero entry, the pivot columns strictly ascend, the pivots are positive,
        and the entries **below** a pivot are reduced into :math:`[0, \\text{pivot})`.  See
        :data:`_CONVENTION`, which also records what desr used to report instead and why.

        Args:
            matrix_ (sympy.Matrix): The matrix in question.

        Returns:
            bool

        >>> is_hnf_row(sympy.eye(4))
        True
        >>> is_hnf_row(sympy.ones(2))
        False

        The pivot of a row is its last nonzero entry, so this is a valid form -- the 2 sits
        to the *left* of row 0's pivot, where nothing constrains it:

        >>> is_hnf_row(sympy.Matrix([[2, 1, 0], [0, 0, 1]]))
        True

        while an entry below a pivot must be reduced below it, and a pivot must be positive:

        >>> is_hnf_row(sympy.Matrix([[0, 1, 0], [0, 2, 1]]))
        False
        >>> is_hnf_row(sympy.Matrix([[0, -1, 0], [0, 0, 1]]))
        False
    '''
    try:
        r = _num_nonzero_rows(matrix_)
    except ValueError:
        return False      # a zero row above a nonzero one: not in any Hermite normal form
    if r == 0:
        return True

    # Reversing the rows and the columns of the nonzero block turns this convention into
    # the Havas-Majewski-Matthews one, and vice versa: a last nonzero entry becomes a
    # first, entries below a pivot become entries above it, and the pivot indices stay
    # ascending because both orders reverse together.  Only the nonzero block is reversed,
    # since zero rows belong at the bottom under either convention.
    return _is_hnf_row_hmm(matrix_[:r, :][::-1, ::-1])


def _is_hnf_row_hmm(matrix_):
    ''' Decide whether a matrix is in row Hermite normal form in the
        Havas-Majewski-Matthews convention: the pivot of a row is its **first** nonzero
        entry and the entries **above** a pivot are reduced.

        Not what desr reports -- see :data:`_CONVENTION` and :func:`is_hnf_row` -- but what
        the FLINT-backed machinery below produces, so the checks in this module's internals
        are written against it.

        >>> _is_hnf_row_hmm(sympy.eye(4))
        True
        >>> _is_hnf_row_hmm(sympy.ones(2))
        False
        >>> _is_hnf_row_hmm(sympy.Matrix([[1, 2, 0], [0, 1, 0]]))
        False
        >>> _is_hnf_row_hmm(sympy.Matrix([[1, 0, 0], [0, 2, 1]]))
        True
        >>> _is_hnf_row_hmm(sympy.Matrix([[1, 0, 0], [0, -2, 1]]))
        False
    '''
    m, n = matrix_.shape

    # The first r rows are non-zero. Determine r and check rows r through m are 0.
    row_is_zero = [all([i == 0 for i in row]) for row in matrix_.tolist()]
    r = m - sum(map(int, row_is_zero))
    for i in range(r, m):
        if not row_is_zero[i]:
            return False

    # If b_{i,j_i} are the indices of the first non-zero entry of row i then j_i are strictly ascending
    indices = _first_nonzero_indices(matrix_[:r, :])
    for i in range(len(indices) - 1):
        if indices[i] >= indices[i + 1]:
            return False

    # Check each of the pivots are positive
    for i in range(r):
        if matrix_[i, indices[i]] <= 0:
            return False

    # Check we are 0 to the right of the pivot
    for i in range(r):
        for k in range(i):
            if not (0 <= matrix_[k, indices[i]] < matrix_[i, indices[i]]):
                return False
    return True

def _canonicalise_hnf_row_hmm(hnf, unimodular_matrix):
    '''Canonicalise a row-style Hermite normal form in place, in the
    Havas-Majewski-Matthews convention -- the form the FLINT machinery produces, not the
    form desr reports.  See :data:`_CONVENTION`.

    Makes every pivot positive and reduces every entry *above* a pivot into the range
    ``[0, pivot)``, applying the same row operations to ``unimodular_matrix`` so that the
    relation ``unimodular_matrix * A == hnf`` is preserved.

    Some HNF algorithms return a valid unimodular triangularisation whose above-pivot
    entries are not fully reduced (e.g. an LLL-based Hermite method whose row sign-
    normalisation can un-reduce a previously reduced entry, with no final reduction sweep).
    Canonicalising here makes the output independent of the backend.  It is a no-op on output
    that is already in canonical HNF, such as FLINT's.

    Args:
        hnf (sympy.Matrix): A row-style HNF candidate (any zero rows are at the bottom).
        unimodular_matrix (sympy.Matrix): The accompanying unimodular row-operation matrix.
    '''
    num_nonzero_rows = sum(1 for row in hnf.tolist() if any(val != 0 for val in row))
    for piv_row, piv_col in enumerate(_first_nonzero_indices(hnf[:num_nonzero_rows, :])):
        if hnf[piv_row, piv_col] < 0:
            hnf[piv_row, :] *= -1
            unimodular_matrix[piv_row, :] *= -1
        pivot = hnf[piv_row, piv_col]
        for above in range(piv_row):
            quotient = hnf[above, piv_col] // pivot
            if quotient != 0:
                hnf[above, :] -= quotient * hnf[piv_row, :]
                unimodular_matrix[above, :] -= quotient * unimodular_matrix[piv_row, :]


def hnf_row_lll(matrix_):
    '''
    Compute the Hermite normal form, acts on the ROWS of a matrix.
    The Havas, Majewski, Matthews LLL method is used.
    We usually take :math:`\\alpha = \\frac{m_1}{n_1}`, with :math:`(m_1,n_1)=(1,1)` to get best results.

    Args:
        matrix_ (sympy.Matrix): Integer mxn matrix, nonzero, at least two rows

    Returns:
        tuple:
            hnf (sympy.Matrix): The row Hermite normal form of matrix\\_.

            unimod (sympy.Matrix): Unimodular matrix representing the row actions.

    :rtype: (sympy.Matrix, sympy.Matrix)

    >>> matrix_ = sympy.Matrix([[2, 0],
    ...                         [3, 3],
    ...                         [0, 0]])
    >>> result = hnf_row_lll(matrix_)
    >>> result[0]
    Matrix([
    [2, 0],
    [1, 3],
    [0, 0]])
    >>> result[1] * matrix_ == result[0]
    True

    >>> hnf_row_lll(sympy.Matrix([[0, -2, 0]]))
    (Matrix([[0, 2, 0]]), Matrix([[-1]]))

    >>> matrix_ = sympy.Matrix([[0, 1, 0, 1], [-1, 0, 1, 0]])
    >>> hnf_row_lll(matrix_)[0]
    Matrix([
    [-1, 0, 1, 0],
    [ 0, 1, 0, 1]])
    '''
    assert len(matrix_.shape) == 2

    m, n = matrix_.shape
    reverse_columns = _reversal(n)

    # Reversing the columns, computing in the Havas-Majewski-Matthews convention and then
    # reversing the nonzero rows back is what turns that convention into the one desr
    # reports; see :data:`_CONVENTION`.  The FLINT machinery underneath is untouched, so the
    # speed it brought is unaffected.  U' A J = H' gives (P U') A = P H' J, since J is its
    # own inverse, and P -- a permutation of the nonzero rows -- keeps U unimodular.
    hnf, unimodular_matrix = _hnf_row_hmm(matrix_ * reverse_columns)

    r = _num_nonzero_rows(hnf)
    reverse_rows = (sympy.diag(_reversal(r), sympy.eye(m - r)) if m > r
                    else _reversal(r)) if r else sympy.eye(m)
    hnf = reverse_rows * hnf * reverse_columns
    unimodular_matrix = reverse_rows * unimodular_matrix

    if not is_hnf_row(hnf):
        raise ValueError('{} not able to be put into row HNF. Output is:\n{}'.format(matrix_, hnf))
    assert (unimodular_matrix * matrix_) == hnf
    return hnf, unimodular_matrix

def _hnf_row_hmm(matrix_):
    '''
    The row Hermite normal form in the Havas-Majewski-Matthews convention, from FLINT.

    Not what desr reports; :func:`hnf_row_lll` reverses this into Hubert & Labahn's form.
    Kept separate because the reversal needs a convention to reverse *from*, and because
    :func:`normal_hnf_col` builds its multiplier by a construction that only works in this
    one.

    Returns:
        tuple: (hnf, unimodular) with ``unimodular * matrix_ == hnf``.
    '''
    h, u = _to_fmpz(matrix_).hnf(transform=True)   # H == U * A
    hnf = _from_fmpz(h)
    unimodular_matrix = _from_fmpz(u)

    if not abs(unimodular_matrix.det()) == 1:
        raise RuntimeError('Row operation matrix {} has determinant {}, not +-1'.format(unimodular_matrix, unimodular_matrix.det()))

    # Canonicalise into Hermite normal form (positive pivots, above-pivot entries in [0, pivot)).
    # This is what makes the result unique and independent of the HNF backend; see
    # _canonicalise_hnf_row_hmm for why some backends return a non-canonical triangularisation.
    _canonicalise_hnf_row_hmm(hnf, unimodular_matrix)

    assert _is_hnf_row_hmm(hnf)
    assert (unimodular_matrix * matrix_) == hnf
    return hnf, unimodular_matrix

def hnf_col_lll(matrix_):
    '''
    Compute the Hermite normal form, acts on the COLUMNS of a matrix.
    The Havas, Majewski, Matthews LLL method is used.
    We usually take :math:`\\alpha = \\frac{m_1}{n_1}`, with :math:`(m_1,n_1)=(1,1)` to get best results.

    Args:
        matrix_ (sympy.Matrix): Integer mxn matrix, nonzero, at least two rows

    Returns:
        tuple:
            hnf (sympy.Matrix): The column Hermite normal form of matrix\\_.
            unimod (sympy.Matrix): Unimodular matrix representing the column actions.

    :rtype: (sympy.Matrix, sympy.Matrix)

    >>> A = sympy.Matrix([[8, 2, 15, 9, 11],
    ...                   [6, 0, 6, 2, 3]])
    >>> h, v = hnf_col(A)
    >>> h
    Matrix([
    [1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0]])
    >>> v
    Matrix([
    [ 1,  -7, -2,  -28,  5],
    [ 4, -28, -7, -104, 20],
    [-1,   7,  2,   27, -5],
    [ 0,   2,  0,    3, -3],
    [ 0,  -1,  0,    0,  2]])
     >>> A * v == h
     True
    '''
    hnf, unimod = hnf_row_lll(matrix_.T)
    return hnf.T, unimod.T


## Normal Hermite multiplier from Hubert-Labahn
def is_normal_hermite_multiplier(hermite_multiplier, matrix_):
    '''
    Determine whether hermite_multiplier is the normal Hermite multiplier of matrix.
    This is the COLUMN version.

    Args:
        hermite_multiplier (sympy.Matrix): Candidate column normal Hermite multiplier
        matrix_ (sympy.Matrix): Matrix

    Returns:
        bool: matrix\\_ * hermite_multiplier is in Hermite normal form and hermite_multiplier is in normal form.
    '''
    if matrix_.rank() != matrix_.shape[0]:
        raise ValueError('Matrix must have full row rank')

    r, n = matrix_.shape
    if hermite_multiplier.shape != (n, n):
        raise ValueError("Matrix dimensions don't match")

    # Check matrix . hermite_multiplier = [H 0]
    prod = matrix_ * hermite_multiplier
    hermite_form, residue = prod[:, :r], prod[:, r:]
    # Condition a)
    if not residue.is_zero_matrix:
        return False
    if not is_hnf_col(hermite_form):
        return False
    # Condition b)
    herm_mul_i, herm_mul_n = hermite_multiplier[:, :r], hermite_multiplier[:, r:]
    if not is_hnf_col(herm_mul_n):
        return False
    # Condition c): V_i is reduced with respect to the pivot rows of V_n.  Hubert & Labahn
    # Prop 2.3(ii)(c): if i_1 < ... < i_{n-r} are those pivot rows then
    # 0 <= [V_i]_{i_j,k} < [V_n]_{i_j,j} for every 1 <= k <= r.
    #
    # This used to iterate `range(r)` -- the first r rows of the matrix -- which are not in
    # general the pivot rows.  So it skipped pivot rows the condition does constrain, and
    # tested rows it says nothing about.  It also bounded by `max(V_n[i, :])` and allowed
    # any sign.  The row maximum is in fact the pivot, but only *on a pivot row*, where the
    # entries left of the pivot are reduced below it and those to the right are zero; off a
    # pivot row it means nothing.  Asking for the pivot directly says what is meant.
    #
    # Together with (a) and (b) this makes the normal Hermite multiplier unique: by
    # Prop 2.3(i) the freedom is V_n -> V_n*U and V_i -> V_i + V_n*M, and (b) pins the
    # first while (c) pins the second.
    for j, pivot_row in enumerate(get_pivot_col_indices(herm_mul_n)):
        pivot = herm_mul_n[pivot_row, j]
        for k in range(r):
            if not 0 <= herm_mul_i[pivot_row, k] < pivot:
                return False

    return True

def normal_hnf_col(matrix_):
    '''
    Return the column HNF and the unique normal Hermite multiplier.

    Args:
        matrix_ (sympy.Matrix): Input matrix.

    Returns:
        tuple: Tuple containing:
            hermite_normal_form (sympy.Matrix): The column Hermite normal form of matrix\\_.
            normal_multiplier (sympy.Matrix): The normal Hermite multiplier.

    This is Hubert & Labahn 2013, example 2.4, p. 486, and the multiplier below is the one
    they print.  Proposition 2.3(ii), pp. 484-485, makes it unique: a Hermite multiplier is
    determined up to :math:`V_n \\to V_n U` and :math:`V_i \\to V_i + V_n M`, and conditions
    (b) and (c) pin those two freedoms in turn.

    >>> A = sympy.Matrix([[8, 2, 15, 9, 11],
    ...                   [6, 0, 6, 2, 3]])
    >>> h, v = normal_hnf_col(A)
    >>> h
    Matrix([
    [1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0]])
    >>> v
    Matrix([
    [-1,  -2, -2,  -2, -1],
    [-3, -14, -7, -13, -7],
    [ 1,   1,  2,   1,  0],
    [ 0,   2,  0,   3,  0],
    [ 0,   1,  0,   0,  2]])
    >>> (A * v) == h
    True
    >>> is_normal_hermite_multiplier(v, A)
    True
    '''
    r, n = matrix_.shape

    # The stacked construction below produces the *normal* multiplier rather than merely a
    # valid one, because the Hermite conditions on [A; I] constrain the identity block too.
    # It only works in the Havas-Majewski-Matthews convention: with pivots read from the
    # bottom instead, every pivot of [A; I] falls in the identity block, and a unimodular
    # matrix in Hermite normal form is the identity -- so the construction would return
    # V = I and H = A.
    #
    # Reversing the variable order first and undoing it afterwards gives the multiplier in
    # the convention desr reports.  V_n's rows reverse, which turns a last nonzero entry
    # into a first and entries below a pivot into entries above; its columns reverse too, so
    # the pivot rows still ascend.  A V = [H, 0] survives untouched because the invariant
    # columns are only permuted among themselves.  See :data:`_CONVENTION`.
    reverse_variables = _reversal(n)
    out, _ = _hnf_col_hmm(sympy.Matrix.vstack(matrix_ * reverse_variables, sympy.eye(n)))
    hermite_normal_form, normal_multiplier = out[:r, :], out[r:, :]

    reverse_invariants = (sympy.diag(sympy.eye(r), _reversal(n - r)) if n > r
                          else sympy.eye(n))
    normal_multiplier = reverse_variables * normal_multiplier * reverse_invariants
    hermite_normal_form = hermite_normal_form * reverse_invariants

    assert (matrix_ * normal_multiplier) == hermite_normal_form
    return hermite_normal_form, normal_multiplier

def _hnf_col_hmm(matrix_):
    '''The column Hermite normal form in the Havas-Majewski-Matthews convention.'''
    hnf, unimod = _hnf_row_hmm(matrix_.T)
    return hnf.T, unimod.T

def normal_hnf_row(matrix_):
    '''
    Return the row HNF and the unique normal Hermite multiplier.

    Args:
        matrix_ (sympy.Matrix): Input matrix.

    Returns:
        tuple: Tuple containing:

            hermite_normal_form (sympy.Matrix): The row Hermite normal form of matrix\\_.
            normal_multiplier (sympy.Matrix): The normal Hermite multiplier.
    '''
    hermite_normal_form, normal_multiplier = normal_hnf_col(matrix_=matrix_.T)
    return hermite_normal_form.T, normal_multiplier.T

## Set defaults
hnf_row = hnf_row_lll
"""
Default function for calculating row Hermite normal forms.

Args:
    matrix_ (sympy.Matrix): Input matrix.

Returns:
    tuple:
        hermite_normal_form (sympy.Matrix): The column Hermite normal form of matrix\\_.

        normal_multiplier (sympy.Matrix): The normal Hermite multiplier.

:rtype: (sympy.Matrix, sympy.Matrix)
"""

hnf_col = hnf_col_lll
"""
Default function for calculating column Hermite normal forms.

Args:
    matrix_ (sympy.Matrix): Input matrix.

Returns:
    tuple: Tuple containing:
        hermite_normal_form (sympy.Matrix): The column Hermite normal form of matrix\\_.

        normal_multiplier (sympy.Matrix): The normal Hermite multiplier.
:rtype: (sympy.Matrix, sympy.Matrix)
"""

## Smith normal form

def expand_matrix(matrix_):
    """
    Given a rectangular :math:`n \\times m` integer matrix, return an :math:`(n+1) \\times (m+1)` matrix where the extra row and
    column are 0 except on the first entry which is 1.

    Parameters
    ----------
    matrix_ : sympy.Matrix
        The rectangular matrix to be expanded

    Returns
    -------
    sympy.Matrix
        An :math:`(n+1) \\times (m+1)` matrix


    >>> matrix_ = sympy.diag(*[1, 1, 2])
    >>> expand_matrix(matrix_)
    Matrix([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 2]])

    >>> matrix_ = sympy.Matrix([1])
    >>> expand_matrix(matrix_)
    Matrix([
    [1, 0],
    [0, 1]])

    >>> matrix_ = sympy.Matrix([[0]])
    >>> expand_matrix(matrix_)
    Matrix([
    [1, 0],
    [0, 0]])

    >>> matrix_ = sympy.Matrix([[]])
    >>> expand_matrix(matrix_)
    Matrix([[1]])

    >>> matrix_ = sympy.Matrix([[1, 2, 3]])
    >>> expand_matrix(matrix_)
    Matrix([
    [1, 0, 0, 0],
    [0, 1, 2, 3]])
    >>> expand_matrix(matrix_.T)
    Matrix([
    [1, 0],
    [0, 1],
    [0, 2],
    [0, 3]])
    """
    if len(matrix_.shape) != 2:
        raise ValueError('expand_matrix called for non-square matrix of dimension {}.'.format(len(matrix_.shape)))

    if len(matrix_) == 0:
        return sympy.Matrix([[1]])

    # Add top row
    matrix_ = sympy.Matrix.vstack(sympy.zeros(1, matrix_.shape[1]), matrix_.copy())
    # Add left column
    matrix_ = sympy.Matrix.hstack(sympy.zeros(matrix_.shape[0], 1), matrix_)
    matrix_[0, 0] = 1
    return matrix_

def _swap_ij_rows(matrices, i, j):
    """
        Given an iterable of matrices, swap all of the ith and jth rows.

        Parameters
        ----------
        matrices : iter
            The matrices to be acted on.
        i : int
            Index of first row
        j : int
            Index of second row

        Returns
        -------
        matrices: tuple
            The resulting matrices.

    >>> matrices = [sympy.eye(3) for _ in range(2)]
    >>> matrices = _swap_ij_rows(matrices, 0, 1)
    >>> matrices[0]
    Matrix([
    [0, 1, 0],
    [1, 0, 0],
    [0, 0, 1]])
    >>> matrices[1]
    Matrix([
    [0, 1, 0],
    [1, 0, 0],
    [0, 0, 1]])
    """
    matrices = [matrix_.copy() for matrix_ in matrices]
    for matrix_ in matrices:
        matrix_.row_swap(i, j)
    return tuple(matrices)

def _swap_ij_cols(matrices, i, j):
    """
        Given an iterable of matrices, swap all of the ith and jth columns. Acts in place.

        Parameters
        ----------
        matrices : iter
            The matrices to be acted on.
        i : int
            Index of first column
        j : int
            Index of second column

        Returns
        -------
        matrices : tuple
            The resulting matrices.

    >>> matrices = [sympy.eye(3) for _ in range(2)]
    >>> matrices = _swap_ij_cols(matrices, 0, 1)
    >>> matrices[0]
    Matrix([
    [0, 1, 0],
    [1, 0, 0],
    [0, 0, 1]])
    >>> matrices[1]
    Matrix([
    [0, 1, 0],
    [1, 0, 0],
    [0, 0, 1]])
    """
    return tuple(matrix_.T for matrix_ in _swap_ij_rows([_matrix.T for _matrix in matrices], i, j))

 # this one is less performant, because it has to convert the entire matrix to bools before checking.  so i have the other one.
def is_non_negative(matrix_):
    """
    Check that all entries in a sympy matrix are at least 0.

    >>> x = sympy.eye(2) * 3
    >>> is_non_negative(x)
    True

    >>> y = sympy.Matrix([[0,1], [12.2, -1]])
    >>> is_non_negative(y)
    False
    """
    for m in matrix_:
        if m<0:
            return False
    return True

def element_wise_lt(matrix_, other):
    """
    Given a rectangular :math:`n \\times m` integer matrix, return a matrix of bools with (i, j)th entry equal to
    matrix_[i,j] < other or other[i,j], depending on the type of other.

    Parameters
    ----------
    matrix_ : sympy.Matrix
        The input matrix
    other : sympy.Matrix, int
        Value(s) to compare it to.

    Returns
    -------
    sympy.Matrix
        :math:`n \\times m` Boolean matrix


    >>> x = sympy.eye(2) * 3
    >>> element_wise_lt(x, 2)
    Matrix([
    [0, 1],
    [1, 0]])
    >>> element_wise_lt(x, sympy.Matrix([[4, -1], [1, 1]]))
    Matrix([
    [1, 0],
    [1, 0]])
    """
    if isinstance(other, sympy.Matrix):
        if other.shape != matrix_.shape:
            raise ValueError('Incompatible shapes: {} != {}'.format(matrix_.shape, other.shape))
        return sympy.Matrix(matrix_.shape[0], matrix_.shape[1], lambda i, j: 1 if matrix_[i, j] < other[i, j] else 0)
    else:
        return sympy.Matrix(matrix_.shape[0], matrix_.shape[1], lambda i, j: 1 if matrix_[i, j] < other else 0)

def is_smf(matrix_):
    """
    Given a rectangular :math:`n \\times m` integer matrix, determine whether it is in Smith normal form or not.

    Parameters
    ----------
    matrix_ : sympy.Matrix
        The rectangular matrix to be decomposed

    Returns
    -------
    bool
        True if in Smith normal form, False otherwise.


    >>> matrix_ = sympy.diag(1, 1, 2)
    >>> is_smf(matrix_)
    True
    >>> matrix_ = sympy.diag(-1, 1, 2)
    >>> is_smf(matrix_)
    False
    >>> matrix_ = sympy.diag(2, 1, 1)
    >>> is_smf(matrix_)
    False
    >>> matrix_ = sympy.diag(1, 2, 0)
    >>> is_smf(matrix_)
    True
    >>> matrix_ = sympy.diag(2, 6, 0)
    >>> is_smf(matrix_)
    True
    >>> matrix_ = sympy.diag(2, 5, 0)
    >>> is_smf(matrix_)
    False
    >>> matrix_ = sympy.diag(0, 1, 1)
    >>> is_smf(matrix_)
    False
    >>> matrix_ = sympy.diag(0)
    >>> is_smf(sympy.diag(0)), is_smf(sympy.diag(1)), is_smf(sympy.Matrix()),
    (True, True, True)

    Check a real example:

    >>> matrix_ = sympy.Matrix([[2, 4, 4],
    ...                         [-6, 6, 12],
    ...                         [10, -4, -16]])
    >>> is_smf(matrix_)
    False

    >>> matrix_ = sympy.diag(2, 6, 12)
    >>> is_smf(matrix_)
    True

    Check it works for non-square matrices:

    >>> matrix_ = sympy.Matrix(4, 5, range(20))
    >>> is_smf(matrix_)
    False

    >>> matrix_ = sympy.Matrix([[1, 0], [1, 2]])
    >>> is_smf(matrix_)
    False
    """
    if len(matrix_) == 0:
        return True

    # Check its a diagonal matrix
    if not matrix_.is_diagonal():
        return False

    # Check all entries are non-negative
    # if any(element_wise_lt(matrix_, 0)): # this one is less performant, because it has to convert the entire matrix to bools before checking.  so i have the other one.
    if not is_non_negative(matrix_):
        return False

    current_int = matrix_[0, 0]
    for i in range(1, min(matrix_.shape)):
        next_int = matrix_[i, i]
        # All zeros go at the end
        if (current_int == 0):
            if (next_int != 0):
                return False

        if (current_int > 0):
            if next_int % current_int:
                return False
        current_int = next_int
    return True

def smf(matrix_):
    """
    Given a rectangular :math:`n \\times m` integer matrix, calculate the Smith Normal Form :math:`S` and multipliers
    :math:`U`, :math:`V` such that :math:`U \\textrm{matrix_} V = S`.

    Parameters
    ----------
    matrix_ : sympy.Matrix
        The rectangular matrix to be decomposed.

    Returns
    -------
    S : sympy.Matrix
        The Smith normal form of matrix\\_.
    U : sympy.Matrix
        :math:`U` (the matrix representing the row operations of the decomposition).
    V : sympy.Matrix
        :math:`V` (the matrix representing the column operations of the decomposition).


    :rtype: (sympy.Matrix, sympy.Matrix, sympy.Matrix)


    >>> matrix_ = sympy.Matrix([[2, 4, 4],
    ...                         [-6, 6, 12],
    ...                         [10, -4, -16]])
    >>> smf(matrix_)[0]
    Matrix([
    [2, 0,  0],
    [0, 6,  0],
    [0, 0, 12]])

    >>> matrix_ = sympy.diag(2, 1, 0)
    >>> smf(matrix_)[0]
    Matrix([
    [1, 0, 0],
    [0, 2, 0],
    [0, 0, 0]])

    >>> matrix_ = sympy.diag(5, 2, 0)
    >>> smf(matrix_)[0]
    Matrix([
    [1,  0, 0],
    [0, 10, 0],
    [0,  0, 0]])
    """
    matrix_ = matrix_.copy()
    transformed = matrix_.copy()

    if is_smf(transformed):
        row_actions = sympy.eye(matrix_.shape[0])
        col_actions = sympy.eye(matrix_.shape[1])
        return transformed, row_actions, col_actions

    # First make a pivot in the top left corner
    transformed, row_actions = hnf_row(matrix_=transformed)
    transformed, col_actions = hnf_col(matrix_=transformed)

    # Now put the lower block into Smith normal form
    block, _row_actions, _col_actions = smf(transformed[1:, 1:])

    # Compose the row and column actions found
    row_actions = row_actions * expand_matrix(_row_actions)
    col_actions = expand_matrix(_col_actions) * col_actions
    transformed[1:, 1:] = block
    assert row_actions * matrix_ * col_actions == transformed

    # Now make sure the pivot is in the right place
    for i in range(min(*transformed.shape) - 1):
        if transformed[i+1,i+1] == 0:
            break
        # Move the pivot entry down if it's 0
        if transformed[i,i] == 0:
            transformed, row_actions = _swap_ij_rows([transformed, row_actions], i, i+1)
            transformed, col_actions = _swap_ij_cols([transformed, col_actions], i, i+1)

    assert row_actions * matrix_ *  col_actions == transformed
    assert transformed[1:, 0].is_zero_matrix
    assert transformed[0, 1:].is_zero_matrix

    # Enforce diagonal entries to divide the next diagonal entry
    for i in range(min(*transformed.shape) - 1):
        assert row_actions * matrix_ * col_actions == transformed
        if transformed[i+1,i+1] == 0:
            break

        # If we don't have division, ensure we do
        if transformed[i+1,i+1] % transformed[i,i]:
            # If we have division in the wrong order, just swap them around
            if not transformed[i,i] % transformed[i+1,i+1]:
                transformed, row_actions = _swap_ij_rows([transformed, row_actions], i, i + 1)
                transformed, col_actions = _swap_ij_cols([transformed, col_actions], i, i + 1)
                # Reset the counter
                i = 0
                continue

            # Add one of column i+1 to column i, and perform a row reduction to get the hcf
            transformed[:, i] += transformed[:, i+1]
            col_actions[:, i] += col_actions[:, i+1]
            assert row_actions * matrix_ * col_actions == transformed

            transformed, _row_actions = hnf_row(transformed)
            row_actions = _row_actions * row_actions
            assert row_actions * matrix_ * col_actions == transformed

            # Now subtract the right multiple of column i from column i+1
            assert (transformed[i, i+1] % transformed[i, i]) == 0
            multiple = int(transformed[i, i+1] / transformed[i, i])
            transformed[:, i+1] -= multiple * transformed[:, i]
            col_actions[:, i+1] -= multiple * col_actions[:, i]
            assert row_actions * matrix_ * col_actions == transformed

            # Reset the loop, so we check from the beginning of the diagonal.
            i = 0

    # _prod = row_actions * matrix_ *  col_actions
    assert row_actions * matrix_ *  col_actions == transformed

    assert is_smf(transformed)
    return transformed, row_actions, col_actions



if __name__ == '__main__':
    import doctest
    doctest.testmod()