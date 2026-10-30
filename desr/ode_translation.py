
import enum

import sympy

from desr.matrix_normal_forms import normal_hnf_col, normal_hnf_row, is_hnf_col
from desr.ode_system import ODESystem
from desr.sympy_helper import sympify_expression
from desr.tex_tools import matrix_to_tex


class ReductionScheme(str, enum.Enum):
    '''
    Which of the three reductions of Hubert & Labahn a system was put through.

    The scheme decides how many variables the reduced system has, whether it carries
    auxiliary variables, and therefore what has to be known to translate values back.  It
    is recorded rather than inferred: inferring it is only correct when
    :meth:`~desr.ode_translation.ODETranslator.translate` did the choosing, and calling
    :meth:`~desr.ode_translation.ODETranslator.translate_general` on a system that
    `translate` would have handed to another scheme is perfectly legal.

    A member is a :class:`str`, so ``scheme == 'general'`` holds and a plain string may be
    passed wherever one of these is expected.

    >>> ReductionScheme.GENERAL
    ReductionScheme.GENERAL
    >>> ReductionScheme.GENERAL == 'general'
    True
    >>> print(ReductionScheme.GENERAL)
    general
    '''

    #: The parameter scheme: the reduced system has :math:`n - r` variables and no
    #: auxiliaries, and it is the only scheme that handles initial conditions.
    PARAMETER = 'parameter'

    #: The dependent-variable scheme, where the scaling action ignores the independent
    #: variable.  The reduced system has :math:`n` variables.
    DEP_VAR = 'dep_var'

    #: The general scheme, which always applies.  :math:`n + 1` variables with the
    #: auxiliaries, :math:`n - r + 1` without.
    GENERAL = 'general'

    def __str__(self):
        # Spelled out so that the printed form is the same on every Python version: the
        # `str` mixin's __str__ and __format__ have changed between releases.
        return self.value

    def __repr__(self):
        # `<ReductionScheme.GENERAL: 'general'>` is the default and reads badly in docs and
        # in error messages, which is where this mostly appears.
        return '{}.{}'.format(type(self).__name__, self.name)

#: The names a reduction may be asked to give the variables it creates.  A naming scheme is a
#: dict over these keys; every key is optional, and an absent one takes the default for
#: whichever reduction scheme runs.  The categories are the ones the invariants fall into:
#: an invariant that is a monomial in the parameters alone, the one that carries the clock,
#: and the ones that involve a dependent variable.  Under the general reduction scheme some
#: of those categories are empty -- there may be no parameter invariants, and the clock may
#: be entangled with a dependent variable, as in Hubert & Labahn 6.6 where `t z1^3` is an
#: invariant -- and the names for them simply go unused.
NAMING_SCHEME_KEYS = ('parameters', 'independent', 'dependents', 'auxiliaries')

#: Defaults for the parameter reduction scheme: desr's own nondimensionalisation names.
PARAMETER_SCHEME_NAMES = {'parameters': 'kappa', 'independent': 'tau', 'dependents': 'nu',
                          'auxiliaries': 'x'}

#: Defaults for the general and dependent-variable schemes: Hubert & Labahn's notation, where
#: the auxiliaries are :math:`x` and the invariants :math:`y`, so that output can be compared
#: with the paper directly.  `independent` is :const:`None`, meaning the independent variable
#: keeps the name it already has.
GENERAL_SCHEME_NAMES = {'parameters': 'kappa', 'independent': None, 'dependents': 'y',
                        'auxiliaries': 'x'}


def _normalise_naming_scheme(naming_scheme):
    '''
    Validate a naming scheme and return it as a dict, which may be partial.

    Args:
        naming_scheme (dict): A dict over some of :data:`NAMING_SCHEME_KEYS`.  :const:`None`
            means "no preferences", and gives an empty dict.

    Returns:
        dict: The scheme, copied, with only the keys the caller actually set.

    A tuple is refused.  Until 2.0 the scheme was a positional triple whose order followed
    the variable order; the variable order changed, which would have silently re-pointed
    every entry, so the positional form is gone rather than quietly re-read:

    >>> _normalise_naming_scheme(('tau', 'nu', 'kappa'))
    Traceback (most recent call last):
        ...
    ValueError: A naming scheme is a dict, not a tuple.  Write {'independent': 'tau', 'dependents': 'nu', 'parameters': 'kappa'} -- the keys are parameters, independent, dependents, auxiliaries, and any you leave out take their default.

    >>> _normalise_naming_scheme({'independent': 'tau', 'dependents': ['u', 'v']}) == {
    ...     'independent': 'tau', 'dependents': ['u', 'v']}
    True

    >>> _normalise_naming_scheme(None)
    {}

    >>> _normalise_naming_scheme({'indep': 'tau'})
    Traceback (most recent call last):
        ...
    ValueError: Unknown naming scheme key 'indep'.  The keys are parameters, independent, dependents, auxiliaries.
    '''
    if naming_scheme is None:
        return {}

    if isinstance(naming_scheme, (tuple, list)):
        raise ValueError(
            "A naming scheme is a dict, not a {}.  Write {{'independent': 'tau', "
            "'dependents': 'nu', 'parameters': 'kappa'}} -- the keys are {}, and any you "
            "leave out take their default.".format(type(naming_scheme).__name__,
                                                   ', '.join(NAMING_SCHEME_KEYS)))

    if not isinstance(naming_scheme, dict):
        raise ValueError('A naming scheme is a dict over the keys {}.  Got a {}.'.format(
            ', '.join(NAMING_SCHEME_KEYS), type(naming_scheme).__name__))

    for key in naming_scheme:
        if key not in NAMING_SCHEME_KEYS:
            raise ValueError("Unknown naming scheme key {!r}.  The keys are {}.".format(
                key, ', '.join(NAMING_SCHEME_KEYS)))

    for key in ('parameters', 'auxiliaries'):
        value = naming_scheme.get(key)
        if value is not None and not isinstance(value, str):
            raise ValueError(
                "The naming scheme's {!r} must be a single string, used as a pattern: the "
                "number of them is not known until the reduction runs.  Got a {}.".format(
                    key, type(value).__name__))

    return dict(naming_scheme)


def _new_symbols(pattern, count, start_at):
    '''
    `count` symbols named `pattern` with consecutive indices appended, counting from `start_at`.
    '''
    return [sympy.Symbol('{}{}'.format(pattern, i + start_at)) for i in range(count)]


def _names_for(pattern, count, start_at, what):
    '''
    Turn one entry of a naming scheme into exactly `count` symbols.

    A string is a pattern and gets indices appended.  A list or tuple is taken literally, and
    must have exactly the right length -- that is the one thing a caller can get wrong here
    without the mistake showing up until much later.

    Args:
        pattern: A string pattern, or an explicit sequence of names.
        count (int): How many symbols are needed.
        start_at (int): Where appended indices begin.
        what (str): What is being named, for the error message.

    Returns:
        list of sympy.Symbol
    '''
    if isinstance(pattern, str):
        return _new_symbols(pattern, count, start_at)

    names = list(pattern)
    if len(names) != count:
        raise ValueError(
            'The naming scheme gives {} names for the {}, but this reduction has {} of '
            'them: {}.'.format(len(names), what, count, names))
    # sympy.Symbol rather than sympify: sympifying a name like 'S' or 'E' would resolve to a
    # sympy singleton instead of a variable of that name.
    return [name if isinstance(name, sympy.Symbol) else sympy.Symbol(str(name))
            for name in names]


def _classify_reduced_derivatives(names, monomials, rates, dynamic):
    '''
    Replace the rate of every parameter-monomial with :const:`None`, so the reduced system
    records it as a parameter rather than as a dependent whose equation says zero.

    A column whose monomial involves none of the original system's *variables* -- its
    independent variable and its dependents -- is a monomial in the parameters alone.  It
    cannot vary under any trajectory, so it is a parameter of the reduced system, and
    storing it as `dy/dt = 0` would assert the different and false thing: that it happens
    to stay constant along the solution.

    This applies to both blocks of the Hermite multiplier.  An auxiliary that is a monomial
    in the parameters is a parameter too -- auxiliary-ness is about being free to choose a
    value, not about whether the value varies.

    Args:
        names (iter of sympy.Symbol): The reduced system's symbols, in order.
        monomials (iter of sympy.Expr): What each one denotes, in the same order.
        rates (iter of sympy.Expr): The computed rate of each, in the same order.
        dynamic (set of sympy.Symbol): The original system's variables.

    Returns:
        list: The rates, with :const:`None` wherever the symbol is a parameter.
    '''
    classified = []
    for name, monomial, rate in zip(names, monomials, rates):
        if monomial.free_symbols & dynamic:
            classified.append(rate)
            continue
        # Proved, and worth failing loudly on: a monomial in the parameters alone has no
        # time dependence, so its computed rate must already be zero.  A nonzero one here
        # would mean the monomial and the rate disagree about what this symbol is.
        if rate != 0:
            raise AssertionError(
                '{} = {} involves no variable of the original system, so it cannot change, '
                'but its computed rate is {}'.format(name, monomial, rate))
        classified.append(None)
    return classified


def _rename_indep_var(system, new_name):
    '''
    Rename a reduced system's independent variable, or leave it alone if `new_name` is None.

    The general and dependent-variable schemes keep the original clock, so their default is
    to leave the name as it is.  A caller who sets `independent` is asking for it to be
    renamed, and gets that -- it is only a relabelling, since under those schemes the clock
    is the system's own independent variable either way.

    Args:
        system (ODESystem): The reduced system.
        new_name: The name to give the independent variable, or :const:`None` to keep it.

    :rtype: ODESystem
    '''
    if new_name is None or sympy.Symbol(str(new_name)) == system.independent_variable:
        return system
    system.rename_indep_var(sympy.Symbol(str(new_name)))
    return system


def _rates_of_change(system):
    '''
    The rate of change of every variable of a system, in variable order, with a parameter's
    absent equation read as zero.

    The scaling formulas want one entry per column of the scaling matrix, and a parameter
    does not change.  :attr:`~desr.ode_system.ODESystem.derivatives` reports :const:`None`
    for a parameter, because a parameter has no equation -- a distinction that matters to
    anyone reading the system and not at all to the quadrature.  It is flattened here, and
    nowhere else, so that "has no equation" and "has the equation zero" stay distinct
    everywhere a human can see them.

    Args:
        system (ODESystem): The system being reduced.

    Returns:
        list: One expression per variable of ``system``.
    '''
    return [sympy.sympify(0) if deriv is None else deriv for deriv in system.derivatives]


def _scaled_rates(variables, factors, divide_by=None):
    '''
    Multiply each new variable by its rate factor, elementwise.

    Args:
        variables (list of sympy.Symbol): The new variables, one per column of ``factors``.
        factors (sympy.Matrix): A row of rate factors.
        divide_by: Optional expression to divide the result through by.

    Returns:
        list: One derivative per variable, empty when there are no variables.

    A reduction can legitimately produce no invariants, or no auxiliaries.  ``dz/dt = z``
    is the smallest example: it reduces completely, leaving only the independent variable
    and one auxiliary.  sympy builds a 0x0 matrix from the empty list and then refuses to
    multiply it elementwise against the 1x0 row of factors, so the empty case is taken
    here rather than left to raise.
    '''
    if not variables:
        return []
    scaled = sympy.Matrix(variables).T.multiply_elementwise(factors)
    if divide_by is not None:
        scaled = scaled / divide_by
    return list(scaled)


def _int_inv(matrix_):
    ''' Given an integer matrix, return the inverse, ensuring we do it right in the integers

    >>> matrix_ = sympy.Matrix([[5, 2, -2],
    ...                         [-7, -3, 3],
    ...                         [17, 8, -7]])
    >>> _int_inv(matrix_)
    Matrix([
    [ 3, 2, 0],
    [-2, 1, 1],
    [ 5, 6, 1]])
    '''
    if abs(matrix_.det()) != 1:
        raise ValueError('Nonunimodular matrix with det {} cannot be inverted over the integers.'.format(matrix_.det()))
    sympy_inverse = sympy.Matrix(matrix_).inv()
    return sympy_inverse

class ODETranslator(object):
    '''
    An object used for translating between systems of ODEs according to a scaling matrix.
    The key data are :attr:`~desr.ode_translation.ODETranslator.scaling_matrix` and
    :attr:`~desr.ode_translation.ODETranslator.herm_mult`, which together contain all the information needed
    (along with a variable order) to reduce an :class:`~desr.ode_system.ODESystem`.

    Args:
        scaling_matrix (sympy.Matrix): Matrix that defines the torus action on the system.
        symbols_domain (iter of sympy.Symbol, optional): An ordering of the variables we expect to act upon.
            If this is not given, we will act on a system according to the position of the variables in
            :attr:`~desr.ode_system.ODESystem.symbols`, as long as there is the correct number of variables.
        hermite_multiplier (sympy.Matrix, optional): User-defined Hermite multiplier, that puts
            :attr:`~desr.ode_translation.ODETranslator.scaling_matrix` into
            column Hermite normal form.
            If not given, the normal Hermite multiplier will be calculated.
        naming_scheme (dict): What to call the variables the reduction creates.  Every key is
            optional; one left out takes the default for whichever reduction scheme runs.

            - `'independent'`: the new name of the independent (time) variable.  The general
              and dependent-variable schemes keep the system's own clock, so their default is
              to leave the name alone; setting this renames it.
            - `'dependents'`: either a pattern, to which indices are appended, or an explicit
              list of names.  A list must have exactly as many entries as the reduction
              produces, which is the number of dependent variables under the parameter
              scheme and the number of invariants under the general one.
            - `'parameters'`: the pattern for the scaling-invariant constants the algorithm
              computes.  How many there will be is not known until the reduction runs, so
              this must be a single string rather than a list.
            - `'auxiliaries'`: the pattern for the auxiliary variables, likewise a single
              string.  Only the general and dependent-variable schemes make any.

            The defaults are :data:`PARAMETER_SCHEME_NAMES` for the parameter scheme and
            :data:`GENERAL_SCHEME_NAMES` for the other two, so that output under the general
            scheme can be compared with Hubert & Labahn directly.  Anything set here is
            honoured by all of them.  A positional tuple is refused; see
            :func:`_normalise_naming_scheme`.
        new_indices_start_at (int): Where the indices for the new variables should start at.  Default is 0.
            Applies to independent variables if not using a user-supplied list, and always to the reduced constants.

    A user-defined Hermite multiplier must be unimodular.  Then every invariant is a product
    of integer powers of its invariant columns, and the reduction and its reverse are
    rational.  A multiplier that is not -- here with :math:`y_0^2` in place of :math:`y_0`,
    so the determinant is :math:`-2` -- is refused at once, rather than failing later when
    its inverse is first needed:

    >>> ODETranslator(sympy.Matrix([[1, 0, 1, 1, -1], [0, 1, 0, -1, 1]]),
    ...                hermite_multiplier=sympy.Matrix([[0, 0, 2,  0, 0],
    ...                                                 [0, 0, 0,  1, 0],
    ...                                                 [1, 1, -2, -1, 0],
    ...                                                 [0, 0, 0,  0, 1],
    ...                                                 [0, 1, 0, -1, 1]]))
    Traceback (most recent call last):
        ...
    ValueError: The Hermite multiplier has determinant -2, so it is not unimodular.  Some invariants are then not products of integer powers of its invariant columns, and reversing the reduction would need roots of degree 2.  Choose invariants from which every invariant can be built with integer powers; see extend_from_invariants.
    '''
    def __init__(self, scaling_matrix, symbols_domain=None, hermite_multiplier=None, naming_scheme=None, new_indices_start_at=0):
        scaling_matrix = scaling_matrix.copy()

        self._scaling_matrix = scaling_matrix
        self._symbols_domain = symbols_domain
        self._new_indices_start_at = new_indices_start_at
        self.set_naming_scheme(naming_scheme)

        if (symbols_domain is not None) and (self.n != len(self.symbols_domain)):
            raise ValueError('{} variables given but we have {} variables (columns) in the scaling matrix'.format(len(self.symbols_domain), self.n))

        if hermite_multiplier is not None:
            determinant = hermite_multiplier.det()
            if abs(determinant) != 1:
                raise ValueError(
                    'The Hermite multiplier has determinant {d}, so it is not unimodular.  '
                    'Some invariants are then not products of integer powers of its invariant '
                    'columns, and reversing the reduction would need roots of degree {index}.  '
                    'Choose invariants from which every invariant can be built with integer '
                    'powers; see extend_from_invariants.'.format(
                        d=determinant, index=abs(determinant)))
            self._scaling_matrix_hnf = scaling_matrix * hermite_multiplier
            if not is_hnf_col(self._scaling_matrix_hnf):
                raise ValueError('{}.{}={} is not in HNF'.format(scaling_matrix,
                                 hermite_multiplier, self._scaling_matrix_hnf))
            self._herm_mult = hermite_multiplier
        else:
            self._scaling_matrix_hnf, self._herm_mult = normal_hnf_col(scaling_matrix)

        self._inv_herm_mult = None
        self._dep_var_inv_herm_mult_cache = {}

    def set_naming_scheme(self, naming_scheme):
        '''
        Set the renaming scheme for the reduced system after translation.  See documentation for the constructor for `ODETranslator`.

        Only the keys the caller sets are stored, so that the remainder can take whichever
        reduction scheme's defaults apply when the reduction actually runs.
        '''
        self._naming_scheme = _normalise_naming_scheme(naming_scheme)

    @property
    def naming_scheme(self):
        '''
        The naming scheme as set, holding only the keys that were given.

        :rtype: dict
        '''
        return dict(self._naming_scheme)

    def _resolved_naming_scheme(self, defaults):
        '''
        The naming scheme with unset keys filled in from `defaults`.

        What the caller asked for always wins; `defaults` only supplies the rest.  This is
        the one place the two sets of defaults meet, so that a scheme passed in is honoured
        by whichever reduction runs, while the names a reduction falls back on stay its own
        -- Hubert & Labahn's `x` and `y` for the general scheme, desr's `tau`, `nu` and
        `kappa` for the parameter scheme.

        Args:
            defaults (dict): :data:`PARAMETER_SCHEME_NAMES` or :data:`GENERAL_SCHEME_NAMES`.

        :rtype: dict
        '''
        resolved = dict(defaults)
        resolved.update(self._naming_scheme)
        return resolved

    def __repr__(self):
        return 'A=\n{}\nV=\n{}\nW=\n{}'.format(self.scaling_matrix.__repr__(),
                                               self.herm_mult.__repr__(),
                                               self.inv_herm_mult.__repr__())

    def to_tex(self):
        '''
        Returns:
            str: The scaling matrix :math:`A`, the Hermite multiplier :math:`V` and :math:`W = V^{-1}`, in beautiful LaTeX.

        >>> print(ODETranslator(sympy.Matrix(range(12)).reshape(3, 4)).to_tex())
        A=
        0 & 1 & 2 & 3 \\\\
        4 & 5 & 6 & 7 \\\\
        8 & 9 & 10 & 11 \\\\
        V=
        -1 & 1 & 2 & 1 \\\\
        1 & 0 & -3 & -2 \\\\
        0 & 0 & 0 & 1 \\\\
        0 & 0 & 1 & 0 \\\\
        W=
        0 & 1 & 2 & 3 \\\\
        1 & 1 & 1 & 1 \\\\
        0 & 0 & 0 & 1 \\\\
        0 & 0 & 1 & 0 \\\\

        '''
        to_print = (self.scaling_matrix, self.herm_mult, self.inv_herm_mult)
        to_print = map(matrix_to_tex, to_print)
        return 'A=\n{}\nV=\n{}\nW=\n{}'.format(*to_print)

    @property
    def scaling_matrix(self):
        """
        Returns:
            sympy.Matrix: The scaling matrix that this translator corresponds to, often denoted :math:`A`.
        """
        return self._scaling_matrix

    @property
    def r(self):
        '''
        Returns: int: The dimension of the scaling action: :math:`r`.  In particular it is the number of rows of :attr:`~scaling_matrix`.
        '''
        return self._scaling_matrix.shape[0]

    @property
    def n(self):
        '''
        Returns:
            int: The number of original variables that the scaling action is acting on: :math:`n`.  In particular, it is the number of columns of :attr:`~scaling_matrix`.
        '''
        return self._scaling_matrix.shape[1]

    @property
    def herm_mult(self):
        '''
        Returns:
            sympy.Matrix: A column Hermite multiplier :math:`V` that puts the :attr:`~scaling_matrix` in column Hermite normal form.  That is: :math:`AV = H` is in column Hermite normal form.
        '''
        return self._herm_mult.copy()

    @property
    def herm_form(self):
        '''
        Returns:
            sympy.Matrix: The Hermite normal form of the scaling matrix: :math:`H = AV`.
        '''
        return self._scaling_matrix_hnf.copy()

    @property
    def herm_mult_i(self):
        '''
        Returns:
            sympy.Matrix: :math:`V_{\\mathfrak{i}}`: the first :math:`r` columns of :math:`V`. The columns represent the auxiliary variables of the reduction.
        '''
        return self.herm_mult[:, :self.r]

    @property
    def herm_mult_n(self):
        '''
        Returns:
            sympy.Matrix: :math:`V_{\\mathfrak{n}}`: the last :math:`n-r` columns of the Hermite multiplier :math:`V``.  The columns represent the invariants of the scaling action.
        '''
        return self.herm_mult[:, self.r:]


    @property
    def inv_herm_mult(self):
        '''
        Returns:
            sympy.Matrix: The inverse of the Hermite multiplier :math:`W=V^{-1}`.
        '''
        if self._inv_herm_mult is None:
            self._inv_herm_mult = _int_inv(self._herm_mult)
        return self._inv_herm_mult.copy()

    @property
    def inv_herm_mult_u(self):
        '''
        Returns:
            sympy.Matrix: :math:`W_{\\mathfrak{u}}`: the first :math:`r` rows of :math:`W`.
        '''
        return self.inv_herm_mult[:self.r, :]

    @property
    def inv_herm_mult_d(self):
        """
        Returns:
            sympy.Matrix: :math:`W_{\\mathfrak{d}}`: the last :math:`n-r` rows of :math:`W`.
        """
        return self.inv_herm_mult[self.r:, :]


    def dep_var_herm_mult(self, indep_var_index=0):
        '''
        Args:
            indep_var_index (int): The index of the independent variable.

        Returns:
            sympy.Matrix: The Hermite multiplier :math:`V`, ignoring the independent variable.

        The independent variable's row of :math:`V` must have exactly one nonzero entry, as
        it does for the normal Hermite multiplier of an actual system -- that is what leaves
        :math:`n - 1` columns and so a square result.  Given a row with more, the result is
        rectangular and inverting it fails with a `NonSquareMatrixError` from somewhere
        further on, which says nothing useful.

        >>> translator = ODETranslator(sympy.Matrix(range(12)).reshape(3, 4))
        >>> translator.herm_mult
        Matrix([
        [-1, 1,  2,  1],
        [ 1, 0, -3, -2],
        [ 0, 0,  0,  1],
        [ 0, 0,  1,  0]])
        >>> translator.dep_var_herm_mult(2)
        Matrix([
        [-1, 1,  2],
        [ 1, 0, -3],
        [ 0, 0,  1]])
        '''
        new_herm_mult = self.herm_mult.copy()
        col_to_delete = new_herm_mult[indep_var_index, :]
        new_herm_mult.row_del(indep_var_index)
        cols_to_keep = [ind for ind, val in enumerate(col_to_delete) if val == 0]
        new_herm_mult = new_herm_mult.extract(range(new_herm_mult.rows), cols_to_keep)
        # new_herm_mult = new_herm_mult[:, ~col_to_delete.astype(bool)]
        return new_herm_mult

    def dep_var_inv_herm_mult(self, indep_var_index=0):
        '''
        Args:
            indep_var_index (int): The index of the independent variable.

        Returns:
            sympy.Matrix:
                The inverse Hermite multiplier :math:`W`, ignoring the independent variable.

        >>> translator = ODETranslator(sympy.Matrix(range(12)).reshape(3, 4))
        >>> translator.inv_herm_mult
        Matrix([
        [0, 1, 2, 3],
        [1, 1, 1, 1],
        [0, 0, 0, 1],
        [0, 0, 1, 0]])
        >>> translator.dep_var_inv_herm_mult(2)
        Matrix([
        [0, 1, 3],
        [1, 1, 1],
        [0, 0, 1]])
        '''
        if indep_var_index not in self._dep_var_inv_herm_mult_cache:
            self._dep_var_inv_herm_mult_cache[indep_var_index] = \
                _int_inv(self.dep_var_herm_mult(indep_var_index=indep_var_index))
        return self._dep_var_inv_herm_mult_cache[indep_var_index].copy()

    @property
    def symbols_domain(self):
        '''
        Returns:
            tuple, None: The variables that the scaling action acts on.
        '''
        return self._symbols_domain

    @classmethod
    def from_ode_system(cls, ode_system, naming_scheme=None, new_indices_start_at = 0):
        '''
        Create a :class:`ODETranslator` given an :class:`ode_system.ODESystem` instance, by taking the maximal scaling matrix.

        Args:
            ode_system (ODESystem)

        :rtype: ODETranslator
        '''
        return cls(scaling_matrix=ode_system.maximal_scaling_matrix(), symbols_domain=ode_system.symbols, naming_scheme=naming_scheme, new_indices_start_at = new_indices_start_at)

    def column_monomials(self):
        '''
        What every column of :attr:`~desr.ode_translation.ODETranslator.herm_mult` stands
        for, as a monomial in the original symbols.

        The auxiliaries come first -- the :math:`r` columns of :math:`V_i` -- then the
        :math:`n - r` invariants of :math:`V_n`, which is the column order of the multiplier
        itself.  So the index into this tuple *is* the column index the
        ``multiplier_*_column`` methods take.

        Returns:
            tuple of sympy.Expr: One monomial per column.

        >>> equations = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
        ...              'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.column_monomials()
        (1/(e_0*k_1), e_0, k_2/(e_0*k_1), k_m1/(e_0*k_1), e_0*k_1*t, c/e_0, s/e_0)
        '''
        variables = self._validate_variables(None, self.n, 'y', True)
        return tuple(scale_action(variables, self.herm_mult))

    def column_names(self, system=None, scheme=None):
        '''
        A name for every column of :attr:`~desr.ode_translation.ODETranslator.herm_mult`,
        in the same order as :meth:`column_monomials`.

        The auxiliary columns are named from the naming scheme's ``auxiliaries`` pattern
        even under the parameter scheme, which keeps none of them: the *columns* exist in
        the multiplier whether or not the reduction retains them, and naming them is what
        lets a column operation refer to one.  :meth:`auxiliary_names` is the other
        question -- what the *reduced system* calls them -- and is empty there.

        Args:
            system (ODESystem): The system being reduced.
            scheme (desr.ode_translation.ReductionScheme, optional): As for
                :meth:`named_symbols`.

        Returns:
            tuple of sympy.Symbol

        >>> equations = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
        ...              'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.column_names(system)
        (x0, x1, kappa0, kappa1, tau, nu0, nu1)
        '''
        naming = self._resolved_naming_scheme(
            PARAMETER_SCHEME_NAMES if scheme == ReductionScheme.PARAMETER
            else GENERAL_SCHEME_NAMES)
        auxiliaries = _new_symbols(naming['auxiliaries'], self.herm_mult_i.shape[1],
                                   self._new_indices_start_at)
        return tuple(auxiliaries) + self.invariant_names(system, scheme)

    def column_of(self, monomial_or_name, system=None, scheme=None):
        '''
        The column index of an invariant or auxiliary, by its monomial or by its name.

        The point is to **demystify** the number, not to hide it.  The index is still an
        :class:`int`, still what the ``multiplier_*_column`` methods take, and still there
        to be quoted in a write-up -- but where it came from is now on the page::

            # before -- why 6 and 2?
            translator.multiplier_add_columns(6, 2)

            # after -- and `clock` and `epsilon` can be printed
            clock = translator.column_of('k_1*s_0*t')
            epsilon = translator.column_of('e_0/s_0')
            translator.multiplier_add_columns(clock, epsilon)

        The ``multiplier_*`` methods deliberately keep taking integers, so that the number a
        reduction depends on stays visible in the code that performs it.

        Args:
            monomial_or_name: A monomial in the original symbols -- `'e_0/s_0'`, or the
                expression itself -- or the name a column's invariant gets, such as `tau`.
                A string is parsed with the library's own convention.
            system (ODESystem): Needed only when looking up by name.
            scheme (desr.ode_translation.ReductionScheme, optional): Likewise.

        Returns:
            int: The column index, in the same numbering the multiplier operations use.

        Raises:
            ValueError: If nothing matches, listing what the columns are.

        **Monomials are tried first**, then names.  A monomial is what the mathematics
        points at, and it needs neither a system nor a scheme to resolve.

        >>> equations = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
        ...              'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.column_of('e_0*k_1*t')
        4
        >>> translator.column_of('c/e_0')
        5

        By name, which needs the system so that the scheme can be worked out:

        >>> translator.column_of('tau', system)
        4
        >>> translator.column_of('nu0', system)
        5

        And an unknown one says what is on offer:

        >>> translator.column_of('nonsense')
        Traceback (most recent call last):
            ...
        ValueError: nonsense is not a column of this Hermite multiplier.  The monomials are 1/(e_0*k_1), e_0, k_2/(e_0*k_1), k_m1/(e_0*k_1), e_0*k_1*t, c/e_0, s/e_0.

        Both companions are indexed *by column*, so the whole table -- which is what a
        write-up wants, column number beside what the column is -- is one zip:

        >>> for column, (name, monomial) in enumerate(zip(translator.column_names(system),
        ...                                               translator.column_monomials())):
        ...     print(column, name, '=', monomial)
        0 x0 = 1/(e_0*k_1)
        1 x1 = e_0
        2 kappa0 = k_2/(e_0*k_1)
        3 kappa1 = k_m1/(e_0*k_1)
        4 tau = e_0*k_1*t
        5 nu0 = c/e_0
        6 nu1 = s/e_0
        '''
        if isinstance(monomial_or_name, str):
            monomial_or_name = sympify_expression(monomial_or_name)

        monomials = self.column_monomials()
        for column, monomial in enumerate(monomials):
            if monomial == monomial_or_name:
                return column

        # A name lookup needs the system, to know which scheme's names to generate, so it is
        # only attempted when one was given.  A monomial lookup never needs it.
        if monomial_or_name.is_Symbol and system is not None:
            for column, name in enumerate(self.column_names(system, scheme)):
                if name == monomial_or_name:
                    return column

        raise ValueError(
            '{} is not a column of this Hermite multiplier.  The monomials are {}.'.format(
                monomial_or_name, ', '.join(map(str, monomials))))

    def multiplier_swap_columns(self, i, j):
        '''
        Swap columns i and j of the Hermite multiplier, checking that we are only swapping 0 columns of the scaling
        matrix.

        Args:
            i (int): Column index
            j (int): Column index

        >>> translator = ODETranslator(sympy.Matrix([[1, 0, 1, 1, -1],
        ...                                            [0, 1, 0, -1, 1]]))
        >>> translator.herm_form
        Matrix([
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0]])
        >>> translator.herm_mult
        Matrix([
        [1, 0, -1, -1,  1],
        [0, 1,  0,  1, -1],
        [0, 0,  1,  0,  0],
        [0, 0,  0,  1,  0],
        [0, 0,  0,  0,  1]])
        >>> translator.inv_herm_mult
        Matrix([
        [1, 0, 1,  1, -1],
        [0, 1, 0, -1,  1],
        [0, 0, 1,  0,  0],
        [0, 0, 0,  1,  0],
        [0, 0, 0,  0,  1]])

        >>> translator.multiplier_swap_columns(2, 3)

        >>> translator.herm_mult
        Matrix([
        [1, 0, -1, -1,  1],
        [0, 1,  1,  0, -1],
        [0, 0,  0,  1,  0],
        [0, 0,  1,  0,  0],
        [0, 0,  0,  0,  1]])

        Check we have recalculated the inverse after a column operation.

        >>> translator.inv_herm_mult
        Matrix([
        [1, 0, 1,  1, -1],
        [0, 1, 0, -1,  1],
        [0, 0, 0,  1,  0],
        [0, 0, 1,  0,  0],
        [0, 0, 0,  0,  1]])

        >>> translator.multiplier_swap_columns(1, 3)
        Traceback (most recent call last):
            ...
        ValueError: Cannot swap non-zero column 1

        >>> translator.multiplier_swap_columns(2, 6)
        Traceback (most recent call last):
            ...
        IndexError: Index out of range: a[6]

        Sympy's natural column accessing is unhappy with negative indices, so make sure we don't pick up any bad habits

        >>> translator = ODETranslator(sympy.Matrix([[1, 0, 1, 0], [0, 1, 0, 1]]))
        >>> translator.herm_mult
        Matrix([
        [1, 0, -1,  0],
        [0, 1,  0, -1],
        [0, 0,  1,  0],
        [0, 0,  0,  1]])
        >>> translator.multiplier_swap_columns(2, -1)
        >>> translator.herm_mult
        Matrix([
        [1, 0,  0, -1],
        [0, 1, -1,  0],
        [0, 0,  0,  1],
        [0, 0,  1,  0]])
        '''
        if i < 0: i += self.herm_mult.cols
        if j < 0: j += self.herm_mult.cols
        if i == j:
            return
        if not self.herm_form.col(i).is_zero_matrix:
            raise ValueError('Cannot swap non-zero column {}'.format(i))
        if not self.herm_form.col(j).is_zero_matrix:
            raise ValueError('Cannot swap non-zero column {}'.format(j))
        self._herm_mult.col_swap(i, j)
        # Blow the cache of the inverse
        self._inv_herm_mult = None
        self._dep_var_inv_herm_mult_cache = {}

    def multiplier_add_columns(self, i, j, power=1):
        '''
        Add column j to column i of the Hermite multiplier, `power` times:
        :math:`C_i \\leftarrow C_i + \\textrm{power} \\cdot C_j`.

        **What `power` means.** The columns are exponent vectors, so adding them multiplies
        the monomials they stand for.  In terms of the invariants rather than the matrix,
        this does

        .. math::

            \\textrm{invariant}_i \\leftarrow \\textrm{invariant}_i
                \\cdot \\textrm{invariant}_j^{\\,\\textrm{power}}

        -- which is why it is a *power* and not a factor: the factor is
        :math:`\\textrm{invariant}_j^{\\textrm{power}}`.  So ``power=1`` multiplies by
        invariant *j*, ``power=-1`` divides by it, and ``power=2`` multiplies by its square.
        One is the default, that being the plain "add this column once".

        Column j must be an invariant column -- one that :attr:`~scaling_matrix` sends to
        zero -- since then :math:`AV` is unchanged by the operation and :math:`V` stays
        unimodular.  Column i may be anything: adding an invariant to another invariant
        changes the basis of invariants, and adding an invariant to an auxiliary column
        changes the choice of auxiliary, which is just as free.  Adding an *auxiliary*
        column anywhere would change :math:`AV`, and is refused.

        :meth:`column_of` is how to get `i` and `j` without writing a bare index.

        Args:
            i (int): Column index of the column to change.
            j (int): Column index of the invariant column to add.
            power (int): How many times to add it, which is the power of invariant *j*
                that invariant *i* is multiplied by.  Defaults to 1.

        >>> translator = ODETranslator(sympy.Matrix([[1, 0, 1, 1, -1],
        ...                                            [0, 1, 0, -1, 1]]))
        >>> translator.herm_form
        Matrix([
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0]])
        >>> translator.herm_mult
        Matrix([
        [1, 0, -1, -1,  1],
        [0, 1,  0,  1, -1],
        [0, 0,  1,  0,  0],
        [0, 0,  0,  1,  0],
        [0, 0,  0,  0,  1]])
        >>> translator.inv_herm_mult
        Matrix([
        [1, 0, 1,  1, -1],
        [0, 1, 0, -1,  1],
        [0, 0, 1,  0,  0],
        [0, 0, 0,  1,  0],
        [0, 0, 0,  0,  1]])
        >>> translator.multiplier_add_columns(2, 3)

        >>> translator.herm_mult
        Matrix([
        [1, 0, -2, -1,  1],
        [0, 1,  1,  1, -1],
        [0, 0,  1,  0,  0],
        [0, 0,  1,  1,  0],
        [0, 0,  0,  0,  1]])

        Check we have recalculated the inverse after a column operation.

        >>> translator.inv_herm_mult
        Matrix([
        [1, 0,  1,  1, -1],
        [0, 1,  0, -1,  1],
        [0, 0,  1,  0,  0],
        [0, 0, -1,  1,  0],
        [0, 0,  0,  0,  1]])

        An invariant column may be added to an auxiliary column too; that changes which
        auxiliary is used, and leaves :math:`AV` alone.

        >>> translator.multiplier_add_columns(1, 3)
        >>> translator.herm_mult.col(1).T
        Matrix([[-1, 2, 0, 1, 0]])
        >>> translator.scaling_matrix * translator.herm_mult == translator.herm_form
        True

        But an auxiliary column cannot be added to anything, since that would change
        :math:`AV`.

        >>> translator.multiplier_add_columns(3, 1)
        Traceback (most recent call last):
            ...
        ValueError: Cannot add column 1: it is an auxiliary column, not an invariant column, so adding it would change the Hermite form.

        >>> translator.multiplier_add_columns(3, 3)
        Traceback (most recent call last):
            ...
        ValueError: Cannot add column 3 to itself

        This is how a choice of invariants other than the normal one is reached.  The normal
        Hermite multiplier of example 6.6 of :cite:`Hubert2013c` gives the paper's own
        answer, so there is nothing to correct here any more -- these operations move *away*
        from it, to an equally valid choice differing by an invariant.

        >>> equations = ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
        ...              'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)']
        >>> system = ODESystem.from_equations(equations)
        >>> system.reorder_variables(['t', 'z1', 'z2'])
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.auxiliary_monomials(), translator.invariant_monomials()
        (Matrix([[t*z1**2]]), Matrix([[t*z1**3, z1**2*z2/t]]))
        >>> translator.multiplier_add_columns(2, 1, power=-1)
        >>> translator.multiplier_add_columns(0, 2, power=-1)
        >>> translator.auxiliary_monomials(), translator.invariant_monomials()
        (Matrix([[t**3*z1**3/z2]]), Matrix([[t*z1**3, z2/(t**2*z1)]]))

        Sympy's natural column accessing is unhappy with negative indices, so make sure we don't pick up any bad habits

        >>> translator = ODETranslator(sympy.Matrix([[1, 0, 1, 0], [0, 1, 0, 1]]))
        >>> translator.herm_mult
        Matrix([
        [1, 0, -1,  0],
        [0, 1,  0, -1],
        [0, 0,  1,  0],
        [0, 0,  0,  1]])
        >>> translator.multiplier_add_columns(2, -1)
        >>> translator.herm_mult
        Matrix([
        [1, 0, -1,  0],
        [0, 1, -1, -1],
        [0, 0,  1,  0],
        [0, 0,  1,  1]])
        '''
        if i < 0: i += self.herm_mult.cols
        if j < 0: j += self.herm_mult.cols
        if i == j:
            raise ValueError('Cannot add column {} to itself'.format(i))
        if not self.herm_form.col(j).is_zero_matrix:
            raise ValueError('Cannot add column {}: it is an auxiliary column, not an '
                             'invariant column, so adding it would change the Hermite '
                             'form.'.format(j))
        self._herm_mult.col_op(i, lambda v, index: v + power * self._herm_mult[index, j])
        # Blow the cache of the inverse
        self._inv_herm_mult = None
        self._dep_var_inv_herm_mult_cache = {}

    def multiplier_negate_column(self, i):
        '''
        Negate column i:
        C_i <- -C_i
        Check that we are only operating with 0 columns of the scaling matrix.

        Args:
            i (int): Column index

        >>> translator = ODETranslator(sympy.Matrix([[1, 0, 1, 1, -1],
        ...                                            [0, 1, 0, -1, 1]]))
        >>> translator.herm_form
        Matrix([
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0]])
        >>> translator.herm_mult
        Matrix([
        [1, 0, -1, -1,  1],
        [0, 1,  0,  1, -1],
        [0, 0,  1,  0,  0],
        [0, 0,  0,  1,  0],
        [0, 0,  0,  0,  1]])
        >>> translator.inv_herm_mult
        Matrix([
        [1, 0, 1,  1, -1],
        [0, 1, 0, -1,  1],
        [0, 0, 1,  0,  0],
        [0, 0, 0,  1,  0],
        [0, 0, 0,  0,  1]])
        >>> translator.multiplier_negate_column(3)

        >>> translator.herm_mult
        Matrix([
        [1, 0, -1,  1,  1],
        [0, 1,  0, -1, -1],
        [0, 0,  1,  0,  0],
        [0, 0,  0, -1,  0],
        [0, 0,  0,  0,  1]])

        Check we have recalculated the inverse after a column operation.

        >>> translator.inv_herm_mult
        Matrix([
        [1, 0, 1,  1, -1],
        [0, 1, 0, -1,  1],
        [0, 0, 1,  0,  0],
        [0, 0, 0, -1,  0],
        [0, 0, 0,  0,  1]])

        >>> translator.multiplier_negate_column(1)
        Traceback (most recent call last):
            ...
        ValueError: Cannot negate non-zero column 1

        Sympy's natural column accessing is unhappy with negative indices, so make sure we don't pick up any bad habits

        >>> translator = ODETranslator(sympy.Matrix([[1, 0, 1], [0, 1, 0]]))
        >>> translator.herm_mult
        Matrix([
        [1, 0, -1],
        [0, 1,  0],
        [0, 0,  1]])
        >>> translator.multiplier_negate_column(-1)
        >>> translator.herm_mult
        Matrix([
        [1, 0,  1],
        [0, 1,  0],
        [0, 0, -1]])
        '''
        if i < 0: i += self.herm_mult.cols
        if not self.herm_form.col(i).is_zero_matrix:
            raise ValueError('Cannot negate non-zero column {}'.format(i))
        self._herm_mult.col_op(i, lambda v, index: -v)
        # Blow the cache of the inverse
        self._inv_herm_mult = None
        self._dep_var_inv_herm_mult_cache = {}

    def translate(self, system, naming_scheme=None, include_aux_vars=True):
        '''
        Translate an :class:`ode_system.ODESystem` into a reduced system.

        First, try the simplest parameter reduction method, then the dependent variable translation (where the scaling action ignores the independent variable) and finally the general reduction scheme.

        Args:
            system (ODESystem): System to reduce.
            naming_scheme (dict): the scheme to use, over the keys in
                :data:`NAMING_SCHEME_KEYS`.  See

        :rtype: ODESystem
        '''

        if naming_scheme is not None:
            self.set_naming_scheme(naming_scheme)

        scheme = self.scheme(system)

        if scheme == ReductionScheme.PARAMETER:
            return self.translate_parameter(system=system)

        if scheme == ReductionScheme.DEP_VAR:
            return self.translate_dep_var(system=system)

        return self.translate_general(system=system, include_aux_vars=include_aux_vars)

    def scheme(self, system):
        '''
        Which reduction :meth:`translate` would use for `system`.

        The three are tried in order of how much they leave intact: the parameter scheme
        first, then the dependent-variable scheme, then the general one, which always
        applies.  :meth:`translate` asks this and dispatches on the answer, so there is one
        copy of the decision rather than two that can disagree.

        Args:
            system (ODESystem): The system to be reduced.

        Returns:
            ReductionScheme

        Raises:
            ValueError: if the system's variables do not match the scaling matrix at all,
                in which case no scheme applies.

        >>> equations = 'dn/dt = n*( r*(1 - n/K) - k*p/(n+d) );dp/dt = s*p*(1 - h*p / n)'.split(';')
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.scheme(system)
        ReductionScheme.PARAMETER

        Hubert & Labahn 6.6 has no parameters at all, and its clock is entangled with a
        dependent variable, so only the general scheme applies:

        >>> hl66 = ODESystem.from_equations(['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
        ...                                  'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])
        >>> ODETranslator.from_ode_system(hl66).scheme(hl66)
        ReductionScheme.GENERAL
        '''
        if self._is_translate_parameter_compatible(system=system):
            return ReductionScheme.PARAMETER

        if ((len(system.symbols) == self.scaling_matrix.shape[1] + 1)
                or (self.scaling_matrix[:, system.independent_variable_index].is_zero_matrix)):
            return ReductionScheme.DEP_VAR

        if len(system.symbols) == self.scaling_matrix.shape[1]:
            return ReductionScheme.GENERAL

        raise ValueError("System doesn't have the right number of variables for translation")

    def expected_variable_count(self, system, scheme=None, include_aux_vars=True):
        '''
        How many variables the reduction of `system` will have.

        The count is scheme-dependent, with :math:`n` the number of variables of `system`
        -- the independent variable included -- and :math:`r` the number of scaling
        symmetries:

        =========================================  =================
        scheme                                     variables
        =========================================  =================
        parameter                                  :math:`n - r`
        dependent-variable                         :math:`n`
        general, with auxiliaries                  :math:`n + 1`
        general, ``include_aux_vars=False``        :math:`n - r + 1`
        =========================================  =================

        **Where the general scheme's extra one comes from.**  It keeps the original
        independent variable, and that is the :math:`+1`.  The remaining :math:`n` are a
        full change of variables -- :math:`r` auxiliaries and :math:`n - r` invariants --
        so the count is :math:`1 + r + (n - r)`.  Dropping the auxiliaries removes exactly
        :math:`r`, giving :math:`n - r + 1`.

        The parameter scheme has no :math:`+1` because nothing is kept: its :math:`r`
        auxiliaries are normalised to 1 and discarded, and the invariant that is a rescaled
        clock *becomes* the new independent variable rather than joining one.  So
        :math:`(p - r)` parameters, one clock and :math:`q` dependents, which is
        :math:`n - r`.

        Under the general scheme a rescaled clock can therefore appear *alongside* the
        original one -- prey-predator keeps :math:`t` and also carries :math:`k t / h` as a
        variable with an equation.  It need not: in Hubert & Labahn 6.6 both invariants
        involve a dependent variable, so there is no second clock, and the :math:`+1` is
        the retained :math:`t` all the same.  :meth:`indep_var_invariant` is what tells the
        two situations apart.

        Args:
            system (ODESystem): The system to be reduced.
            scheme (desr.ode_translation.ReductionScheme, optional): The scheme that will be used.  **Pass this
                whenever you choose the scheme yourself.**  Left out, it is inferred by
                :meth:`scheme`, which is only the right answer when :meth:`translate` does
                the choosing: calling :meth:`translate_general` on a system that
                :meth:`translate` would have given to another scheme is legal, and gives a
                different count.  Inferring where the caller had already chosen is a
                mistake that produces a plausible wrong number rather than an error, and it
                has been made twice in this library's history.
            include_aux_vars (bool): Whether the auxiliaries will be kept.  Only the general
                scheme has any, so this is ignored by the other two.

        Returns:
            int

        >>> equations = 'dn/dt = n*( r*(1 - n/K) - k*p/(n+d) );dp/dt = s*p*(1 - h*p / n)'.split(';')
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> len(system.symbols), translator.r
        (9, 3)
        >>> translator.expected_variable_count(system)
        6
        >>> len(translator.translate(system).symbols)
        6

        The general scheme on the same system gives a different count, which is why it has
        to be said rather than inferred:

        >>> translator.expected_variable_count(system, scheme=ReductionScheme.GENERAL)
        10
        >>> len(translator.translate_general(system).symbols)
        10
        >>> translator.expected_variable_count(system, scheme='general', include_aux_vars=False)
        7
        >>> len(translator.translate_general(system, include_aux_vars=False).symbols)
        7
        '''
        if scheme is None:
            scheme = self.scheme(system)
        else:
            scheme = ReductionScheme(scheme)

        n = len(system.symbols)

        if scheme == ReductionScheme.PARAMETER:
            return n - self.r
        if scheme == ReductionScheme.DEP_VAR:
            return n
        return n + 1 if include_aux_vars else n - self.r + 1



    def translate_parameter(self, system):
        ''' Translate according to parameter scheme '''

        # make the substitutions
        to_sub = self.translate_parameter_substitutions(system=system)

        variable_map = self.variable_map(system)
        new_indep = variable_map[system.independent_variable]

        new_deriv_dict = {}
        # derivative_dict holds only the variables that have an equation, so the parameters
        # are already absent and need no skipping.
        for var, deriv in system.derivative_dict.items():
            Var = variable_map[var]
            # the chain rule makes substitution more complicated
            derivative_factor = (to_sub[var] / Var) * (new_indep / to_sub[system.independent_variable])

            new_deriv_dict[Var] = (deriv.subs(to_sub) / derivative_factor).expand()

        # Build the reduced system from the variables the reduction is known to produce,
        # rather than from the symbols that happen to survive into the right-hand sides.
        # An invariant constant that cancels out of every equation is still a variable of
        # the reduced system -- it is one of the n - r coordinates, and the substitutions
        # returned by translate_parameter_substitutions refer to it -- but inferring the
        # variables from the expressions, as ODESystem.from_dict does, would lose it.
        dependents = [variable_map[var] for var in system.dependents]
        new_constants = self._new_constants(system)
        new_variables = [new_indep] + dependents + new_constants
        new_derivatives = ([new_deriv_dict[new_indep]]
                           + [new_deriv_dict[var] for var in dependents]
                           + [None] * len(new_constants))

        # The parameter scheme normalises its auxiliaries to 1 and discards them, so there
        # are none to record.
        reduced_system = ODESystem(new_variables, new_derivatives, indep_var=new_indep,
                                   is_reduced=True,
                                   scheme=ReductionScheme.PARAMETER,
                                   num_scaling_symmetries=self.r)
        reduced_system.default_order_variables()

        # Since z(0) / z_0 is a rational invariant, we can rewrite it using to_sub substitutions.
        if system.initial_conditions:
            # For parameter substitution, the dependent variables keep their name
            new_initial_conditions = {variable_map[var]: (init_cond.subs(to_sub) * variable_map[var] / to_sub[var]).expand()
                                      for var, init_cond in system.initial_conditions.items()}
            reduced_system.update_initial_conditions(new_initial_conditions)

        if system.constraints:
            for eqn in system.constraints:
                reduced_system.add_constraint(eqn.lhs.subs(to_sub), eqn.rhs.subs(to_sub))

        # rename the time variable
        reduced_system.rename_indep_var(new_indep)

        return reduced_system



    def translate_dep_var(self, system):
        '''
        Given a system of ODEs, translate the system into a simplified version. Assume we are only working on
        dependent variables, not the independent variable.

        Args:
            system (ODESystem): System to reduce.

        :rtype: ODESystem

        >>> equations = 'dz1/dt = z1*(1+z1*z2);dz2/dt = z2*(1/t - z1*z2)'.split(';')
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.translate_dep_var(system=system)
        dt/dt = 1
        dx0/dt = x0*(-y0 - 1)
        dy0/dt = y0*(1 + 1/t)
        '''
        if system.initial_conditions:
            raise NotImplementedError('`translate_dep_var` not yet implemented for systems with initial conditions')
        # First check that our scaling action doesn't act on the independent variable
        if self.n == len(system.symbols):
            if not self.scaling_matrix[:, system.independent_variable_index].is_zero_matrix:
                raise ValueError('The independent variable of\n{}\nis acted on by\n{}.'.format(system, self))
            scaling_matrix = self.scaling_matrix.copy()
            scaling_matrix.col_del(system.independent_variable_index)
            new_herm_mult = self.dep_var_herm_mult(indep_var_index=system.independent_variable_index)
        else:
            assert self.scaling_matrix.shape[1] == len(system.symbols) - 1
            scaling_matrix = self.scaling_matrix.copy()
            new_herm_mult = self.dep_var_herm_mult(indep_var_index=system.independent_variable_index)
            assert new_herm_mult[-1, :].is_zero_matrix
            new_herm_mult.col_del(system.independent_variable_index)
            new_herm_mult.row_del(system.independent_variable_index)

        if self.symbols_domain is not None:
            symbols_domain = list(self.symbols_domain)
            symbols_domain.pop(system.independent_variable_index)
        else:
            symbols_domain = None

        # silviana says: i think it is weird that an ODETranslator is made here.  we're already in one.  why make another one?
        reduced_scaling = ODETranslator(scaling_matrix=scaling_matrix,
                                         symbols_domain=symbols_domain,
                                         hermite_multiplier=new_herm_mult)

        # y = sympy.Matrix(scale_action(system.symbols, self.herm_mult_n))
        #print 'y = ', sympy.Matrix(scale_action(system.symbols, self.herm_mult_n))
        scheme = self._resolved_naming_scheme(GENERAL_SCHEME_NAMES)
        num_inv_var = reduced_scaling.herm_mult_n.shape[1]
        invariant_names = _names_for(scheme['dependents'], num_inv_var,
                                         self._new_indices_start_at, 'invariants')

        # x = sympy.Matrix(scale_action(system.symbols, self.herm_mult_i))
        #print 'x = ', sympy.Matrix(scale_action(system.symbols, self.herm_mult_i))
        num_aux_var = reduced_scaling.herm_mult_i.shape[1]
        auxiliary_names = _new_symbols(scheme['auxiliaries'], num_aux_var,
                                           self._new_indices_start_at)

        system_var_no_indep = list(system.symbols)
        system_var_no_indep.pop(system.independent_variable_index)
        to_sub = dict(zip(system_var_no_indep, scale_action(invariant_names, reduced_scaling.inv_herm_mult_d)))
        system_derivatives = _rates_of_change(system)
        # fywd = F(y^(W_d)) in Hubert Labahn
        fywd = sympy.Matrix([(f / v).subs(to_sub, simultaneous=True).expand() for v, f in
                            zip(system.symbols, system_derivatives)]).T
        fywd.col_del(system.independent_variable_index)
        dydt = _scaled_rates(invariant_names, fywd * reduced_scaling.herm_mult_n)
        dxdt = _scaled_rates(auxiliary_names, fywd * reduced_scaling.herm_mult_i)

        new_variables = [system.independent_variable] + list(auxiliary_names) + list(invariant_names)

        # The dependent-variable scheme's columns are indexed by the symbols other than the
        # independent variable, so that is what the monomials are built from -- and why none
        # of them can involve the clock.
        dynamic = set(system.dependents)
        monomials = (list(scale_action(system_var_no_indep, reduced_scaling.herm_mult_i))
                     + list(scale_action(system_var_no_indep, reduced_scaling.herm_mult_n)))
        rates = _classify_reduced_derivatives(
            list(auxiliary_names) + list(invariant_names), monomials,
            list(dxdt) + list(dydt), dynamic)
        new_derivatives = [sympy.sympify(1)] + rates

        reduced = ODESystem(new_variables, new_derivatives,
                            indep_var=system.independent_variable, is_reduced=True,
                            scheme=ReductionScheme.DEP_VAR,
                            num_scaling_symmetries=self.r,
                            auxiliary_names=auxiliary_names)
        return _rename_indep_var(reduced, scheme['independent'])

    def translate_general(self, system, include_aux_vars = True):
        '''
        The most general reduction scheme.
        If there are :math:`n` variables (including the independent variable) then there will be a system of
        :math:`n-r+1` invariants and :math:`r` auxiliary variables.

        Args:
            system (ODESystem): System to reduce.

        :rtype: ODESystem

        >>> equations = 'dn/dt = n*( r*(1 - n/K) - k*p/(n+d) );dp/dt = s*p*(1 - h*p / n)'.split(';')
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)

        Six of the ten symbols are monomials in the original parameters alone, so they are
        reported as parameters and only three equations need integrating:

        >>> translator.translate_general(system=system)
        dt/dt = 1
        dy3/dt = y3/t
        dy4/dt = y4*(-y1*y3*y4 + y1*y3 - y3*y5/(y0 + y4))/t
        dy5/dt = y5*(y2*y3 - y2*y3*y5/y4)/t
        parameters: x0, x1, x2, y0, y1, y2

        >>> translator.translate_general(system=system, include_aux_vars=False)
        dt/dt = 1
        dy3/dt = y3/t
        dy4/dt = y4*(-y1*y3*y4 + y1*y3 - y3*y5/(y0 + y4))/t
        dy5/dt = y5*(y2*y3 - y2*y3*y5/y4)/t
        parameters: y0, y1, y2

        The auxiliaries are recorded whatever kind they turn out to be:

        >>> reduced = translator.translate_general(system=system)
        >>> reduced.auxiliary_names
        (x0, x1, x2)
        >>> reduced.scheme, reduced.num_scaling_symmetries
        (ReductionScheme.GENERAL, 3)
        '''
        if system.initial_conditions:
            raise NotImplementedError('General translation not yet implemented for systems with initial conditions')
        # y = sympy.Matrix(scale_action(system.symbols, self.herm_mult_n))
        #print 'y = ', sympy.Matrix(scale_action(system.symbols, self.herm_mult_n))
        scheme = self._resolved_naming_scheme(GENERAL_SCHEME_NAMES)
        num_inv_var = self.herm_mult_n.shape[1]
        invariant_names = _names_for(scheme['dependents'], num_inv_var,
                                         self._new_indices_start_at, 'invariants')

        if include_aux_vars:
            # x = sympy.Matrix(scale_action(system.symbols, self.herm_mult_i))
            #print 'x = ', sympy.Matrix(scale_action(system.symbols, self.herm_mult_i))
            num_aux_var = self.herm_mult_i.shape[1]
            auxiliary_names = _new_symbols(scheme['auxiliaries'], num_aux_var,
                                               self._new_indices_start_at)

        to_sub = dict(zip(system.symbols, scale_action(invariant_names, self.inv_herm_mult_d)))
        system_derivatives = _rates_of_change(system)

        # fywd = F(y^(W_d)) in Hubert Labahn
        fywd = sympy.Matrix([(system.independent_variable * f / v).subs(to_sub, simultaneous=True).expand() for v, f in
                            zip(system.symbols, system_derivatives)]).T
        dydt = _scaled_rates(invariant_names, fywd * self.herm_mult_n,
                             divide_by=system.independent_variable)

        if include_aux_vars:
            dxdt = _scaled_rates(auxiliary_names, fywd * self.herm_mult_i,
                                 divide_by=system.independent_variable)

        dynamic = set(system.dependents) | {system.independent_variable}

        # scale_action against system.symbols, not the *_monomials() accessors:
        # those go through _validate_variables and fall back to placeholder x_i and y_i
        # names when symbols_domain is unset, which would pair each name with somebody
        # else's monomial.
        if include_aux_vars:
            named = list(auxiliary_names) + list(invariant_names)
            monomials = (list(scale_action(system.symbols, self.herm_mult_i))
                         + list(scale_action(system.symbols, self.herm_mult_n)))
            rates = list(dxdt) + list(dydt)
            new_variables = [system.independent_variable] + named
            assert len(new_variables) == len(rates) + 1 == len(system.symbols) + 1
        else:
            # don't include the aux vars.
            named = list(invariant_names)
            monomials = list(scale_action(system.symbols, self.herm_mult_n))
            rates = list(dydt)
            new_variables = [system.independent_variable] + named
            # see remark 3.2 in the second arxiv version of the paper, where we discuss this number
            assert len(new_variables) == len(rates) + 1 == len(system.symbols) - self.r + 1

        new_derivatives = ([sympy.sympify(1)]
                           + _classify_reduced_derivatives(named, monomials, rates, dynamic))

        reduced = ODESystem(new_variables, new_derivatives,
                            indep_var=system.independent_variable, is_reduced=True,
                            scheme=ReductionScheme.GENERAL,
                            num_scaling_symmetries=self.r,
                            auxiliary_names=auxiliary_names if include_aux_vars else ())
        return _rename_indep_var(reduced, scheme['independent'])

    def _is_translate_parameter_compatible(self, system):
        '''
        Check whether a system satisfies the conditions of the parameter scheme of translation.

        The scheme needs each of the independent variable and the dependent variables to
        contribute exactly one invariant of its own, and the surplus parameters to become
        the reduced constants.  That is what the block form below says, reading
        :math:`V = [V_i, V_n]` with :math:`p` parameters and :math:`m = 1 + q` the rest::

                       r aux        p - r          m
                     +-----------+-----------+-----------+
            p rows   |     *     |     *     |     *     |
                     +-----------+-----------+-----------+
            m rows   |     0     |     0     |    I_m    |
                     +-----------+-----------+-----------+

        The zero block says the scaling acts only through the parameters, leaving the
        independent and dependent variables unrescaled; the identity says each of them
        appears in exactly one invariant, with exponent one and in no other.

        Those rows are the **last** :math:`m`, which is why the variables must be ordered
        :math:`(c, t, z)` -- see
        :meth:`~desr.ode_system.ODESystem.default_order_variables`.  The order is not a
        free choice beside the Hermite normal form convention but the other half of it:
        under the mirrored convention the same blocks sit at the top and the order is
        :math:`(t, z, c)`.  Ordering the parameters last *here* would satisfy every block
        condition while putting the identity on parameters, so each parameter would get its
        own invariant and the independent variable would be smeared across several.  The
        test therefore pins which variables occupy those rows, not merely the shape.
        '''
        p = system.num_parameters
        m = len(system.symbols) - p

        # First check our system's variables are in the required order:
        # parameters, then the independent variable, then the dependents.
        if system.independent_variable_index != p:
            return False
        # Ask which variables are the parameters, rather than which derivatives are zero: a
        # conserved variable's derivative is also zero, so testing the derivative would call
        # a system eligible as soon as one was ordered among the parameters.
        if p and set(system.symbols[:p]) != set(system.parameters):
            return False

        # Now check our transformation is valid: that V and W have the required forms.
        if not self.herm_mult_i[p:, :].is_zero_matrix:
            return False
        if not self.herm_mult_n[p:, p - self.r:] == sympy.eye(m):
            return False
        if not self.herm_mult_n[p:, :p - self.r].is_zero_matrix:
            return False

        if not (self.inv_herm_mult[p:, :] ==
                sympy.Matrix.hstack(sympy.zeros(m, p), sympy.eye(m))):
            return False

        return True

    def variable_map(self, system):
        """
        Construct a dictionary mapping old variable names to new ones, using the stored renaming scheme.

        The naming scheme is passed in at construct time, not here.

        todo: use sympify so that erroneous naming schemes can be detected, but this requires using the _clash1 namespace or something.
        """
        scheme = self._resolved_naming_scheme(PARAMETER_SCHEME_NAMES)
        num_variables_wo_time = len(system.symbols) - system.num_parameters - 1  # Excluding indep

        new_vars = _names_for(scheme['dependents'], num_variables_wo_time,
                              self._new_indices_start_at, 'dependent variables')

        independent = scheme['independent']
        the_map = {system.independent_variable: system.independent_variable if independent is None
                                     else sympy.Symbol(independent)}

        for v,V in zip(system.dependents, new_vars):
            the_map[v] = V

        return the_map

    def _new_constants(self, system):
        '''
        The invariant constants of the reduced system, named by the naming scheme.

        There are :math:`k - r` of them, where :math:`k` is the number of constants of the
        original system: the reduction normalises :math:`r` of them away.  This is the one
        place the list is built, so that :meth:`translate_parameter` and
        :meth:`translate_parameter_substitutions` cannot disagree about what exists.

        Args:
            system (ODESystem): The system being reduced.

        :rtype: list of sympy.Symbol
        '''
        scheme = self._resolved_naming_scheme(PARAMETER_SCHEME_NAMES)
        return _new_symbols(scheme['parameters'], system.num_parameters - self.r,
                            self._new_indices_start_at)

    def translate_parameter_substitutions(self, system):
        '''
        Given a system, determine the substitutions made in the parameter reduction.

        Args:
            system (ODESystem):
                The system in question.

        The substitutions come back in the system's own variable order -- parameters, then
        the independent variable, then the dependents -- so they can be read straight down
        against :attr:`~desr.ode_system.ODESystem.symbols`.

        :rtype: dict

        >>> equations = 'dn/dt = n*( r*(1 - n/K) - k*p/(n+d) );dp/dt = s*p*(1 - h*p / n)'.split(';')
        >>> system = ODESystem.from_equations(equations)
        >>> system.symbols
        (K, d, h, k, r, s, t, n, p)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.translate_parameter_substitutions(system=system)
        {K: 1, d: kappa0, h: 1, k: 1, r: kappa1, s: kappa2, t: tau, n: nu0, p: nu1}

        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system, naming_scheme={'parameters': 'c', 'independent': 't', 'dependents': ['n','p']})
        >>> translator.translate_parameter_substitutions(system=system)
        {K: 1, d: c0, h: 1, k: 1, r: c1, s: c2, t: t, n: n, p: p}

        >>> translator = ODETranslator.from_ode_system(system, new_indices_start_at=1)
        >>> translator.translate_parameter_substitutions(system=system)
        {K: 1, d: kappa1, h: 1, k: 1, r: kappa2, s: kappa3, t: tau, n: nu1, p: nu2}
        '''

        p = system.num_parameters
        m = len(system.symbols) - p  # The independent variable and the dependents

        if not self._is_translate_parameter_compatible(system):
            err_str = ['System is not compatible for parameter translation.']
            err_str.append('System may not be ordered properly. Must be parameters, independent, dependent. Order is: {}'.format(system.symbols))
            err_str.append('Transformation may not be appropriate. Should be 0, 0, I. Transformation is:')
            err_str.append(repr(self.herm_mult_i[p:, :]))
            err_str.append(repr(self.herm_mult_n[p:, :p - self.r]))
            err_str.append(repr(self.herm_mult_n[p:, p - self.r:]))
            raise ValueError('\n'.join(err_str))

        # Extract the right bits of W.  The rows of W_d that belong to the reduced constants
        # are the first p - r of them, since those invariants are the leading columns of
        # V_n under this ordering; the columns follow the variable order (c, t, z).
        inv_herm_mult_d = self.inv_herm_mult_d[:p - self.r, :]
        W_c = inv_herm_mult_d[:, :p]
        W_t = inv_herm_mult_d[:, p:p + 1]
        W_v = inv_herm_mult_d[:, p + 1:]

        variable_map = self.variable_map(system)
        # Form new constants
        new_consts = self._new_constants(system)

        # Built in the system's own variable order -- parameters, independent, dependents --
        # so that reading the substitutions against system.symbols needs no reshuffling.
        to_sub = {}

        # Scale parameters
        scale = scale_action(new_consts, W_c)
        assert len(system.parameters) == len(scale.T)
        for const, _const_scale in zip(system.parameters, scale.T):
            to_sub[const] = _const_scale

        # Scale t
        to_sub[system.independent_variable] = scale_action(new_consts, W_t)[0] * variable_map[system.independent_variable]

        # Scale dependents
        scale = scale_action(new_consts, W_v)
        assert len(system.dependents) == len(scale.T)
        for dep_var, _const_scale in zip(system.dependents, scale.T):
            to_sub[dep_var] = _const_scale * variable_map[dep_var]

        return to_sub


    def reverse_translate(self, variables, indep_var_index=0):
        """
        Given the solutions of a reduced system, reverse translate them into solutions of the original system.

        Note:
            This *doesn't* reverse translate the parameter reduction scheme.

            It *only* guesses between :meth:`~desr.ode_translation.ODETranslator.reverse_translate_general`
            and :meth:`~desr.ode_translation.ODETranslator.reverse_translate_dep_var`,
            by the number of variables given: :math:`n + 1` of them, the independent variable
            first, means the general scheme; :math:`n` or :math:`n - 1` means the dependent
            variable scheme.

        Args:
            variables (iter of sympy.Expression): The solution auxiliary variables :math:`x(t)` and solution invariants
                :math:`y(t)` of the reduced system.
                Variables should be ordered auxiliary variables followed by invariants.
            indep_var_index (int): The location of the independent variable.

        :rtype: tuple

        The dependent variable scheme, on example 6.4 of :cite:`Hubert2013c`:

        >>> equations = 'dz1/dt = z1*(1+z1*z2);dz2/dt = z2*(1/t - z1*z2)'.split(';')
        >>> system = ODESystem.from_equations(equations)
        >>> system.reorder_variables(['t', 'z1', 'z2'])
        >>> translator = ODETranslator(scaling_matrix=sympy.Matrix([[0, 1, -1]]),
        ...                              hermite_multiplier=sympy.Matrix([[0, 1, 0],
        ...                                                               [1, 0, 1],
        ...                                                               [0, 0, 1]]))
        >>> x0, y0 = sympy.var('x0 y0')
        >>> translator.reverse_translate((x0, y0), system.independent_variable_index)
        (x0, y0/x0)

        The general scheme, on example 6.6, dispatched on the independent variable coming
        first.  This branch needs genuine solutions rather than bare symbols, since it
        checks that the reconstructed independent variable is a constant multiple of the
        reduced one; the paper's exact solution serves.

        >>> translator = ODETranslator(scaling_matrix=sympy.Matrix([[3, -1, 5]]),
        ...                              hermite_multiplier=sympy.Matrix([[1, 1, -1],
        ...                                                               [2, 3, 2],
        ...                                                               [0, 0, 1]]))
        >>> exact = (sympy.var('t'),
        ...          sympy.sympify('c3/(t**(1/3)*(ln(t-c1)-ln(t)+c2)**(2/3))'),
        ...          sympy.sympify('c1/(t*(ln(t-c1)-ln(t)+c2))'),
        ...          sympy.sympify('t/(c1 - t)'))
        >>> translator.reverse_translate(exact) == translator.reverse_translate_general(exact)
        True
        """
        if len(variables) == self.scaling_matrix.shape[1]:
            return self.reverse_translate_dep_var(variables=variables, indep_var_index=indep_var_index)
        elif len(variables) == self.scaling_matrix.shape[1] - 1:
            return self.reverse_translate_dep_var(variables=variables, indep_var_index=indep_var_index)
        elif len(variables) == self.scaling_matrix.shape[1] + 1:
            return self.reverse_translate_general(variables=variables, system_indep_var_index=indep_var_index)
        else:
            raise ValueError('Incorrect number of variables for reverse translation')

    def reverse_translate_parameter(self, variables):
        '''
        Given the solutions of a reduced system, reverse translate them into solutions of the original system.

        Args:
            variables (iter of sympy.Expression):
                :math:`r` constants (auxiliary variables) followed by the independent, dependent and invariant constants.

        Example 7.5 (Prey-predator model) from :cite:`Hubert2013c`.

        .. math::
            \\frac{dn}{dt} = n  \\left( r \\left(1 - \\frac{n}{K} \\right) - \\frac{k p}{n + d} \\right)

            \\frac{dp}{dt} = s p \\left(1 - \\frac{h p}{n} \\right)

        The paper normalises :math:`K`, :math:`h` and :math:`r` away and keeps :math:`d`,
        :math:`k` and :math:`s` as the reduced constants.  Ordering those three first is all
        it takes to reach that choice, since the leading parameters are the ones normalised;
        this example used to hand-craft a nine by nine Hermite multiplier instead, because
        desr reported the mirror of the paper's normal form and no ordering could get there.

        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )','dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> system.symbols
        (K, d, h, k, r, s, t, n, p)
        >>> system.reorder_variables(['K', 'h', 'r', 'd', 'k', 's', 't', 'n', 'p'])
        >>> translator = ODETranslator.from_ode_system(system,
        ...                                            naming_scheme={'parameters': 'c', 'independent': 't', 'dependents': ['n','p']})
        >>> translator.translate_parameter(system)
        dt/dt = 1
        dn/dt = -c1*n*p/(c0 + n) - n**2 + n
        dp/dt = c2*p - c2*p**2/n
        parameters: c0, c1, c2
        >>> translator.translate_parameter_substitutions(system)
        {K: 1, h: 1, r: 1, d: c0, k: c1, s: c2, t: t, n: n, p: p}
        >>> soln_reduced = sympy.var('x1, x2, x3, t, n, p, c0, c1, c2')
        >>> translator.reverse_translate_parameter(soln_reduced)
        Matrix([[x2, x2/x3, 1/x1, t*x2, n*x2/(x1*x3), p/x1, c0*x1, c1*x2, c2*x3]])
        '''
        return scale_action(variables, self.inv_herm_mult)

    def reverse_translate_dep_var(self, variables, indep_var_index):
        '''
        Given the solutions of a (dep_var) reduced system, reverse translate them into solutions of the original system.

        Args:
            variables (iter of sympy.Expression): The solution auxiliary variables :math:`x(t)` and solution invariants
                :math:`y(t)` of the reduced system.
                Variables should be ordered auxiliary variables followed by invariants.
            indep_var_index (int): The location of the independent variable.

        :rtype: tuple

        Example 6.4 from :cite:`Hubert2013c`.
        We will just take the matrices from the paper, rather than use our own
        (which are constructed using different conventions).

        .. math::
            \\frac{dz_1}{dt} = z_1 \\left( 1+ z_1 z_2 \\right)

            \\frac{dz_2}{dt} = z_2 \\left(\\frac{1}{t} - z_1 z_2 \\right)

        >>> equations = 'dz1/dt = z1*(1+z1*z2);dz2/dt = z2*(1/t - z1*z2)'.split(';')
        >>> system = ODESystem.from_equations(equations)
        >>> var_order = ['t', 'z1', 'z2']
        >>> system.reorder_variables(var_order)
        >>> scaling_matrix = sympy.Matrix([[0, 1, -1]])
        >>> hermite_multiplier = sympy.Matrix([[0, 1, 0],
        ...                                    [1, 0, 1],
        ...                                    [0, 0, 1]])
        >>> translator = ODETranslator(scaling_matrix=scaling_matrix,
        ...                              hermite_multiplier=hermite_multiplier)
        >>> reduced = translator.translate_dep_var(system)
        >>> # Check reduction
        >>> answer = ODESystem.from_equations('dx0/dt = x0*(y0 + 1);dy0/dt = y0*(1 + 1/t)'.split(';'))
        >>> reduced == answer
        True
        >>> reduced
        dt/dt = 1
        dx0/dt = x0*(y0 + 1)
        dy0/dt = y0*(1 + 1/t)
        >>> t_var, c1_var, c2_var = sympy.var('t c1 c2')
        >>> reduced_soln = (c2_var*sympy.exp(t_var+c1_var*(1-t_var)*sympy.exp(t_var)),
        ...                 c1_var * t_var * sympy.exp(t_var))
        >>> original_soln = translator.reverse_translate_dep_var(reduced_soln, system.independent_variable_index)
        >>> original_soln == (reduced_soln[0], reduced_soln[1] / reduced_soln[0])
        True
        >>> original_soln
        (c2*exp(c1*(1 - t)*exp(t) + t), c1*t*exp(t)*exp(-c1*(1 - t)*exp(t) - t)/c2)

        Giving all three reduced variables -- the auxiliary, then the invariants, of which
        the first is :math:`t` itself -- rebuilds all three originals:

        >>> x0, y0 = sympy.var('x0 y0')
        >>> translator.reverse_translate_dep_var((x0, t_var, y0), system.independent_variable_index)
        (t, x0, y0/x0)
        '''
        if len(variables) == self.scaling_matrix.shape[1]:
            # All n reduced variables, the independent variable's invariant among them: the
            # full inverse multiplier rebuilds all n originals, independent variable included.
            return type(variables)(scale_action(variables, self.inv_herm_mult))
        elif len(variables) == self.scaling_matrix.shape[1] - 1:
            return type(variables)(scale_action(variables, self.dep_var_inv_herm_mult(indep_var_index=indep_var_index)))
        else:
            raise ValueError('Incorrect number of variables for reverse translation')


    def reverse_translate_general(self, variables, system_indep_var_index=0):
        '''
        Given an iterable of variables, or exprs, reverse translate into the original variables.
        Here we expect t as the first variable, since we need to divide by it and substitute


        Given the solutions of a (general) reduced system, reverse translate them into solutions of the original system.

        Args:
            variables (iter of sympy.Expression): The independent variable :math:`t`,
                the solution auxiliary variables :math:`x(t)` and
                the solution invariants :math:`y(t)` of the reduced system.
                Variables should be ordered: independent variable, auxiliary variables then invariants.
            indep_var_index (int): The location of the independent variable.

        :rtype: tuple

        Example 6.6 from :cite:`Hubert2013c`.
        We will just take the matrices from the paper, rather than use our own
        (which are constructed using different conventions).

        .. math::
            \\frac{dz_1}{dt} = \\frac{z_1 \\left(z_1^5 z_2 - 2 \\right)}{3t}

            \\frac{dz_2}{dt} = \\frac{z_2 \\left(10 - 2 z_1^5 z_2 + \\frac{3 z_1^2 z_2}{t} \\right)}{3t}

        >>> equations = ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
        ...              'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)']
        >>> system = ODESystem.from_equations(equations)
        >>> var_order = ['t', 'z1', 'z2']
        >>> system.reorder_variables(var_order)
        >>> scaling_matrix = sympy.Matrix([[3, -1, 5]])
        >>> hermite_multiplier = sympy.Matrix([[1, 1, -1],
        ...                                    [2, 3, 2],
        ...                                    [0, 0, 1]])
        >>> translator = ODETranslator(scaling_matrix=scaling_matrix,
        ...                              hermite_multiplier=hermite_multiplier)
        >>> reduced = translator.translate_general(system)
        >>> # Quickly check the reduction
        >>> answer = ODESystem.from_equations(['dx0/dt = x0*(2*y0*y1/3 - 1/3)/t',
        ...                                    'dy0/dt = y0*(y0*y1 - 1)/t',
        ...                                    'dy1/dt = y1*(y1 + 1)/t'])
        >>> reduced == answer
        True
        >>> reduced
        dt/dt = 1
        dx0/dt = x0*(2*y0*y1/3 - 1/3)/t
        dy0/dt = y0*(y0*y1 - 1)/t
        dy1/dt = y1*(y1 + 1)/t

        Check reverse translation:

        >>> reduced_soln = (sympy.var('t'),
        ...                 sympy.sympify('c3/(t**(1/3)*(ln(t-c1)-ln(t)+c2)**(2/3))'),  # x
        ...                 sympy.sympify('c1/(t*(ln(t-c1)-ln(t)+c2))'),  # y1
        ...                 sympy.sympify('t/(c1 - t)'))  # y2
        >>> original_soln = translator.reverse_translate_general(reduced_soln, system.independent_variable_index)
        >>> original_soln_answer = [#reduced_soln[1] ** 3 / (reduced_soln[2] ** 2)  # x^3 / y1^2
        ...                         reduced_soln[2] / reduced_soln[1],  # y1 / x
        ...                         reduced_soln[1] ** 5 * reduced_soln[3] / reduced_soln[2] ** 4]  # x^5 y2 / y1^4
        >>> original_soln_answer = tuple([soln.subs({sympy.var('t'): sympy.sympify('t / (c3**3 / c1**2)')})
        ...                               for soln in original_soln_answer])
        >>> original_soln == original_soln_answer
        True
        >>> original_soln[0]
        c1/(c3*(c1**2*t/c3**3)**(2/3)*(c2 - log(c1**2*t/c3**3) + log(c1**2*t/c3**3 - c1))**(1/3))
        >>> original_soln[1]
        c3**5*(c1**2*t/c3**3)**(10/3)*(c2 - log(c1**2*t/c3**3) + log(c1**2*t/c3**3 - c1))**(2/3)/(c1**4*(-c1**2*t/c3**3 + c1))
        '''
        if len(variables) != self.scaling_matrix.shape[1] + 1:
            raise ValueError('Incorrect number of variables for reverse translation')

        indep_var, variables = variables[0], variables[1:]  # Indep var is always first in reverse translation only
        raw_solutions = scale_action(variables, self.inv_herm_mult)
        indep_const = raw_solutions[system_indep_var_index] / indep_var  # Shift up by 1 as our first variable is the independent variable
        raw_solutions.col_del(system_indep_var_index)

        # Check indep_const is independent of t
        indep_const_atoms = indep_const.atoms(sympy.Symbol)
        if indep_var in indep_const_atoms:
            raise ValueError('{} is not independent of the independent variable'.format(indep_const))
        # And if we can, the original variables
        if (self.symbols_domain is not None) and len(set(self.symbols_domain).intersection(indep_const_atoms)):
            raise ValueError('{} is not independent of {}'.format(indep_const,
                                                                  set(self.symbols_domain).intersection(indep_const_atoms)))


        to_sub = {indep_var: indep_var / indep_const}
        solns = [soln.subs(to_sub) for soln in raw_solutions]
        return tuple(solns)

    def _validate_variables(self, variables, num_var, stem, use_domain_var):
        '''
        Make sure we have the right number of variables if given, else make them up.
        If use_domain_var is True, then use the domain variables before making up new ones.
        '''
        if variables is None:
            if use_domain_var and (self.symbols_domain is not None):
                variables = self.symbols_domain
            else:
                variables = sympy.var(', '.join('{}{}'.format(stem, i) for i in range(num_var)))

        if len(variables) != num_var:
            raise ValueError('Expecting {} variables not {}'.format(num_var, variables))
        return variables

    def invariant_monomials(self, variables=None):
        '''
        The scaling invariants, as monomials in the original symbols, read off
        :attr:`~desr.ode_translation.ODETranslator.herm_mult_n`.

        **Monomials, and the name says so**, because that is what a scaling reduction
        produces and the distinction matters when reading a reduced system.  Each column of
        :math:`V_n` is an exponent vector, so each invariant is a Laurent monomial --
        :math:`d/K`, :math:`k t / h` -- never a sum.  That is not a limitation of this
        implementation but what makes the reduction work: Hubert & Labahn's result is that
        these monomials *generate*, so every rational invariant of the scaling action is a
        rational function of them.  A quantity like :math:`E + C` may well be invariant, but
        it is expressed **through** these generators rather than appearing among them.

        The companion is :meth:`invariant_names`, the names the reduction gives these
        monomials -- :math:`\\kappa_0`, :math:`\\tau`, :math:`\\nu_0`.  Monomial and name are
        the two faces of one invariant, and :meth:`named_symbols` pairs them up.

        Args:
            variables (iter of sympy.Symbol, optional): The symbols to write the monomials
                in.  Defaults to :attr:`~desr.ode_translation.ODETranslator.symbols_domain`,
                and failing that to placeholder :math:`y_i`.

        Returns:
            sympy.Matrix: A :math:`1 \\times (n - r)` row of monomials, in the column order
                of :math:`V_n` -- so it reads directly against the scaling matrix.

        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )','dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.invariant_monomials()
        Matrix([[d/K, h*r/k, h*s/k, k*t/h, n/K, h*p/K]])
        '''
        variables = self._validate_variables(variables, self.n, 'y', True)
        return scale_action(variables, self.herm_mult_n)

    def invariant_names(self, system=None, scheme=None):
        '''
        The names the reduction gives its invariants, in the same order as
        :meth:`invariant_monomials`.

        Args:
            system (ODESystem): The system being reduced.
            scheme (desr.ode_translation.ReductionScheme, optional): As for
                :meth:`named_symbols`.

        Returns:
            tuple of sympy.Symbol: Names -- :math:`\\kappa_0`, :math:`\\tau`,
                :math:`\\nu_0` -- as sympy symbols, not strings.

        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )','dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.invariant_names(system)
        (kappa0, kappa1, kappa2, tau, nu0, nu1)
        >>> translator.invariant_monomials()
        Matrix([[d/K, h*r/k, h*s/k, k*t/h, n/K, h*p/K]])
        '''
        if scheme is None:
            scheme = self.scheme(system)
        else:
            scheme = ReductionScheme(scheme)

        count = self.herm_mult_n.shape[1]

        if scheme == ReductionScheme.PARAMETER:
            # The columns of V_n run (parameter invariants, the clock, the dependents), and
            # the reduced system's symbols are built in that order.
            variable_map = self.variable_map(system)
            names = (list(self._new_constants(system))
                     + [variable_map[system.independent_variable]]
                     + [variable_map[v] for v in system.dependents])
        else:
            naming = self._resolved_naming_scheme(GENERAL_SCHEME_NAMES)
            names = _names_for(naming['dependents'], count,
                               self._new_indices_start_at, 'invariants')

        if len(names) != count:
            raise ValueError('the {} scheme names {} invariants but there are {}'.format(
                scheme, len(names), count))
        return tuple(names)

    def auxiliary_monomials(self, variables=None):
        '''
        The auxiliaries, as monomials in the original symbols, read off
        :attr:`~desr.ode_translation.ODETranslator.herm_mult_i`.

        Monomials for the same reason as :meth:`invariant_monomials`: each column of
        :math:`V_i` is an exponent vector.  These are the :math:`r` normalised directions --
        the reduction is free to set each to any value, and the parameter scheme sets them
        all to 1 and discards them.

        The companion is :meth:`auxiliary_names`.

        Args:
            variables (iter of sympy.Symbol, optional): The symbols to write the monomials
                in.  Defaults to :attr:`~desr.ode_translation.ODETranslator.symbols_domain`,
                and failing that to placeholder :math:`x_i`.

        Returns:
            sympy.Matrix: A :math:`1 \\times r` row of monomials, in the column order of
                :math:`V_i`.

        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )','dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.auxiliary_monomials()
        Matrix([[h/k, K, K/h]])
        '''
        variables = self._validate_variables(variables, self.n, 'x', True)
        return scale_action(variables, self.herm_mult_i)

    def auxiliary_names(self, system=None, scheme=None):
        '''
        The names the reduction gives its auxiliaries, in the same order as
        :meth:`auxiliary_monomials`.

        Empty under the parameter scheme, which normalises its auxiliaries to 1 and keeps
        none of them, so there is nothing to name.

        Args:
            system (ODESystem): The system being reduced.
            scheme (desr.ode_translation.ReductionScheme, optional): As for
                :meth:`named_symbols`.

        Returns:
            tuple of sympy.Symbol

        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )','dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.auxiliary_names(system, scheme='general')
        (x0, x1, x2)
        >>> translator.auxiliary_names(system, scheme='parameter')
        ()
        '''
        if scheme is None:
            scheme = self.scheme(system)
        else:
            scheme = ReductionScheme(scheme)

        if scheme == ReductionScheme.PARAMETER:
            return ()

        naming = self._resolved_naming_scheme(GENERAL_SCHEME_NAMES)
        return tuple(_new_symbols(naming['auxiliaries'], self.herm_mult_i.shape[1],
                                  self._new_indices_start_at))

    def _dynamic_rows(self, system):
        '''
        The row indices of the scaling matrix that belong to variables which move: the
        independent variable and the dependent variables.

        The rest are parameter rows.  Every classification below is a question about which
        of these rows an invariant's exponents touch.
        '''
        dynamic = set(system.dependents) | {system.independent_variable}
        return [i for i, variable in enumerate(system.symbols) if variable in dynamic]

    def _kind_by_column(self, system):
        '''
        The kind of each invariant, by column index of :attr:`herm_mult_n`.

        A three-way partition on which rows the column's exponents touch:

        - ``'parameter'``   -- no dynamic row at all, so a monomial in the original
          parameters alone.
        - ``'independent'`` -- the independent variable's row but no dependent's, so it is
          a rescaled clock and nothing else.
        - ``'dependent'``   -- at least one dependent variable's row, whether or not it
          also involves the independent variable.

        The kind is structural and so does not depend on which reduction scheme runs; only
        the *name* an invariant ends up with does.  Kept separate for that reason.
        '''
        indep_row = system.independent_variable_index
        dependent_rows = [i for i, variable in enumerate(system.symbols)
                          if variable in set(system.dependents)]

        kinds = []
        for column in range(self.herm_mult_n.shape[1]):
            if any(self.herm_mult_n[row, column] != 0 for row in dependent_rows):
                kinds.append('dependent')
            elif self.herm_mult_n[indep_row, column] != 0:
                kinds.append('independent')
            else:
                kinds.append('parameter')
        return kinds

    def symbol_kinds(self, system, scheme=None, include_aux_vars=True):
        '''
        What kind of symbol each variable of the reduced system is.

        A dict over **every** variable the reduction produces, mapping it to one of four
        kinds.  The keys are the same as :meth:`named_symbols`', so the two join: one
        says what a symbol *means*, the other says what it *is*.

        - ``'independent'`` -- the variable to integrate against.  There is always exactly
          one.  Under the parameter scheme it is the rescaled clock, :math:`\\tau`; under
          the general and dependent-variable schemes the reduction keeps the original
          independent variable, so it is :math:`t` itself.
        - ``'dependent'``   -- a variable of the reduced system with an equation.
        - ``'parameter'``   -- an invariant that is a monomial in the original parameters
          alone, so it does not change.
        - ``'auxiliary'``   -- a normalised direction, from
          :attr:`~desr.ode_translation.ODETranslator.herm_mult_i`.  Only the general and
          dependent-variable schemes produce any, and only when they are kept.

        The distinction the kinds exist for is between the first three: a parameter and a
        conserved dependent both have rate zero, and only the structure tells them apart.

        A kind is the variable's **role in the reduced system**, which is not quite the same
        as the structure of its monomial.  Under the general scheme an invariant can be a
        rescaled clock -- :math:`k t / h` -- and still be a ``'dependent'`` here, because
        it has an equation and evolves while :math:`t` is what the system is integrated
        against.  :meth:`indep_var_invariant` answers the other question, whether any
        invariant is a rescaled clock and nothing else, and returns :const:`None` when none
        is.

        Args:
            system (ODESystem): The system being reduced.
            scheme (desr.ode_translation.ReductionScheme, optional): The scheme that will be
                used, since the names depend on it.  The kinds of the invariants do not --
                those are structural.
            include_aux_vars (bool): Whether the auxiliaries will be kept.

        Returns:
            dict: Reduced variable to its kind, in the reduced system's variable order.

        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )','dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.symbol_kinds(system)
        {kappa0: 'parameter', kappa1: 'parameter', kappa2: 'parameter', tau: 'independent', nu0: 'dependent', nu1: 'dependent'}

        Under the general scheme the auxiliaries appear, and `t` is kept as the clock:

        >>> translator.symbol_kinds(system, scheme='general')
        {t: 'independent', x0: 'auxiliary', x1: 'auxiliary', x2: 'auxiliary', y0: 'parameter', y1: 'parameter', y2: 'parameter', y3: 'dependent', y4: 'dependent', y5: 'dependent'}

        Which is what the reduced system turns out to have:

        >>> translator.translate_general(system).symbols
        (t, x0, x1, x2, y0, y1, y2, y3, y4, y5)

        Hubert & Labahn 6.6 has no invariant that is a rescaled clock alone -- both involve
        a dependent variable -- but the reduction still integrates against `t`:

        >>> hl66 = ODESystem.from_equations(['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
        ...                                  'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])
        >>> hl = ODETranslator.from_ode_system(hl66)
        >>> hl.symbol_kinds(hl66)
        {t: 'independent', x0: 'auxiliary', y0: 'dependent', y1: 'dependent'}
        >>> hl.indep_var_invariant(hl66) is None
        True
        '''
        names = list(self.named_symbols(system, scheme, include_aux_vars))
        return dict(zip(names, self._kind_by_variable(system, scheme, include_aux_vars)))

    def _kind_by_variable(self, system, scheme=None, include_aux_vars=True):
        '''
        The kind of each variable of the reduced system, in its variable order.
        '''
        if scheme is None:
            scheme = self.scheme(system)
        else:
            scheme = ReductionScheme(scheme)

        invariant_kinds = self._kind_by_column(system)

        if scheme == ReductionScheme.PARAMETER:
            # The clock invariant *becomes* the reduced system's independent variable, so
            # the structural kinds are already the roles.
            return invariant_kinds

        # The general and dependent-variable schemes keep the original independent variable
        # and prepend it, then the auxiliaries, then the invariants.  So an invariant that
        # is structurally a rescaled clock -- `k*t/h`, say -- is nonetheless a *dependent*
        # variable of this reduced system: it has an equation and it evolves, while `t`
        # remains what everything is integrated against.  Hence exactly one 'independent'
        # either way, and it is the reduced system's own independent variable.
        auxiliaries = ['auxiliary'] * (self.herm_mult_i.shape[1] if include_aux_vars else 0)
        roles = ['dependent' if kind == 'independent' else kind
                 for kind in invariant_kinds]
        return ['independent'] + auxiliaries + roles

    def _invariants_by_kind(self, system, variables=None):
        '''
        The invariants grouped by what kind of symbol each is.

        The rule is **structural**, read off the exponents in
        :attr:`~desr.ode_translation.ODETranslator.herm_mult_n`, and never off the rate of
        change:

        - :attr:`ReductionScheme.PARAMETER`'s value, ``'parameter'``, if the invariant's
          exponents are zero on every *dynamic* row -- the independent variable's and the
          dependent variables' -- so it is a monomial in the original parameters alone.
        - ``'dependent'`` otherwise: it has at least one nonzero exponent on a dynamic row.

        **Why not look at the rate.** A reduced parameter always has rate zero, and a
        reduced dependent *may* have rate zero because it is conserved, so the rate cannot
        tell them apart while the exponents always can.  The case is not hypothetical: a
        catalyst in a reaction network gives an invariant like ``k_0_1*w/k_2_3``, which is
        conserved because the catalyst is, and which must stay a dependent variable --
        calling it a parameter would assert that it is a constant of the model rather than
        a quantity conserved along the trajectory.

        Args:
            system (ODESystem): The system being reduced.
            variables (iter of sympy.Symbol, optional): As for :meth:`invariants`.

        Returns:
            dict: ``{'parameter': (...), 'dependent': (...)}``, each in invariant order.
                Both keys are always present, so an empty group is an empty tuple rather
                than a :exc:`KeyError`.

        Private: :meth:`parameter_invariants` and :meth:`dependent_invariants` are the two
        groups, and :meth:`symbol_kinds` is the per-symbol answer, so a third public way
        to ask the same question would be sprawl.
        '''
        invariants = self.invariant_monomials(variables)
        grouped = {'parameter': [], 'independent': [], 'dependent': []}
        for column, kind in enumerate(self._kind_by_column(system)):
            grouped[kind].append(invariants[column])
        return {kind: tuple(found) for kind, found in grouped.items()}

    def parameter_invariants(self, system, variables=None):
        '''
        The invariants that are monomials in the original parameters alone.

        These become the reduced system's parameters.  There are :math:`p - r` of them when
        the parameter scheme applies, and there may be none at all -- Hubert & Labahn 6.6
        has no parameters to begin with.

        Args:
            system (ODESystem): The system being reduced.
            variables (iter of sympy.Symbol, optional): As for :meth:`invariants`.

        Returns:
            tuple of sympy.Expr

        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )','dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.parameter_invariants(system)
        (d/K, h*r/k, h*s/k)
        '''
        return self._invariants_classified_as('parameter', system, variables)

    def dependent_invariants(self, system, variables=None):
        '''
        The invariants that involve the independent variable or a dependent variable.

        These become the reduced system's dependent variables.  There are :math:`q` of
        them, one per dependent variable of the original system -- the clock is reported
        separately by :meth:`indep_var_invariant`, and the constants by
        :meth:`parameter_invariants`.

        Args:
            system (ODESystem): The system being reduced.
            variables (iter of sympy.Symbol, optional): As for :meth:`invariants`.

        Returns:
            tuple of sympy.Expr

        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )','dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.dependent_invariants(system)
        (n/K, h*p/K)

        There are exactly as many as the original system has dependent variables, which is
        the precise sense in which a reduction loses symbols only among the parameters:

        >>> len(translator.dependent_invariants(system)) == system.num_dependents
        True
        '''
        return self._invariants_classified_as('dependent', system, variables)

    def indep_var_invariant(self, system, variables=None):
        '''
        The invariant that carries the clock, or :const:`None` if there is not one.

        It is the invariant whose exponents are nonzero on the independent variable's row
        and zero on every dependent variable's row.  Under the parameter scheme there is
        exactly one.  Under the general scheme there may be none, because the clock can be
        entangled with a dependent variable: in Hubert & Labahn 6.6 the invariant is
        :math:`t z_1^3`, and there is no invariant that is a function of :math:`t` and the
        parameters alone.

        :const:`None` states that rather than papering over it.  A caller that needs a new
        independent variable and gets :const:`None` is being told something true about the
        reduction, which is better than being handed an invariant that also moves.

        Args:
            system (ODESystem): The system being reduced.
            variables (iter of sympy.Symbol, optional): As for :meth:`invariants`.

        Returns:
            sympy.Expr or None

        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )','dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.indep_var_invariant(system)
        k*t/h

        Hubert & Labahn 6.6 has no separable clock:

        >>> hl66 = ODESystem.from_equations(['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
        ...                                  'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])
        >>> translator = ODETranslator.from_ode_system(hl66)
        >>> translator.invariant_monomials()
        Matrix([[t*z1**3, z1**2*z2/t]])
        >>> translator.indep_var_invariant(hl66) is None
        True
        '''
        invariants = self.invariant_monomials(variables)
        indep_row = system.independent_variable_index
        dependent_rows = [i for i, variable in enumerate(system.symbols)
                          if variable in set(system.dependents)]

        for column in range(self.herm_mult_n.shape[1]):
            if self.herm_mult_n[indep_row, column] == 0:
                continue
            if any(self.herm_mult_n[row, column] != 0 for row in dependent_rows):
                continue
            return invariants[column]
        return None

    def _invariants_classified_as(self, kind, system, variables=None):
        '''
        The `kind` group of :meth:`_invariants_by_kind`.
        '''
        return self._invariants_by_kind(system, variables)[kind]

    def named_symbols(self, system, scheme=None, include_aux_vars=True):
        '''
        What each variable of the reduced system stands for.

        A dict from the reduced system's variable to the monomial in the original variables
        that it denotes -- ``{kappa0: d/K, ..., tau: k*t/h, nu0: n/K, ...}``.  The keys are
        plain symbols, so this is the map to reach for when you have a reduced system in
        hand and want to know what one of its symbols means: an invariant is both a new
        name and the monomial that name stands for, and this is the two of them together.

        :meth:`symbol_kinds` is keyed the same way and says what each symbol *is*.

        Args:
            system (ODESystem): The system being reduced.
            scheme (desr.ode_translation.ReductionScheme, optional): The scheme that will be
                used, since the names depend on it -- the parameter scheme gives `kappa`,
                `tau` and `nu` while the general scheme gives `x` and `y`.  Inferred by :meth:`scheme`
                when left out, with the caveat in :meth:`expected_variable_count`.
            include_aux_vars (bool): Whether the auxiliaries will be kept.

        Returns:
            dict: Reduced variable to the monomial it denotes, in the reduced system's
                variable order.  Under the general and dependent-variable schemes the
                original independent variable is kept, and maps to itself.

        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )','dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.named_symbols(system)
        {kappa0: d/K, kappa1: h*r/k, kappa2: h*s/k, tau: k*t/h, nu0: n/K, nu1: h*p/K}

        So reading a reduced system needs no reconstruction of monomials:

        >>> reduced = translator.translate(system)
        >>> reduced.parameters
        (kappa0, kappa1, kappa2)
        >>> translator.named_symbols(system)[reduced.parameters[0]]
        d/K

        Under the general scheme the auxiliaries are included, and `t` maps to itself:

        >>> translator.named_symbols(system, scheme='general')
        {t: t, x0: h/k, x1: K, x2: K/h, y0: d/K, y1: h*r/k, y2: h*s/k, y3: k*t/h, y4: n/K, y5: h*p/K}
        '''
        if scheme is None:
            scheme = self.scheme(system)
        else:
            scheme = ReductionScheme(scheme)

        # Built from the per-block accessors rather than repeating their naming logic, so
        # that the pairing here cannot disagree with what invariant_names and
        # auxiliary_names say on their own.
        invariant_names = list(self.invariant_names(system, scheme))
        invariant_monomials = list(self.invariant_monomials())

        if scheme == ReductionScheme.PARAMETER:
            names, monomials = invariant_names, invariant_monomials
        else:
            # These two schemes keep the original independent variable, which maps to
            # itself, then the auxiliaries, then the invariants.
            if include_aux_vars:
                auxiliary_names = list(self.auxiliary_names(system, scheme))
                auxiliary_monomials = list(self.auxiliary_monomials())
            else:
                auxiliary_names, auxiliary_monomials = [], []
            names = ([system.independent_variable] + auxiliary_names + invariant_names)
            monomials = ([system.independent_variable] + auxiliary_monomials
                         + invariant_monomials)

        if len(names) != len(monomials):
            raise ValueError(
                'the {} scheme names {} symbols but there are {}'.format(
                    scheme, len(names), len(monomials)))

        return dict(zip(names, monomials))

    def rewrite_rules(self, variables=None):
        '''
        Given a set of variables, print the rewrite rules.
        These are the rules that allow you to write any rational invariant in terms of the
        :attr:`~desr.ode_translation.ODETranslator.invariants`.

        Args:
            variables (iter of sympy.Symbol):
                The names of the variables appearing in the rational invariant.

        Returns:
            dict:
                A :code:`dict` whose keys are the given variables and whose values are the values to be substituted.

        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )','dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.rewrite_rules()
        {K: 1, d: d/K, h: 1, k: 1, r: h*r/k, s: h*s/k, t: k*t/h, n: n/K, p: h*p/K}

        For example, :math:`r t` is an invariant.
        Substituting in the above mapping gives us a way to write it in terms of our generating set of invariants:

        .. math::
            r t \\mapsto \\left( \\frac{r}{s} \\right) \\left( s t \\right) = r t
        '''
        variables = self._validate_variables(variables, self.n, 'z', True)

        rules = scale_action(variables, self.herm_mult_n * self.inv_herm_mult_d).T
        assert len(rules) == len(variables)

        # Create a dict of rules
        rules = dict(zip(variables, rules))

        return rules

    def moving_frame(self, variables=None):
        '''
        Given a finite-dimensional Lie group :math:`G` acting on a manifold :math:`M`,
        a moving frame is defined as a :math:`G`-equivariant map :math:`\\rho : M \\rightarrow G`.

        We can construct a moving frame given a rational section,
        by sending a point :math:`x \\in M` to the Lie group element
        :math:`\\rho ( x )` such that :math:`\\rho ( x ) \\cdot x` lies on :math:`K`.


        >>> equations = 'dz1/dt = z1*(1+z1*z2);dz2/dt = z2*(1/t - z1*z2)'.split(';')
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> moving_frame = translator.moving_frame()
        >>> moving_frame
        (z1,)

        Now we can check that the moving frame moves a point to the cross-section.

        >>> translator.rational_section()
        Matrix([[1 - 1/z1]])
        >>> scale_action(moving_frame, translator.scaling_matrix)  # \\rho(x)
        Matrix([[1, 1/z1, z1]])
        >>> # :math:`\\rho(x).x` should satisfy the equations of the rational section.
        >>> scale_action(moving_frame, translator.scaling_matrix).multiply_elementwise(sympy.Matrix(system.symbols).T)
        Matrix([[t, 1, z1*z2]])


        Another example:

        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )','dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> moving_frame = translator.moving_frame()
        >>> moving_frame
        (k/h, 1/K, h/K)

        Now we can check that the moving frame moves a point to the cross-section.

        >>> translator.rational_section()
        Matrix([[h - 1/k, K - 1, K - 1/h]])
        >>> scale_action(moving_frame, translator.scaling_matrix)  # \\rho(x)
        Matrix([[1/K, 1/K, 1/h, 1/k, h/k, h/k, k/h, 1/K, h/K]])
        >>> # :math:`\\rho(x).x` should satisfy the equations of the rational section.
        >>> scale_action(moving_frame, translator.scaling_matrix).multiply_elementwise(sympy.Matrix(system.symbols).T)
        Matrix([[1, d/K, 1, 1, h*r/k, h*s/k, k*t/h, n/K, h*p/K]])

        See :cite:`Fels1999` for more details.
        '''
        variables = self._validate_variables(variables, self.n, 'z', True)

        rules = scale_action(variables, - self.herm_mult_i)
        assert len(rules) == self.r
        return tuple(rules)

    def rational_section(self, variables=None):
        '''
        Given a finite-dimensional Lie group :math:`G` acting on a manifold :math:`M`,
        a cross-section :math:`K \\subset M` is a subset
        that intersects each orbit of :math:`M` in exactly one place.

        Args:
            variables (iter of sympy.Symbol, optional):
                Variables of the system.

        Returns:
            sympy.Matrix:
                A set of equations defining a variety that is a cross-section.

        >>> equations = 'dz1/dt = z1*(1+z1*z2);dz2/dt = z2*(1/t - z1*z2)'.split(';')
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.rational_section()
        Matrix([[1 - 1/z1]])


        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )','dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> translator = ODETranslator.from_ode_system(system)
        >>> translator.rational_section()
        Matrix([[h - 1/k, K - 1, K - 1/h]])


        See :cite:`Fels1999` for more details.
        '''
        variables = self._validate_variables(variables, self.n, 'z', True)

        # Positive values of Vi
        vip = self.herm_mult_i
        vip = vip.applyfunc(lambda x: x if x > 0 else 0)
        # vip[vip < 0] = 0

        # Negative values of Vi
        vin = self.herm_mult_i
        vin = vin.applyfunc(lambda x: x if x < 0 else 0)

        return scale_action(variables, vip) - scale_action(variables, vin)

    ## Choosing invariants
    def extend_from_invariants(self, invariant_choice):
        '''
        Extend a given set of invariants, expressed as a matrix of exponents, to find a Hermite multiplier that will
        rewrite the system in terms of invariant_choice.


        Parameters
        ----------
        invariant_choice : sympy.Matrix
            The :math:`n \\times k` matrix representing the invariants, where :math:`n` is the number of variables and :math:`k` is the
            number of given invariants.


        Returns
        -------
        ODETranslator
            An ODETranslator representing the rewrite rules in terms of the given invariants.


        >>> variables = sympy.symbols(' '.join(['y{}'.format(i) for i in range(6)]))
        >>> ode_translator = ODETranslator(sympy.Matrix([[1, 0, 3, 0, 2, 2],
        ...                                                [0, 2, 0, 1, 0, 1],
        ...                                                [2, 0, 0, 3, 0, 0]]))
        >>> ode_translator.invariant_monomials(variables=variables)
        Matrix([[y2*y3**2/(y0**3*y1), y4**3/y2**2, y4*y5**2/(y1*y2**2)]])

        Now we can express two new invariants, which we think are more interesting, as a matrix.
        We pick the product of the first two invariants, and the product of the last two invariants:
        :code:`y0**3 * y1 * y2/(y3**2 * y4**3)` and :code:`y1 *  y2**2/(y4 * y5**2)`

        >>> new_inv = sympy.Matrix([[3, 1, 1, -2, -3, 0],
        ...                         [0, 1, 2, 0, -1, -2]]).T

        >>> new_translator = ode_translator.extend_from_invariants(new_inv)
        >>> new_translator.invariant_monomials(variables=variables)
        Matrix([[y0**3*y1*y2/(y3**2*y4**3), y1*y2**2/(y4*y5**2), y4**3/y2**2]])
        '''
        ## Step 1: Check we have invariants
        choice_actions = self.scaling_matrix * invariant_choice
        if not choice_actions.is_zero_matrix:
            raise ValueError('The desired combinations {} are not invariants of the scaling action.'.format(invariant_choice))

        ## Step 2: Try and extend the choices by a basis of invariants
        ## Step 2a: Extend (W_d . invariant_choice) to a unimodular matrix
        column_operations = extend_rectangular_matrix(self.inv_herm_mult_d * invariant_choice)

        # Step 2b: Apply these column operations to the current Hermite multiplier
        hermite_multiplier = self.herm_mult.copy()
        hermite_multiplier[:, self.r:] = self.herm_mult_n * column_operations
        max_scal = ODETranslator(scaling_matrix=self.scaling_matrix, symbols_domain=self.symbols_domain,
                                  hermite_multiplier=hermite_multiplier)

        return max_scal

    ## Determine whether an expression is an invariant or not
    def is_invariant_expr(self, expr, variable_order=None):
        '''
        Determine whether an expression is an invariant or not

        Parameters
        ----------
        expr : sympy.Expr
            A symbolic expression which may or may not be an invariant under the scaling action.
        variable_order : list, tuple
            An ordered list that determines how the matrix acts on the variables.

        Returns
        -------
        bool
            True if the expression is an invariant.
        '''
        if variable_order is not None:
            variable_order = tuple(variable_order)
            if len(variable_order) != self.n:
                raise ValueError('{} variables given but we have {} actions'.format(len(variable_order), self.n))
        else:
            if self.symbols_domain is None:
                raise ValueError('No variable order found to determine invariance.')
            variable_order = self.symbols_domain

        missing_var = set(variable_order).difference(expr.atoms(sympy.Symbol))
        if len(missing_var):
            raise ValueError('Unknown action on variables: {}'.format())

        raise NotImplemented()



def scale_action(vect, scaling_matrix):
    '''
    Given a vector of sympy expressions, determine the action defined by scaling_matrix
    I.e. Given vect, calculate vect^scaling_matrix (in notation of Hubert Labahn).

    Args:
        vect (iter of sympy.Expression):
            Element of the algebraic torus that is going to act on the system.
        scaling_matrix (sympy.Matrix):
            Matrix that defines the action on the system.

    Returns:
        sympy.Matrix:
            Matrix that, when multiplied element-wise with an element of the manifold, gives the result of the action.


    Example 3.3

    >>> x = sympy.var('x')
    >>> A = sympy.Matrix([[2, 3]])
    >>> scale_action([x], A)
    Matrix([[x**2, x**3]])

    Example 3.4

    >>> v = sympy.var('m v')
    >>> A = sympy.Matrix([[6, 0, -4, 1, 3], [0, 3, 1, -4, 3]])
    >>> scale_action(v, A)
    Matrix([[m**6, v**3, v/m**4, m/v**4, m**3*v**3]])
    '''
    assert len(vect) == scaling_matrix.shape[0]
    out = []
    for col in scaling_matrix.T.tolist():
        mon = 1
        for var, power in zip(vect, col):
            mon *= var ** power
        out.append(mon)
    return sympy.Matrix(out).T

def extend_rectangular_matrix(matrix_, check_unimodular=True):
    """
    Given a rectangular :math:`n \\times m` integer matrix, extend it to a unimodular one by appending columns.

    Parameters
    ----------
    matrix_ : sympy.Matrix
        The rectangular matrix to be extended.


    Returns
    -------
    sympy.Matrix
        Square matrix of determinant 1.


    >>> matrix_ = sympy.Matrix([[1, 0],
    ...                         [0, 1],
    ...                         [0, 0],])
    >>> extend_rectangular_matrix(matrix_)
    Matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]])

    By default, we will throw if we cannot extend a matrix to a unimodular matrix.

    >>> matrix_[0, 0] = 2
    >>> extend_rectangular_matrix(matrix_=matrix_)
    Traceback (most recent call last):
        ...
    ValueError: Unable to extend the matrix
    Matrix([[2, 0], [0, 1], [0, 0]])
    to a unimodular matrix.

    We can extend some more interesting matrices.

    >>> matrix_ = sympy.Matrix([[3, 2],
    ...                         [-2, 1],
    ...                         [5, 6],])
    >>> extend_rectangular_matrix(matrix_)
    Matrix([
    [ 3, 2, -3],
    [-2, 1,  1],
    [ 5, 6, -6]])
    """
    matrix_ = matrix_.copy()
    if len(matrix_.shape) != 2:
        raise ValueError('Can only extend arrays of dimension 2, not {}'.format(len(matrix_.shape)))

    n, m = matrix_.shape

    # Assume we have more rows than columns
    if n < m:
        return extend_rectangular_matrix(matrix_=matrix_.T, check_unimodular=check_unimodular).T
    elif n == m:
        return matrix_

    # Compute the unique normal row Hermite multiplier, so that multiplier . matrix_ = [I; 0].
    # matrix_ extends to a unimodular matrix exactly when this row HNF is [I; 0]; then
    # [matrix_ | extension] = multiplier^-1, so the canonical extension is the last n - m
    # columns of multiplier^-1.  Using the *normal* (unique) Hermite multiplier makes the
    # extension unique and independent of the underlying HNF backend, unlike the non-unique
    # transformation matrices returned by the Smith normal form.
    hermite_normal_form, multiplier = normal_hnf_row(matrix_)

    if check_unimodular:
        # Extendable to a determinant +-1 matrix exactly when the row HNF is the identity on top of zeros.
        if (hermite_normal_form[:m, :] != sympy.eye(m)) or (not hermite_normal_form[m:, :].is_zero_matrix):
            raise ValueError('Unable to extend the matrix\n{}\nto a unimodular matrix.'.format(matrix_))

    extension = _int_inv(multiplier)[:, m:]

    extended = sympy.Matrix.hstack(matrix_, extension)

    assert extended.shape[0] == extended.shape[1] == n
    if check_unimodular:
        if abs(extended.det()) != 1:
            raise RuntimeError('Extended matrix\n{}\nhas determinant {}, not +-1'.format(extended, extended.det()))

    return extended

__all__ = ['ODETranslator', 'ReductionScheme', 'scale_action', 'extend_rectangular_matrix',
           'NAMING_SCHEME_KEYS', 'PARAMETER_SCHEME_NAMES', 'GENERAL_SCHEME_NAMES']

if __name__ == '__main__':
    import doctest
    doctest.testmod()