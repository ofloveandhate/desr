"""
Numerical translation between an :class:`~desr.ode_system.ODESystem` and its reduction.

:mod:`desr` reduces a system symbolically.  This module carries *numbers* across the same
reduction: it turns parameter values and initial conditions of the original system into
their counterparts in the reduced system, and turns a numerical solution of the reduced
system back into a solution of the original one.

This module needs :mod:`numpy`.  It is deliberately not imported by ``desr/__init__.py``,
so :mod:`desr` itself keeps no numerical dependencies.

All three reduction schemes are supported, and they differ in what you must supply to come
back.  The reduction hides :math:`r` scaling symmetries in :math:`r` *auxiliary* variables.

- Under the general and dependent-variable schemes the auxiliaries are variables of the
  reduced system, so a solution of it already carries them and nothing else is needed.
- Under the parameter scheme they are not, and neither they are when the general scheme is
  used with ``include_aux_vars=False``.  Then :meth:`NumericTranslation.reverse` needs
  ``known_values``: :math:`r` values of the original system that you already know.  Without
  them the original system is not determined, since an entire :math:`r`-parameter family of
  them share the same reduction.
"""

import sympy

try:
    import numpy as np
except ImportError:  # pragma: no cover
    raise ImportError('desr.numerics requires numpy.  Install it with `pip install numpy`.')

from .ode_translation import ODETranslation
from .ode_system import ODESystem

__all__ = ['NumericTranslation', 'InsufficientKnownValues', 'ConflictingKnownValues',
           'CONSISTENCY_TOLERANCE', 'PARAMETER_SCHEME', 'DEP_VAR_SCHEME', 'GENERAL_SCHEME']

PARAMETER_SCHEME = 'parameter'
DEP_VAR_SCHEME = 'dep_var'
GENERAL_SCHEME = 'general'


class InsufficientKnownValues(ValueError):
    '''
    Raised when the values supplied to :meth:`NumericTranslation.reverse` do not pin down a
    unique original system.  The message says which values would.
    '''
    pass


class ConflictingKnownValues(ValueError):
    '''
    Raised when the values supplied to :meth:`NumericTranslation.reverse` cannot all hold at
    once given the reduced values: they determine the same scale, and disagree about it.
    '''
    pass


# Relative tolerance within which known values must agree with each other and with the
# reduced values.  Known values are measurements and reduced values are numerical solutions,
# so this is looser than machine precision but far tighter than any real disagreement.
CONSISTENCY_TOLERANCE = 1e-6


def _first(value):
    '''The first sample of an array, or the value itself if it is a scalar.'''
    value = np.asarray(value)
    return value.reshape(-1)[0] if value.ndim else value[()]


def _monomial(values, exponents, column):
    '''
    Evaluate one monomial: the product of ``values`` raised to the powers in one column.

    Every exponent produced by the reduction is an integer, so negative values are raised to
    integer powers and come through exactly.  There is no need to work in logarithms, and no
    restriction to positive data.
    '''
    result = 1.0
    for i in range(exponents.rows):
        power = exponents[i, column]
        if power != 0:
            result = result * values[i] ** int(power)
    return result


def scheme_of(system, translation):
    '''
    Report which reduction scheme :meth:`~desr.ode_translation.ODETranslation.translate`
    applies to this system, following the same order of preference.

    Args:
        system (ODESystem): The original system.
        translation (ODETranslation): Its reduction.

    Returns:
        str: One of :data:`PARAMETER_SCHEME`, :data:`DEP_VAR_SCHEME`, :data:`GENERAL_SCHEME`.
    '''
    if translation._is_translate_parameter_compatible(system):
        return PARAMETER_SCHEME
    if ((len(system.variables) == translation.scaling_matrix.shape[1] + 1)
            or translation.scaling_matrix[:, system.indep_var_index].is_zero_matrix):
        return DEP_VAR_SCHEME
    if len(system.variables) == translation.scaling_matrix.shape[1]:
        return GENERAL_SCHEME
    raise ValueError("System doesn't have the right number of variables for translation")


def _dep_var_translation(system, translation):
    '''
    Rebuild the reduction that :meth:`~desr.ode_translation.ODETranslation.translate_dep_var`
    works with internally: the same scaling action, with the independent variable removed.
    '''
    index = system.indep_var_index
    scaling_matrix = translation.scaling_matrix.copy()
    hermite_multiplier = translation.dep_var_herm_mult(indep_var_index=index)
    if translation.n == len(system.variables):
        scaling_matrix.col_del(index)
    else:
        hermite_multiplier.col_del(index)
        hermite_multiplier.row_del(index)
    return ODETranslation(scaling_matrix=scaling_matrix,
                          hermite_multiplier=hermite_multiplier)


class NumericTranslation(object):
    '''
    Carry numbers back and forth across a reduction.

    Args:
        system (ODESystem): The original system.
        translation (ODETranslation): The reduction of that system.
        reduced (ODESystem, optional): The reduced system.  Computed with
            :meth:`~desr.ode_translation.ODETranslation.translate` if not given.  Pass it
            explicitly when it was produced with non-default options, such as
            ``include_aux_vars=False``.

    Dictionaries handed in and out are keyed by the actual symbols of each system.  Values
    may be scalars or :class:`numpy.ndarray`, so a whole time series translates in one call.

    Attributes:
        scheme (str): Which reduction scheme applies.
        reduced (ODESystem): The reduced system.

    Take the prey-predator system, which reduces by the parameter scheme.

    >>> import numpy
    >>> from desr.numerics import NumericTranslation
    >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)']
    >>> system = ODESystem.from_equations(equations)
    >>> translation = ODETranslation.from_ode_system(system)
    >>> numeric = NumericTranslation(system, translation)
    >>> numeric.scheme
    'parameter'

    Three of the eight variables are scaling symmetries, so three values must be known to
    come back.

    >>> numeric.r
    3
    >>> numeric.carries_auxiliaries
    False

    Translate a state of the original system into the reduced one.

    >>> t, n, p, K, d, h, k, r, s = system.variables
    >>> original = {t: 2.0, n: 1.0, p: 0.5, K: 5.0, d: 0.7, h: 1.3, k: 0.9, r: 1.1, s: 0.4}
    >>> reduced_values = numeric.forward(original)
    >>> sorted(str(x) for x in reduced_values)
    ['kappa0', 'kappa1', 'kappa2', 'nu0', 'nu1', 'tau']

    And back again, given three values we already know.

    >>> recovered = numeric.reverse(reduced_values, known_values={d: 0.7, s: 0.4, k: 0.9})
    >>> [round(float(recovered[x]), 10) for x in (t, n, p)]
    [2.0, 1.0, 0.5]
    >>> [round(float(recovered[x]), 10) for x in (K, h, r)]
    [5.0, 1.3, 1.1]
    '''

    def __init__(self, system, translation, reduced=None):
        self.system = system
        self.translation = translation
        self.reduced = translation.translate(system) if reduced is None else reduced

        # `translate` prefers the parameter scheme, but a caller may have asked for another
        # one explicitly, so choose the scheme whose bookkeeping actually fits the reduced
        # system we were handed rather than the one `translate` would have picked.
        layout = None
        for candidate in self._applicable_schemes():
            layout = self._layout_for(candidate)
            if layout is not None:
                self.scheme = candidate
                break
        if layout is None:
            raise ValueError(
                'Reduced system {} does not match any reduction of {} by this translation.  '
                'Check that `reduced` came from this system and translation.'.format(
                    tuple(self.reduced.variables), tuple(system.variables)))

        (self._acted, source, self._auxiliaries, self._invariants,
         self._shares_indep_var) = layout
        self._herm_mult_i = source.herm_mult_i
        self._herm_mult_n = source.herm_mult_n
        self._inv_herm_mult = source.inv_herm_mult
        self._growth_rates = None

    def _applicable_schemes(self):
        '''The schemes this system admits, in the order `translate` prefers them.'''
        schemes = []
        if self.translation._is_translate_parameter_compatible(self.system):
            schemes.append(PARAMETER_SCHEME)
        columns = self.translation.scaling_matrix.shape[1]
        if ((len(self.system.variables) == columns + 1)
                or self.translation.scaling_matrix[
                    :, self.system.indep_var_index].is_zero_matrix):
            schemes.append(DEP_VAR_SCHEME)
        if len(self.system.variables) == columns:
            schemes.append(GENERAL_SCHEME)
        if not schemes:
            raise ValueError(
                "System doesn't have the right number of variables for translation")
        return schemes

    def _layout_for(self, scheme):
        '''
        Work out how the reduced system's variables correspond to the reduction's matrices
        under one scheme, or return None if they cannot.
        '''
        if scheme == PARAMETER_SCHEME:
            source = self.translation
            # The reduced system carries no auxiliaries, and its variables line up
            # one-for-one with the columns of V_n -- the new independent variable among
            # them, since it is itself an invariant.
            if len(self.reduced.variables) != source.herm_mult_n.cols:
                return None
            return (list(self.system.variables), source, [],
                    list(self.reduced.variables), False)

        if scheme == DEP_VAR_SCHEME:
            source = _dep_var_translation(self.system, self.translation)
            acted = [v for v in self.system.variables if v != self.system.indep_var]
        else:
            source = self.translation
            acted = list(self.system.variables)

        # Here the reduced system keeps the original independent variable, and its
        # remaining variables are the auxiliaries followed by the invariants.
        rest = [v for v in self.reduced.variables if v != self.reduced.indep_var]
        if len(rest) == source.herm_mult_i.cols + source.herm_mult_n.cols:
            return (acted, source, rest[:source.herm_mult_i.cols],
                    rest[source.herm_mult_i.cols:], True)
        if len(rest) == source.herm_mult_n.cols:
            # Reduced with include_aux_vars=False; the auxiliaries were dropped.
            return (acted, source, [], rest, True)
        return None

    @property
    def r(self):
        '''
        int: The number of scaling symmetries, and so the number of auxiliary variables the
        reduction hides.  When the reduced system does not carry them, this is how many
        values :meth:`reverse` needs in ``known_values``.
        '''
        return self._herm_mult_i.cols

    @property
    def carries_auxiliaries(self):
        '''
        bool: Whether the reduced system carries the auxiliary variables itself.  When it
        does, :meth:`reverse` needs no ``known_values``.
        '''
        return bool(self._auxiliaries)

    @property
    def auxiliary_variables(self):
        '''tuple: The auxiliary variables of the reduced system, empty if it has none.'''
        return tuple(self._auxiliaries)

    @property
    def invariant_variables(self):
        '''tuple: The variables of the reduced system that are invariants.'''
        return tuple(self._invariants)

    def forward(self, values):
        '''
        Translate values of the original system into values of the reduced system.

        Args:
            values (dict): Values of variables of the original system.  Anything the given
                values determine is returned, so a partial dictionary is fine -- passing
                only the constants and the independent variable translates the time axis.

        Returns:
            dict: Keyed by the variables of the reduced system.

        Raises:
            ValueError: If ``values`` names something that is not a variable of the original
                system, or determines nothing at all.
        '''
        values = self._check(values, self.system.variables, 'the original system')
        ordered = [values.get(v) for v in self._acted]

        result = {}
        result.update(self._evaluate(ordered, self._herm_mult_i, self._auxiliaries))
        result.update(self._evaluate(ordered, self._herm_mult_n, self._invariants))
        if self._shares_indep_var and values.get(self.system.indep_var) is not None:
            result[self.reduced.indep_var] = values[self.system.indep_var]

        if not result:
            raise ValueError(
                'The values given determine nothing.  Supply values for more of: '
                '{}.'.format(', '.join(map(str, self._acted))))
        return result

    def reverse(self, values, known_values=None, auxiliaries=None, invariants_at=None,
                tolerance=CONSISTENCY_TOLERANCE, **solver_options):
        '''
        Translate values of the reduced system back into values of the original system.

        There are two situations, told apart by the shape of ``values``.

        **A single point** -- every value a scalar.  Any variable of the original system may
        be given in ``known_values``, since at one point a dependent variable is just a
        number.  Purely arithmetic.

        **A series** -- some value an array, such as a solution over time.  Then
        ``known_values`` may name only *constants* of the original system.  A dependent
        variable is a function of time, not a number, and the series already says how it
        varies; giving it a single value would mix point translation with series
        translation.  If the reduced system does not carry its auxiliary variables and
        they vary along the solution, they are recovered by integrating their equations
        -- see :meth:`recover_auxiliaries` -- which needs :mod:`scipy`.

        Args:
            values (dict): Values of variables of the reduced system.  A partial dictionary
                is fine; anything it determines is returned.
            known_values (dict, optional): Values of variables of the original system that
                you already know.  Needed when the reduced system does not carry the
                auxiliary variables -- see :attr:`carries_auxiliaries` -- and rejected when
                it does, since then the original system is already determined.  In a series,
                constants only.
            auxiliaries (iter, optional): The :attr:`r` auxiliary values themselves, in the
                order of :attr:`auxiliary_variables`, as an alternative to ``known_values``.
                Obtain them from :meth:`auxiliaries_from` at a single point.  This is how a
                point value of a dependent variable is used to translate a series.
            invariants_at (callable, optional): Series only, and only when the auxiliaries
                have to be integrated: ``f(t)`` returning the invariants at time ``t``.  See
                :meth:`recover_auxiliaries`.
            tolerance (float): Relative tolerance for the known values to agree with each
                other and with the reduced values.
            **solver_options: Series only: passed to :func:`scipy.integrate.solve_ivp` if
                the auxiliaries have to be integrated.

        Returns:
            dict: Keyed by the variables of the original system.

        Raises:
            InsufficientKnownValues: If the auxiliaries are needed and neither
                ``known_values`` nor ``auxiliaries`` pins them down.
            ConflictingKnownValues: If the known values contradict each other.
        '''
        values = self._check(values, self.reduced.variables, 'the reduced system')
        if known_values is not None and auxiliaries is not None:
            raise ValueError('Give known_values or auxiliaries, not both.')

        if any(np.ndim(v) for v in values.values()):
            found = self._reverse_series(values, known_values, auxiliaries, invariants_at,
                                         tolerance, solver_options)
        else:
            found = self._reverse_point(values, known_values, auxiliaries, tolerance)

        ordered = list(found) + [values.get(v) for v in self._invariants]
        result = self._evaluate(ordered, self._inv_herm_mult, self._acted)
        if self._shares_indep_var and values.get(self.reduced.indep_var) is not None:
            # Under these schemes the reduced system keeps the original independent
            # variable, so take it as given rather than reconstructing it -- the monomial
            # would agree, but only to the accuracy of everything that went into it.
            result[self.system.indep_var] = values[self.reduced.indep_var]

        if not result:
            raise ValueError(
                'The values given determine nothing.  Supply values for more of: '
                '{}.'.format(', '.join(map(str, self._invariants))))
        return result

    def _carried_auxiliaries(self, values, known_values, auxiliaries):
        '''Read the auxiliaries from the reduced values when the reduced system has them.'''
        if known_values or auxiliaries is not None:
            raise ValueError(
                'Neither known_values nor auxiliaries is needed here: the reduced system '
                'carries the auxiliary variables {} itself, so the original system is '
                'already determined.'.format(', '.join(map(str, self._auxiliaries))))
        missing = [x for x in self._auxiliaries if values.get(x) is None]
        if missing:
            raise InsufficientKnownValues(
                'No value given for the auxiliary variable{} {}, which the original system '
                'cannot be recovered without.'.format(
                    '' if len(missing) == 1 else 's', ', '.join(map(str, missing))))
        return [values[x] for x in self._auxiliaries]

    def _given_auxiliaries(self, auxiliaries):
        auxiliaries = list(auxiliaries)
        if len(auxiliaries) != self.r:
            raise ValueError('Expected {} auxiliary value{}, not {}.'.format(
                self.r, '' if self.r == 1 else 's', len(auxiliaries)))
        if any(np.ndim(x) for x in auxiliaries):
            raise ValueError('The auxiliaries must be single numbers: their values at one '
                             'point, from which a series is rebuilt.')
        return auxiliaries

    def _reverse_point(self, values, known_values, auxiliaries, tolerance):
        '''
        One point.  The auxiliaries need not be constant along a solution to have a value
        here, so this works under every scheme, and any variable may be a known value.
        '''
        if self.carries_auxiliaries:
            return self._carried_auxiliaries(values, known_values, auxiliaries)
        if auxiliaries is not None:
            return self._given_auxiliaries(auxiliaries)
        return self.auxiliaries_from(values, known_values or {}, tolerance)

    def _reverse_series(self, values, known_values, auxiliaries, invariants_at, tolerance,
                        solver_options):
        '''
        A series.  Known values must be constants; a dependent variable given a single
        value is the caller mixing the two situations up.
        '''
        if self.carries_auxiliaries:
            return self._carried_auxiliaries(values, known_values, auxiliaries)
        return self.recover_auxiliaries(values, known_values, auxiliaries,
                                        invariants_at=invariants_at, tolerance=tolerance,
                                        **solver_options)

    def _reject_time_varying(self, known_values):
        '''In a series, a known value of anything but a constant is a category error.'''
        constants = set(self.system.constant_variables)
        for variable in known_values:
            if variable in constants:
                continue
            hint = self.system.initial_conditions.get(variable)
            if hint is not None:
                remedy = ('If you know its initial value, supply {}, which stands for '
                          '{}(0).'.format(hint, variable))
            else:
                remedy = ('To use a value of it at one point, compute the auxiliaries there '
                          'with auxiliaries_from() and pass them as auxiliaries=.')
            raise ValueError(
                '{} is a function of time, but a series is being translated, and the series '
                'already says how {} varies.  A single value for it mixes point translation '
                'with series translation.  {}'.format(variable, variable, remedy))

    def auxiliaries_from(self, values, known_values, tolerance=CONSISTENCY_TOLERANCE,
                         allowed=None):
        '''
        Solve for the auxiliary variables at one point, from values you already know there.

        Each original variable is a monomial in the auxiliaries and the invariants, so every
        known value gives one equation in the :attr:`r` auxiliaries.  Any number of known
        values may be supplied: too few, or a set that happens to constrain the same scale
        twice, raises :class:`InsufficientKnownValues` with a message saying what would
        complete it; more than :attr:`r` is fine provided they agree, and
        :class:`ConflictingKnownValues` is raised if they do not.

        This is the point situation of :meth:`reverse`: ``values`` should be scalars, the
        reduced values at the one point where the known values hold, and any variable of the
        original system may be known there.  The result can then be handed to
        :meth:`reverse` as ``auxiliaries=`` to translate a whole series.

        Args:
            values (dict): Values of variables of the reduced system at one point.
            known_values (dict): Values of variables of the original system at that point,
                each a single number.
            tolerance (float): Relative tolerance for the known values to agree with each
                other and with the reduced values.
            allowed (iter, optional): Restrict which variables may be known.  Used by the
                series situation, where only constants may.

        Returns:
            list: The auxiliary values, in the order of the columns of
            :attr:`~desr.ode_translation.ODETranslation.herm_mult_i`.
        '''
        known_values = {sympy.sympify(k): v for k, v in known_values.items()}
        allowed = list(self._acted) if allowed is None else list(allowed)
        unusable = set(known_values) - set(allowed)
        if unusable:
            raise InsufficientKnownValues(
                'Cannot use {} to determine the original system.  Expected some of: '
                '{}.'.format(', '.join(sorted(map(str, unusable))),
                             ', '.join(map(str, allowed))))
        for variable, known in known_values.items():
            if np.ndim(known):
                raise ValueError(
                    'The known value of {} must be a single number, not an array: its value '
                    'at the one point being translated.'.format(variable))

        W = self._inv_herm_mult
        columns = [self._acted.index(k) for k in known_values]
        # One row per known value: its exponents in the auxiliaries.
        coefficients = sympy.Matrix(len(columns), self.r,
                                    lambda row, j: W[j, columns[row]])

        # Divide out of each known value the part the invariants contribute.  What is left
        # is a monomial in the auxiliaries alone.
        invariants = [None if values.get(y) is None else _first(values[y])
                      for y in self._invariants]
        residuals = []
        for (variable, known), i in zip(known_values.items(), columns):
            needed = [self._invariants[j] for j in range(len(self._invariants))
                      if W[self.r + j, i] != 0]
            absent = [y for y in needed if values.get(y) is None]
            if absent:
                raise InsufficientKnownValues(
                    'Cannot use the known value of {} without values for {} from the '
                    'reduced system.'.format(variable, ', '.join(map(str, absent))))
            if known == 0:
                raise ValueError(
                    'The known value of {} is zero, and a zero cannot determine a scale.  '
                    'Every variable of the original system is a product of powers of the '
                    'reduced ones, so no rescaling turns a nonzero value into zero or back.  '
                    'Supply the value of a variable that is nonzero here.'.format(variable))
            invariant_part = _monomial(invariants, W[self.r:, :], i)
            if invariant_part == 0:
                zero = [y for y in needed if _first(values[y]) == 0]
                raise ValueError(
                    'The known value of {} cannot be used: at this point the reduced '
                    'value{} of {} {} zero, and a zero cannot determine a scale.  Supply the '
                    'value of a variable that is nonzero here instead.'.format(
                        variable, '' if len(zero) == 1 else 's', ', '.join(map(str, zero)),
                        'is' if len(zero) == 1 else 'are'))
            residuals.append(known / invariant_part)

        rank = coefficients.rank()
        if rank < self.r:
            # Before saying more are needed, make sure the ones given do not already
            # contradict each other -- that is the more useful thing to report.
            self._check_consistent(coefficients, residuals, known_values, tolerance)
            raise InsufficientKnownValues(
                self._shortfall(coefficients, known_values, invariants, allowed))

        # Solve on an independent subset, then hold every known value to the answer.
        independent = list(coefficients.T.rref()[1])
        square = coefficients.extract(independent, range(self.r))
        exponents = square.inv()
        if not all(entry.is_integer for entry in exponents):
            raise NotImplementedError(
                'Determining the auxiliaries from {} needs fractional powers, which are not '
                'yet supported.  Choose a different set of known values.'.format(
                    ', '.join(map(str, known_values))))
        auxiliaries = [_monomial([residuals[row] for row in independent], exponents.T, j)
                       for j in range(self.r)]

        disagreements = []
        for row, ((variable, known), residual) in enumerate(zip(known_values.items(),
                                                                residuals)):
            implied = _monomial(auxiliaries, coefficients.T, row)
            if not np.isclose(implied, residual, rtol=tolerance, atol=0):
                disagreements.append((variable, known, known * implied / residual))
        if disagreements:
            raise ConflictingKnownValues(self._conflict(known_values, disagreements))
        return auxiliaries

    @property
    def auxiliaries_are_constant(self):
        '''
        bool: Whether the auxiliary variables are constant along a solution.

        When they are, recovering the original system is arithmetic.  When they are not,
        they satisfy :math:`dx_j/dt = x_j H_j`, which :meth:`recover_auxiliaries` integrates
        when a series is translated.

        The auxiliaries are always constant under the parameter scheme, and they are
        constant under the general scheme whenever the original system is autonomous -- so
        the integration is rarer than it looks.
        '''
        return all(rate == 0 for rate in self.auxiliary_growth_rates())

    def auxiliary_growth_rates(self):
        '''
        The relative growth rates of the auxiliary variables.

        The reduction gives each auxiliary an equation of the form
        :math:`dx_j / dt = x_j H_j`, where :math:`H_j` involves the invariants and the
        independent variable but never the auxiliaries themselves.  The auxiliaries are
        therefore a quadrature on top of the invariants, not coupled to them, which is what
        makes ``include_aux_vars=False`` a legitimate thing to do.

        Returns:
            list: One :math:`H_j` per auxiliary, as a :mod:`sympy` expression.

        An autonomous system has auxiliaries that do not move, even under the general
        scheme, so recovering the original system needs no integration.

        >>> from desr.numerics import NumericTranslation
        >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)']
        >>> system = ODESystem.from_equations(equations)
        >>> translation = ODETranslation.from_ode_system(system)
        >>> autonomous = NumericTranslation(system, translation,
        ...                                 translation.translate_general(system))
        >>> autonomous.scheme
        'general'
        >>> autonomous.auxiliary_growth_rates()
        [0, 0, 0]
        >>> autonomous.auxiliaries_are_constant
        True

        A system with the independent variable in its right-hand side need not be so kind.

        >>> equations = ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
        ...              'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)']
        >>> system = ODESystem.from_equations(equations)
        >>> system.reorder_variables(['t', 'z1', 'z2'])
        >>> translation = ODETranslation.from_ode_system(system)
        >>> moving = NumericTranslation(
        ...     system, translation,
        ...     translation.translate_general(system, include_aux_vars=False))
        >>> moving.carries_auxiliaries
        False
        >>> moving.auxiliaries_are_constant
        False
        >>> moving.auxiliary_growth_rates()
        [(2*y0*(y1 + 1)/3 + y1)/(t*y0)]
        '''
        if self._growth_rates is None:
            if self.scheme == PARAMETER_SCHEME:
                # The parameter scheme normalises only constants, so its auxiliaries never
                # move, and the reduced system does not carry them to ask.
                self._growth_rates = [sympy.Integer(0)] * self.r
            else:
                if self.carries_auxiliaries:
                    carrier, auxiliaries = self.reduced, self._auxiliaries
                elif self.scheme == GENERAL_SCHEME:
                    carrier = self.translation.translate_general(self.system,
                                                                 include_aux_vars=True)
                    auxiliaries = [v for v in carrier.variables
                                   if v != carrier.indep_var][:self.r]
                else:
                    carrier = self.translation.translate_dep_var(self.system)
                    auxiliaries = [v for v in carrier.variables
                                   if v != carrier.indep_var][:self.r]
                rates = [sympy.simplify(carrier.derivative_dict[x] / x)
                         for x in auxiliaries]
                for rate, x in zip(rates, auxiliaries):
                    if rate.has(*auxiliaries):
                        raise ValueError(
                            'The equation for {} depends on an auxiliary variable, so the '
                            'auxiliaries are not a quadrature.'.format(x))
                self._growth_rates = rates
        return list(self._growth_rates)

    def recover_auxiliaries(self, values, known_values=None, auxiliaries=None,
                            invariants_at=None, tolerance=CONSISTENCY_TOLERANCE,
                            **solver_options):
        '''
        Recover the auxiliary variables along a series, when the reduced system does not
        carry them.

        This is what :meth:`reverse` does for a series; call it directly to see the
        auxiliaries themselves rather than the original system rebuilt from them.

        The auxiliaries obey :math:`dx_j / dt = x_j H_j`, driven by the invariants and
        uncoupled from each other, so they are recovered by integrating
        :math:`d(\\log|x_j|)/dt = H_j` once along the solution from their values at the
        first sample.  Their sign cannot change, since :math:`x_j = 0` is invariant, so it
        is fixed at the start and only the magnitude is integrated.  When the auxiliaries
        turn out to be constant no integration happens at all.

        The values at the first sample come from one of two places: from ``known_values``,
        which in a series may name only *constants* of the original system; or from
        ``auxiliaries`` directly, as computed by :meth:`auxiliaries_from` at a point -- the
        way to use a point value of a dependent variable, such as an initial condition for
        a system that has no constant standing for it.

        Args:
            values (dict): Values of variables of the reduced system: the independent
                variable as an array of times, and every invariant as an array of the same
                length or a scalar.
            known_values (dict, optional): Values of constants of the original system.
            auxiliaries (iter, optional): The :attr:`r` auxiliary values at the first
                sample, as an alternative to ``known_values``.
            invariants_at (callable, optional): ``f(t)`` returning the invariants at time
                ``t``, in the order of :attr:`invariant_variables`.  Pass the ``sol``
                attribute of a :func:`scipy.integrate.solve_ivp` result computed with
                ``dense_output=True`` for the most accurate integration.  Without it the
                sampled values are interpolated.
            tolerance (float): Relative tolerance for the known values to agree with each
                other and with the reduced values.
            **solver_options: Passed to :func:`scipy.integrate.solve_ivp`.

        Returns:
            list: One value per auxiliary, in the order of the columns of
            :attr:`~desr.ode_translation.ODETranslation.herm_mult_i`: an array along the
            series if it moves, a scalar if it is constant.
        '''
        values = self._check(values, self.reduced.variables, 'the reduced system')
        if self.carries_auxiliaries:
            raise ValueError(
                'The reduced system carries the auxiliary variables {} itself, so there is '
                'nothing to recover.  Read them from the solution.'.format(
                    ', '.join(map(str, self._auxiliaries))))
        if known_values is not None and auxiliaries is not None:
            raise ValueError('Give known_values or auxiliaries, not both.')

        times = values.get(self.reduced.indep_var)
        if times is None:
            raise ValueError('No values given for {}, the independent variable of the '
                             'reduced system.'.format(self.reduced.indep_var))
        times = np.asarray(times, dtype=float)
        absent = [y for y in self._invariants if values.get(y) is None]
        if absent:
            raise ValueError('No values given for the invariant{} {}.'.format(
                '' if len(absent) == 1 else 's', ', '.join(map(str, absent))))
        sampled = np.array([np.broadcast_to(np.asarray(values[y], dtype=float), times.shape)
                            for y in self._invariants])

        if auxiliaries is not None:
            start = self._given_auxiliaries(auxiliaries)
        else:
            known_values = {sympy.sympify(k): v for k, v in (known_values or {}).items()}
            self._reject_time_varying(known_values)
            at_first = {self.reduced.indep_var: times[0]}
            at_first.update(dict(zip(self._invariants, sampled[:, 0])))
            constants = [v for v in self._acted if v in self.system.constant_variables]
            start = self.auxiliaries_from(at_first, known_values, tolerance,
                                          allowed=constants)

        rates = self.auxiliary_growth_rates()
        if all(rate == 0 for rate in rates):
            # Constant along the series, so they stay scalars; anything rebuilt from them
            # and a constant invariant is then a scalar too, as a constant should be.
            return list(start)
        return self._integrate_auxiliaries(times, sampled, start, rates, invariants_at,
                                           solver_options)

    def _integrate_auxiliaries(self, times, sampled, start, rates, invariants_at, options):
        '''Integrate d(log|x|)/dt = H along the solution.  Needs scipy.'''
        try:
            from scipy.integrate import solve_ivp
        except ImportError:  # pragma: no cover
            raise ImportError(
                'Recovering auxiliary variables that vary along the solution needs scipy.  '
                'Install it with `pip install scipy`.')

        if invariants_at is None:
            from scipy.interpolate import CubicSpline
            invariants_at = CubicSpline(times, sampled, axis=1)

        rate_of = sympy.lambdify([self.reduced.indep_var, list(self._invariants)],
                                 list(rates), modules='numpy')
        signs = [np.sign(float(x)) for x in start]
        if any(sign == 0 for sign in signs):
            raise ValueError(
                'An auxiliary variable is zero at the first time point, where the change of '
                'variables is singular, so the solution cannot be translated back.')

        integrated = solve_ivp(
            lambda t, u: np.asarray(rate_of(t, np.asarray(invariants_at(t))), dtype=float),
            (times[0], times[-1]), np.log(np.abs([float(x) for x in start])),
            t_eval=times, **dict({'rtol': 1e-10, 'atol': 1e-12}, **options))
        if not integrated.success:  # pragma: no cover
            raise RuntimeError('Could not integrate the auxiliary variables: {}'.format(
                integrated.message))
        return [sign * np.exp(row) for sign, row in zip(signs, integrated.y)]

    def _evaluate(self, ordered, exponents, outputs):
        result = {}
        for j, out in enumerate(outputs):
            required = [i for i in range(exponents.rows) if exponents[i, j] != 0]
            if all(ordered[i] is not None for i in required):
                result[out] = _monomial(ordered, exponents, j)
        return result

    def _check(self, values, variables, description):
        values = {sympy.sympify(k): v for k, v in values.items()}
        unknown = set(values) - set(variables)
        if unknown:
            raise ValueError(
                'Not variables of {}: {}.  Expected some of: {}.'.format(
                    description, ', '.join(sorted(map(str, unknown))),
                    ', '.join(map(str, variables))))
        return values

    def _check_consistent(self, coefficients, residuals, known_values, tolerance):
        '''
        With too few independent equations to solve, the known values can still contradict
        each other -- two that constrain the same scale but disagree about it.  Check in
        log-magnitude space, where the equations are linear.
        '''
        if coefficients.rows < 2:
            return
        magnitudes = np.log(np.abs(np.array(residuals, dtype=float)))
        matrix = np.array(coefficients.tolist(), dtype=float)
        fit = np.linalg.lstsq(matrix, magnitudes, rcond=None)[0]
        misfit = np.abs(matrix @ fit - magnitudes)
        if misfit.max() > tolerance:
            worst = int(misfit.argmax())
            variables = list(known_values)
            raise ConflictingKnownValues(
                'The known values {} cannot all hold at once given the reduced values: they '
                'constrain the same scale and disagree about it.  The disagreement is largest '
                'at {}.'.format(', '.join(map(str, variables)), variables[worst]))

    def _conflict(self, known_values, disagreements):
        lines = ['The known values disagree with each other given the reduced values.']
        for variable, given, implied in disagreements:
            lines.append('    {} was given as {:g}, but the others imply {:g}.'.format(
                variable, float(given), float(implied)))
        lines.append('Check the values, or supply fewer of them: {} independent value{} '
                     'determine the original system.'.format(self.r, '' if self.r == 1 else 's'))
        return '\n'.join(lines)

    def _shortfall(self, coefficients, known_values, invariants, allowed):
        '''
        Explain what is still free after the given known values, and what would fix it.

        The kernel of the coefficient matrix is the set of rescalings of the auxiliaries that
        leave every known value unchanged.  Mapped through W, each kernel direction is a
        rescaling of the original system; the variables it moves are exactly the ones whose
        values would help -- provided they may be known here, and can be used, which a
        variable whose reduced value is zero at this point cannot.
        '''
        W = self._inv_herm_mult
        kernel = coefficients.nullspace() if coefficients.rows else [
            sympy.eye(self.r).col(j) for j in range(self.r)]
        missing = len(kernel)

        candidates = []
        for i, variable in enumerate(self._acted):
            if variable in known_values or variable not in allowed:
                continue
            moved = any(sum(W[j, i] * direction[j] for j in range(self.r)) != 0
                        for direction in kernel)
            needed = [j for j in range(len(self._invariants)) if W[self.r + j, i] != 0]
            usable = all(invariants[j] is not None and invariants[j] != 0 for j in needed)
            if moved and usable:
                candidates.append(variable)

        given = ', '.join(map(str, known_values))
        if not known_values:
            opening = 'No known values were supplied, and {} {} needed.'.format(
                self.r, 'is' if self.r == 1 else 'are')
        else:
            opening = '{} {} not determine the original system: {} more value{} {} needed.'.format(
                given, 'does' if len(known_values) == 1 else 'do',
                missing, '' if missing == 1 else 's', 'is' if missing == 1 else 'are')

        lines = [opening]
        if missing == 1 and known_values:
            lines.append('The rescaling {} leaves {} unchanged.'.format(
                self._describe_rescaling(kernel[0]), given))
        lines.append('The reduced system is shared by an entire {}-parameter family of '
                     'original systems, and this does not choose between them.'.format(
                         self.r))
        if missing == 1:
            lines.append('Supply the value of one of: {}.'.format(
                ', '.join(map(str, candidates))))
        else:
            lines.append('Supply values for {} of: {}.  Not every choice works; if yours '
                         'does not, this message will say what is still free.'.format(
                             missing, ', '.join(map(str, candidates))))
        return '\n'.join(lines)

    def _describe_rescaling(self, direction):
        '''Render one kernel direction as the rescaling of original variables it induces.'''
        W = self._inv_herm_mult
        by_power = {}
        for i, variable in enumerate(self._acted):
            power = sum(W[j, i] * direction[j] for j in range(self.r))
            if power != 0:
                by_power.setdefault(power, []).append(str(variable))
        parts = []
        for power in sorted(by_power, key=lambda p: (-sympy.sign(p), abs(p))):
            factor = 'lambda' if abs(power) == 1 else 'lambda^{}'.format(abs(power))
            parts.append('{} by {}'.format(', '.join(by_power[power]),
                                           factor if power > 0 else '1/' + factor))
        return ' and '.join(parts)
