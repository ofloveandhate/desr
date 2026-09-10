"""
Numerical translation between an :class:`~desr.ode_system.ODESystem` and its reduction.

:mod:`desr` reduces a system symbolically.  This module carries *numbers* across the same
reduction: it turns parameter values and initial conditions of the original system into
their counterparts in the reduced system, and turns a numerical solution of the reduced
system back into a solution of the original one.

This module needs :mod:`numpy`.  It is deliberately not imported by ``desr/__init__.py``,
so :mod:`desr` itself keeps no numerical dependencies.

Reduction discards information, and no amount of arithmetic brings it back.  A reduced
system with :math:`r` scaling symmetries is shared by an entire :math:`r`-parameter family
of original systems.  To come back you must supply :math:`r` values you already know --
see :meth:`NumericTranslation.reverse`.
"""

import sympy

try:
    import numpy as np
except ImportError:  # pragma: no cover
    raise ImportError('desr.numerics requires numpy.  Install it with `pip install numpy`.')

from .ode_translation import ODETranslation
from .ode_system import ODESystem

__all__ = ['NumericTranslation', 'InsufficientKnownValues']


class InsufficientKnownValues(ValueError):
    '''
    Raised when the values supplied to :meth:`NumericTranslation.reverse` do not pin down a
    unique original system.
    '''
    pass


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


class NumericTranslation(object):
    '''
    Carry numbers back and forth across a reduction.

    Args:
        system (ODESystem): The original system.
        translation (ODETranslation): The reduction of that system.
        reduced (ODESystem, optional): The reduced system.  Computed with
            :meth:`~desr.ode_translation.ODETranslation.translate` if not given.

    Dictionaries handed in and out are keyed by the actual symbols of each system.  Values
    may be scalars or :class:`numpy.ndarray`, so a whole time series translates in one call.
    '''

    def __init__(self, system, translation, reduced=None):
        self.system = system
        self.translation = translation
        self.reduced = translation.translate(system) if reduced is None else reduced

        if not translation._is_translate_parameter_compatible(system):
            raise NotImplementedError(
                'desr.numerics currently handles the parameter reduction scheme only.  The '
                'dependent-variable and general schemes are not yet supported.')

        # In the parameter scheme the reduced system has no auxiliary variables, and its
        # variables line up one-for-one with the columns of the Hermite multiplier V_n.
        self._invariants = list(self.reduced.variables)
        if len(self._invariants) != translation.herm_mult_n.cols:
            raise ValueError(
                'Reduced system has {} variables but the reduction produces {} '
                'invariants.'.format(len(self._invariants), translation.herm_mult_n.cols))

    @property
    def r(self):
        '''
        int: The number of scaling symmetries, and so the number of values that
        :meth:`reverse` needs supplied in ``known_values``.
        '''
        return self.translation.r

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
        ordered = [values.get(v) for v in self.system.variables]
        return self._evaluate(ordered, self.translation.herm_mult_n, self._invariants,
                              values, 'forward')

    def reverse(self, values, known_values):
        '''
        Translate values of the reduced system back into values of the original system.

        Args:
            values (dict): Values of variables of the reduced system.  A partial dictionary
                is fine; anything it determines is returned.
            known_values (dict): Values of :attr:`r` variables of the *original* system that
                you already know.  These choose one system out of the family that share this
                reduction.

        Returns:
            dict: Keyed by the variables of the original system.

        Raises:
            InsufficientKnownValues: If the wrong number of values is supplied, or if the
                ones supplied do not pin down a unique original system.
        '''
        values = self._check(values, self._invariants, 'the reduced system')
        auxiliaries = self.auxiliaries_from(values, known_values)
        ordered = list(auxiliaries) + [values.get(v) for v in self._invariants]
        available = dict(values)
        available.update({sympy.Symbol('_aux{}'.format(i)): a
                          for i, a in enumerate(auxiliaries)})
        names = [sympy.Symbol('_aux{}'.format(i)) for i in range(self.r)] + self._invariants
        return self._evaluate(ordered, self.translation.inv_herm_mult, self.system.variables,
                              dict(zip(names, ordered)), 'reverse', names=names)

    def auxiliaries_from(self, values, known_values):
        '''
        Solve for the auxiliary variables implied by ``known_values``.

        Each original variable is a monomial in the auxiliaries and the invariants.  The
        :attr:`r` variables named in ``known_values`` therefore give :attr:`r` equations in
        the :attr:`r` auxiliaries, solved here exactly over the rationals.

        Args:
            values (dict): Values of variables of the reduced system.
            known_values (dict): Values of :attr:`r` variables of the original system.

        Returns:
            list: The auxiliary values, in the order of the columns of
            :attr:`~desr.ode_translation.ODETranslation.herm_mult_i`.
        '''
        known_values = {sympy.sympify(k): v for k, v in known_values.items()}
        unknown = set(known_values) - set(self.system.variables)
        if unknown:
            raise InsufficientKnownValues(
                'Not variables of the original system: {}.  Expected some of: {}.'.format(
                    ', '.join(sorted(map(str, unknown))),
                    ', '.join(map(str, self.system.variables))))
        if len(known_values) != self.r:
            raise InsufficientKnownValues(self._shortfall(known_values))

        W = self.translation.inv_herm_mult
        columns = [self.system.variables.index(k) for k in known_values]

        # log v = V log x + U log y, with V square.  Invert V exactly over the rationals.
        coefficients = sympy.Matrix([[W[j, i] for j in range(self.r)] for i in columns])
        if coefficients.det() == 0:
            raise InsufficientKnownValues(self._shortfall(known_values, degenerate=True))

        # Divide out the part of each known value contributed by the invariants.
        invariants = [values.get(v) for v in self._invariants]
        residuals = []
        for known, i in zip(known_values.values(), columns):
            needed = [self._invariants[j] for j in range(len(self._invariants))
                      if W[self.r + j, i] != 0]
            missing = [n for n in needed if values.get(n) is None]
            if missing:
                raise InsufficientKnownValues(
                    'Cannot use the known value of {} without values for {} from the reduced '
                    'system.'.format(self.system.variables[i], ', '.join(map(str, missing))))
            residuals.append(known / _monomial(invariants, W[self.r:, :], i))

        exponents = coefficients.inv().T
        return [_monomial(residuals, exponents, j) for j in range(self.r)]

    def _evaluate(self, ordered, exponents, outputs, available, direction, names=None):
        names = list(available) if names is None else names
        result = {}
        for j, out in enumerate(outputs):
            required = [i for i in range(exponents.rows) if exponents[i, j] != 0]
            if all(ordered[i] is not None for i in required):
                result[out] = _monomial(ordered, exponents, j)
        if not result:
            raise ValueError(
                'The values given determine nothing when translating {}.  Supply values for '
                'more of: {}.'.format(direction, ', '.join(map(str, names))))
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

    def _shortfall(self, known_values, degenerate=False):
        given = ', '.join(map(str, known_values)) or 'nothing'
        if degenerate:
            opening = ('The {n} known value{s} supplied ({given}) do not determine the '
                       'original system: they overlap, and leave some of it free.')
        else:
            opening = ('Need {r} known value{plural} of the original system, but {n} '
                       '({given}) {were} supplied.')
        return (opening + '\nThe reduced system is shared by an entire {r}-parameter family '
                'of original systems, so this is not enough to choose between them.\n'
                'Supply values for {r} of: {candidates}.').format(
            r=self.r, n=len(known_values), given=given,
            s='' if len(known_values) == 1 else 's',
            plural='' if self.r == 1 else 's',
            were='was' if len(known_values) == 1 else 'were',
            candidates=', '.join(map(str, self.system.variables[1:])))
