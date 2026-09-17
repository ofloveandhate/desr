Recovering auxiliary variables by quadrature
============================================

The point of reduction is to carry fewer symbols.  An application that solves a system
thousands of times -- fitting parameters, say -- wants the smallest system it can get, and
:meth:`~desr.ode_translation.ODETranslation.translate_general` offers
``include_aux_vars=False`` for exactly that: it returns the invariants alone and throws the
auxiliary variables away.

This page shows why that is safe, and how to get the original system back afterwards
anyway.  The previous example recovered a solution by arithmetic; here there is an
integration to do first.


When can the auxiliaries be dropped?
------------------------------------

The general reduction produces two blocks of equations: the invariants :math:`y`, which are
what the reduction is for, and the auxiliaries :math:`x`, which carry the scaling
information that the invariants discard.  Dropping the second block is only legitimate if
the first does not depend on it.

It does not.  Every auxiliary satisfies an equation of the form

.. math::
    :nowrap:

    \begin{align}
    \frac{dx_j}{dt} &= x_j \, H_j
    \end{align}

where :math:`H_j` involves the invariants and the independent variable but never the
auxiliaries themselves.  The auxiliaries are therefore a *quadrature* sitting on top of the
invariants rather than part of the system, and solving without them loses nothing that
cannot be recomputed.  :mod:`desr` checks this rather than assuming it.


A system whose auxiliaries move
-------------------------------

We use example 6.6 of :cite:`Hubert2013c`, which has the independent variable in its
right-hand side.

.. math::
    :nowrap:

    \begin{align}
    \frac{dz_1}{dt} &= \frac{z_1 \left(z_1^5 z_2 - 2 \right)}{3t} \\
    \frac{dz_2}{dt} &= \frac{z_2 \left(10 - 2 z_1^5 z_2 + \frac{3 z_1^2 z_2}{t} \right)}{3t}
    \end{align}

Build it, and reduce.

    >>> import numpy as np
    >>> from scipy.integrate import solve_ivp
    >>> from desr.numerics import NumericTranslation

    >>> equations = ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
    ...              'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)']
    >>> system = ODESystem.from_equations(equations)
    >>> system.reorder_variables(['t', 'z1', 'z2'])
    >>> translation = ODETranslation.from_ode_system(system)
    >>> translation.invariants()
    Matrix([[t*z1**3, z1**5*z2]])
    >>> translation.auxiliaries()
    Matrix([[z1**4*z2]])

Reducing the usual way keeps the auxiliary :math:`x_0` as a variable of the reduced system:

    >>> translation.translate_general(system)
    dt/dt = 1
    dx0/dt = x0*(2*y1/3 + 2/3 + y1/y0)/t
    dy0/dt = y0*(y1 - 1)/t
    dy1/dt = y1*(y1 + y1/y0)/t

Note that :math:`x_0` appears nowhere but in its own equation.  Asking for the invariants
alone drops it, leaving one fewer equation to integrate:

    >>> reduced = translation.translate_general(system, include_aux_vars=False)
    >>> reduced
    dt/dt = 1
    dy0/dt = y0*(y1 - 1)/t
    dy1/dt = y1*(y1 + y1/y0)/t


What was dropped
----------------

:class:`~desr.numerics.NumericTranslation` recognises a reduction with its auxiliaries
missing.  Pass the reduced system explicitly, since it was not made with the default
options.

    >>> numeric = NumericTranslation(system, translation, reduced)
    >>> numeric.scheme
    'general'
    >>> numeric.carries_auxiliaries
    False
    >>> numeric.r
    1

It can still say what the dropped equation was.  This is the :math:`H_j` from above,
and it mentions :math:`y_0`, :math:`y_1` and :math:`t` -- but no :math:`x`.

    >>> numeric.auxiliary_growth_rates()
    [(2*y0*(y1 + 1)/3 + y1)/(t*y0)]
    >>> numeric.auxiliaries_are_constant
    False

Because they are not constant, they will have to be integrated.


Solving the reduced system
--------------------------

    >>> t, z1, z2 = system.variables
    >>> y0, y1 = numeric.invariant_variables

    >>> def as_function(a_system):
    ...     variables = list(a_system.non_constant_variables)
    ...     rhs = [a_system.derivative_dict[x] for x in variables]
    ...     return sympy.lambdify([a_system.indep_var, variables], rhs, modules='numpy')

    >>> start = numeric.forward({t: 1.0, z1: 0.5, z2: 0.2})
    >>> [float(start[y]) for y in (y0, y1)]
    [0.125, 0.00625]

    >>> times = np.linspace(1.0, 2.0, 200)
    >>> solution = solve_ivp(as_function(reduced), (1.0, 2.0),
    ...                      [float(start[y0]), float(start[y1])],
    ...                      t_eval=times, rtol=1e-11, atol=1e-13, dense_output=True)


Translating back
----------------

Translating a series back needs the auxiliary at its first sample, from which its own
equation carries it forward.  For a system with constants, the known values that fix it
would be constants.  This system has none -- its variables are :math:`t`, :math:`z_1` and
:math:`z_2`, and that is all -- so the only thing we can know is a *point* value: where the
solution started.

A point value of a dependent variable is not something a series can be told directly.  The
series already says how :math:`z_1` varies, so a single number for it mixes two different
kinds of statement, and :meth:`~desr.numerics.NumericTranslation.reverse` says so.

    >>> reduced_values = {reduced.indep_var: solution.t,
    ...                   y0: solution.y[0], y1: solution.y[1]}
    >>> numeric.reverse(reduced_values, known_values={z1: 0.5})
    Traceback (most recent call last):
        ...
    ValueError: z1 is a function of time, but a series is being translated, and the series already says how z1 varies.  A single value for it mixes point translation with series translation.  To use a value of it at one point, compute the auxiliaries there with auxiliaries_from() and pass them as auxiliaries=.

So it is done in two steps that say out loud they are two different facts.  First, at the
one point where :math:`z_1 = 0.5` is known -- the start, whose reduced values are the
scalars in ``start`` -- determine the auxiliary:

    >>> at_start = numeric.auxiliaries_from(start, known_values={z1: 0.5})
    >>> [round(float(x), 6) for x in at_start]
    [0.0125]

That is :math:`z_1^4 z_2` at :math:`t = 1`, as it should be.  Then hand it to
:meth:`~desr.numerics.NumericTranslation.reverse` for the series, which integrates it
along the solution.  Passing the solver's dense output lets that integration evaluate the
invariants wherever it likes.

    >>> recovered = numeric.reverse(reduced_values, auxiliaries=at_start,
    ...                             invariants_at=solution.sol)

Compare against solving the original system directly.

    >>> reference = solve_ivp(as_function(system), (1.0, 2.0), [0.5, 0.2],
    ...                       t_eval=times, rtol=1e-11, atol=1e-13)
    >>> bool(np.max(np.abs(recovered[t] - times)) < 1e-12)
    True
    >>> bool(np.max(np.abs(recovered[z1] - reference.y[0])) < 1e-6)
    True
    >>> bool(np.max(np.abs(recovered[z2] - reference.y[1])) < 1e-6)
    True

The agreement here is looser than in the arithmetic case, and it should be.  Reconstructing
:math:`z_2` raises the recovered auxiliary to the fifth power, so whatever error the
quadrature made is amplified along with it.  That is the price of the smaller system, and it
is paid once at the end rather than on every solve.

Internally the integration is done on :math:`\log |x_j|`, where the equation
:math:`dx_j/dt = x_j H_j` is exactly linear.  The auxiliaries are monomials and can span a
wide range of magnitudes, which logarithms handle comfortably.  Their sign cannot change,
since :math:`x_j = 0` is invariant, so it is read once from the starting value and only the
magnitude is integrated.


The whole workflow, drawn
-------------------------

``examples/quadrature_recovery.py`` runs all of the above and draws it.

.. plot:: ../../examples/quadrature_recovery.py

The middle panel is the one this page is about.  :math:`x_0` was never solved for -- it is
absent from the reduced system entirely -- and the line is what
:meth:`~desr.numerics.NumericTranslation.recover_auxiliaries` produces by integrating its
equation along the invariants.  The circles are its true value :math:`z_1^4 z_2`, read off
the direct solution, which the quadrature never sees.  They agree to about
:math:`3 \times 10^{-11}`.

Given that curve, the right-hand panel is arithmetic: each original variable is a product of
integer powers of :math:`x_0` and the invariants.


When there is nothing to integrate
----------------------------------

The quadrature above is the exception rather than the rule.  Auxiliaries only move when the
scaling acts on the dependent variables in a way that a plain rescaling of parameters cannot
reproduce, which in practice means the independent variable appearing in the right-hand
side.  An autonomous system -- which covers most models one would want to fit -- has
constant auxiliaries even under the general scheme.

    >>> equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)']
    >>> autonomous = ODESystem.from_equations(equations)
    >>> translation = ODETranslation.from_ode_system(autonomous)
    >>> dropped = translation.translate_general(autonomous, include_aux_vars=False)
    >>> numeric = NumericTranslation(autonomous, translation, dropped)
    >>> numeric.auxiliary_growth_rates()
    [0, 0, 0]
    >>> numeric.auxiliaries_are_constant
    True

:meth:`~desr.numerics.NumericTranslation.reverse` notices and skips the integration
entirely, so nothing is lost to it and :mod:`scipy` is never called.  Michaelis-Menten,
Lotka-Volterra and the chemical reaction networks all land here.
