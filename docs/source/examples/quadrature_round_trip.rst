Recovering auxiliary variables by quadrature
============================================

The point of reduction is to carry fewer symbols.  An application that solves a system
thousands of times -- fitting parameters, say -- wants the smallest system it can get, and
:meth:`~desr.ode_translation.ODETranslation.translate_general` offers
``include_aux_vars=False`` for exactly that: it returns the invariants alone and throws the
auxiliary variables away.

This page discusses this reduction, and how to get the original system back afterwards
anyway.  Here there is an, integration to do first.


When can the auxiliaries be dropped?
------------------------------------

The general reduction produces two blocks of equations: the invariants :math:`y`, which are
what the reduction is for, and the auxiliaries :math:`x`, which carry the scaling
information that the invariants discard (we're using notation compatible with :cite:`Hubert2013c`).  Dropping the :math:`x` block is only legitimate if
the :math:`y` block does not depend on it.

Fortunately, it does not.  Every auxiliary variable satisfies an equation of the form

.. math::
    :nowrap:

    \begin{align}
    \frac{dx_j}{dt} &= x_j \, H_j
    \end{align}

where :math:`H_j` involves the invariants and the independent variable but never the
auxiliaries themselves.  The auxiliaries are therefore a *quadrature* sitting on top of the
invariants rather than part of the system, and solving without them loses nothing that
cannot be recomputed.  :mod:`desr` checks this rather than assuming it.


Example: a system whose auxiliaries move
-----------------------------------------

Check out example 6.6 of :cite:`Hubert2013c` (starts very bottom of page 503).  The system is originally given as

.. math::
    :nowrap:

    \begin{align}
    t \frac{dz_1}{dt} &= z_1 \left(-\frac{2}{3} + \frac{1}{3} z_1^5 z_2 \right) \\
    t \frac{dz_2}{dt} &= z_2 \left( \frac{10}{3} - \frac{2}{3} z_1^5 z_2 + \frac{z_1^2 z_2}{t} \right)
    \end{align}

But desr expects the derivative only on the right hand side:

.. math::
    :nowrap:

    \begin{align}
    \frac{dz_1}{dt} &= \frac{z_1}{t} \left(-\frac{2}{3} + \frac{1}{3} z_1^5 z_2 \right) \\
    \frac{dz_2}{dt} &= \frac{z_2}{t} \left( \frac{10}{3} - \frac{2}{3} z_1^5 z_2 + \frac{z_1^2 z_2}{t} \right)
    \end{align}

which has the independent variable in its
right-hand side.

Build it, and reduce.

    >>> import numpy as np
    >>> from scipy.integrate import solve_ivp
    >>> from desr.numerics import NumericTranslation

    >>> eq1 = 'dz1/dt = z1/t * (-2/3 + 1/3 *z1**5 *z2 )'
    >>> eq2 = 'dz2/dt = z2/t * ( 10/3 - 2/3 *z1**5 *z2 + z1**2*z2/t )'
    >>> eqns = [eq1, eq2]
    >>> system = ODESystem.from_equations(eqns)
    >>> system.reorder_variables(['t', 'z1', 'z2'])
    >>> translation = ODETranslation.from_ode_system(system)
    >>> translation.invariants()
    Matrix([[t*z1**3, z1**5*z2]])
    >>> translation.auxiliaries()
    Matrix([[z1**4*z2]])

The paper has :math:`y_1 = t z_1^3`, :math:`y_2 = z_1^2 z_2 / t` and the auxiliary
:math:`t z_1^2`.  Those span the same lattice as desr's -- the paper's :math:`y_2` is desr's
:math:`z_1^5 z_2 / (t z_1^3)` -- and two column operations on the Hermite multiplier move
desr onto the paper's basis, so that everything below can be compared with p. 504 directly.

    >>> translation.multiplier_add_columns(2, 1, -1)   # y2  <-  y2 / y1
    >>> translation.multiplier_add_columns(0, 2, -1)   # x   <-  x / y2
    >>> translation.invariants()
    Matrix([[t*z1**3, z1**2*z2/t]])
    >>> translation.auxiliaries()
    Matrix([[t*z1**2]])

Reducing the usual way keeps the auxiliary :math:`x_0` as a variable of the reduced system,
and gives exactly the paper's equations (6.6) for :math:`x`, :math:`y_1` and :math:`y_2`:

    >>> translation.translate_general(system)
    dt/dt = 1
    dx0/dt = x0*(2*y0*y1/3 - 1/3)/t
    dy0/dt = y0*(y0*y1 - 1)/t
    dy1/dt = y1*(y1 + 1)/t

Note that :math:`x_0` appears nowhere but in its own equation.  Asking for the invariants
alone drops it, leaving one fewer equation to integrate:

    >>> reduced = translation.translate_general(system, include_aux_vars=False)
    >>> reduced
    dt/dt = 1
    dy0/dt = y0*(y0*y1 - 1)/t
    dy1/dt = y1*(y1 + 1)/t


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
    [(2*y0*y1 - 1)/(3*t)]
    >>> numeric.auxiliaries_are_constant
    False

Because they are not constant, they will have to be integrated.


The paper's exact solution
--------------------------

Hubert and Labahn solve the reduced system in closed form (p. 504).  With their
:math:`\ln(t - c_1)` written as :math:`\ln(c_1 - t)` -- the same derivative, and real for
the constants used below --

.. math::
    :nowrap:

    \begin{align}
    x &= \frac{c_3}{t^{1/3} \, L^{2/3}}, \qquad
    y_1 = \frac{c_1}{t \, L}, \qquad
    y_2 = \frac{t}{c_1 - t}, \qquad
    L = \ln(c_1 - t) - \ln t + c_2 .
    \end{align}

Write these down and check that they satisfy the system desr produced, which is the most
direct confirmation that the two reductions are the same.

    >>> t, z1, z2 = system.variables
    >>> y0, y1 = numeric.invariant_variables
    >>> with_aux = translation.translate_general(system)
    >>> x0 = with_aux.variables[1]
    >>> c1, c2, c3 = sympy.symbols('c1 c2 c3')
    >>> L = sympy.log(c1 - t) - sympy.log(t) + c2
    >>> exact = {x0: c3 / (t**sympy.Rational(1, 3) * L**sympy.Rational(2, 3)),
    ...          y0: c1 / (t * L),
    ...          y1: t / (c1 - t)}
    >>> [sympy.simplify(sympy.diff(exact[v], t) - with_aux.derivative_dict[v].subs(exact)) == 0
    ...  for v in (x0, y0, y1)]
    [True, True, True]

Fix the constants by the starting point we will use, :math:`t = 1`, :math:`z_1 = 0.5`,
:math:`z_2 = 0.2`, whose reduced values are

    >>> start = numeric.forward({t: 1.0, z1: 0.5, z2: 0.2})
    >>> [float(start[y]) for y in (y0, y1)]
    [0.125, 0.05]

and whose auxiliary is :math:`t z_1^2 = 0.25`.

    >>> constants = sympy.solve([exact[y1].subs(t, 1) - sympy.Rational(1, 20),
    ...                          exact[y0].subs(t, 1) - sympy.Rational(1, 8),
    ...                          exact[x0].subs(t, 1) - sympy.Rational(1, 4)],
    ...                         [c1, c2, c3], dict=True)[0]
    >>> constants
    {c1: 21, c2: 168 - log(20), c3: 21**(2/3)}

    >>> y0_exact = sympy.lambdify(t, exact[y0].subs(constants), modules='numpy')
    >>> y1_exact = sympy.lambdify(t, exact[y1].subs(constants), modules='numpy')


Solving the reduced system
--------------------------

Now solve the two-equation system numerically and hold it to the exact solution.

    >>> def as_function(a_system):
    ...     variables = list(a_system.non_constant_variables)
    ...     rhs = [a_system.derivative_dict[x] for x in variables]
    ...     return sympy.lambdify([a_system.indep_var, variables], rhs, modules='numpy')

    >>> times = np.linspace(1.0, 2.0, 200)
    >>> solution = solve_ivp(as_function(reduced), (1.0, 2.0),
    ...                      [float(start[y0]), float(start[y1])],
    ...                      t_eval=times, rtol=1e-11, atol=1e-13, dense_output=True)
    >>> bool(np.max(np.abs(solution.y[0] - y0_exact(times))) < 1e-9)
    True
    >>> bool(np.max(np.abs(solution.y[1] - y1_exact(times))) < 1e-9)
    True


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
    [0.25]

That is :math:`t z_1^2` at :math:`t = 1`, as it should be.  Then hand it to
:meth:`~desr.numerics.NumericTranslation.reverse` for the series, which integrates it
along the solution.  Passing the solver's dense output lets that integration evaluate the
invariants wherever it likes.

    >>> recovered = numeric.reverse(reduced_values, auxiliaries=at_start,
    ...                             invariants_at=solution.sol)

The paper also reverse-translates its exact solution, and desr's symbolic
:meth:`~desr.ode_translation.ODETranslation.reverse_translate_general` does the same
computation: it rebuilds :math:`z_1` and :math:`z_2` from :math:`x`, :math:`y_1` and
:math:`y_2`, then rescales time by a constant that the arbitrary :math:`c_i` leave free --
here :math:`c_3^3 / c_1^2`.  With the constants fixed by our starting point that factor is
exactly :math:`1`, which is why the numerical route needs no time rescaling at all.

    >>> sympy.simplify(constants[c3]**3 / constants[c1]**2)
    1
    >>> z1_exact, z2_exact = translation.reverse_translate_general(
    ...     (t, exact[x0], exact[y0], exact[y1]))
    >>> z1_exact = sympy.lambdify(t, z1_exact.subs(constants), modules='numpy')
    >>> z2_exact = sympy.lambdify(t, z2_exact.subs(constants), modules='numpy')
    >>> [round(float(f(1.0)), 12) for f in (z1_exact, z2_exact)]
    [0.5, 0.2]

So the recovered solution can be held to the exact one, not merely to another numerical
solve.

    >>> bool(np.max(np.abs(recovered[t] - times)) < 1e-12)
    True
    >>> bool(np.max(np.abs(recovered[z1] - z1_exact(times))) < 1e-8)
    True
    >>> bool(np.max(np.abs(recovered[z2] - z2_exact(times))) < 1e-7)
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

.. plot:: _examples/quadrature_recovery.py

The middle panel is the one this page is about.  :math:`x_0` was never solved for -- it is
absent from the reduced system entirely -- and the line is what
:meth:`~desr.numerics.NumericTranslation.recover_auxiliaries` produces by integrating its
equation along the invariants.  The circles are its exact value, the paper's :math:`x`
with the constants above, which the quadrature never sees.  They agree to about
:math:`10^{-10}`.

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
