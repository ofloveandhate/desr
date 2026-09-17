Recovering auxiliary variables by quadrature
============================================

One reason to do a scaling symmetry reduction is to produce a system fewer symbols.  An application that solves a system
thousands of times -- fitting parameters, say -- wants the smallest system it can get, and
:meth:`~desr.ode_translation.ODETranslation.translate_general` offers
``include_aux_vars=False`` for exactly that: it returns the invariants and reduced equations, and throws
auxiliary variables away.

This page discusses this reduction, and how to get the original system back afterwards
anyway.


What are auxiliaries?
------------------------------------

The general reduction produces two blocks of equations: invariants :math:`y`, which are
what the reduction is for, and the auxiliaries :math:`x` (we're using notation compatible with :cite:`Hubert2013c`).  

Dropping the :math:`x` block is only legitimate if
the :math:`y` block does not depend on it.
Fortunately, :math:`x` does not.  Every auxiliary variable satisfies an equation of the form

.. math::
    :nowrap:

    \begin{align}
    \frac{dx_j}{dt} &= x_j \, H_j
    \end{align}

where :math:`H_j` involves the invariants and the independent variable but never the
auxiliaries themselves.  The auxiliaries are therefore a *quadrature* sitting on top of the
invariants, rather than part of the system.  So, solving the :math:`y` equations without them loses nothing that
cannot be recomputed. 


Example
-----------------------------------------

Check out Example 6.6 of :cite:`Hubert2013c` (starts very bottom of page 503).  The system is originally given as

.. math::
    :nowrap:

    \begin{align}
    t \frac{dz_1}{dt} &= z_1 \left(-\frac{2}{3} + \frac{1}{3} z_1^5 z_2 \right) \\
    t \frac{dz_2}{dt} &= z_2 \left( \frac{10}{3} - \frac{2}{3} z_1^5 z_2 + \frac{z_1^2 z_2}{t} \right)
    \end{align}

But desr expects only the formal name of the derivative on the left hand side:

.. math::
    :nowrap:

    \begin{align}
    \frac{dz_1}{dt} &= \frac{z_1}{t} \left(-\frac{2}{3} + \frac{1}{3} z_1^5 z_2 \right) \\
    \frac{dz_2}{dt} &= \frac{z_2}{t} \left( \frac{10}{3} - \frac{2}{3} z_1^5 z_2 + \frac{z_1^2 z_2}{t} \right)
    \end{align}

Build it in Python with desr.

    >>> import numpy as np
    >>> from scipy.integrate import solve_ivp
    >>> from desr.numerics import NumericTranslation

    >>> eq1 = 'dz1/dt = z1/t * (-2/3 + 1/3 *z1**5 *z2 )'
    >>> eq2 = 'dz2/dt = z2/t * ( 10/3 - 2/3 *z1**5 *z2 + z1**2*z2/t )'
    >>> eqns = [eq1, eq2]
    >>> system = ODESystem.from_equations(eqns)

And nondimensionalize.

    >>> system.reorder_variables(['t', 'z1', 'z2'])
    >>> translation = ODETranslation.from_ode_system(system)
    >>> translation.invariants()
    Matrix([[t*z1**3, z1**5*z2]])
    >>> translation.auxiliaries()
    Matrix([[z1**4*z2]])

Hubert-Labahn has :math:`y_1 = t z_1^3`, :math:`y_2 = z_1^2 z_2 / t` and the auxiliary
:math:`t z_1^2`.  Those are products of integer powers of desr's, and desr's of theirs --
the paper's :math:`y_2` is desr's :math:`z_1^5 z_2 / (t z_1^3)`; two column operations on the Hermite multiplier move
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

Notably, :math:`x_0` appears nowhere but in its own equation.  Asking for the invariants
alone drops it, leaving one fewer equation to integrate:

    >>> reduced_system = translation.translate_general(system, include_aux_vars=False)
    >>> reduced_system
    dt/dt = 1
    dy0/dt = y0*(y0*y1 - 1)/t
    dy1/dt = y1*(y1 + 1)/t


What was dropped
----------------

:class:`~desr.numerics.NumericTranslation` is down with a reduction with its auxiliaries
dropped.  Pass the reduced system explicitly, since it was made with the option to drop auxiliaries.

    >>> numeric_translation = NumericTranslation(system, translation, reduced_system)
    >>> numeric_translation.scheme
    'general'
    >>> numeric_translation.carries_auxiliaries
    False
    >>> numeric_translation.r
    1

It can still say what the dropped equation was.  This is the :math:`H_j` from above,
and it mentions :math:`y_0`, :math:`y_1` and :math:`t` -- but no :math:`x`.

    >>> numeric_translation.auxiliary_growth_rates()
    [(2*y0*y1 - 1)/(3*t)]
    >>> numeric_translation.auxiliaries_are_constant
    False

Because such auxiliaries are not constant, they will have to be integrated.


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
    >>> y0, y1 = numeric_translation.invariant_variables
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

    >>> initial_conditions = numeric_translation.forward({t: 1.0, z1: 0.5, z2: 0.2})
    >>> [float(initial_conditions[y]) for y in (y0, y1)]
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

Now solve the two-equation system numerically.

    >>> def as_function(a_system):
    ...     variables = list(a_system.non_constant_variables)
    ...     rhs = [a_system.derivative_dict[x] for x in variables]
    ...     return sympy.lambdify([a_system.indep_var, variables], rhs, modules='numpy')

    >>> times = np.linspace(1.0, 2.0, 200)
    >>> solution = solve_ivp(as_function(reduced_system), (1.0, 2.0),
    ...                      [float(initial_conditions[y0]), float(initial_conditions[y1])],
    ...                      t_eval=times, rtol=1e-11, atol=1e-13, dense_output=True)


The numerical solution is very close to the exact:

    >>> bool(np.max(np.abs(solution.y[0] - y0_exact(times))) < 1e-9)
    True
    >>> bool(np.max(np.abs(solution.y[1] - y1_exact(times))) < 1e-9)
    True


Reverse translation
---------------------

Translating a series back to original variables needs the auxiliary at its first sample, from which its own
equation carries it forward.  For a system with constants, the known values that fix it
would be constants.  This system has none -- its variables are :math:`t`, :math:`z_1` and
:math:`z_2`, and that is all -- so the only thing we can know is a *point* value: where the
solution started.

A point value of a dependent variable is not something a series can be told directly.  The
series already says how :math:`z_1` varies, so a single number for it mixes two different
kinds of statement, and :meth:`~desr.numerics.NumericTranslation.reverse` says so.

    >>> reduced_values = {reduced_system.indep_var: solution.t,
    ...                   y0: solution.y[0], y1: solution.y[1]}
    >>> numeric_translation.reverse(reduced_values, known_values={z1: 0.5})
    Traceback (most recent call last):
        ...
    ValueError: z1 is a function of time, but a series is being translated, and the series already says how z1 varies.  A single value for it mixes point translation with series translation.  To use a value of it at one point, compute the auxiliaries there with auxiliaries_from() and pass them as auxiliaries=.

So reverse translation of a system like this is done in two steps.  

First, at the one point where :math:`z_1 = 0.5` is known -- the start, whose reduced values are the
scalars in ``initial_conditions`` -- determine the auxiliary:

    >>> at_start = numeric_translation.auxiliaries_from(initial_conditions, known_values={z1: 0.5})
    >>> [round(float(x), 6) for x in at_start]
    [0.25]

That value is :math:`t z_1^2` at :math:`t = 1`, as it should be.  Hand it to
:meth:`~desr.numerics.NumericTranslation.reverse`, which integrates it
along the solution.  (Passing the solver's dense output lets that integration evaluate the
invariants wherever it likes.)

    >>> recovered = numeric_translation.reverse(reduced_values, auxiliaries=at_start,
    ...                             invariants_at=solution.sol)

Hubert-Labahn reverse-translates the exact solution, and desr's symbolic
:meth:`~desr.ode_translation.ODETranslation.reverse_translate_general` does the same
computation: it rebuilds :math:`z_1` and :math:`z_2` from :math:`x`, :math:`y_1` and
:math:`y_2`, then rescales time by a constant that the arbitrary :math:`c_i` leave free --
here :math:`c_3^3 / c_1^2`.  With the constants fixed by our starting point that factor is
exactly :math:`1`.

    >>> sympy.simplify(constants[c3]**3 / constants[c1]**2)
    1
    >>> z1_exact, z2_exact = translation.reverse_translate_general(
    ...     (t, exact[x0], exact[y0], exact[y1]))
    >>> z1_exact = sympy.lambdify(t, z1_exact.subs(constants), modules='numpy')
    >>> z2_exact = sympy.lambdify(t, z2_exact.subs(constants), modules='numpy')
    >>> [round(float(f(1.0)), 12) for f in (z1_exact, z2_exact)]
    [0.5, 0.2]

Internally the integration is done on :math:`\log |x_j|`, where the equation
:math:`dx_j/dt = x_j H_j` is exactly linear.  The auxiliaries are monomials and can span a
wide range of magnitudes, which logarithms handle comfortably.  Their sign cannot change,
since :math:`x_j = 0` is invariant, so it is read once from the starting value and only the
magnitude is integrated.


Visualization
-------------------------

Plotting the numerical results:

.. plot:: _examples/quadrature_recovery.py

The left depicts the numerical solution of the nondimensionalized system, sans auxiliary variables.  
The middle panel depicts the auxiliary variable, obtained by quadrature.  This curve is what
:meth:`~desr.numerics.NumericTranslation.recover_auxiliaries` produces by integrating their
equations along the invariants.  The circles are this example's exact values.  

The right-hand panel is each original variable, obtained from reverse translation.


When there is nothing to integrate
----------------------------------

The quadrature above is needed only when the general reduction scheme was needed and auxiliaries were omitted when computing the reduced system.  Auxiliaries only appear when the
scaling acts on the dependent variables.  
