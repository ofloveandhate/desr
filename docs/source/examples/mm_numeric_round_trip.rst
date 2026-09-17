Un-translating a numerical solution
========================================================

The nondimensionalization algorithm implemented in ``desr`` reduces the number of parameters or variables in a system, making systems simpler or perhaps easier to reason about or simulate.  This example discusses how to think about numerical solutions (or experimental data!) in the context of a reduced model.


Having solved a nondimensionalized system, we might want the answer back in the original units.  Let's walk
through that round trip using the Michaelis-Menten equations: we solve the *reduced* system
numerically, and recover the solution of the *original* one, then check the two agree.

:mod:`desr` provides some tools in :mod:`desr.numerics`, which needs :mod:`numpy`.  Due to this dependence, is not imported by
``desr``.  The numerics submodule optional, so users who only work symbolically are free of this dependency.  This example 
uses :mod:`scipy` to do the integrating, but a user is of course free to use any solver, or provide numeric data from any source.

The system and its reduction
----------------------------

For self-containedness we include this here.  Take the two-variable Michaelis-Menten system with an initial condition on :math:`s` (and without imposing the constant :math:`K_m`).

.. math::
    :nowrap:

    \begin{align}
    \frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\
    \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c \\
    s(0) &= s_0
    \end{align}

We build the system in Python from its TeX, and add the initial condition. 
The variable order matters: the reduction algorithm tends to reduce away the variables that come last.

    >>> import numpy as np
    >>> from scipy.integrate import solve_ivp
    >>> from desr.numerics import NumericTranslation

    >>> system_tex = '''\frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\\\
    ...          \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c'''
    >>> system = ODESystem.from_tex(system_tex)
    >>> system.update_initial_conditions({'s': 's_0'})
    >>> system.reorder_variables(['t', 's', 'c', 'k_m1', 'k_2', 'k_1', 'e_0', 's_0'])

Reduce it, naming the new variables :math:`\tau, u, v` and the new constants :math:`c_i`.

    >>> translation = ODETranslation.from_ode_system(system,
    ...                                              naming_scheme=('tau', ['u', 'v'], 'c'))
    >>> translation.invariants()
    Matrix([[k_1*s_0*t, s/s_0, c/s_0, k_m1/(k_1*s_0), k_2/(k_1*s_0), e_0/s_0]])
    >>> reduced_system = translation.translate(system)
    >>> reduced_system
    dtau/dtau = 1
    du/dtau = c0*v - c2*u + u*v
    dv/dtau = -c0*v - c1*v + c2*u - u*v
    dc0/dtau = 0
    dc1/dtau = 0
    dc2/dtau = 0
    u(0) = 1

Five constants have become three -- revealing that there are :math:`r = 2` scaling symmetries in this system:

    >>> translation.r
    2

The number of scaling symmetries is also the number of parameters we must provide later, when reverse translating from nondimensionalized to the original symbols. 

For programming convenience, let's give the variables names we can use below.

    >>> t, s, c, k_m1, k_2, k_1, e_0, s_0 = system.variables
    >>> tau, u, v, c0, c1, c2 = reduced_system.variables


Numerically solve the original system
-------------------------------------------

Choose some arbitrary parameter values and solve directly.  This is the answer we want to recover later via reverse translation.

    >>> parameters = {k_1: 1.5, k_m1: 0.9, k_2: 0.4, e_0: 0.3, s_0: 2.0}
    >>> def as_function(a_system, constants):
    ...     variables = list(a_system.non_constant_variables)
    ...     rhs = [a_system.derivative_dict[x].subs(constants) for x in variables]
    ...     return sympy.lambdify([a_system.indep_var, variables], rhs, modules='numpy')

    >>> times = np.linspace(0.0, 8.0, 9)
    >>> reference_soln = solve_ivp(as_function(system, parameters), (0.0, 8.0), [2.0, 0.0],
    ...                       t_eval=times, rtol=1e-10, atol=1e-12)


Translate parameter values
---------------------------

Let's translate the values of our original parameters, into the dimensionless parameters.  To this end,  
:class:`~desr.numerics.NumericTranslation` translates numeric values for parameters and variables across the same reduction that
:class:`~desr.ode_translation.ODETranslation` applies to symbols.

    >>> numeric_translation = NumericTranslation(system, translation, reduced_system)

The :meth:`~desr.numerics.NumericTranslation.forward`  function takes values of the original variables and
returns values of the reduced ones.  Here it converts the starting state and the parameters.

    >>> translated_values = numeric_translation.forward({t: 0.0, s: 2.0, c: 0.0, **parameters})
    >>> [float(translated_values[x]) for x in (tau, u, v)]
    [0.0, 1.0, 0.0]
    >>> [float(translated_values[x]) for x in (c0, c1, c2)]
    [0.3, 0.13333333333333333, 0.15]

You don't have to translate all the values at once, the function accepts a partial dictionary. For example, passing
a dict containing only parameters and time converts the time axis, since :math:`\tau = k_1 s_0 t`.

    >>> reduced_times = numeric_translation.forward({t: times, **parameters})[tau]
    >>> [float(x) for x in reduced_times[:3]]
    [0.0, 3.0, 6.0]

Do note, however, that if a dict does NOT contain enough values to do the reduction, then an error occurs.  So if I only supply the value of ``k_1``, then the reduction is missing both ``s_0`` and ``t``:
    
    >>> numeric_translation.forward({t: times, k_1: 1.5})
    Traceback (most recent call last):
    ValueError: The values given determine nothing.  Supply values for more of: t, s, c, k_m1, k_2, k_1, e_0, s_0.


Numerically solve the reduced system
-------------------------------------

Now that we have translated the parameters we're using, let's solve the reduced system numerically.

    >>> reduced_parameters = {c0: translated_values[c0], c1: translated_values[c1], c2: translated_values[c2]}
    >>> reduced_solution = solve_ivp(as_function(reduced_system, reduced_parameters),
    ...                      (reduced_times[0], reduced_times[-1]),
    ...                      [float(translated_values[u]), float(translated_values[v])],
    ...                      t_eval=reduced_times, rtol=1e-10, atol=1e-12)


Solving-then-reducing is equivalent to reducing-then-solving
-------------------------------------------------------------

It is worth checking the two routes agree.  There are two ways to get
from the original system to its solution, and they should give the same answer.  That is, the following diagram should commute.

.. image:: ../figures/commuting_square.svg
   :alt: A commuting square: forward and reverse across the top and bottom, solve down each side.
   :align: center
   :width: 100%

Going right-then-down is what we just did: reduce, then solve.  Going down-then-right is the
other route: solve the original system, then translate its solution.  We have the original
solution already, in ``reference_soln``.  Since the invariants are :math:`\tau = k_1 s_0 t`,
:math:`u = s / s_0` and :math:`v = c / s_0`, translating it is exactly what
:meth:`~desr.numerics.NumericTranslation.forward` does.  It takes the whole trajectory
in a single call.

    >>> translated_solution = numeric_translation.forward({
    ...                               t: times, 
    ...                               s: reference_soln.y[0], 
    ...                               c: reference_soln.y[1],
    ...                               **parameters})

The reduced solution we computed above and the translated original solution are the same
curve, on the same time axis.  Here we verify that the numerical error between the two is small.

    >>> bool(np.max(np.abs(translated_solution[tau] - reduced_solution.t)) < 1e-10)
    True
    >>> bool(np.max(np.abs(translated_solution[u] - reduced_solution.y[0])) < 1e-8)
    True
    >>> bool(np.max(np.abs(translated_solution[v] - reduced_solution.y[1])) < 1e-8)
    True

So the square commutes: reducing and solving may be done in either order.  That is the
property that makes the reduced system worth solving at all, and the correspondence is pretty much exact.


Reverse translation
---------------------------------------

The reduced system cannot tell us which original system it came from.  Remember, :math:`r = 2`, so the reduced system's parameters and solutions are shared by an
entire two-parameter family, all of which fit equally well.  So
:meth:`~desr.numerics.NumericTranslation.reverse` asks for
:attr:`~desr.numerics.NumericTranslation.r` values that you already know.  We'll use the rate :math:`k_1` and the initial substrate :math:`s_0`, both of which an experimenter would probably have.

Values may be arrays, so the whole time series converts in one call.

    >>> reduced_soln_as_dict = {tau: reduced_solution.t, u: reduced_solution.y[0], v: reduced_solution.y[1],
    ...           c0: translated_values[c0], 
    ...           c1: translated_values[c1], 
    ...           c2: translated_values[c2]}
    >>> recovered_soln = numeric_translation.reverse(reduced_soln_as_dict, known_values={k_1: 1.5, s_0: 2.0})

We get a dict back.  The independent variable comes back as the original :math:`t`.

    >>> bool(np.max(np.abs(recovered_soln[t] - times)) < 1e-10)
    True

And the trajectories agree with the direct solve.

    >>> bool(np.max(np.abs(recovered_soln[s] - reference_soln.y[0])) < 1e-8)
    True
    >>> bool(np.max(np.abs(recovered_soln[c] - reference_soln.y[1])) < 1e-8)
    True

Neat fact: Constants come back reverse-translated, too, which is what a fitting workflow might be after.  So if you fit
:math:`c_0, c_1, c_2` in the reduced system, then un-translate, 
you read off the original rates.

    >>> [round(float(recovered_soln[x]), 10) for x in (k_m1, k_2, e_0)]
    [0.9, 0.4, 0.3]


Visualization
-------------------------

The reduced system on the left and the original on the right are the same
curves on rescaled axes, which is what the commuting square above says they must be.  The
circles are the recovered solution sitting on the directly computed one.

.. plot:: _examples/michaelis_menten_numeric.py

The right-hand panel is worth dwelling on.  The disagreement between the two routes
never exceeds :math:`10^{-11}`, two orders of magnitude below the tolerance the integrator
was asked for.  Each original
variable is a product of integer powers of the reduced ones -- so translation in either direction introduces practically no error of
its own.  Thus, the right hand panel is the solver's noise, not that of the translation.


Without enough known parameter values, reverse-translation is impossible
-------------------------------------------------------------------------

Asking for the original system without enough known parameter values supplied raises an exception.
The message says what is still undetermined -- the rescaling that the values given do not
rule out -- and which variables would settle it.

    >>> numeric_translation.reverse(reduced_soln_as_dict, known_values={k_1: 1.5})
    Traceback (most recent call last):
        ...
    desr.numerics.InsufficientKnownValues: k_1 does not determine the original system: 1 more value is needed.
    The rescaling s, c, k_m1, k_2, e_0, s_0 by lambda and t by 1/lambda leaves k_1 unchanged.
    The reduced system is shared by an entire 2-parameter family of original systems, and this does not choose between them.
    Supply the value of one of: s, k_m1, k_2, e_0, s_0.

Even though :math:`r = 2`, two values are not always sufficient.  For example, knowing :math:`k_{-1}` and :math:`k_2` tells us
only about :math:`c_0` and :math:`c_1`, which between them pin down the single combination
:math:`k_1 s_0` rather than both factors.  The message names the rescaling that both values
are blind to, and offers only variables that would actually break it -- :math:`c` is not
among them, since :math:`c(0) = 0` and a zero cannot fix a scale.

    >>> numeric_translation.reverse(reduced_soln_as_dict, known_values={k_m1: 0.9, k_2: 0.4})
    Traceback (most recent call last):
        ...
    desr.numerics.InsufficientKnownValues: k_m1, k_2 do not determine the original system: 1 more value is needed.
    The rescaling s, c, e_0, s_0 by lambda and k_1 by 1/lambda leaves k_m1, k_2 unchanged.
    The reduced system is shared by an entire 2-parameter family of original systems, and this does not choose between them.
    Supply the value of one of: s, k_1, e_0, s_0.

Only *constants* may be known values when a series is translated.  :math:`s` is a function
of time, and the series -- :math:`u(\tau)`, which is :math:`s/s_0` -- already says how it
varies; a single number for it is a different kind of statement, and mixing the two is
refused rather than resolved by some convention.  The constant that stands for its initial
value is :math:`s_0`, and the message says so.

    >>> numeric_translation.reverse(reduced_soln_as_dict, known_values={k_1: 1.5, s: 2.0})
    Traceback (most recent call last):
        ...
    ValueError: s is a function of time, but a series is being translated, and the series already says how s varies.  A single value for it mixes point translation with series translation.  If you know its initial value, supply s_0, which stands for s(0).

At a single point, by contrast, :math:`s` *is* a number, and may be known.  Here the point
is the start, where :math:`u = 1` and so :math:`s = s_0`:

    >>> at_start = {x: reduced_soln_as_dict[x][0] for x in (tau, u, v)}
    >>> at_start.update({x: reduced_soln_as_dict[x] for x in (c0, c1, c2)})
    >>> point = numeric_translation.reverse(at_start, known_values={k_1: 1.5, s: 2.0})
    >>> [round(float(point[x]), 10) for x in (s_0, k_m1, k_2, e_0)]
    [2.0, 0.9, 0.4, 0.3]

Values that contradict each other are an error too, rather than a silent choice of one of
them.

    >>> numeric_translation.reverse(reduced_soln_as_dict,
    ...                             known_values={k_1: 1.5, s_0: 2.0, k_2: 99.0})
    Traceback (most recent call last):
        ...
    desr.numerics.ConflictingKnownValues: The known values disagree with each other given the reduced values.
        k_2 was given as 99, but the others imply 0.4.
    Check the values, or supply fewer of them: 2 independent values determine the original system.

More values than needed are fine, so long as they agree.

    >>> extra = numeric_translation.reverse(reduced_soln_as_dict,
    ...                                     known_values={k_1: 1.5, s_0: 2.0, k_2: 0.4})
    >>> [round(float(extra[x]), 10) for x in (k_m1, k_2, e_0)]
    [0.9, 0.4, 0.3]

