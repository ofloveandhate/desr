


Michaelis-Menten with initial conditions and a constant
=============================================================

We derive the analysis of the Michaelis-Menten equations found in :cite:`Segel1989`: in a systematic manner.



Michaelis-Menten model, two variables
-----------------------------------------

There are two chemical species :math:`E` and :math:`S`, which reversibly combine to form :math:`ES`, which decomposes into :math:`E` and :math:`P`.

.. math::
    :nowrap:

    \[
    E+S\leftrightarrow ES \rightarrow E+P
    \]

After using the law of mass action, we are able to translate to the ODE system

.. math::
    :nowrap:

    \begin{align}
    \frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\
    \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c \\
    s(0) &= s_0.
    \end{align}



Model reduction in `desr`
-----------------------------------

First we must create the system in Python.  We start with the same system, represented in a string in TeX.

    >>> system_tex = '''\frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\\\
    ...          \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c'''
    >>> system = ODESystem.from_tex(system_tex)

Add the constraint with the Michaelis-Menten constant :math:`K_m = k_2 + k_{-1}`.
    
    >>> system.add_constraint('K_m', '(k_m1 + k_2) / k_1')

To match the analysis in Murray, We need to add in our initial condition for :math:`s`.

    >>> system.update_initial_conditions({'s': 's_0'})


Choose the variable order.

    >>> system.reorder_variables(['s_0', 'e_0', 'k_1', 'k_2', 'k_m1', 'K_m', 't', 's', 'c'])

Now the exponent matrix can be read against that order: one row per variable, in the order
just chosen, and one column per monomial appearing in the system.

    >>> system.exponent_matrix()
    Matrix([
    [0, 0,  0, 0, 0, 0,  0, -1,  0,  0],
    [1, 0,  0, 0, 0, 0,  1,  0,  0,  0],
    [1, 1,  0, 1, 0, 0,  1,  0,  1,  0],
    [0, 0,  0, 0, 1, 0,  0,  0,  0,  1],
    [0, 0,  1, 0, 0, 1,  0,  0, -1, -1],
    [0, 0,  0, 0, 0, 0,  0,  0,  1,  0],
    [1, 1,  1, 1, 1, 1,  1,  0,  0,  0],
    [0, 0, -1, 1, 0, 0,  1,  1,  0,  0],
    [0, 1,  1, 0, 0, 0, -1,  0,  0,  0]])

Compute the reduced system by making an `ODETranslator` from the original system,

    >>> translator = ODETranslator.from_ode_system(system, naming_scheme={'parameters': 'c', 'independent': 't', 'dependents': ['s','c']})

Optionally, observe some properties.  First, the scaling matrix :math:`A`:

    >>> translator.scaling_matrix
    Matrix([
    [0, 0, -1, -1, -1, 0, 1, 0, 0],
    [1, 1, -1,  0,  0, 1, 0, 1, 1]])

The corresponding scaling invariants:

    >>> translator.invariant_monomials()
    Matrix([[e_0/s_0, k_2/(k_1*s_0), k_m1/(k_1*s_0), K_m/s_0, k_1*s_0*t, s/s_0, c/s_0]])

These don't match those from Murray, but they are very close.   We'll recover it next.


Just for curiosity, here's what we get if we translate the system at this point:

    >>> translator.translate(system)
    dt/dt = 1
    dc/dt = -c*c1 - c*c2 - c*s + c0*s
    ds/dt = c*c2 + c*s - c0*s
    s(0) = 1
    c3 == c1 + c2
    parameters: c0, c1, c2, c3

Column operations to match the canonical system
*************************************************

Returning back to the original system with initial conditions,

    >>> system
    dt/dt = 1
    ds/dt = c*k_1*s + c*k_m1 - e_0*k_1*s
    dc/dt = -c*k_1*s - c*k_2 - c*k_m1 + e_0*k_1*s
    s(0) = s_0
    K_m == (k_2 + k_m1)/k_1
    parameters: s_0, e_0, k_1, k_2, k_m1, K_m
    >>> translator = ODETranslator.from_ode_system(system, naming_scheme={'parameters': 'c', 'independent': 'tau', 'dependents': ['u','v']})

Do some column operations on the Hermite multiplier

Rather than writing the column numbers in directly, ask which column holds which invariant.
The numbers are the same, and still worth knowing -- but now they are derived on the page
instead of appearing from nowhere.

    >>> clock = translator.column_of('k_1*s_0*t')
    >>> scale = translator.column_of('e_0/s_0')
    >>> concentration = translator.column_of('c/s_0')
    >>> clock, scale, concentration
    (6, 2, 8)

    >>> translator.multiplier_add_columns(clock, scale, 1)  # Scale time by e_0 not s_0
    >>> translator.multiplier_add_columns(concentration, scale, -1)  # Scale c by e_0

Optionally, see the invariants and substitutions

    >>> translator.invariant_monomials()
    Matrix([[e_0/s_0, k_2/(k_1*s_0), k_m1/(k_1*s_0), K_m/s_0, e_0*k_1*t, s/s_0, c/e_0]])

We have substitutions to make.  Note that many get normalized to :math:`1`.

    >>> translator.translate_parameter_substitutions(system)
    {s_0: 1, e_0: c0, k_1: 1, k_2: c1, k_m1: c2, K_m: c3, t: tau/c0, s: u, c: c0*v}

Translate

    >>> reduced_system = translator.translate(system)
    >>> reduced_system
    dtau/dtau = 1
    du/dtau = c2*v + u*v - u
    dv/dtau = -c1*v/c0 - c2*v/c0 - u*v/c0 + u/c0
    u(0) = 1
    c3 == c1 + c2
    parameters: c0, c1, c2, c3

Now we have a system, but the names of symbols don't match.  We'll do substitutions.
    
    >>> to_sub =  {'c0': 'epsilon', 'c1': 'L', 'c2':'K-L', 'c3':'K'}
    >>> reduced_system.subs(to_sub, subs_constraints=True)
    dtau/dtau = 1
    du/dtau = K*v - L*v + u*v - u
    dv/dtau = -K*v/epsilon - u*v/epsilon + u/epsilon
    u(0) = 1
    parameters: K, L, epsilon

The only step that `desr` cannot do is multiply :math:`\frac{dv}{d\tau}` by :math:`\epsilon`.  That has to be done by hand.  And finally we have the result from :cite:`murray` (6.13).

.. math::
    :nowrap:

    \begin{align*}
    \frac{du}{d\tau} & = -u + (u+K-\lambda)v & u(0)&=1 \\
    \epsilon \frac{dv}{d\tau} &= u - (u+K)v    &  v(0)&=0.
    \end{align*}


