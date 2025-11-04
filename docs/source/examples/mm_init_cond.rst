


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


    >>> system.power_matrix()
    Matrix([
    [1, 1, 1,  1, 1, 1,  1,  0,  0,  0],
    [0, 0, 0, -1, 0, 1,  1,  0,  0,  0],
    [1, 0, 0,  1, 0, 0, -1,  1,  0,  0],
    [0, 0, 0,  1, 1, 0,  0,  0,  0,  0],
    [1, 0, 0,  1, 1, 1,  0,  0,  1,  0],
    [0, 1, 0,  0, 0, 0,  0,  0,  0,  1],
    [0, 0, 1,  0, 0, 0,  1,  0, -1, -1],
    [0, 0, 0,  0, 0, 0,  0,  0,  1,  0],
    [0, 0, 0,  0, 0, 0,  0, -1,  0,  0]])




Choose the variable order.

    >>> system.reorder_variables(['t', 's', 'c', 'K_m', 'k_m1', 'k_2', 'k_1', 'e_0', 's_0'])

Compute the reduced system by making an `ODETranslation` from the original system,

    >>> translation = ODETranslation.from_ode_system(system, renaming_scheme=('t',['s','c'], 'c'))

Optionally, observe some properties.  First, the scaling matrix :math:`A`:

    >>> translation.scaling_matrix
    Matrix([
    [1, 0, 0, 0, -1, -1, -1, 0, 0],
    [0, 1, 1, 1,  0,  0, -1, 1, 1]])

The corresponding scaling invariants:

    >>> translation.invariants()
    Matrix([[k_1*s_0*t, s/s_0, c/s_0, K_m/s_0, k_m1/(k_1*s_0), k_2/(k_1*s_0), e_0/s_0]])

These don't match those from Murray, but they are very close.   We'll recover it next.


Just for curiosity, here's what we get if we translate the system at this point:

    >>> translation.translate(system)
    dt/dt = 1
    dc/dt = -c*c1 - c*c2 - c*s + c3*s
    ds/dt = c*c1 + c*s - c3*s
    dc1/dt = 0
    dc2/dt = 0
    dc3/dt = 0
    dc0/dt = 0
    s(0) = 1
    c0 == c1 + c2

Column operations to match the canonical system
*************************************************

Returning back to the original system with initial conditions,

    >>> system
    dt/dt = 1
    ds/dt = c*k_1*s + c*k_m1 - e_0*k_1*s
    dc/dt = -c*k_1*s - c*k_2 - c*k_m1 + e_0*k_1*s
    dK_m/dt = 0
    dk_m1/dt = 0
    dk_2/dt = 0
    dk_1/dt = 0
    de_0/dt = 0
    ds_0/dt = 0
    s(0) = s_0
    K_m == (k_2 + k_m1)/k_1
    >>> translation = ODETranslation.from_ode_system(system, renaming_scheme=('tau',['u','v'], 'c'))

Do some column operations on the Hermite multiplier

    >>> translation.multiplier_add_columns(2, 8, 1)  # Scale time by e_0 not s_0
    >>> translation.multiplier_add_columns(4, 8, -1)  # Scale c by e_0

Optionally, see the invariants and substitutions

    >>> translation.invariants()
    Matrix([[e_0*k_1*t, s/s_0, c/e_0, K_m/s_0, k_m1/(k_1*s_0), k_2/(k_1*s_0), e_0/s_0]])

We have substitutions to make.  Note that many get normalized to :math:`1`.

    >>> translation.translate_parameter_substitutions(system)
    {t: tau/c3, s: u, c: c3*v, K_m: c0, k_m1: c1, k_2: c2, k_1: 1, e_0: c3, s_0: 1}

Translate

    >>> reduced_system = translation.translate(system)
    >>> reduced_system
    dtau/dtau = 1
    du/dtau = c1*v + u*v - u
    dv/dtau = -c1*v/c3 - c2*v/c3 - u*v/c3 + u/c3
    dc1/dtau = 0
    dc2/dtau = 0
    dc3/dtau = 0
    dc0/dtau = 0
    u(0) = 1
    c0 == c1 + c2

Now we have a system, but the names of symbols don't match.  We'll do substitutions.
    
    >>> to_sub =  {'c0': 'K', 'c1': 'K-L', 'c2':'L', 'c3':'epsilon'}
    >>> reduced_system.diff_subs(to_sub, subs_constraints=True)
    dtau/dtau = 1
    du/dtau = K*v - L*v + u*v - u
    dv/dtau = -K*v/epsilon - u*v/epsilon + u/epsilon
    dK/dtau = 0
    dL/dtau = 0
    depsilon/dtau = 0
    u(0) = 1

The only step that `desr` cannot do is multiply :math:`\frac{dv}{d\tau}` by :math:`\epsilon`.  That has to be done by hand.  And finally we have the result from :cite:`murray` (6.13).

.. math::
    :nowrap:

    \begin{align*}
    \frac{du}{d\tau} & = -u + (u+K-\lambda)v & u(0)&=1 \\
    \epsilon \frac{dv}{d\tau} &= u - (u+K)v    &  v(0)&=0.
    \end{align*}


