Multiple timescales: Segel-Slemrod analysis
============================================


This is our third example using the two-variable formulation of the Michaelis-Menten system.  In this example, our goal is to matching Segel and Slemrod's analysis :cite:`Segel1989`, and in Chapter 6 of Murray :cite:`murray`.


The model
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
    \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c
    \end{align}

This is the starting point for our reduction.

    >>> system_tex = '''\frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\\\
    ... \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c \\\\'''


"Outer" equations from Segel 1989
-------------------------------------

Let's recover :cite:`Segel1989` "outer" equations (8a-b).
    
Start with a fresh system.

    >>> system = ODESystem.from_tex(system_tex)

Add initial conditions

    >>> system.update_initial_conditions({'s': 's_0'})

the Michaelis constant

    >>> system.add_constraint('K_m', '(k_2 + k_m1) / k_1')

Choose a variable order, with initial conditions at end so the algorithm will normalize by them.

    >>> system.reorder_variables(['t', 's', 'c', 'k_m1', 'k_2', 'k_1', 'K_m', 'e_0', 's_0'])

    >>> translation = ODETranslation.from_ode_system(system, renaming_scheme=('H',['y','z'],'c'))
    >>> translation.multiplier_add_columns(2, -1, 1)  # Multiply time by epsilon
    >>> translation.multiplier_add_columns(4, -1, -1)  # divide column 4 by epsilon

The invariants match :cite:`Segel1989` (9):

    >>> translation.invariants()
    Matrix([[e_0*k_1*t, s/s_0, c/e_0, k_m1/(k_1*s_0), k_2/(k_1*s_0), K_m/s_0, e_0/s_0]])

Here are the substitutions the software will perform:

    >>> translation.translate_parameter_substitutions(system)
    {t: H/c3, s: y, c: c3*z, k_m1: c0, k_2: c1, k_1: 1, K_m: c2, e_0: c3, s_0: 1}

Translate and obtain a system

    >>> reduced_system = translation.translate(system)

Do a few substitutions, and we have :cite:`Segel1989` (8a-b)

    >>> reduced_system.diff_subs({'c1': 'lam', 'c0': 'mu-lam', 'c3':'epsilon', 'c2':'mu'},
    ...                    subs_constraints=True,
    ...                    expand_after=True,
    ...                    factor_after=True)
    dH/dH = 1
    dy/dH = -lam*z + mu*z + y*z - y
    dz/dH = -(mu*z + y*z - y)/epsilon
    depsilon/dH = 0
    dlam/dH = 0
    dmu/dH = 0
    y(0) = 1






Reduction to Segel (21a-d)
------------------------------


First create the original system.

    >>> system = ODESystem.from_tex(system_tex)

Add initial conditions

    >>> system.update_initial_conditions({'s': 's_0'})

the Michaelis constant

    >>> system.add_constraint('K_m', '(k_2 + k_m1) / k_1')

The constant :math:`\epsilon = \frac{e_0}{s_0 + K_m}`, :cite:`murray` (6.18), and :cite:`Segel1989` p. 451

    >>> system.add_constraint('epsilon', 'e_0 / (s_0 + K_m)')

Choose a variable order, with initial conditions at end so the algorithm will normalize by them.

    >>> system.reorder_variables(['t', 's', 'c', 'epsilon', 'k_m1', 'k_2', 'k_1', 'K_m', 'e_0', 's_0'])

A sanity check

    >>> system.variables
    (t, s, c, epsilon, k_m1, k_2, k_1, K_m, e_0, s_0)

Now we can construct the `ODETranslation`

    >>> translation = ODETranslation.from_ode_system(system, renaming_scheme=('tau',['s','c'],'c'))
    >>> translation.scaling_matrix
    Matrix([
    [1, 0, 0, 0, -1, -1, -1, 0, 0, 0],
    [0, 1, 1, 0,  0,  0, -1, 1, 1, 1]])

and translate the system

    >>> translation.translate(system)
    dtau/dtau = 1
    dc/dtau = -c*c1 - c*c2 - c*s + c4*s
    ds/dtau = c*c1 + c*s - c4*s
    dc1/dtau = 0
    dc2/dtau = 0
    dc4/dtau = 0
    dc3/dtau = 0
    dc0/dtau = 0
    s(0) = 1
    c3 == c1 + c2
    c0 == c4/(c3 + 1)

The scaling invariants:

    >>> translation.invariants()
    Matrix([[k_1*s_0*t, s/s_0, c/s_0, epsilon, k_m1/(k_1*s_0), k_2/(k_1*s_0), K_m/s_0, e_0/s_0]])

This is not the desired form of the system, we need to do column operations on the Hermite multiplier of the translation.  

To get :math:`\tau = \frac{e_0 k_1 t}{\epsilon}`, modify Hermite multiplier as :math:`V_3' = V_3 + V_{10} - V_6`



Note that indices in Python start at 0, but in math / the paper start at 1, so they're off-by-one.

    >>> # Scale t correctly to t/t_C = k_1 L t = e_0 k_1 t / epsilon
    >>> translation.multiplier_add_columns(2, 5, -1)
    >>> translation.multiplier_add_columns(2, 9, 1)

To get :math:`v = \frac{c}{s_0 \epsilon}`, we need to modify the Hermite multiplier as :math:`V_5' = V_5 - V_6`.  In code,

    >>> # Scale s correctly to s / s_0
    >>> # Scale c correctly to c / (e_0 s_0 / L) = c / (s_0 epsilon)
    >>> translation.multiplier_add_columns(4, 5, -1)


We also want :math:`\sigma = \frac{s_0}{K_m}`, so :math:`V_9' = -V_9`.

    >>> # sigma = s_0 / K_m
    >>> translation.multiplier_negate_column(8)


To get :math:`\kappa = \frac{k_{-1}}{k_2}` in :cite:`Segel1989` (17), or :math:`\rho = \frac{k_{-1}}{k_2}` in :cite:`murray` (6.20), we need :math:`V_7' = V_7 - V_8`

    >>> # Find kappa = rho = k_{-1} / k_2 = (K_m k_1 / k_2) - 1
    >>> translation.multiplier_add_columns(6, 7, -1)


We have the invariants up to this point of

    >>> translation.invariants()
    Matrix([[e_0*k_1*t/epsilon, s/s_0, c/(epsilon*s_0), epsilon, k_m1/k_2, k_2/(k_1*s_0), s_0/K_m, e_0/s_0]])
    

Making some substitutions leads to :cite:`Segel1989` (21a-d)
    
    >>> reduced_system = translation.translate(system)
    >>> reduced_system = reduced_system.diff_subs({'c2': '1 / (sigma * (kappa + 1))',
    ...                                          'c4': 'epsilon * (1 + 1 / sigma)'},subs_constraints=True,
    ...                    expand_after=True,
    ...                    factor_after=True)
    >>> reduced_system = reduced_system.diff_subs({'c0': 'epsilon',
    ...                          'c3': 'sigma',
    ...                          'c1': 'kappa'},
    ...                         subs_constraints=True,
    ...                         expand_after=True,
    ...                         factor_after=True)
    >>> reduced_system
    dtau/dtau = 1
    dc/dtau = -(c*s*sigma + c - s*sigma - s)/(sigma + 1)
    ds/dtau = epsilon*(c*kappa*s*sigma + c*kappa + c*s*sigma - kappa*s*sigma - kappa*s - s*sigma - s)/((kappa + 1)*(sigma + 1))
    depsilon/dtau = 0
    dkappa/dtau = 0
    dsigma/dtau = 0
    s(0) = 1
    1/sigma == kappa/(sigma*(kappa + 1)) + 1/(sigma*(kappa + 1))


To get :cite:`murray` (6.21), we need to rename some variables

    >>> reduced_system_murray = reduced_system.diff_subs({'c': 'v', 's':'u'},
    ...                         subs_constraints=True,
    ...                         expand_after=True,
    ...                         factor_after=True)
    >>> reduced_system_murray
    dtau/dtau = 1
    dv/dtau = -(sigma*u*v - sigma*u - u + v)/(sigma + 1)
    du/dtau = epsilon*(kappa*sigma*u*v - kappa*sigma*u - kappa*u + kappa*v + sigma*u*v - sigma*u - u)/((kappa + 1)*(sigma + 1))
    depsilon/dtau = 0
    dkappa/dtau = 0
    dsigma/dtau = 0
    u(0) = 1
    1/sigma == kappa/(sigma*(kappa + 1)) + 1/(sigma*(kappa + 1))


After the pre-steady state
-------------------------------

Computing equations :cite:`Segel1989` (24a-b)

    >>> # Scale t correctly to t/t_S = k_2 epsilon t
    >>> import copy
    >>> translation_24 = copy.deepcopy(translation)
    >>> translation_24.multiplier_add_columns(2, 7, 1)
    >>> translation_24.multiplier_add_columns(2, -1, -1)
    >>> translation_24.multiplier_add_columns(2, 5, 2)
    >>> translation_24.invariants()
    Matrix([[epsilon*k_2*t, s/s_0, c/(epsilon*s_0), epsilon, k_m1/k_2, k_2/(k_1*s_0), s_0/K_m, e_0/s_0]])
    >>> reduced_system = translation_24.translate(system)
    >>> reduced_system = reduced_system.diff_subs({'c2': '1 / (sigma * (kappa + 1))',
    ...                                          'c4': 'epsilon * (1 + 1 / sigma)'
    ...                                          },subs_constraints=True,
    ...                    expand_after=True,
    ...                    factor_after=True)
    >>> reduced_system.diff_subs({'c0': 'epsilon',
    ...                          'c3': 'sigma',
    ...                          'c1': 'kappa',
    ...                          's': 'u', 'c': 'v'},
    ...                         subs_constraints=True,
    ...                         expand_after=True,
    ...                         factor_after=True)
    dtau/dtau = 1
    dv/dtau = -(kappa + 1)*(sigma*u*v - sigma*u - u + v)/epsilon
    du/dtau = kappa*sigma*u*v - kappa*sigma*u - kappa*u + kappa*v + sigma*u*v - sigma*u - u
    depsilon/dtau = 0
    dkappa/dtau = 0
    dsigma/dtau = 0
    u(0) = 1
    1/sigma == kappa/(sigma*(kappa + 1)) + 1/(sigma*(kappa + 1))




Another method
---------------------------------

We can work with the expression :math:`s_0 + K_m` systematically by adding a new variable :math:`L = s_0 + K_m` via a constraint.

    >>> # Substitute K_m into the equations
    >>> system_tex_l = system_tex.replace('k_{-1}', '(K - k_2)').replace('K', 'K_m k_1')
    >>> system_l = ODESystem.from_tex(system_tex_l)
    >>> system_l
    dt/dt = 1
    dc/dt = -c*k_1*s - c*k_2 - c*(K_m*k_1 - k_2) + e_0*k_1*s
    ds/dt = c*k_1*s + c*(K_m*k_1 - k_2) - e_0*k_1*s
    dK_m/dt = 0
    de_0/dt = 0
    dk_1/dt = 0
    dk_2/dt = 0
    >>> system_l.update_initial_conditions({'s': 's_0'})
    >>> system_l.add_constraint('L', 's_0 + K_m')

Check that if we keep :math:`L` at the end, we have the same reduced system as before

    >>> system_l.reorder_variables(['t', 's', 'c', 'K_m', 'k_2', 'k_1', 'e_0', 'L', 's_0'])
    >>> max_scal = ODETranslation.from_ode_system(system_l, renaming_scheme=('t',['s','c'], 'c'))
    >>> max_scal.scaling_matrix
    Matrix([
    [1, 0, 0, 0, -1, -1, 0, 0, 0],
    [0, 1, 1, 1,  0, -1, 1, 1, 1]])
    >>> max_scal.invariants()
    Matrix([[k_1*s_0*t, s/s_0, c/s_0, K_m/s_0, k_2/(k_1*s_0), e_0/s_0, L/s_0]])
    >>> max_scal.translate(system_l)
    dt/dt = 1
    dc/dt = -c*c0 - c*s + c2*s
    ds/dt = c*c0 - c*c1 + c*s - c2*s
    dc0/dt = 0
    dc1/dt = 0
    dc2/dt = 0
    dc3/dt = 0
    s(0) = 1
    c3 == c0 + 1

Now we put :math:`L` into the mix:

    >>> system_l.reorder_variables(['t', 's', 'c', 'k_2', 'k_1', 'e_0', 's_0', 'L', 'K_m'])
    >>> translation = ODETranslation.from_ode_system(system_l, renaming_scheme=('t',['s','c'], 'c'))
    >>> # Scale t correctly to t/t_C = k_1 L t
    >>> translation.multiplier_add_columns(2, -1, 1)
    >>> # Scale s correctly to s / s_0
    >>> translation.multiplier_add_columns(3, -2, -1)
    >>> # Scale c correctly to c / (e_0 s_0 / L)
    >>> translation.multiplier_add_columns(4, 6, -1)
    >>> translation.multiplier_add_columns(4, 7, -1)
    >>> translation.multiplier_add_columns(4, -1, 1)
    >>> # Find kappa = k_{-1} / k_2 = (K_m k_1 / k_2) - 1
    >>> translation.multiplier_negate_column(5)
    >>> # Find epsilon = e_0 / L
    >>> translation.multiplier_add_columns(6, -1, -1)
    >>> # Find sigma = s_0 / K_m
    >>> translation.invariants()
    Matrix([[L*k_1*t, s/s_0, L*c/(e_0*s_0), K_m*k_1/k_2, e_0/L, s_0/K_m, L/K_m]])

We now have:

.. math::
    :nowrap:

    \begin{align}
    c_0 &= \kappa + 1  = K_m k_1/k_2 \\
    c_1 &= \epsilon = \frac{e_0}{L} \\
    c_2 &= \sigma  = \frac{s_0}{K_m}\\
    c_3 &= \frac{L}{K_m}
    \end{align}


Which, after some trivial rearrangement, gives us exactly equations (24) from Segel and (6.21) from Murray:

    >>> translation.translate(system_l)
    dt/dt = 1
    dc/dt = -c*c2*s/c3 - c/c3 + s
    ds/dt = c*c1*c2*s/c3 + c*c1/c3 - c*c1/(c0*c3) - c1*s
    dc0/dt = 0
    dc1/dt = 0
    dc2/dt = 0
    dc3/dt = 0
    s(0) = 1
    c3 == c2 + 1

To get the equations on the long or slow timescale :math:`t_c` from Murray, multiply :math:`\tau = Lk_1t` by the factors
:math:`\frac{e_0}{L} \frac{k_2}{K_m*k_1}\frac{K_m}{L} = \frac{c_1}{c_0c_3}`

    >>> translation.multiplier_add_columns(2, 6, 1)
    >>> translation.multiplier_add_columns(2, 5, -1)
    >>> translation.multiplier_add_columns(2, -1, -1)
    >>> translation.invariants()
    Matrix([[e_0*k_2*t/L, s/s_0, L*c/(e_0*s_0), K_m*k_1/k_2, e_0/L, s_0/K_m, L/K_m]])
    >>> reduced_system = translation.translate(system_l)
    >>> reduced_system
    dt/dt = 1
    dc/dt = -c*c0*c2*s/c1 - c*c0/c1 + c0*c3*s/c1
    ds/dt = c*c0*c2*s + c*c0 - c - c0*c3*s
    dc0/dt = 0
    dc1/dt = 0
    dc2/dt = 0
    dc3/dt = 0
    s(0) = 1
    c3 == c2 + 1


If we do a few substitutions, we get (6.23) from Murray.

    >>> reduced_system.diff_subs({'c':'v', 's':'u', 'c2':'sigma', 'c1':'epsilon', 'c3':'sigma+1', 'c0':'something'},
    ...                         subs_constraints=True,
    ...                         expand_after=True,
    ...                         factor_after=True)



