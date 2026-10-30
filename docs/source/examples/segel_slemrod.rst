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


Primary method
-----------------

"Outer" equations from Segel 1989
***********************************

Let's recover :cite:`Segel1989` "outer" equations (8a-b).
    
Start with a fresh system.

    >>> system = ODESystem.from_tex(system_tex)

Add initial conditions

    >>> system.update_initial_conditions({'s': 's_0'})

the Michaelis constant

    >>> system.add_constraint('K_m', '(k_2 + k_m1) / k_1')

Choose a variable order, with initial conditions at the front so the algorithm will normalize by them.

    >>> system.reorder_variables(['s_0', 'e_0', 'K_m', 'k_1', 'k_2', 'k_m1', 't', 's', 'c'])

    >>> translator = ODETranslator.from_ode_system(system, naming_scheme={'parameters': 'c', 'independent': 'H', 'dependents': ['y','z']})

The column numbers come from :meth:`~desr.ode_translation.ODETranslator.column_of`, so the
page says which invariant each operation acts on rather than carrying a bare index.

    >>> clock = translator.column_of('k_1*s_0*t')
    >>> scale = translator.column_of('e_0/s_0')
    >>> concentration = translator.column_of('c/s_0')
    >>> clock, scale, concentration
    (6, 2, 8)

    >>> translator.multiplier_add_columns(clock, scale, 1)  # Multiply time by e_0/s_0
    >>> translator.multiplier_add_columns(concentration, scale, -1)  # Scale c by e_0 rather than s_0

The invariants match :cite:`Segel1989` (9):

    >>> translator.invariant_monomials()
    Matrix([[e_0/s_0, K_m/s_0, k_2/(k_1*s_0), k_m1/(k_1*s_0), e_0*k_1*t, s/s_0, c/e_0]])

Here are the substitutions the software will perform:

    >>> translator.translate_parameter_substitutions(system)
    {s_0: 1, e_0: c0, K_m: c1, k_1: 1, k_2: c2, k_m1: c3, t: H/c0, s: y, c: c0*z}

Translate and obtain a system

    >>> reduced_system = translator.translate(system)

Do a few substitutions, and we have :cite:`Segel1989` (8a-b).  Note the spelling
:math:`\lambda` has to take: ``lambda`` is a reserved word in Python, so :mod:`sympy` cannot
parse it, and its own convention -- which :mod:`desr` follows -- is ``lamda``.  That spelling
typesets as :math:`\lambda`, so :func:`desr.tex_tools.expr_to_tex` renders it correctly.

    >>> reduced_system.subs({'c0': 'epsilon', 'c1': 'mu', 'c2': 'lamda', 'c3': 'mu-lamda'},
    ...                    subs_constraints=True,
    ...                    expand_after=True,
    ...                    factor_after=True)
    dH/dH = 1
    dy/dH = -lamda*z + mu*z + y*z - y
    dz/dH = -(mu*z + y*z - y)/epsilon
    y(0) = 1
    parameters: epsilon, lamda, mu






Reduction to Segel (21a-d)
***********************************


First create the original system.

    >>> system = ODESystem.from_tex(system_tex)

Add initial conditions

    >>> system.update_initial_conditions({'s': 's_0'})

the Michaelis constant

    >>> system.add_constraint('K_m', '(k_2 + k_m1) / k_1')

The constant :math:`\epsilon = \frac{e_0}{s_0 + K_m}`, :cite:`murray` (6.18), and :cite:`Segel1989` p. 451

    >>> system.add_constraint('epsilon', 'e_0 / (s_0 + K_m)')

Choose a variable order, with initial conditions at the front so the algorithm will normalize by them.

    >>> system.reorder_variables(['s_0', 'e_0', 'K_m', 'k_1', 'k_2', 'k_m1', 'epsilon', 't', 's', 'c'])

A sanity check

    >>> system.symbols
    (s_0, e_0, K_m, k_1, k_2, k_m1, epsilon, t, s, c)

Now we can construct the `ODETranslator`

    >>> translator = ODETranslator.from_ode_system(system, naming_scheme={'parameters': 'c', 'independent': 'tau', 'dependents': ['s','c']})
    >>> translator.scaling_matrix
    Matrix([
    [0, 0, 0, -1, -1, -1, 0, 1, 0, 0],
    [1, 1, 1, -1,  0,  0, 0, 0, 1, 1]])

and translate the system

    >>> translator.translate(system)
    dtau/dtau = 1
    dc/dtau = -c*c2 - c*c3 - c*s + c0*s
    ds/dtau = c*c3 + c*s - c0*s
    s(0) = 1
    c1 == c2 + c3
    c4 == c0/(c1 + 1)
    parameters: c0, c1, c2, c3, c4

The scaling invariants:

    >>> translator.invariant_monomials()
    Matrix([[e_0/s_0, K_m/s_0, k_2/(k_1*s_0), k_m1/(k_1*s_0), epsilon, k_1*s_0*t, s/s_0, c/s_0]])

This is not the desired form of the system, we need to do column operations on the Hermite multiplier of the translator.  

To get :math:`\tau = \frac{e_0 k_1 t}{\epsilon}`, modify Hermite multiplier as :math:`V_8' = V_8 + V_3 - V_7`



Note that indices in Python start at 0, but in math / the paper start at 1, so they're off-by-one.

Name each column once, by the invariant it holds, and the five operations below read as what
they do.  The numbers are printed because they are what the write-up quotes:

    >>> clock = translator.column_of('k_1*s_0*t')
    >>> scale = translator.column_of('e_0/s_0')
    >>> eps = translator.column_of('epsilon')
    >>> concentration = translator.column_of('c/s_0')
    >>> michaelis = translator.column_of('K_m/s_0')
    >>> k_m1_ratio = translator.column_of('k_m1/(k_1*s_0)')
    >>> k_2_ratio = translator.column_of('k_2/(k_1*s_0)')
    >>> clock, scale, eps, concentration, michaelis, k_m1_ratio, k_2_ratio
    (7, 2, 6, 9, 3, 5, 4)

    >>> # Scale t correctly to t/t_C = k_1 L t = e_0 k_1 t / epsilon
    >>> translator.multiplier_add_columns(clock, scale, 1)
    >>> translator.multiplier_add_columns(clock, eps, -1)

To get :math:`v = \frac{c}{s_0 \epsilon}`, we need to modify the Hermite multiplier as :math:`V_{10}' = V_{10} - V_7`.  In code,

    >>> # Scale s correctly to s / s_0
    >>> # Scale c correctly to c / (e_0 s_0 / L) = c / (s_0 epsilon)
    >>> translator.multiplier_add_columns(concentration, eps, -1)


We also want :math:`\sigma = \frac{s_0}{K_m}`, so :math:`V_4' = -V_4`.

    >>> # sigma = s_0 / K_m
    >>> translator.multiplier_negate_column(michaelis)


To get :math:`\kappa = \frac{k_{-1}}{k_2}` in :cite:`Segel1989` (17), or :math:`\rho = \frac{k_{-1}}{k_2}` in :cite:`murray` (6.20), we need :math:`V_6' = V_6 - V_5`

    >>> # Find kappa = rho = k_{-1} / k_2 = (K_m k_1 / k_2) - 1
    >>> translator.multiplier_add_columns(k_m1_ratio, k_2_ratio, -1)


We have the invariants up to this point of

    >>> translator.invariant_monomials()
    Matrix([[e_0/s_0, s_0/K_m, k_2/(k_1*s_0), k_m1/k_2, epsilon, e_0*k_1*t/epsilon, s/s_0, c/(epsilon*s_0)]])
    

Making some substitutions leads to :cite:`Segel1989` (21a-d)
    
    >>> reduced_system = translator.translate(system)
    >>> reduced_system = reduced_system.subs({'c0': 'epsilon * (1 + 1 / sigma)',
    ...                                          'c2': '1 / (sigma * (kappa + 1))'},subs_constraints=True,
    ...                    expand_after=True,
    ...                    factor_after=True)
    >>> reduced_system = reduced_system.subs({'c1': 'sigma',
    ...                          'c3': 'kappa',
    ...                          'c4': 'epsilon'},
    ...                         subs_constraints=True,
    ...                         expand_after=True,
    ...                         factor_after=True)
    >>> reduced_system
    dtau/dtau = 1
    dc/dtau = -(c*s*sigma + c - s*sigma - s)/(sigma + 1)
    ds/dtau = epsilon*(c*kappa*s*sigma + c*kappa + c*s*sigma - kappa*s*sigma - kappa*s - s*sigma - s)/((kappa + 1)*(sigma + 1))
    s(0) = 1
    1/sigma == kappa/(sigma*(kappa + 1)) + 1/(sigma*(kappa + 1))
    parameters: epsilon, kappa, sigma


To get :cite:`murray` (6.21), we need to rename some variables

    >>> reduced_system_murray = reduced_system.subs({'c': 'v', 's':'u'},
    ...                         subs_constraints=True,
    ...                         expand_after=True,
    ...                         factor_after=True)
    >>> reduced_system_murray
    dtau/dtau = 1
    dv/dtau = -(sigma*u*v - sigma*u - u + v)/(sigma + 1)
    du/dtau = epsilon*(kappa*sigma*u*v - kappa*sigma*u - kappa*u + kappa*v + sigma*u*v - sigma*u - u)/((kappa + 1)*(sigma + 1))
    u(0) = 1
    1/sigma == kappa/(sigma*(kappa + 1)) + 1/(sigma*(kappa + 1))
    parameters: epsilon, kappa, sigma


After the pre-steady state
***********************************

Computing equations :cite:`Segel1989` (24a-b).  Note that we start from a copy of the translator
above, rather than a fresh one: :cite:`Segel1989` (24a-b) uses the same state variables as (21a-d),
and only the clock changes.  The three column operations act on the time column alone.

    >>> # Scale t correctly to t/t_S = k_2 epsilon t
    >>> import copy
    >>> translator_24 = copy.deepcopy(translator)
    >>> # the clock is now e_0*k_1*t/epsilon, the other three columns are unchanged
    >>> clock_24 = translator_24.column_of('e_0*k_1*t/epsilon')
    >>> clock_24, k_2_ratio, scale, eps
    (7, 4, 2, 6)
    >>> translator_24.multiplier_add_columns(clock_24, k_2_ratio, 1)
    >>> translator_24.multiplier_add_columns(clock_24, scale, -1)
    >>> translator_24.multiplier_add_columns(clock_24, eps, 2)
    >>> translator_24.invariant_monomials()
    Matrix([[e_0/s_0, s_0/K_m, k_2/(k_1*s_0), k_m1/k_2, epsilon, epsilon*k_2*t, s/s_0, c/(epsilon*s_0)]])
    >>> reduced_system = translator_24.translate(system)
    >>> reduced_system = reduced_system.subs({'c0': 'epsilon * (1 + 1 / sigma)',
    ...                                          'c2': '1 / (sigma * (kappa + 1))'
    ...                                          },subs_constraints=True,
    ...                    expand_after=True,
    ...                    factor_after=True)
    >>> reduced_system.subs({'c': 'v',
    ...                          'c1': 'sigma',
    ...                          'c3': 'kappa',
    ...                          'c4': 'epsilon',
    ...                          's': 'u'},
    ...                         subs_constraints=True,
    ...                         expand_after=True,
    ...                         factor_after=True)
    dtau/dtau = 1
    dv/dtau = -(kappa + 1)*(sigma*u*v - sigma*u - u + v)/epsilon
    du/dtau = kappa*sigma*u*v - kappa*sigma*u - kappa*u + kappa*v + sigma*u*v - sigma*u - u
    u(0) = 1
    1/sigma == kappa/(sigma*(kappa + 1)) + 1/(sigma*(kappa + 1))
    parameters: epsilon, kappa, sigma




Another method: using :math:`L = s_0 + K_m`
-------------------------------------------------

Starting from the same system, Michaelis-Menten, with :math:`K_m = \frac{k_{-1} + k_2}{k_1}` and initial condition :math:`s(0) = s_0`.

    >>> # Substitute K_m into the equations
    >>> system_tex_l = system_tex.replace('k_{-1}', '(K - k_2)').replace('K', 'K_m k_1')
    >>> system_l = ODESystem.from_tex(system_tex_l)
    >>> system_l
    dt/dt = 1
    dc/dt = -c*k_1*s - c*k_2 - c*(K_m*k_1 - k_2) + e_0*k_1*s
    ds/dt = c*k_1*s + c*(K_m*k_1 - k_2) - e_0*k_1*s
    parameters: K_m, e_0, k_1, k_2
    >>> system_l.update_initial_conditions({'s': 's_0'})


Fast transient time scale
****************************

Let's first work in the fast transient time scale (:cite:`murray` eqn 6.15), we have :math:`t_c = \frac{1}{k_1 (s_0 + K_m)}`.  Our time variable is :math:`\tau = \frac{t}{t_c} = k_1 (s_0 + K_m) \, t`.

We can work with the expression :math:`s_0 + K_m` systematically by adding a new variable :math:`L = s_0 + K_m` via a constraint.

    >>> system_l.add_constraint('L', 's_0 + K_m')

Check that if we keep :math:`L` out of the normalising positions -- here :math:`s_0` and :math:`k_1`
are the two parameters that get normalised -- we have the same reduced system as before

    >>> system_l.reorder_variables(['s_0', 'L', 'e_0', 'k_1', 'k_2', 'K_m', 't', 's', 'c'])
    >>> translator = ODETranslator.from_ode_system(system_l, naming_scheme={'parameters': 'c', 'independent': 't', 'dependents': ['s','c']})
    >>> translator.scaling_matrix
    Matrix([
    [0, 0, 0, -1, -1, 0, 1, 0, 0],
    [1, 1, 1, -1,  0, 1, 0, 1, 1]])
    >>> translator.invariant_monomials()
    Matrix([[L/s_0, e_0/s_0, k_2/(k_1*s_0), K_m/s_0, k_1*s_0*t, s/s_0, c/s_0]])
    >>> translator.translate(system_l)
    dt/dt = 1
    dc/dt = -c*c3 - c*s + c1*s
    ds/dt = -c*c2 + c*c3 + c*s - c1*s
    s(0) = 1
    c0 == c3 + 1
    parameters: c0, c1, c2, c3

Putting :math:`K_m` at the front, we have

    >>> system_l.reorder_variables(['K_m', 'L', 's_0', 'e_0', 'k_1', 'k_2', 't', 's', 'c'])
    >>> translator = ODETranslator.from_ode_system(system_l, naming_scheme={'parameters': 'c', 'independent': 'tau', 'dependents': ['s','c']})
    >>> clock = translator.column_of('K_m*k_1*t')
    >>> total = translator.column_of('L/K_m')
    >>> substrate = translator.column_of('s/K_m')
    >>> sigma = translator.column_of('s_0/K_m')
    >>> complexed = translator.column_of('c/K_m')
    >>> enzyme = translator.column_of('e_0/K_m')
    >>> k_2_ratio = translator.column_of('k_2/(K_m*k_1)')
    >>> clock, total, substrate, sigma, complexed, enzyme, k_2_ratio
    (6, 2, 7, 3, 8, 4, 5)

    >>> # Scale t correctly to t/t_C = k_1 L t
    >>> translator.multiplier_add_columns(clock, total, 1)
    >>> # Scale s correctly to s / s_0
    >>> translator.multiplier_add_columns(substrate, sigma, -1)
    >>> # Scale c correctly to c / (e_0 s_0 / L)
    >>> translator.multiplier_add_columns(complexed, total, 1)
    >>> translator.multiplier_add_columns(complexed, enzyme, -1)
    >>> translator.multiplier_add_columns(complexed, sigma, -1)
    >>> # Find epsilon = e_0 / L
    >>> translator.multiplier_add_columns(enzyme, total, -1)
    >>> # Find kappa = k_{-1} / k_2 = (K_m k_1 / k_2) - 1
    >>> translator.multiplier_negate_column(k_2_ratio)
    >>> # sigma = s_0 / K_m needs no operation
    >>> translator.invariant_monomials()
    Matrix([[L/K_m, s_0/K_m, e_0/L, K_m*k_1/k_2, L*k_1*t, s/s_0, L*c/(e_0*s_0)]])


In the invariants, first come the scale-invariant coefficients, then the time variable, then the two
"space" variables.  The coefficients are:

.. math::
    :nowrap:

    \begin{align}
    c_0 &= \frac{L}{K_m} = \sigma + 1 \\
    c_1 &= \frac{s_0}{K_m} = \sigma \textnormal{ (in both Segel and Murray)} \\
    c_2 &= \frac{e_0}{L} = \frac{e_0}{s_0 + K_m} = \epsilon  \\
    c_3 &= \frac{K_m k_1}{k_2} = \kappa + 1 \textnormal{ (in Segel)} = \rho + 1 \textnormal{ (in Murray)}
    \end{align}

The substitutions that will be made:

    >>> translator.translate_parameter_substitutions(system_l)
    {K_m: 1, L: c0, s_0: c1, e_0: c0*c2, k_1: 1, k_2: 1/c3, t: tau/c0, s: c1*s, c: c*c1*c2}

Here's the reduced system at this point:

    >>> reduced_system = translator.translate(system_l)
    >>> reduced_system
    dtau/dtau = 1
    dc/dtau = -c*c1*s/c0 - c/c0 + s
    ds/dtau = c*c1*c2*s/c0 + c*c2/c0 - c*c2/(c0*c3) - c2*s
    s(0) = 1
    c0 == c1 + 1
    parameters: c0, c1, c2, c3

To get the system in (21a-d) Segel, let's do the substitutions into known symbol names

    >>> system_segel = reduced_system.subs({'c0':'sigma+1', 'c1':'sigma', 'c2':'epsilon', 'c3':'kappa+1'},expand_after=True, factor_after=False)
    >>> system_segel
    dtau/dtau = 1
    dc/dtau = -c*s*sigma/(sigma + 1) - c/(sigma + 1) + s
    ds/dtau = c*epsilon*s*sigma/(sigma + 1) - c*epsilon/(kappa*sigma + kappa + sigma + 1) + c*epsilon/(sigma + 1) - epsilon*s
    s(0) = 1
    parameters: epsilon, kappa, sigma

..
   _Silviana checked these on November 4, 2025, by hand in research notebook 10, page 75.

The equation for :math:`\frac{dc}{d \tau}` is exactly the same as in Segel (21b).  The formula for :math:`\frac{ds}{d \tau}` looks different from Segel (21a), but algebraically it's exactly the same.  It's just difficult to get Sympy to factor it into exactly the same form.

These are exactly the same equations as in Murray (6.21) upon replacement of :math:`(s,c) = (u,v)`, and :math:`\rho = \kappa`.






Long/slow time scale
****************************

This is the second route on this page to :cite:`Segel1989` (23) and (24a-b); the first was
"After the pre-steady state" above, which reached the same timescale from the :math:`K_m`
and :math:`\epsilon` constraint system.

To get the equations on the long/slow timescale :math:`t_c` from Murray, multiply :math:`\tau = Lk_1t` by the factors
:math:`\frac{e_0}{L} \frac{k_2}{K_m*k_1}\frac{K_m}{L} = \frac{c_2}{c_3c_0}`

    >>> translator.herm_mult.shape
    (9, 9)
    >>> # the clock is now L*k_1*t; epsilon and kappa+1 sit where e_0/K_m and k_2/(K_m*k_1) did
    >>> clock = translator.column_of('L*k_1*t')
    >>> epsilon = translator.column_of('e_0/L')
    >>> kappa_plus_one = translator.column_of('K_m*k_1/k_2')
    >>> clock, epsilon, kappa_plus_one, total
    (6, 4, 5, 2)
    >>> translator.multiplier_add_columns(clock, epsilon, 1)
    >>> translator.multiplier_add_columns(clock, kappa_plus_one, -1)
    >>> translator.multiplier_add_columns(clock, total, -1)

Inspect

    >>> translator.invariant_monomials()
    Matrix([[L/K_m, s_0/K_m, e_0/L, K_m*k_1/k_2, e_0*k_2*t/L, s/s_0, L*c/(e_0*s_0)]])

    >>> translator.translate_parameter_substitutions(system_l)
    {K_m: 1, L: c0, s_0: c1, e_0: c0*c2, k_1: 1, k_2: 1/c3, t: c3*tau/c2, s: c1*s, c: c*c1*c2}

Translate and reduce the system

    >>> reduced_system = translator.translate(system_l, naming_scheme={'parameters': 'c', 'independent': 'T', 'dependents': ['s','c']})
    >>> reduced_system
    dT/dT = 1
    dc/dT = -c*c1*c3*s/c2 - c*c3/c2 + c0*c3*s/c2
    ds/dT = c*c1*c3*s + c*c3 - c - c0*c3*s
    s(0) = 1
    c0 == c1 + 1
    parameters: c0, c1, c2, c3

These are not the symbols we are looking for 🤖.  Looking at the invariants suggests a path forward.  


First, the time and state variables:

.. math::
    :nowrap:

    \begin{align}
    e_0*k_2*t/L &= T & & \textnormal{ Segel (23)}\\
    s/s_0 &= \frac{S}{\bar{S_0}} & & \textnormal{ Segel (20a)}\\
    L*c/(e_0*s_0) &= \frac{C}{\bar{C}} & & \textnormal{ Segel (20b)}
    \end{align}

Then, the coefficients

.. math::
    :nowrap:

    \begin{align}
    c_0 &= L/K_m & &= \sigma+1 \\
    c_1 &= s_0/K_m & &= \sigma \\
    c_2 &= e_0/L & &= \epsilon \\
    c_3 &= K_m*k_1/k_2 & &=  \kappa+1
    \end{align}


we do a few substitutions, we get (6.23) from Murray.

    >>> system_murray = reduced_system.subs({#     ...                         'c':'v', 's':'u',
    ...                         'c0':'sigma+1', 'c1':'sigma',
    ...                         'c2':'epsilon', 'c3':'kappa+1'},
    ...                         subs_constraints=True,
    ...                         expand_after=True,
    ...                         factor_after=True)
    >>> system_murray
    dT/dT = 1
    dc/dT = -(kappa + 1)*(c*s*sigma + c - s*sigma - s)/epsilon
    ds/dT = c*kappa*s*sigma + c*kappa + c*s*sigma - kappa*s*sigma - kappa*s - s*sigma - s
    s(0) = 1
    parameters: epsilon, kappa, sigma


The equation for :math:`\frac{dc}{dT}` is exactly the same as Segel (24b), after some minor algebra.  
The equation for :math:`\frac{ds}{dT}` also requires a bit of basic algebra, and is equal to Segal (24a).

We recover Murray (6.23) with a bit more algebra, not shown here.

..
   _Silviana checked these on November 4, 2025, by hand in research notebook 10, page 77.


