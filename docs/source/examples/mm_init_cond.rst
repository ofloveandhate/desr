


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
    s(0) &= s_0, & c(0) &= 0.
    \end{align}

Both initial conditions are part of the model: the complex starts at nothing.



Model reduction in `desr`
-----------------------------------

First we must create the system in Python.  We start with the same system, represented in a string in TeX.

    >>> system_tex = '''\frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\\\
    ...          \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c'''
    >>> system = ODESystem.from_tex(system_tex)

Add the constraint with the Michaelis-Menten constant :math:`K_m = k_2 + k_{-1}`.
    
    >>> system.add_constraint('K_m', '(k_m1 + k_2) / k_1')

To match the analysis in Murray, we need our initial conditions.  Both of them: the
substrate starts at :math:`s_0`, and the complex starts at nothing.

    >>> system.update_initial_conditions(['s(0) = s_0', 'c(0) = 0'])
    >>> system.initial_conditions
    {s: s_0, c: c_0}

A number still carries units in a system that has not been reduced -- a zero concentration
is a concentration -- so it is recorded as the *value* of a symbol, and the symbol is what
carries the units into the scaling.  The name is inferred from the variable's:

    >>> system.initial_condition_values
    {c: 0}

Reading the two together: :math:`c(0)` is :math:`c_0`, and :math:`c_0` is zero.

Choose the variable order.  :math:`c_0` goes towards the back, which is the end that is
*kept*: the front is the end that gets scaled by, and a quantity known to be zero must not
be scaled by.

    >>> system.reorder_variables(['s_0', 'e_0', 'k_1', 'k_2', 'k_m1', 'K_m', 'c_0', 't', 's', 'c'])

Now the exponent matrix can be read against that order: one row per variable, in the order
just chosen, and one column per monomial appearing in the system.

    >>> system.exponent_matrix()
    Matrix([
    [0, 0,  0, 0, 0, 0,  0, -1,  0,  0,  0],
    [1, 0,  0, 0, 0, 0,  1,  0,  0,  0,  0],
    [1, 1,  0, 1, 0, 0,  1,  0,  0,  1,  0],
    [0, 0,  0, 0, 1, 0,  0,  0,  0,  0,  1],
    [0, 0,  1, 0, 0, 1,  0,  0,  0, -1, -1],
    [0, 0,  0, 0, 0, 0,  0,  0,  0,  1,  0],
    [0, 0,  0, 0, 0, 0,  0,  0, -1,  0,  0],
    [1, 1,  1, 1, 1, 1,  1,  0,  0,  0,  0],
    [0, 0, -1, 1, 0, 0,  1,  1,  0,  0,  0],
    [0, 1,  1, 0, 0, 0, -1,  0,  1,  0,  0]])

Compute the reduced system by making an `ODETranslator` from the original system,

The reduced parameters are named :math:`\kappa_i`.  Here :math:`c` is a dependent variable
and its initial value is :math:`c_0`, which a parameter named :math:`c_0` would be hard to
read beside.

    >>> translator = ODETranslator.from_ode_system(system, naming_scheme={'parameters': 'kappa', 'independent': 't', 'dependents': ['s','c']})

Optionally, observe some properties.  First, the scaling matrix :math:`A`:

    >>> translator.scaling_matrix
    Matrix([
    [0, 0, -1, -1, -1, 0, 0, 1, 0, 0],
    [1, 1, -1,  0,  0, 1, 1, 0, 1, 1]])

The corresponding scaling invariants:

    >>> translator.invariant_monomials()
    Matrix([[e_0/s_0, k_2/(k_1*s_0), k_m1/(k_1*s_0), K_m/s_0, c_0/s_0, k_1*s_0*t, s/s_0, c/s_0]])

These don't match those from Murray, but they are very close.   We'll recover it next.

The reduction divides by these symbols, so none of them may be zero.  :math:`c_0` is not
among them.

    >>> translator.symbols_that_must_not_be_zero()
    (s_0, k_1)

Just for curiosity, here's what we get if we translate the system at this point:

    >>> translator.translate(system)
    dt/dt = 1
    dc/dt = -c*kappa1 - c*kappa2 - c*s + kappa0*s
    ds/dt = c*kappa2 + c*s - kappa0*s
    c(0) = kappa4
    s(0) = 1
    kappa3 == kappa1 + kappa2
    parameters: kappa0, kappa1, kappa2, kappa3, kappa4
    parameter values: kappa4 = 0

:math:`\kappa_4` is :math:`c_0/s_0`, so it is zero as well: the complex starts at nothing
in the reduced system too.

Column operations to match the canonical system
*************************************************

Returning back to the original system with initial conditions,

    >>> system
    dt/dt = 1
    ds/dt = c*k_1*s + c*k_m1 - e_0*k_1*s
    dc/dt = -c*k_1*s - c*k_2 - c*k_m1 + e_0*k_1*s
    s(0) = s_0
    c(0) = c_0
    K_m == (k_2 + k_m1)/k_1
    parameters: s_0, e_0, k_1, k_2, k_m1, K_m, c_0
    parameter values: c_0 = 0
    >>> translator = ODETranslator.from_ode_system(system, naming_scheme={'parameters': 'kappa', 'independent': 'tau', 'dependents': ['u','v']})

Do some column operations on the Hermite multiplier

Rather than writing the column numbers in directly, ask which column holds which invariant.
The numbers are the same, and still worth knowing -- but now they are derived on the page
instead of appearing from nowhere.

    >>> clock = translator.column_of('k_1*s_0*t')
    >>> scale = translator.column_of('e_0/s_0')
    >>> concentration = translator.column_of('c/s_0')
    >>> clock, scale, concentration
    (7, 2, 9)

    >>> translator.multiplier_add_columns(clock, scale)  # Scale time by e_0 not s_0
    >>> translator.multiplier_add_columns(concentration, scale, power=-1)  # Scale c by e_0

Optionally, see the invariants and substitutions

    >>> translator.invariant_monomials()
    Matrix([[e_0/s_0, k_2/(k_1*s_0), k_m1/(k_1*s_0), K_m/s_0, c_0/s_0, e_0*k_1*t, s/s_0, c/e_0]])

We have substitutions to make.  Note that many get normalized to :math:`1`.

    >>> translator.translate_parameter_substitutions(system)
    {s_0: 1, e_0: kappa0, k_1: 1, k_2: kappa1, k_m1: kappa2, K_m: kappa3, c_0: kappa4, t: tau/kappa0, s: u, c: kappa0*v}

Translate

    >>> reduced_system = translator.translate(system)
    >>> reduced_system
    dtau/dtau = 1
    du/dtau = kappa2*v + u*v - u
    dv/dtau = -kappa1*v/kappa0 - kappa2*v/kappa0 - u*v/kappa0 + u/kappa0
    u(0) = 1
    v(0) = kappa4/kappa0
    kappa3 == kappa1 + kappa2
    parameters: kappa0, kappa1, kappa2, kappa3, kappa4
    parameter values: kappa4 = 0

The complex is now measured in units of the enzyme, so :math:`v(0)` is a ratio of two
reduced parameters.  It is zero for any :math:`\kappa_0`, since :math:`\kappa_4` is zero.
An initial condition can be read as its expression or as its value:

    >>> reduced_system.initial_condition_values
    {u: 1, v: 0}

Reading the reduction as a graph
*************************************************

The Hermite multiplier is the weighted biadjacency matrix of a bipartite graph: the
original symbols on one side, the reduced ones on the other, and an edge carrying the
exponent with which its source appears in its target.  So the reduction can be drawn as
well as printed, and the zero can be followed through it.

    >>> graph = translator.symbol_graph(system)
    >>> graph.number_of_nodes(), graph.number_of_edges()
    (20, 22)

:math:`c_0` has exactly one edge.  It is positive, and it runs to an invariant rather than
to an auxiliary -- which is the whole of why a zero initial condition survives this
reduction.  Nothing divides by :math:`c_0`, so the invariant it feeds is simply zero.

    >>> list(graph.out_edges(('original', sympy.Symbol('c_0')), data='exponent'))
    [(('original', c_0), ('reduced', kappa4), 1)]
    >>> graph.nodes[('reduced', sympy.Symbol('kappa4'))]['monomial']
    c_0/s_0
    >>> graph.nodes[('reduced', sympy.Symbol('kappa4'))]['value']
    0

The symbols the reduction *does* divide by are those with an edge into an auxiliary, which
is set to 1, or a negative edge into an invariant.  Coming back is a different question
over the other set of symbols, and gets a different answer:

    >>> translator.symbols_that_must_not_be_zero(system)
    (s_0, e_0, k_1)
    >>> translator.reduced_symbols_that_must_not_be_zero(system)
    (x0, x1, kappa0)

The zero returns the way it came.  Every node says what it equals in terms of the other
side, so the reverse graph gives :math:`c_0` as a monomial in the reduced symbols, with
:math:`\kappa_4` a factor of it:

    >>> back = translator.symbol_graph(system, direction='reverse')
    >>> back.nodes[('original', sympy.Symbol('c_0'))]['monomial']
    kappa4*x1

Drawing it is one further call --
:meth:`~desr.ode_translation.ODETranslator.draw_symbol_graph` for a figure straight from
the translator, or :func:`desr.symbol_graph.draw_symbol_graph` on a graph you already
hold.

.. plot:: _examples/michaelis_menten_graph.py

The source layer is on the left in both panels, so the reduced column ends the first and
begins the second at the same heights, and the round trip
:math:`c_0 \rightarrow \kappa_4 \rightarrow c_0` reads as one path across the pair.  Shape
says what a symbol is and a heavy border marks what must not be zero, so the three symbols
listed above are the three outlined ones, and :math:`c_0` is plainly not among them.

The hatched hexagons are the auxiliaries, and they are the one thing the two panels
disagree about.  Going out they are normalised to 1, which is how :math:`s_0` and
:math:`k_1` are absorbed.  Coming back they are unknowns: recovering the original values
needs two of them supplied from somewhere, which is what
:class:`desr.numerics.InsufficientKnownValues` is raised about.

Drawing needs :mod:`networkx` and :mod:`matplotlib`, which the rest of `desr` does not --
``pip install desr[graphs]``.

Matching the names in the literature
*************************************************

Now we have a system, but the names of symbols don't match.  We'll do substitutions.
Murray has no name for :math:`\kappa_4`, because he sets it to zero; :math:`\gamma` until
we do the same.

    >>> to_sub =  {'kappa0': 'epsilon', 'kappa1': 'L', 'kappa2':'K-L', 'kappa3':'K', 'kappa4': 'gamma'}
    >>> reduced_system = reduced_system.subs(to_sub, subs_constraints=True)
    >>> reduced_system
    dtau/dtau = 1
    du/dtau = K*v - L*v + u*v - u
    dv/dtau = -K*v/epsilon - u*v/epsilon + u/epsilon
    u(0) = 1
    v(0) = gamma/epsilon
    parameters: K, L, epsilon, gamma
    parameter values: gamma = 0

Murray writes three groups, not four, because he sets the initial complex to zero.
Substituting the value the system holds for :math:`\gamma` does the same, and the group
goes with it.  Substituting any other value is refused, since it would contradict
:math:`c(0) = 0`.

    >>> reduced_system.subs({'gamma': 0})
    dtau/dtau = 1
    du/dtau = K*v - L*v + u*v - u
    dv/dtau = -K*v/epsilon - u*v/epsilon + u/epsilon
    u(0) = 1
    v(0) = 0
    parameters: K, L, epsilon

The only step that `desr` cannot do is multiply :math:`\frac{dv}{d\tau}` by :math:`\epsilon`.  That has to be done by hand.  And finally we have the result from :cite:`murray` (6.13).

.. math::
    :nowrap:

    \begin{align*}
    \frac{du}{d\tau} & = -u + (u+K-\lambda)v & u(0)&=1 \\
    \epsilon \frac{dv}{d\tau} &= u - (u+K)v    &  v(0)&=0.
    \end{align*}


