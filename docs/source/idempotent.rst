.. _idempotent:

Reducing a reduced system
===========================

A reduction removes the scaling symmetries a system has.  Applied a second time, to its
own output, it has nothing further to remove.  This page shows that on two examples.

In both, the reduced system is rebuilt from its symbols and equations alone before it
is reduced again.  A system that `translate` returns carries a record of the reduction
that made it; rebuilding it sets that record aside, so the second reduction sees only
the equations.


When nothing is carried
-------------------------

The parameter scheme normalises every auxiliary, so its reduced system contains only
invariants.  Michaelis-Menten is such a system:

    >>> system = ODESystem.from_equations(
    ...     ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...      'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
    >>> once = ODETranslator.from_ode_system(system).translate(system)
    >>> once
    dtau/dtau = 1
    dnu0/dtau = -kappa0*nu0 - kappa1*nu0 - nu0*nu1 + nu1
    dnu1/dtau = kappa1*nu0 + nu0*nu1 - nu1
    parameters: kappa0, kappa1

Rebuilt as a plain system, it has no scaling symmetry.  Its exponent matrix has full
row rank, so its scaling matrix has no rows (see :ref:`no-symmetry`):

    >>> plain = ODESystem(once.symbols, once.derivatives,
    ...                   indep_var=once.independent_variable)
    >>> plain.exponent_matrix().rank(), len(plain.symbols)
    (5, 5)
    >>> plain.maximal_scaling_matrix().shape
    (0, 5)

Reducing it returns the same system, under the same names:

    >>> twice = ODETranslator.from_ode_system(plain).translate(plain)
    >>> twice
    dtau/dtau = 1
    dnu0/dtau = -kappa0*nu0 - kappa1*nu0 - nu0*nu1 + nu1
    dnu1/dtau = kappa1*nu0 + nu0*nu1 - nu1
    parameters: kappa0, kappa1
    >>> twice.named_symbols
    {kappa0: kappa0, kappa1: kappa1, tau: tau, nu0: nu0, nu1: nu1}

Every symbol denotes itself.  The names the second reduction generates coincide with the
ones already there, and since each means the same quantity in both systems that is not a
clash.


When an auxiliary is carried
------------------------------

The general scheme keeps an auxiliary that the scaling moves.  The enzyme chain's
reduction keeps one, :math:`x_3`:

    >>> enzyme = ChemicalReactionNetwork.from_diagram(
    ...     'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
    >>> first = ODETranslator.from_ode_system(enzyme).translate(enzyme)
    >>> first
    dtau/dtau = 1
    dx3/dtau = -x3*y1 + x3*y2*y4
    dy1/dtau = 0
    dy2/dtau = 0
    dy3/dtau = y1*y3 + y1 - y2*y3*y4 - y3
    dy4/dtau = y1*y4 - y2*y4**2 - y2*y4 + y3
    >>> first.auxiliary_names
    (x3,)

This system does have a scaling symmetry left, and it is the auxiliary's own.  The
equation for :math:`x_3` has the form :math:`dx_3/d\tau = x_3\,H(y)`, which is unchanged
when :math:`x_3` is multiplied by a constant, and no other equation mentions
:math:`x_3`:

    >>> kept = ODESystem(first.symbols, first.derivatives,
    ...                  indep_var=first.independent_variable)
    >>> kept.maximal_scaling_matrix()
    Matrix([[0, 1, 0, 0, 0, 0]])
    >>> kept.scaling_symmetries()
    [{x3: lamda*x3}]

That freedom is the constant of integration of the quadrature.  It is not a symmetry
the first reduction missed: the invariants :math:`y_1, \ldots, y_4` and the clock are
untouched by it.

Reducing again uses that one symmetry, with :math:`x_3` as its auxiliary, and returns
the same system under new names.  The names are chosen here so that none repeats:

    >>> names = {'parameters': 'q', 'dependents': 'w', 'auxiliaries': 'u'}
    >>> second = ODETranslator.from_ode_system(kept, naming_scheme=names).translate(kept)
    >>> second.named_symbols
    {tau: tau, u0: x3, w0: y1, w1: y2, w2: y3, w3: y4}
    >>> second
    dtau/dtau = 1
    du0/dtau = u0*(-w0 + w1*w3)
    dw0/dtau = 0
    dw1/dtau = 0
    dw2/dtau = w2*(w0 + w0/w2 - w1*w3 - 1)
    dw3/dtau = w3*(w0 - w1*w3 - w1 + w2/w3)

Each symbol of the second reduction denotes one symbol of the first, and the equations
are the same ones written in those names.


Summary
---------

After one reduction a system has either no scaling symmetry, or only the symmetries
that rescale an auxiliary it still carries.  In the first case a second reduction is the
identity.  In the second it is a relabelling.  In neither does the number of symbols or
the content of the equations change.
