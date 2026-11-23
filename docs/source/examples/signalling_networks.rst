.. _signalling_networks:

Five signalling networks
==========================

Mass-action models of phosphorylation cycles and cascades.  They are larger than the
other examples, with up to eleven species and twelve rate constants, and they all
reduce in the same way.  This page works the smallest in full and tabulates the rest.

One-site modification
-----------------------

A substrate :math:`S_0` is phosphorylated to :math:`S_1` by a kinase :math:`E`, and
dephosphorylated by a phosphatase :math:`F`.  :math:`X` is the complex of the kinase
with :math:`S_0`, and :math:`Y` the complex of the phosphatase with :math:`S_1`.

.. math::

    S_0 + E \;\underset{k_{-1}}{\overset{k_1}{\rightleftharpoons}}\; X
    \;\overset{k_2}{\longrightarrow}\; S_1 + E ,
    \qquad
    S_1 + F \;\underset{k_{-3}}{\overset{k_3}{\rightleftharpoons}}\; Y
    \;\overset{k_4}{\longrightarrow}\; S_0 + F .

The kinase is written ``Ev`` below, and :math:`k_{-1}` as ``km1``.  The equations for
the two complexes are those of the free enzymes with the sign changed, since each
enzyme is either free or bound:

    >>> S0, S1, Ev, F, X, Y = sympy.symbols('S0 S1 Ev F X Y')
    >>> rates = {
    ...     S0: sympify_expression('k4*Y + km1*X - k1*S0*Ev'),
    ...     S1: sympify_expression('k2*X + km3*Y - k3*S1*F'),
    ...     Ev: sympify_expression('(km1 + k2)*X - k1*S0*Ev'),
    ...     F: sympify_expression('(km3 + k4)*Y - k3*S1*F'),
    ... }
    >>> rates[X] = -rates[Ev]
    >>> rates[Y] = -rates[F]

The total amount of substrate, free and bound, is conserved:

    >>> sympy.expand(rates[S0] + rates[S1] + rates[X] + rates[Y])
    0

There are six species and six rate constants:

    >>> system = ODESystem.from_dict(rates)
    >>> system.dependents
    (Ev, F, S0, S1, X, Y)
    >>> system.parameters
    (k1, k2, k3, k4, km1, km3)

Its scaling symmetries
------------------------

    >>> translator = ODETranslator.from_ode_system(system)
    >>> translator.num_scaling_symmetries
    2
    >>> system.dimensions_of()
    {k1: 1/([Ev]*[t]), k2: 1/[t], k3: 1/([Ev]*[t]), k4: 1/[t], km1: 1/[t], km3: 1/[t], t: [t], Ev: [Ev], F: [Ev], S0: [Ev], S1: [Ev], X: [Ev], Y: [Ev]}

Two symmetries, and two units.  One is time.  The other is concentration, and every
species has it: a complex forms from two species and breaks up into two, and its
equation has a term in each, so none can be rescaled alone.  The rate constants of the
binding steps, :math:`k_1` and :math:`k_3`, are per concentration per time; the other
four are per time.

The reduction chooses a unit for each from the rate constants:

    >>> translator.auxiliary_monomials()
    Matrix([[1/k2, k2/k1]])

Time is measured in units of :math:`1/k_2`, the time the kinase takes to convert a
bound substrate, and concentration in units of :math:`k_2/k_1`.

The reduced system
--------------------

    >>> reduced = translator.translate(system)
    >>> reduced
    dtau/dtau = 1
    dnu0/dtau = kappa2*nu4 - nu0*nu2 + nu4
    dnu1/dtau = -kappa0*nu1*nu3 + kappa1*nu5 + kappa3*nu5
    dnu2/dtau = kappa1*nu5 + kappa2*nu4 - nu0*nu2
    dnu3/dtau = -kappa0*nu1*nu3 + kappa3*nu5 + nu4
    dnu4/dtau = -kappa2*nu4 + nu0*nu2 - nu4
    dnu5/dtau = kappa0*nu1*nu3 - kappa1*nu5 - kappa3*nu5
    parameters: kappa0, kappa1, kappa2, kappa3
    >>> reduced.named_symbols
    {kappa0: k3/k1, kappa1: k4/k2, kappa2: km1/k2, kappa3: km3/k2, tau: k2*t, nu0: Ev*k1/k2, nu1: F*k1/k2, nu2: S0*k1/k2, nu3: S1*k1/k2, nu4: X*k1/k2, nu5: Y*k1/k2}

Six species and four parameters.  Two rate constants have been used as units, and the
four parameters left are the other four rate constants, each divided by the one of
:math:`k_1` and :math:`k_2` that has the same unit.  The species keep their roles:
each reduced variable is one species, in units of :math:`k_2/k_1`.

The conservation law is still there.  The same four reduced variables sum to a
constant:

    >>> names = {monomial: symbol for symbol, monomial in reduced.named_symbols.items()}
    >>> k1, k2 = sympy.symbols('k1 k2')
    >>> total = sum(reduced.derivative_dict[names[species * k1 / k2]]
    ...             for species in (S0, S1, X, Y))
    >>> sympy.expand(total)
    0

The other four
----------------

The same construction gives four larger networks.  They are built in
``tests/example_networks.py``, each with an assertion of the conservation law it
satisfies.

- **Two-step phosphorelay.**  A signal activates :math:`X_1`, whose activated form
  passes its phosphate group to :math:`X_2`; each has a phosphatase.
- **Two-site modification.**  A substrate with two phosphorylation sites,
  :math:`S_0 \to S_1 \to S_2`, with a kinase and a phosphatase for each step.
- **Two-substrate modification.**  One kinase acting on two substrates, each with its
  own phosphatase.
- **Two-layer cascade.**  The activated form of one substrate is the kinase for a
  second.

.. list-table::
   :header-rows: 1
   :widths: 40 15 15 15 15

   * - network
     - species
     - rate constants
     - symmetries
     - rate constants left
   * - one-site modification
     - 6
     - 6
     - 2
     - 4
   * - two-step phosphorelay
     - 10
     - 10
     - 2
     - 8
   * - two-site modification
     - 11
     - 12
     - 2
     - 10
   * - two-substrate modification
     - 11
     - 12
     - 2
     - 10
   * - two-layer cascade
     - 11
     - 12
     - 2
     - 10

All five have the same two symmetries, and all five reductions choose the same two
units: time in :math:`1/k_2` and concentration in :math:`k_2/k_1`.  In each, the
reduction removes two rate constants, whatever the size of the network.

The script this page replaces identified the two-layer cascade with equations
(35)–(38) of `arXiv:1101.5799 <https://arxiv.org/abs/1101.5799>`_.

The reduction applies the parameter scheme to each of them; see :ref:`three-schemes`.
