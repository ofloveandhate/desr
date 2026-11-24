.. _linear_compartment:

A two-compartment model
=========================

A linear model with two compartments.  Material enters compartment 2 at the rate
:math:`u_2`, leaves it at the rate :math:`a_2`, and passes from it to compartment 1 at
the rate :math:`a_{12}`.  Compartment 1 loses material at the rate :math:`a_1`.

.. math::

    \frac{dz_1}{dt} &= -a_1 z_1 + a_{12} z_2 \\
    \frac{dz_2}{dt} &= -(a_2 + a_{12})\, z_2 + u_2

The model is from Weissburg et al. :cite:`WeissburgEtal`, and is Example 13.6 of
DiStefano :cite:`DiStefano`.

The page reduces the model, and then shows what happens to the reduction when the
total loss rate :math:`a_2 + a_{12}` is replaced by a single symbol.  It is the effect
described in :ref:`dimensional_consistency`, on a model where it has a physical
reading.

The model as written
----------------------

    >>> system = ODESystem.from_tex(r'''
    ... \frac{dz_1}{dt} &= -a_1 z_1 + a_12 z_2\\
    ... \frac{dz_2}{dt} &= -(a_2 + a_12) z_2 + u_2
    ... ''')
    >>> system.parameters
    (a_1, a_12, a_2, u_2)
    >>> translator = ODETranslator.from_ode_system(system)
    >>> translator.num_scaling_symmetries
    2
    >>> system.dimensions_of()
    {a_1: 1/[t], a_12: 1/[t], a_2: 1/[t], u_2: [z_1]/[t], t: [t], z_1: [z_1], z_2: [z_1]}

Two symmetries: the unit of time, and the unit the material is measured in.  There is
one unit of material, shared by both compartments.  The equations say so, because
:math:`a_{12}` is added to :math:`a_2`, which makes it a rate, and the same
:math:`a_{12}` multiplies :math:`z_2` in the equation for :math:`z_1`: what leaves
compartment 2 arrives in compartment 1 in the unit it left in.

Of the four parameters, two remain:

    >>> reduced = translator.translate(system)
    >>> reduced
    dtau/dtau = 1
    dnu0/dtau = kappa0*nu1 - nu0
    dnu1/dtau = -kappa0*nu1 - kappa1*nu1 + 1
    parameters: kappa0, kappa1
    >>> reduced.named_symbols
    {kappa0: a_12/a_1, kappa1: a_2/a_1, tau: a_1*t, nu0: a_1*z_1/u_2, nu1: a_1*z_2/u_2}

Time is measured in units of :math:`1/a_1`, and material in units of
:math:`u_2 / a_1`.  The two parameters left are the other two rates, each as a
multiple of :math:`a_1`.

With one symbol for the total loss rate
-----------------------------------------

Write :math:`b_{22} = a_2 + a_{12}` for the total rate of loss from compartment 2, and
put it in the equation:

.. math::

    \frac{dz_2}{dt} = -b_{22}\, z_2 + u_2

Given to the reduction in that form:

    >>> lumped = ODESystem.from_tex(r'''
    ... \frac{dz_1}{dt} &= -a_1 z_1 + a_12 z_2\\
    ... \frac{dz_2}{dt} &= -b_22 z_2 + u_2
    ... ''')
    >>> translator = ODETranslator.from_ode_system(lumped)
    >>> translator.num_scaling_symmetries
    3
    >>> lumped.dimensions_of()
    {a_1: 1/[t], a_12: [z_1]/([t]*[z_2]), b_22: 1/[t], u_2: [z_2]/[t], t: [t], z_1: [z_1], z_2: [z_2]}

There are now three symmetries.  The two compartments have separate units, and
:math:`a_{12}` has become a conversion factor between them.  Nothing in these two
equations says that :math:`a_{12}` is a rate, because it is no longer added to one.

The reduction has one parameter:

    >>> reduced = translator.translate(lumped)
    >>> reduced
    dtau/dtau = 1
    dnu0/dtau = -nu0 + nu1
    dnu1/dtau = -kappa0*nu1 + 1
    parameters: kappa0
    >>> reduced.named_symbols
    {kappa0: b_22/a_1, tau: a_1*t, nu0: a_1**2*z_1/(a_12*u_2), nu1: a_1*z_2/u_2}

This is a correct reduction of the two equations, and it is true of the model as well:
put :math:`b_{22} = a_2 + a_{12}` into it and it still holds.  It says the dynamics
depend on the rates only through :math:`\kappa_0 = (a_2 + a_{12}) / a_1`.

What it costs is the common unit.  The two compartments are now measured in different
units, :math:`z_1` in :math:`a_{12} u_2 / a_1^2` and :math:`z_2` in :math:`u_2 / a_1`,
so the reduced :math:`\nu_0` and :math:`\nu_1` cannot be added to give a total amount.
And :math:`\kappa_0` is not a product of powers of the model's parameters, because
the extra symmetry rescales :math:`a_{12}` while leaving :math:`b_{22}`, which
contains it, alone.  Whether that is acceptable depends on what the reduced system is
for.

Keeping the relation
----------------------

:meth:`~desr.ode_system.ODESystem.add_constraint` attaches the relation to the lumped
system.  The sum is then part of the system again, and the units it implies come back:

    >>> related = ODESystem.from_tex(r'''
    ... \frac{dz_1}{dt} &= -a_1 z_1 + a_12 z_2\\
    ... \frac{dz_2}{dt} &= -b_22 z_2 + u_2
    ... ''')
    >>> related.add_constraint(lhs='b_22', rhs='a_2 + a_12')
    >>> translator = ODETranslator.from_ode_system(related)
    >>> translator.num_scaling_symmetries
    2
    >>> related.dimensions_of()
    {a_1: 1/[t], a_12: 1/[t], b_22: 1/[t], u_2: [z_1]/[t], a_2: 1/[t], t: [t], z_1: [z_1], z_2: [z_1]}
    >>> translator.translate(related)
    dtau/dtau = 1
    dnu0/dtau = kappa0*nu1 - nu0
    dnu1/dtau = -kappa1*nu1 + 1
    kappa1 == kappa0 + kappa2
    parameters: kappa0, kappa1, kappa2

Two symmetries, one unit of material.  The reduced system has three parameters and a
relation among them, which is the constraint written in the reduced symbols; two of
the three are free, as in the model as written.

If the relation itself is not wanted, saying that the two have the same unit is
enough:

    >>> same_unit = ODESystem.from_tex(r'''
    ... \frac{dz_1}{dt} &= -a_1 z_1 + a_12 z_2\\
    ... \frac{dz_2}{dt} &= -b_22 z_2 + u_2
    ... ''')
    >>> same_unit.require_same_unit({'b_22': 'a_12'})
    >>> translator = ODETranslator.from_ode_system(same_unit)
    >>> translator.num_scaling_symmetries
    2
    >>> translator.translate(same_unit)
    dtau/dtau = 1
    dnu0/dtau = kappa0*nu1 - nu0
    dnu1/dtau = -kappa1*nu1 + 1
    parameters: kappa0, kappa1

Summary
---------

.. list-table::
   :header-rows: 1
   :widths: 46 18 18 18

   * - system
     - symmetries
     - units of material
     - free parameters left
   * - as written, with :math:`a_2 + a_{12}`
     - 2
     - one
     - 2
   * - with :math:`b_{22}` in place of the sum
     - 3
     - two
     - 1
   * - with :math:`b_{22}` and the constraint
     - 2
     - one
     - 2
   * - with :math:`b_{22}` and the same unit as :math:`a_{12}`
     - 2
     - one
     - 2
