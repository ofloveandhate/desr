.. _linear_compartment:

A two-compartment model
=========================

Consider a linear model with two compartments.  Material enters compartment 2 at the rate
:math:`u_2`, is lost from it at the rate :math:`a_2`, and passes from it to
compartment 1 at the rate :math:`a_{12}`.  Compartment 1 loses material at the rate
:math:`a_1`.  Nothing passes from compartment 1 back to compartment 2.

.. math::

    \frac{dz_1}{dt} &= -a_1 z_1 + a_{12} z_2 \\
    \frac{dz_2}{dt} &= -(a_2 + a_{12})\, z_2 + u_2

This is the model of the second experiment in Example 13.6 of DiStefano
:cite:`DiStefano`, an example the book bases on the kinetic study of the HIV envelope
protein gp120 by Weissburg et al. :cite:`WeissburgEtal`.  In that experiment, the
protein is given into muscle (compartment 2), and is taken up in one direction
into plasma (compartment 1).  The book names the rate constants
:math:`k_{12}` for the uptake into plasma, :math:`k_{01}` for elimination from plasma
and :math:`k_{02}` for loss in muscle; they are :math:`a_{12}`, :math:`a_1` and
:math:`a_2` here.

The book describes the experiment in words and draws the model in its Fig. 13.7.
Weissburg et al. give their own model as a diagram of compartments in
Scheme 1.

This page nondimensionalises the model, and then shows what happens to the symmetry reduction when the
total loss rate :math:`a_2 + a_{12}` is replaced by a single symbol.  It is the effect
described in :ref:`dimensional_consistency`, but this time in a model with physical meaning.

The model as written
----------------------

    >>> system = ODESystem.from_tex(r'''
    ... \frac{dz_1}{dt} &= -a_1 z_1 + a_12 z_2\\
    ... \frac{dz_2}{dt} &= -(a_2 + a_12) z_2 + u_2
    ... ''')
    >>> system.parameters
    (a_1, a_12, a_2, u_2)

The reduction uses as units the parameters at the front of the symbol order; see
:meth:`~desr.ode_system.ODESystem.reorder_symbols`.  The input :math:`u_2` and the
loss rate :math:`a_2` of compartment 2 go first:

    >>> system.reorder_symbols(['u_2', 'a_2', 'a_1', 'a_12', 't', 'z_1', 'z_2'])
    >>> translator = ODETranslator(system)
    >>> translator.num_scaling_symmetries
    2
    >>> system.dimensions_of()
    {u_2: [z_1]/[t], a_2: 1/[t], a_1: 1/[t], a_12: 1/[t], t: [t], z_1: [z_1], z_2: [z_1]}

Two symmetries: the unit of time, and the unit the material is measured in.  There is
one unit of material, shared by both compartments.  That is, it's helpful to think of addition and equality as constraining dimensions of those quantities.
:math:`a_{12}` is added to :math:`a_2`, which makes them have the same units: a rate. Similarly,
:math:`a_{12}` multiplies :math:`z_2` in the equation for :math:`z_1`: what leaves
compartment 2 arrives in compartment 1 in the unit it left in.


Of the four parameters, two remain after scaling symmetry reduction:

    >>> reduced = translator.translate()
    >>> reduced
    dtau/dtau = 1
    dnu0/dtau = -kappa0*nu0 + kappa1*nu1
    dnu1/dtau = -kappa1*nu1 - nu1 + 1
    parameters: kappa0, kappa1
    >>> reduced.named_symbols
    {kappa0: a_1/a_2, kappa1: a_12/a_2, tau: a_2*t, nu0: a_2*z_1/u_2, nu1: a_2*z_2/u_2}

Time is measured in units of :math:`1/a_2`, and material in units of
:math:`u_2 / a_2`.  The two parameters left are the other two rates, each as a
multiple of :math:`a_2`.

Naively substituting erases a constraint
-----------------------------------------

Write :math:`b_{22} = a_2 + a_{12}` for the total rate of loss from compartment 2, and
put it in the equation:

.. math::

    \frac{dz_2}{dt} = -b_{22}\, z_2 + u_2

Given to the reduction in that form:

    >>> substituted = ODESystem.from_tex(r'''
    ... \frac{dz_1}{dt} &= -a_1 z_1 + a_12 z_2\\
    ... \frac{dz_2}{dt} &= -b_22 z_2 + u_2
    ... ''')
    >>> substituted.reorder_symbols(['u_2', 'b_22', 'a_12', 'a_1', 't', 'z_1', 'z_2'])
    >>> translator = ODETranslator(substituted)
    >>> translator.num_scaling_symmetries
    3
    >>> substituted.dimensions_of()
    {u_2: [z_2]/[t], b_22: 1/[t], a_12: [z_1]/([t]*[z_2]), a_1: 1/[t], t: [t], z_1: [z_1], z_2: [z_2]}

The symbol order is the same, with :math:`b_{22}` in the place of :math:`a_2`, and
:math:`a_{12}` third.  There now appear to be three symmetries, and so the first three parameters
are used as units.  The two compartments have separate units, and
:math:`a_{12}` has become a conversion factor between them.  Nothing in these two
equations says that :math:`a_{12}` is a rate, because it is no longer added to a rate.
The substitution lost the constraint.

The reduced system has only one parameter:

    >>> reduced = translator.translate()
    >>> reduced
    dtau/dtau = 1
    dnu0/dtau = -kappa0*nu0 + nu1
    dnu1/dtau = 1 - nu1
    parameters: kappa0
    >>> reduced.named_symbols
    {kappa0: a_1/b_22, tau: b_22*t, nu0: b_22**2*z_1/(a_12*u_2), nu1: b_22*z_2/u_2}

This is a correct reduction of the two equations, and it is true of the model as well:
put :math:`b_{22} = a_2 + a_{12}` into it and it still holds.  It says the dynamics
depend on the rates only through :math:`\kappa_0 = a_1 / (a_2 + a_{12})`.

What it costs is the common unit.  The two compartments are now measured in different
units, :math:`z_1` in :math:`a_{12} u_2 / b_{22}^2` and :math:`z_2` in
:math:`u_2 / b_{22}`,
so the reduced :math:`\nu_0` and :math:`\nu_1` cannot be added to give a total amount.
And :math:`\kappa_0` is not a product of powers of the model's parameters, because
the extra symmetry rescales :math:`a_{12}` while leaving :math:`b_{22}`, which
contains it, alone.

Keeping the relation preserves the scaling symmetries
------------------------------------------------------

The function :meth:`~desr.ode_system.ODESystem.add_constraint` attaches the relation to the substituted
system.  The sum is then part of the system again, and the units it implies come back:

    >>> substituted_but_kept = ODESystem.from_tex(r'''
    ... \frac{dz_1}{dt} &= -a_1 z_1 + a_12 z_2\\
    ... \frac{dz_2}{dt} &= -b_22 z_2 + u_2
    ... ''')
    >>> substituted_but_kept.add_constraint(lhs='b_22', rhs='a_2 + a_12')
    >>> substituted_but_kept.reorder_symbols(['a_2', 'u_2', 'a_1', 'a_12', 'b_22', 't', 'z_1', 'z_2'])
    >>> translator = ODETranslator(substituted_but_kept)
    >>> translator.num_scaling_symmetries
    2
    >>> substituted_but_kept.dimensions_of()
    {a_2: 1/[t], u_2: [z_1]/[t], a_1: 1/[t], a_12: 1/[t], b_22: 1/[t], t: [t], z_1: [z_1], z_2: [z_1]}
    >>> translator.translate()
    dtau/dtau = 1
    dnu0/dtau = -kappa0*nu0 + kappa1*nu1
    dnu1/dtau = -kappa2*nu1 + 1
    kappa2 == kappa1 + 1
    parameters: kappa0, kappa1, kappa2

Two symmetries, one unit of material.  The units are those of the model as written,
:math:`1/a_2` for time and :math:`u_2/a_2` for material.  The reduced system has three
parameters and a relation among them, which is the constraint written in the reduced
symbols; two of the three are free, as in the model as written.

If you don't want to keep the relation itself, then merely saying that the two have the same unit is
enough:

    >>> sys_with_same_unit = ODESystem.from_tex(r'''
    ... \frac{dz_1}{dt} &= -a_1 z_1 + a_12 z_2\\
    ... \frac{dz_2}{dt} &= -b_22 z_2 + u_2
    ... ''')
    >>> sys_with_same_unit.require_same_unit({'b_22': 'a_12'})
    >>> sys_with_same_unit.reorder_symbols(['u_2', 'b_22', 'a_1', 'a_12', 't', 'z_1', 'z_2'])
    >>> translator = ODETranslator(sys_with_same_unit)
    >>> translator.num_scaling_symmetries
    2
    >>> translator.translate()
    dtau/dtau = 1
    dnu0/dtau = -kappa0*nu0 + kappa1*nu1
    dnu1/dtau = 1 - nu1
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
