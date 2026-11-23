.. _dimensional_consistency:

A change of variables and the units it implies
================================================

The scaling symmetries of a system are read off its equations: two parameters that are
added together must scale together, and that is how the reduction learns they have the
same unit.  Rewriting the equations can remove such a sum, and the information goes
with it.  This page follows one small system through that.

The system as written
-----------------------

.. math::

    \frac{dx}{dt} = (a + b)\,x + b

The parameter :math:`b` appears twice.  Added to :math:`a`, it is a rate, with the unit
of :math:`1/t`.  Standing alone, it has the unit of :math:`x/t`.  Both can hold only if
:math:`x` has no unit of its own, and that is what the reduction finds:

    >>> system = ODESystem.from_tex(r'\frac{dx}{dt} &= (a + b) x + b')
    >>> translator = ODETranslator.from_ode_system(system)
    >>> translator.num_scaling_symmetries
    1
    >>> system.dimensions_of()
    {a: 1/[t], b: 1/[t], t: [t], x: 1}

There is one scaling symmetry, a change of the unit of time.  Removing it leaves one
parameter:

    >>> translator.invariant_monomials()
    Matrix([[b/a, a*t, x]])
    >>> translator.translate(system)
    dtau/dtau = 1
    dx/dtau = kappa0*x + kappa0 + x
    parameters: kappa0

The same equation with a new name for the sum
------------------------------------------------

Write :math:`c = a + b`, and put :math:`c` in the equation in place of the sum:

.. math::

    \frac{dx}{dt} = c\,x + b

Given to the reduction in that form:

    >>> substituted = ODESystem.from_tex(r'\frac{dx}{dt} &= (c) x + b')
    >>> translator = ODETranslator.from_ode_system(substituted)
    >>> translator.num_scaling_symmetries
    2
    >>> substituted.dimensions_of()
    {b: [x]/[t], c: 1/[t], t: [t], x: [x]}

Now there are two symmetries, and :math:`x` has a unit of its own.  Nothing in this
equation adds :math:`b` to a rate, so nothing says :math:`b` is one; it is free to
carry the unit :math:`x/t` alone, and :math:`x` is free to be rescaled.

The reduction removes both parameters:

    >>> translator.invariant_monomials()
    Matrix([[c*t, c*x/b]])
    >>> translator.translate(substituted)
    dtau/dtau = 1
    dnu0/dtau = nu0 + 1

This is a correct reduction of the equation :math:`dx/dt = c x + b`.  With
:math:`\nu = c x / b` and :math:`\tau = c t` that equation is :math:`d\nu/d\tau = \nu + 1`
whatever :math:`b` and :math:`c` are, provided neither is zero.  It is also a true
statement about the original system, with :math:`c = a + b`.  What has changed is the
meaning of the result.  The new variable :math:`\nu = (a + b)\,x / b` is not a product
of powers of :math:`a`, :math:`b` and :math:`x`, so it is not a choice of units for the
quantities the model was written in.  If :math:`b` and :math:`a` are physical rates
measured in the same unit, the second symmetry rescales :math:`b` and leaves
:math:`c = a + b` alone, which no change of units does.

Two ways to keep the information
----------------------------------

**Say that the two have the same unit.**
:meth:`~desr.ode_system.ODESystem.require_same_unit` records that :math:`c` and
:math:`b` scale together, as the sum did:

    >>> constrained = ODESystem.from_tex(r'\frac{dx}{dt} &= (c) x + b')
    >>> constrained.require_same_unit({'c': 'b'})
    >>> translator = ODETranslator.from_ode_system(constrained)
    >>> translator.num_scaling_symmetries
    1
    >>> constrained.dimensions_of()
    {b: 1/[t], c: 1/[t], t: [t], x: 1}
    >>> translator.translate(constrained)
    dtau/dtau = 1
    dx/dtau = kappa0*x + 1
    parameters: kappa0

One symmetry again, :math:`x` without a unit again, and one parameter left.

**Keep the original equation and state the relation.**
:meth:`~desr.ode_system.ODESystem.add_constraint` attaches :math:`c = a + b` to the
system without substituting it.  The sum is still in the equation, so the units are
still implied, and the constraint is carried through the reduction:

    >>> related = ODESystem.from_tex(r'\frac{dx}{dt} &= (a + b) x + b')
    >>> related.add_constraint('c', 'a+b')
    >>> translator = ODETranslator.from_ode_system(related)
    >>> translator.num_scaling_symmetries
    1
    >>> translator.translate(related)
    dtau/dtau = 1
    dx/dtau = kappa0*x + kappa0 + x
    kappa1 == kappa0 + 1
    parameters: kappa0, kappa1

The four side by side
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 46 18 18 18

   * - system
     - symmetries
     - unit of :math:`x`
     - parameters left
   * - :math:`(a + b)\,x + b`
     - 1
     - none
     - 1
   * - :math:`c\,x + b`
     - 2
     - its own
     - 0
   * - :math:`c\,x + b`, with :math:`c` and :math:`b` in the same unit
     - 1
     - none
     - 1
   * - :math:`(a + b)\,x + b`, with the constraint :math:`c = a + b`
     - 1
     - none
     - 2, related by the constraint

A change of variables that is consistent with the units leaves the number of scaling
symmetries as it was.  One that hides a sum can raise it.  The reduction cannot tell
which was intended; :meth:`~desr.ode_system.ODESystem.dimensions_of` shows what units
it has inferred, and is the place to check.
