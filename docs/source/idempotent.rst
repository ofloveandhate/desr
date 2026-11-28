.. _idempotent:

Reducing a reduced system is idempotent
==========================================

A reduction removes the scaling symmetries a system has.  Applied a second time, to its
own output, it has nothing further to remove.
That is, reduction is an idempotent operation.

This page demonstrates idempotency on two examples.
In both, the reduced system is rebuilt from its symbols and equations alone before it
is reduced again.  A system that `translate` returns carries a record of the reduction
that made it; rebuilding it sets that record aside, so the second reduction sees only
the equations.  The point is, reducing a reduced system gives the same system back.
Reduced systems must be fixed points in the space of ODE systems
under the action of scaling symmetry reduction.


Example 1 - no auxiliaries
---------------------------

The parameter scheme normalises every auxiliary, so its reduced system contains only
invariants.  Michaelis-Menten is such a system:

    >>> system = ODESystem.from_equations(
    ...     ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...      'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
    >>> result_first = ODETranslator(system).translate()
    >>> result_first
    dtau/dtau = 1
    dnu0/dtau = -kappa0*nu0 - kappa1*nu0 - nu0*nu1 + nu1
    dnu1/dtau = kappa1*nu0 + nu0*nu1 - nu1
    parameters: kappa0, kappa1

Rebuilt as a plain system, it has no scaling symmetry.  Its exponent matrix has full
row rank, so its scaling matrix has no rows (see :ref:`no-symmetry`):

    >>> rebuilt_from_result = ODESystem(result_first.derivative_dict,
    ...                                 indep_var=result_first.independent_variable)
    >>> rebuilt_from_result.exponent_matrix().rank(), len(rebuilt_from_result.symbols)
    (5, 5)
    >>> rebuilt_from_result.maximal_scaling_matrix().shape
    (0, 5)

Reducing it returns the same system, under the same names:

    >>> result_second = ODETranslator(rebuilt_from_result).translate()
    >>> result_second
    dtau/dtau = 1
    dnu0/dtau = -kappa0*nu0 - kappa1*nu0 - nu0*nu1 + nu1
    dnu1/dtau = kappa1*nu0 + nu0*nu1 - nu1
    parameters: kappa0, kappa1
    >>> result_second.named_symbols
    {kappa0: kappa0, kappa1: kappa1, tau: tau, nu0: nu0, nu1: nu1}

Every symbol denotes itself.  An invariant that is one of the original symbols keeps its
name, and with no symmetry left that is all of them.


Example 2 - includes an auxiliary
----------------------------------

The general scheme can keep an auxiliary that the scaling moves,
particularly when the auxiliary is a variable.  The enzyme chain's
reduction keeps one, :math:`x_3`:

    >>> enzyme = ChemicalReactionNetwork.from_diagram(
    ...     'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
    >>> reduction_first_enzyme = ODETranslator(enzyme).translate()
    >>> reduction_first_enzyme
    dtau/dtau = 1
    dx3/dtau = -x3*y1 + x3*y2*y4
    dy1/dtau = 0
    dy2/dtau = 0
    dy3/dtau = y1*y3 + y1 - y2*y3*y4 - y3
    dy4/dtau = y1*y4 - y2*y4**2 - y2*y4 + y3
    >>> reduction_first_enzyme.auxiliary_names
    (x3,)

This system does have a scaling symmetry left, and it is the auxiliary's own.  The
equation for :math:`x_3` has the form :math:`dx_3/d\tau = x_3\,H(y)`, which is unchanged
when :math:`x_3` is multiplied by a constant, and no other equation mentions
:math:`x_3`:

    >>> rebuilt_reduced_enzyme = ODESystem(reduction_first_enzyme.derivative_dict,
    ...                                    indep_var=reduction_first_enzyme.independent_variable)
    >>> rebuilt_reduced_enzyme.maximal_scaling_matrix()
    Matrix([[0, 1, 0, 0, 0, 0]])
    >>> rebuilt_reduced_enzyme.scaling_symmetries()
    [{x3: lamda*x3}]

That freedom is the constant of integration of the quadrature.  It is not a symmetry
the first reduction missed: the invariants :math:`y_1, \ldots, y_4` and the clock are
untouched by it.

Reducing a second time sees that one symmetry, with :math:`x_3` as its auxiliary:

    >>> second_result_enzyme = ODETranslator(rebuilt_reduced_enzyme).translate()
    >>> second_result_enzyme.named_symbols
    {tau: tau, x0: x3, y1: y1, y2: y2, y3: y3, y4: y4}
    >>> second_result_enzyme
    dtau/dtau = 1
    dx0/dtau = x0*(-y1 + y2*y4)
    dy1/dtau = 0
    dy2/dtau = 0
    dy3/dtau = y3*(y1 + y1/y3 - y2*y4 - 1)
    dy4/dtau = y4*(y1 - y2*y4 - y2 + y3/y4)

The clock and the four invariants are unchanged and keep their names.  The auxiliary is
given a new name, as an auxiliary always is: the second reduction has one, so it is
:math:`x_0`, and it denotes the :math:`x_3` of the first.  The equations are the same.

A third reduction changes absolutely nothing.  The one auxiliary is named :math:`x_0` already:

    >>> again = ODESystem(second_result_enzyme.derivative_dict,
    ...                   indep_var=second_result_enzyme.independent_variable)
    >>> ODETranslator(again).translate().named_symbols
    {tau: tau, x0: x0, y1: y1, y2: y2, y3: y3, y4: y4}

Silviana admits, her decision to make the subscripts correspond to the partition of symbols regardless of the keeping of auxiliaries or normalization, makes the third reduction the one that literally changes nothing.  C'est la vie.  I think it's a very sensible system as-implemented.

Summary
---------

After one reduction a system has either no scaling symmetry, or only the symmetries
that rescale an auxiliary it still carries.  In the first case a second reduction is the
identity.  In the second it renames the auxiliary, by position, and leaves everything
else as it is.  In neither does the number of symbols or the content of the equations
change.
