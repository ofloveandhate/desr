.. _dimensions:

Dimensions, units and normalisation
=====================================

A scaling symmetry is a fundamentally a statement about dimensions.  If rescaling both :math:`t` and
the rate constants together leaves a system unchanged, then :math:`t` and
those constants are measured in units that trade against each other, and the system
never depended on which units you picked.  `desr` computes these relationships from the "scaling matrix",
and does not require you to manually specify dimensions for your symbols (though you can).
This page explains how desr computes and represents dimensions.


Dimensions
-----------------------

The scaling matrix :math:`A` is computed as part of the Hermite Normal Form of the exponent matrix,
which is described and illustrated explictly elsewhere in the documentation.
:math:`A` has

* one row per scaling symmetry and
* one column per symbol, and

its entries are the weights.  Two symbols with the same column scale together, so they
are measured in the same thing.  That is a dimension, and `desr` names it after a
symbol that carries it alone:

    >>> system = ODESystem.from_equations(['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...                                    'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
    >>> system.symbols
    (e_0, k_1, k_2, k_m1, t, c, s)

So the symbols correspond in this order, in the scaling matrix's columns.

    >>> system.maximal_scaling_matrix()
    Matrix([
    [0, -1, -1, -1, 1, 0, 0],
    [1, -1,  0,  0, 0, 1, 1]])

The scaling matrix tells you how a symbol scales, with a scaling coefficent. The values are the exponents in the scaling actions, that leave the system invariant.  From :math:`A`, we can infer the dimensions of the symbols.  I think it is useful to think of each scaling action as corresponding to a dimension.  The first column of (0,1) says that the first symbol does not scale with the first dimension, whatever it is, but it DOES scale with the second dimension since it has an exponent of 1.

Note that in the next snippet, desr is naming the dimensions of the symbols, after the symbols themselves.

    >>> system.dimensions_of()
    {e_0: [e_0], k_1: 1/([e_0]*[t]), k_2: 1/[t], k_m1: 1/[t], t: [t], c: [e_0], s: [e_0]}

Desr chooses the first column with all zeros and a single 1, to name a dimension.  In the result above, there are three columns (0,1), so any of those three columns names the same dimension, and desr just chose the first one.  The first one is not special.

:math:`[t]` is read "the dimension of :math:`t`".  Nothing has been declared; these
come out of the equations.  The reading is the one a modeller would make anyway —
:math:`k_1` is a second-order rate constant, so it is one over a concentration times a
time — and it is derived rather than assumed.

When there are many symbols, grouping by sharing dimension is useful.  Check this out:

    >>> system.symbols_by_dimension()
    (([e_0], (e_0, c, s)), (1/([e_0]*[t]), (k_1,)), (1/[t], (k_2, k_m1)), ([t], (t,)))

Any given monomial in the symbols of a system has a dimension too.  You can use this to check that whether a monomial you have in mind is dimensionless:

    >>> system.dimension_of('k_1*e_0*t')
    1

The result of `1` means dimensionless.

The symmetries themselves can be printed as actions rather than as a matrix, which is
usually what you want to see:

    >>> for action in system.scaling_symmetries():
    ...     print(action)
    {k_1: k_1/lamda_1, k_2: k_2/lamda_1, k_m1: k_m1/lamda_1, t: lamda_1*t}
    {e_0: e_0*lamda_2, k_1: k_1/lamda_2, c: c*lamda_2, s: lamda_2*s}

``lamda`` is sympy's spelling of the Greek letter (because `lambda` is a keyword in Python).


Specifying physical dimensions yourself
-----------------------------------------

Structural dimensions say which symbols share a unit.  They do not say what that unit
*is*, because nothing in the equations can.  Declaring one or two propagates to the
rest:

    >>> sys_w_declared_dimensions = ODESystem.from_equations(
    ...     ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...      'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
    >>> sys_w_declared_dimensions.set_dimensions({'t': 'time', 'e_0': 'concentration'})
    >>> sys_w_declared_dimensions.physical_dimensions_of()
    {e_0: concentration, k_1: 1/(concentration*time), k_2: 1/time, k_m1: 1/time, t: time, c: concentration, s: concentration}

Two declarations determined seven symbols.  Naming changes no matrix: the reduction is
the same reduction either way, and :func:`~desr.ode_system.ODESystem.set_dimension` is
a label on structure that was already there.

Declare too little and the part no declaration reaches keeps `desr`'s own name:

    >>> partial = ODESystem.from_equations(
    ...     ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...      'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
    >>> partial.set_dimension('t', 'time')
    >>> partial.physical_dimensions_of()
    {e_0: [e_0], k_1: 1/([e_0]*time), k_2: 1/time, k_m1: 1/time, t: time, c: [e_0], s: [e_0]}

Time is named, and the other dimension is still ``[e_0]``, the dimension of whatever
:math:`e_0` is.  That is everything the declaration says and nothing more: :math:`k_1`
is an inverse time per unit of :math:`e_0`, whatever that unit turns out to be.  With
nothing declared at all, every dimension is `desr`'s own, as in
:meth:`~desr.ode_system.ODESystem.dimensions_of`.


Normalisation as a choice of units
------------------------------------

In the language of Hubert-Labahn, the result of the forward-translated system
is split into two pieces: auxiliaries, and nondimensionals.
An auxiliary is a direction which a scaling moves.
When an auxiliary is a monomial purely in the parameters of the system, it cannot
vary with time, so it can be "normalized" to 1.   There's a subtle thing to teach here...

The decoupling of auxiliaries is not an assumption.  No non-dimensional
symbol ever depends on an auxiliary.   Instead, a auxiliary **is** a unit, and a unit measured in itself
is 1.  Setting :math:`1/(e_0 k_1)` to 1 says "measure time in units of
:math:`1/(e_0 k_1)`", which merely is a choice of ruler.

Desr offers ways to compute around this: :meth:`~desr.ode_translation.ODETranslator.unit_table` is a way to inspect the bookkeeping, with one
row per symbol the reduction produces.  Each row carries what the symbol ``denotes`` as
a monomial, its ``role``, what it is ``of``, what it is ``measured_in``, and its
``dimension``.

Here it is with **nothing declared**, which is all `desr` needs:

    >>> bare = ODESystem.from_equations(
    ...     ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...      'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
    >>> def print_table(translator):
    ...     for row in translator.unit_table():
    ...         print('{:<7} {:<14} {:<5} {!s:<5} {!s:<5} {}'.format(
    ...             str(row['symbol']), str(row['denotes']), row['role'],
    ...             row['of'], row['measured_in'], row['dimension']))
    >>> print_table(ODETranslator(bare))
    x0      1/(e_0*k_1)    unit  None  None  [t]
    x1      e_0            unit  None  None  [e_0]
    kappa0  k_2/(e_0*k_1)  value k_2   1/x0  1
    kappa1  k_m1/(e_0*k_1) value k_m1  1/x0  1
    tau     e_0*k_1*t      value t     x0    1
    nu0     c/e_0          value c     x1    1
    nu1     s/e_0          value s     x1    1

The two auxiliaries are the units, and the table already knows one measures whatever
:math:`t` is measured in and the other whatever :math:`e_0` is.  Everything else is a
value — a quantity divided by the unit it is measured in — and every value's dimension
is 1, which is what makes the reduced system dimensionless.
Aside from Silviana: this is SUPER cool.  We never had to name the dimensions or units of any symbol,
just sensibly name our variables, and the system itself contains the structure that
time and concentration are the two dimensions in this system.  Amazing.

Now the **same table with two declarations**:

    >>> print_table(ODETranslator(sys_w_declared_dimensions))
    x0      1/(e_0*k_1)    unit  None  None  time
    x1      e_0            unit  None  None  concentration
    kappa0  k_2/(e_0*k_1)  value k_2   1/x0  1
    kappa1  k_m1/(e_0*k_1) value k_m1  1/x0  1
    tau     e_0*k_1*t      value t     x0    1
    nu0     c/e_0          value c     x1    1
    nu1     s/e_0          value s     x1    1

Only the last column changed, from `desr`'s own names for the dimensions to yours.
That is the overlay in one picture: every monomial, every role and every pairing is the
same, because declaring a dimension changes no matrix.

The counts differ, though, and that is the one thing a declaration buys:

    >>> ODETranslator(bare).num_units_and_dimensions()
    {'units': 2, 'dimensions': None}
    >>> ODETranslator(sys_w_declared_dimensions).num_units_and_dimensions()
    {'units': 2, 'dimensions': 2}

The reduction coins two units either way.  How many *dimensions* those units measure is
:const:`None` until the declarations reach every symbol — the structural names are
placeholders, and counting distinct placeholders would have reported a number nobody
had earned.

A reduction records the units it chose, in the caller's own symbols:

    >>> ODETranslator(sys_w_declared_dimensions).normalisations()
    {e_0: 1, k_1: 1}

Read that as "measure concentration in units of :math:`e_0`, and time in units of
:math:`1/(e_0 k_1)`", not as "assume :math:`e_0 = 1`".  A note for the curious: for a system which used the
parameter reduction scheme, all the auxiliaries are parameters, and the :math:`r` constants of
integration correspond to the :math:`r` units.


Finally, under the general and dependent-variable schemes the reduced system carries the same
record of normalisations, since there normalising is a step that can be declined:

    >>> enzyme = ChemicalReactionNetwork.from_diagram(
    ...     'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
    >>> ODETranslator(enzyme).translate().normalisations
    {k_0_1: 1, k_2_3: 1, k_4_5: 1}

.. note::

   Substituting a normalisation into an invariant gives you the invariant's *value* in
   those units, not an equality between quantities.  The dimensionless reading is the
   one available before substituting.
