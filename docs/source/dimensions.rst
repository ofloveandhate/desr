.. _dimensions:

Dimensions, units and normalisation
=====================================

A scaling symmetry is a statement about dimensions.  If rescaling :math:`t` and
rescaling the rate constants together leaves a system unchanged, then :math:`t` and
those constants are measured in units that trade against each other, and the system
never depended on which units you picked.  `desr` can say all of that from the
scaling matrix alone, before anyone declares what anything means.


Structural dimensions
-----------------------

The scaling matrix :math:`A` has one row per symmetry and one column per symbol, and
its entries are the weights.  Two symbols with the same column scale together, so they
are measured in the same thing.  That is a dimension, and `desr` names it after a
symbol that carries it alone:

    >>> system = ODESystem.from_equations(['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...                                    'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
    >>> system.dimensions_of()
    {e_0: [e_0], k_1: 1/([e_0]*[t]), k_2: 1/[t], k_m1: 1/[t], t: [t], c: [e_0], s: [e_0]}

:math:`[t]` is read "the dimension of :math:`t`".  Nothing has been declared; these
come out of the equations.  The reading is the one a modeller would make anyway —
:math:`k_1` is a second-order rate constant, so it is one over a concentration times a
time — and it is derived rather than assumed.

Symbols sharing a dimension are the more useful grouping when there are many of them:

    >>> system.symbols_by_dimension()
    (([e_0], (e_0, c, s)), (1/([e_0]*[t]), (k_1,)), (1/[t], (k_2, k_m1)), ([t], (t,)))

Any monomial has a dimension too, which is how to check that a group you have in mind
is dimensionless:

    >>> system.dimension_of('k_1*e_0*t')
    1

The symmetries themselves can be printed as actions rather than as a matrix, which is
usually what you want to see:

    >>> for action in system.scaling_symmetries():
    ...     print(action)
    {k_1: k_1/lamda_1, k_2: k_2/lamda_1, k_m1: k_m1/lamda_1, t: lamda_1*t}
    {e_0: e_0*lamda_2, k_1: k_1/lamda_2, c: c*lamda_2, s: lamda_2*s}

``lamda`` is sympy's spelling of the letter.  Each of these is checked against the
equations, not merely read off the matrix it came from.


Physical dimensions are an overlay
------------------------------------

Structural dimensions say which symbols share a unit.  They do not say what that unit
*is*, because nothing in the equations can.  Declaring one or two propagates to the
rest:

    >>> declared = ODESystem.from_equations(
    ...     ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...      'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
    >>> declared.set_dimensions({'t': 'time', 'e_0': 'concentration'})
    >>> declared.physical_dimensions_of()
    {e_0: concentration, k_1: 1/(concentration*time), k_2: 1/time, k_m1: 1/time, t: time, c: concentration, s: concentration}

Two declarations determined seven symbols.  Naming changes no matrix: the reduction is
the same reduction either way, and :func:`~desr.ode_system.ODESystem.set_dimension` is
a label on structure that was already there.

Declare too little and the symbols a declaration cannot reach are :const:`None`, not 1:

    >>> partial = ODESystem.from_equations(
    ...     ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...      'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
    >>> partial.set_dimension('t', 'time')
    >>> partial.physical_dimensions_of()
    {e_0: None, k_1: None, k_2: 1/time, k_m1: 1/time, t: time, c: None, s: None}

:const:`None` is "nobody has said", which is a different claim from "dimensionless".
Reporting 1 here would have been a quiet lie about :math:`k_1`.


Normalisation is a choice of units
------------------------------------

An auxiliary variable is a direction the scaling moves, and the reduction is free to
set it to any value.  When an auxiliary is a monomial in the parameters alone it cannot
vary, so it can be set to 1 once and for all.

That is not an assumption.  The auxiliary **is** a unit, and a unit measured in itself
is 1.  Setting :math:`1/(e_0 k_1)` to 1 says "measure time in units of
:math:`1/(e_0 k_1)`", which is a choice of ruler and changes no physics.

:meth:`~desr.ode_translation.ODETranslator.unit_table` is the whole bookkeeping, one
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
    >>> print_table(ODETranslator.from_ode_system(bare))
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
is 1, which is what makes the reduced system dimensionless.  None of that needed a
declaration.

Now the **same table with two declarations**:

    >>> print_table(ODETranslator.from_ode_system(declared))
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

    >>> ODETranslator.from_ode_system(bare).num_units_and_dimensions()
    {'units': 2, 'dimensions': None}
    >>> ODETranslator.from_ode_system(declared).num_units_and_dimensions()
    {'units': 2, 'dimensions': 2}

The reduction coins two units either way.  How many *dimensions* those units measure is
:const:`None` until the declarations reach every symbol — the structural names are
placeholders, and counting distinct placeholders would have reported a number nobody
had earned.

A reduction records the units it chose, in the caller's own symbols:

    >>> ODETranslator.from_ode_system(declared).normalisations()
    {e_0: 1, k_1: 1}

Under the general and dependent-variable schemes the reduced system carries the same
record, since there normalising is a step that can be declined:

    >>> enzyme = ChemicalReactionNetwork.from_diagram(
    ...     'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
    >>> ODETranslator.from_ode_system(enzyme).translate(enzyme).normalisations
    {k_0_1: 1, k_2_3: 1, k_4_5: 1}

Read that as "measure concentration in units of :math:`e_0`, and time in units of
:math:`1/(e_0 k_1)`", not as "assume :math:`e_0 = 1`".  The :math:`r` constants of
integration are the :math:`r` units.

.. note::

   Substituting a normalisation into an invariant gives you the invariant's *value* in
   those units, not an equality between quantities.  The dimensionless reading is the
   one available before substituting.
