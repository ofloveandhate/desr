.. _symbol-kinds:

What kind of thing each symbol is
===================================

Two questions can be asked of any symbol, and they are independent of each other:

1. **Does it evolve?**  A parameter sits still, a dependent variable moves, and exactly one
   symbol is the independent variable everything is integrated against.
2. **Does the scaling move it?**  The symmetry rescales some quantities and leaves others
   alone.  The ones it leaves alone are the invariants; the ones it moves are, on the
   reduced side, the auxiliaries.

Neither answer implies the other.  Keeping them apart is what lets a reduction be
described exactly, and running them together is what makes the auxiliaries confusing.


Does it evolve?
-----------------

The symbols partition into three groups, and `desr` names the collective as well as the
parts:

    >>> system = ODESystem.from_equations(['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...                                    'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
    >>> system.symbols
    (e_0, k_1, k_2, k_m1, t, c, s)
    >>> system.parameters
    (e_0, k_1, k_2, k_m1)
    >>> system.variables
    (t, c, s)
    >>> system.independent_variable
    t
    >>> system.dependents
    (c, s)

A parameter is a symbol with no equation.  A dependent is a symbol with one, **including
when that equation says zero** -- a conserved quantity is a dependent variable whose rate
happens to vanish, not a parameter.


Does the scaling move it?
----------------------------

The reduction splits the columns of the Hermite multiplier in two.  The invariants are the
monomials the scaling leaves fixed; the auxiliaries are the directions it moves along, and
there are :math:`r` of them, one per scaling symmetry.

    >>> translator = ODETranslator.from_ode_system(system)
    >>> translator.auxiliary_monomials()
    Matrix([[1/(e_0*k_1), e_0]])
    >>> translator.invariant_monomials()
    Matrix([[k_2/(e_0*k_1), k_m1/(e_0*k_1), e_0*k_1*t, c/e_0, s/e_0]])

This says nothing about whether a symbol evolves.  :math:`1/(e_0 k_1)` is moved by the
scaling and is also perfectly constant in time, because :math:`e_0` and :math:`k_1` are.

The question can be asked of an original symbol too, and there it is whether that symbol's
column of the scaling matrix is zero.  A zero column means the symbol is already
dimensionless: the scaling leaves it exactly where it is, and it survives the reduction as
itself rather than as part of some combination.

    >>> already = ODESystem.from_equations(['dn/dt = r*n*(1 - n/K)',
    ...                                     'dp/dt = p*(q - 1)'])
    >>> untouched = ODETranslator.from_ode_system(already)
    >>> [symbol for index, symbol in enumerate(already.symbols)
    ...  if untouched.scaling_matrix[:, index].is_zero_matrix]
    [q, r, t]
    >>> untouched.invariant_monomials()
    Matrix([[q, r, t, n/K]])

:math:`q` is added to a pure number, so the equations give it no room to be rescaled, and
the reduction cannot remove it -- it is a parameter of the reduced system for the same
reason it was dimensionless to begin with.  Most models have no such symbol, and then
every original symbol is moved by the scaling.


The two together
------------------

A symbol is a point in the product of the two, not a value of a single attribute.  Laying
the reduced symbols out on both axes shows which combinations actually occur.

For Michaelis-Menten, where the parameter scheme applies:

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * -
     - moved by the scaling
     - fixed by it
   * - parameter
     - :math:`1/(e_0 k_1)`, :math:`e_0`
     - :math:`\kappa_0`, :math:`\kappa_1`
   * - dependent
     - *empty*
     - :math:`\nu_0`, :math:`\nu_1`
   * - independent
     - *empty*
     - :math:`\tau`

The two empty cells are the whole of what the parameter scheme is.  Every auxiliary is a
monomial in the parameters alone, so every auxiliary is constant, so all :math:`r` of them
can be set to 1 once and discarded:

    >>> translator.translate_parameter_substitutions(system)
    {e_0: 1, k_1: 1, k_2: kappa0, k_m1: kappa1, t: tau, c: nu0, s: nu1}
    >>> translator.translate(system)
    dtau/dtau = 1
    dnu0/dtau = -kappa0*nu0 - kappa1*nu0 - nu0*nu1 + nu1
    dnu1/dtau = kappa1*nu0 + nu0*nu1 - nu1
    parameters: kappa0, kappa1

The reduced system keeps no auxiliaries at all, which is why there are no names for them
to have:

    >>> translator.auxiliary_names(system)
    ()


When the cell is not empty
-----------------------------

Fill in that cell and the reduction behaves differently.  An enzyme-catalysed chain, built
from its reaction network:

    >>> enzyme = ChemicalReactionNetwork.from_diagram(
    ...     'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
    >>> enzyme.parameters
    (k_0_1, k_2_3, k_4_5)
    >>> enzyme.dependents
    (e, f, x, y, z)

Its four auxiliaries do not agree with each other.  Three are monomials in rate constants;
the fourth is a chemical species:

    >>> transl = ODETranslator.from_ode_system(enzyme)
    >>> transl.scheme(enzyme)
    ReductionScheme.GENERAL
    >>> transl.auxiliary_monomials()
    Matrix([[1/k_2_3, k_2_3/k_0_1, k_2_3/k_4_5, x]])

So the table has an occupant where Michaelis-Menten had a hole:

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * -
     - moved by the scaling
     - fixed by it
   * - parameter
     - :math:`x_0`, :math:`x_1`, :math:`x_2`
     - *empty*
   * - dependent
     - :math:`x_3`
     - :math:`y_0 \ldots y_4`
   * - independent
     - *empty*
     - *empty*

and a reduction told to keep every auxiliary sorts them the same way, three into the
parameters and one into the dependents:

    >>> reduced = transl.translate_general(enzyme, normalise_parameters=False)
    >>> reduced.auxiliary_names
    (x0, x1, x2, x3)
    >>> reduced.parameters
    (x0, x1, x2)
    >>> reduced.dependents
    (x3, y0, y1, y2, y3, y4)

:math:`x_3` is the one with an equation, and it is the only auxiliary that has one:

    >>> reduced.derivative_dict[sympy.Symbol('x3')]
    x3*(-y0*y1 + y0*y2*y4)/t

The independent variable's row is empty on both sides here because this scheme keeps the
system's own :math:`t` rather than making a new one.


Which row an auxiliary is in decides what becomes of it
----------------------------------------------------------

The table is not only bookkeeping.  The three auxiliaries in the **parameter** row are
constant, so a reduction can set each to 1 and be done with it; :math:`x_3` is in the
dependent row, moves, and cannot be.  That is what a reduction does by default, so the one
above had to be asked not to:

    >>> reduced = transl.translate(enzyme)
    >>> reduced.auxiliary_names
    (x3,)
    >>> reduced.normalisations
    {k_0_1: 1, k_2_3: 1, k_4_5: 1}

**So an auxiliary a reduced system still carries is always one that moves**, and the
parameter row of the table is always empty in a default reduction -- not because such
auxiliaries do not arise, but because they are exactly the ones that are normalised away.
Both kinds exist; only one is kept.  To see the other, ask as above.

What that buys is a reduced system whose parameters its equations use.  Compare
:math:`x_3`'s equation with the one above: every :math:`1/t` has gone, because normalising
:math:`k_{23}` makes the invariant :math:`y_0 = k_{23} t` into :math:`t` itself, which the
reduced system already carries.

    >>> reduced.derivative_dict[sympy.Symbol('x3')]
    -x3*y1 + x3*y2*y4
    >>> reduced.parameters
    ()


Constants or quadratures
---------------------------

The difference between those two tables is the difference between recovering the original
system with constants and recovering it with integrals.

:cite:`Hubert2013c` writes the auxiliary system as :math:`\frac{dx}{dt} = x \star
F(y^{W_\partial}) \cdot V_i` and calls it a quadrature: once the reduced system is solved,
each auxiliary follows by one integration.  The integrand is built from the logarithmic
derivatives of the original symbols, and **a parameter contributes zero to it**, having no
derivative at all.

So an auxiliary that is a monomial in the parameters alone has :math:`d(\log x)/dt = 0`.
It is constant, and it supplies a constant of integration rather than an integral.  When
every auxiliary is of that kind -- the two empty cells above -- the whole recovery needs
only :math:`r` constants, and the reduction can normalise the auxiliaries to 1 and forget
them.  When one is not, as with :math:`x_3`, that one has to be integrated and the reduced
system has to carry it.

This is why the distinction is worth drawing rather than glossing: it decides whether a
reduction gives back the original system for the price of :math:`r` arbitrary constants,
or for the price of :math:`r` integrals.


Reading it off a figure
--------------------------

:meth:`~desr.ode_translation.ODETranslator.symbol_graph` carries both axes as separate
node attributes -- ``kind`` for the first and ``auxiliary`` for the second -- and
:func:`~desr.symbol_graph.draw_symbol_graph` gives each its own visual channel, so a
hatched square is an auxiliary that is a constant and a hatched circle is one that must be
integrated.  There is a worked figure in :doc:`examples/mm_init_cond`.

    >>> graph = transl.symbol_graph(enzyme)
    >>> sorted((str(data['symbol']), data['kind'])
    ...        for _, data in graph.nodes(data=True) if data['auxiliary'])
    [('x0', 'parameter'), ('x1', 'parameter'), ('x2', 'parameter'), ('x3', 'dependent')]
