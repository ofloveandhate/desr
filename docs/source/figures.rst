.. _figures:

Drawing a reduction
=====================

`desr` draws two figures of a reduction, and both are library calls rather than
scripts to adapt.  Give either one a system and it draws that system:

.. code-block:: python

    from desr import ODESystem
    from desr.figures import flowchart_figure, symbol_map_figure

    system = ODESystem.from_equations(
        ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
         'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])

    flowchart_figure(system, 'Michaelis-Menten').savefig('flowchart.pdf')
    symbol_map_figure(system, 'Michaelis-Menten').savefig('symbol_map.pdf')

Both need the ``graphs`` extra::

    pip install desr[graphs]

which brings matplotlib and networkx.  ``import desr`` needs neither.


Nothing in either figure computes
-----------------------------------

Every matrix, monomial, dimension, name and substitution on either page comes from a
getter.  Neither drawing derives anything of its own, and neither carries narrative
written for one particular model.

That is the point of them rather than a stylistic preference.  A figure that works
things out for itself is only as true as the figure's own arithmetic, and it is true
only for the system it was drawn for.  A figure that is a pure consumer of the library
is as right as the library is, for every system — including one `desr` has never seen,
which is the case you care about if you are drawing your own model.

It also means the figures are a test of the interface.  Each time one needed a quantity
the library could not supply, the answer was to add the getter rather than to compute
it in the drawing; several methods exist because a figure asked for them.


The flowchart
---------------

:func:`~desr.figures.flowchart.flowchart_figure` draws the reduction as nine numbered
boxes, from the original system through the exponent matrix :math:`K`, its row Hermite
normal form, the scaling matrix :math:`A`, the normal Hermite multiplier :math:`V`, its
inverse :math:`W`, and out to the reduced system.  The matrices are shown with their
entries, with the rows and columns labelled by the symbols they belong to.

It is page-text wide and as tall as it needs to be.  Past about twenty symbols a
page-wide figure of every matrix entry does not exist, and it raises rather than
shrinking the type past legibility — the honest answer, if not a useful one.

.. code-block:: python

    figure = flowchart_figure(system, 'Michaelis-Menten')

Pass ``width`` for a journal column, and ``font`` to set the base point size that
everything else is sized against.


The symbol map
----------------

:func:`~desr.figures.symbol_map.symbol_map_figure` draws the same reduction as a chain
of five systems — original, reduced, normalised, nondimensionalised, original again —
with the four maps between them.  The outer two maps are exponent graphs, drawn from
the Hermite multiplier; the middle two are substitutions, because by then the change of
coordinates is already made and what remains is a choice of units acting on the reduced
symbols.

Between them sits the unit table from :ref:`dimensions`, which spans both middle maps,
since choosing the units is what makes the clock identifiable at all.

.. code-block:: python

    figure = symbol_map_figure(system, 'Michaelis-Menten')

``steps`` replaces the words on the four arrows and ``notes`` adds a line of your own
under any stage, for the things `desr` cannot work out — what a particular model's
units mean, or why a stage is worth a second look.  Left out, the arrows describe the
reduction that actually ran.


The graph underneath
----------------------

Both figures are built on :mod:`desr.symbol_graph`, which is the exponent graph on its
own: a bipartite graph from the original symbols to the reduced ones, with an edge per
nonzero exponent, signed and labelled.

    >>> translator = ODETranslator.from_ode_system(ODESystem.from_equations(
    ...     ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...      'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']))
    >>> graph = translator.symbol_graph()
    >>> graph.number_of_nodes()
    14

:meth:`~desr.ode_translation.ODETranslator.draw_symbol_graph` draws it by itself, which
is what you want when the maps are the point and the matrices are not.
