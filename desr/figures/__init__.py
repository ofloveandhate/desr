"""
Figures that explain a reduction, drawn from a system desr has never seen.

Two of them, and they are halves of the same explanation:

:func:`~desr.figures.symbol_map.symbol_map_figure`
    the reduction as a chain of systems, with the map between each pair drawn out --
    the exponent graphs on the outer maps, and a table of units and values in the
    middle.  It holds no matrices.

:func:`~desr.figures.flowchart.flowchart_figure`
    the same reduction with the matrices on the page and a step named on every arrow,
    in nine numbered boxes.  It holds no equations until the last one.

Neither computes anything.  Every matrix, monomial, dimension and substitution on
either page comes from a method of :class:`~desr.ode_system.ODESystem` or
:class:`~desr.ode_translation.ODETranslator`, so a figure cannot draw a derivation the
library does not perform, and a figure drawn for a model desr has never met is as
right as the library is.

Both need :mod:`matplotlib`, and the symbol map needs :mod:`networkx` through
:mod:`desr.symbol_graph`; they are the ``graphs`` extra.  Neither is imported by
``import desr``.

    >>> from desr import ODESystem                        # doctest: +SKIP
    >>> from desr.figures import flowchart_figure         # doctest: +SKIP
    >>> figure = flowchart_figure(ODESystem.from_equations([...]))   # doctest: +SKIP
    >>> figure.savefig('reduction.pdf')                   # doctest: +SKIP
"""

from desr.figures.flowchart import flowchart_figure
from desr.figures.symbol_map import symbol_map_figure

__all__ = ['flowchart_figure', 'symbol_map_figure']
