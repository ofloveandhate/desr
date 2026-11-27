"""A reduction drawn as a chain of systems, with each map written out above it.

    original -> reduced -> normalised -> nondimensionalised -> original
      forward   adopt units   rescale clock     reverse

Five systems in a row, an arrow between each pair naming what that step does, and above
them the maps themselves.  The systems sit at the seams, because the column a map ends
on is the column the next one begins from, and that shared column *is* the system they
share.

The two outer maps are changes of coordinates, drawn as the exponent graphs that realise
them.  The two middle maps are one table, read downwards: both act on the reduced
system's own symbols, the first on its auxiliaries and the second on its clock.

What each stage is
-------------------

**original** -- what you wrote down, with its dimensions.

**reduced** -- the reduction, every column kept.  Its invariants are
scale-invariant and so dimensionless already; what it also carries is the auxiliaries,
which is where the dimensions went.  For a general-scheme system this stage reports
parameters its own equations never use, which is the thing that made them worth
normalising.

**normalised** -- the auxiliaries that are monomials in the parameters alone are
constant, so each is set to 1.  That is an adoption of units rather than an assumption:
:math:`x_0 = 1/k_{23}` becomes the unit of time, after which the model's own
characteristic time measures 1 by definition.

**nondimensionalised** -- an invariant that is the clock times parameters,
:math:`\\tau = k t`, is taken as the clock in place of :math:`t`.  That is a change of
independent variable, and it does not depend on the first step.

This system is nondimensionalised whichever scheme ran, because invariants are exactly
the dimensionless quantities and it is built from them.  What varies is whether anything
dimensional stands beside it: an auxiliary that survives is one the scaling moves, so it
carries dimension and has to be integrated.  That is a quadrature, and the box counts
them.  A reduction with none of them is the parameter scheme.

**original again** -- the round trip closes, given back the units it chose.  It will not
guess at them: a normalised auxiliary is 1 in the units the reduction chose, and taking
that for its value in the caller's units would hand back a rescaled system wearing the
original's symbols.
"""

import matplotlib
import matplotlib.pyplot as plt
import sympy
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.transforms import Bbox

from desr.figures._palette import (FONT_RC, INERT, PLAIN, QUIET, UNIT_COLOUR,
                                   VALUE_COLOUR)
from desr.ode_translation import ODETranslator
from desr.symbol_graph import (symbol_graph, draw_symbol_graph, node_positions,
                               ORIGINAL, REDUCED)
# expr_to_tex rather than sympy.latex: desr writes a negative subscript `k_m1`, and
# sympy sets that as a literal m.
from desr.tex_tools import dimension_to_tex, display_tex, expr_to_tex


def is_autonomous(system):
    """Whether any right-hand side mentions the independent variable."""
    clock = system.independent_variable
    return not any(rate.has(clock)
                   for variable, rate in system.derivative_dict.items()
                   if variable != clock and rate is not None)


def auxiliary_split(system):
    """
    The carried auxiliaries, split by what they cost you.

    An auxiliary built from the original parameters alone cannot vary, so it is a *unit*
    the normalisation can set to 1, after which it is gone.  One the scaling moves
    carries dimension and has to be integrated, which is a *quadrature*.
    Counting them together would say a reduction owes four integrals when it owes one.

    Returns:
        tuple: ``(units, quadratures)``.
    """
    monomials, origin = system.named_symbols, system.original_system
    if not system.auxiliary_names or monomials is None or origin is None:
        return 0, 0
    moving = set(origin.dependents) | {origin.independent_variable}
    integrals = sum(1 for name in system.auxiliary_names
                    if monomials.get(name) is not None
                    and monomials[name].free_symbols & moving)
    return len(system.auxiliary_names) - integrals, integrals


def summary_lines(system):
    """
    What a box says about a system.

    The shape of the system, which is what changes along the chain: how many symbols and
    of which kinds, whether it is autonomous, and how many units and quadratures its
    auxiliaries come to.  No symbol names -- a reduction's are generated, and the graph
    beside the box is where a symbol is named.
    """
    def plural(count, noun):
        return '{} {}{}'.format(count, noun, '' if count == 1 else 's')

    # The three kinds add up to the symbol count, so the list is the count broken out
    # rather than four facts side by side.  What a reduction does is move symbols
    # between these three lines, and a reader can watch it happen along the chain.
    dependents = len(system.dependents)
    parameters = len(system.parameters)
    lines = [plural(len(system.symbols), 'symbol'),
             '1 independent',
             plural(dependents, 'dependent')]
    lines.append(plural(parameters, 'parameter') if parameters else 'no parameters')

    units, integrals = auxiliary_split(system)
    if units:
        lines.append(plural(units, 'unit'))
    if integrals:
        lines.append(plural(integrals, 'quadrature'))
    # Last, under the counts.  It is the one line that is not a count of anything.
    if is_autonomous(system):
        lines.append('autonomous')
    return lines



def normalise_pairs(translator, system, reduced, normalised):
    """
    The first of the two middle maps: a parameter auxiliary goes to 1.

    An auxiliary built from the original parameters alone cannot vary, so it can be set
    to 1 once and for all.  That is a choice of units, not an assumption, which is why
    the arrow is a ``|->`` and never an ``=``.
    """
    import sympy

    monomials = reduced.named_symbols
    survivors = set(normalised.symbols)
    return [(symbol, sympy.Integer(1), monomials.get(symbol), 'normalise')
            for symbol in reduced.symbols if symbol not in survivors]


def rescale_pairs(translator, system, reduced, normalised, final):
    """
    The second of the two middle maps: an invariant that is the clock times parameters
    is taken as the clock.

    **Nothing is asserted about a rate constant.**  The invariant :math:`k_{23} t` is
    not :math:`t` "because :math:`k_{23} = 1`"; :math:`k_{23}` is whatever the model
    says it is.  The invariant is the clock measured in a unit: where :math:`t` divided
    by it is one of the auxiliaries :math:`x_0`,

    .. math::  y_c = t / x_0

    and :math:`x_0`, which has the dimensions of time, is the system's
    **characteristic time**.  The entry is written that way when it can be.  The unit
    need not be a single auxiliary -- it can be an auxiliary times a dimensionless
    parameter -- and then the entry is the invariant's monomial in the caller's own
    symbols.

    Returns:
        list: ``(symbol, image, what it equals, 'rescale')``, where "what it equals" is
        written in the reduced symbols where that is possible.
    """
    import sympy

    clock = final.independent_variable
    monomials = normalised.named_symbols
    survivors = set(final.symbols)
    units = {monomial: name
             for name, monomial in (reduced.named_symbols or {}).items()
             if name in reduced.auxiliary_names}

    pairs = []
    for symbol in normalised.symbols:
        if symbol in survivors:
            continue
        monomial = monomials.get(symbol)
        characteristic = (sympy.simplify(clock / monomial)
                          if monomial is not None else None)
        # Say it in the units, where it is a statement about measurement.  Only a
        # characteristic time that is not one of the auxiliaries falls back to the
        # monomial in the caller's own symbols.
        if characteristic in units:
            reads_as = clock / units[characteristic]
        else:
            reads_as = monomial
        pairs.append((symbol, clock, reads_as, 'rescale'))
    return pairs



#: An arrow here is grey rather than coloured, because the two middle steps do not
#: divide the way the symbols do: normalising touches every rust symbol, while
#: rescaling the clock touches one blue one and makes it the independent variable,
#: which is the one thing blue does not mean.  The palest grey is a step that does
#: nothing.  See :mod:`desr.figures._palette` for the rest.


def dimensions_shown(system, scaling_matrix):
    """
    The dimensions to put on a figure: in the declared names where declarations reach,
    in desr's own elsewhere.

    The same call :meth:`~desr.ode_translation.ODETranslator.unit_table` makes, so the
    two halves of the panel are never reading from different tables.
    """
    return system.physical_dimensions_of(scaling_matrix)


def draw_unit_table(ax, rows, font, heights, normalisations, title,
                    dimensions=None, origins=None):
    """
    Both tables, in one box: what the system is measured in, and what the reduction
    coins to measure it with.

    On the left, every symbol the caller wrote down and its dimension.  On the right,
    one row per symbol of the reduced system: what it denotes as a monomial in the
    caller's symbols, the map the chain's middle two steps apply to it, and in words
    what it is.  A unit says which dimension it measures; a value says which quantity it
    measures and in which unit, and is set as that quantity against the rest of its
    monomial, so the unit stands as a factor.

    The map is where the monomial goes once the units are adopted, worked out by making
    the substitution and shown wherever anything changed.  One rule for every row, so
    :math:`x_0 \\mapsto 1` and :math:`y_1 \\mapsto e` are the same kind of statement --
    an arrow and never an equals, because neither claims that :math:`1/k_{23}` *is* 1 or
    that :math:`y_1` *is* :math:`e`.  A row whose monomial survives the substitution
    carries no arrow, which is what marks a quadrature and the values measured in it.

    The columns are placed by measuring what is in them.  The expressions are monomials
    of no fixed width, so each column starts a fixed gap after the widest entry of the
    one before, and the arrows share a column of their own.  The horizontal unit is the
    inch: the axes are left exactly as wide as the table came out, and
    :func:`box_unit_table` draws the outline once the caller has put them where they go.

    Args:
        rows (list): :meth:`~desr.ode_translation.ODETranslator.unit_table`.
        normalisations (dict):
            :meth:`~desr.ode_translation.ODETranslator.normalisations`.
        heights, origins (dict): the heights the graph alongside gives the reduced and
            the original symbols, so either can be followed across the figure.
        dimensions (dict, optional): ``{symbol: dimension}`` for the original system,
            from :func:`dimensions_shown`.  Left out, only the right half is drawn.

    Returns:
        tuple: the width of the table in inches, and the artists
        :func:`box_unit_table` is to enclose.
    """
    ax.set_ylim(-0.18, 1.42)
    ax.axis('off')
    if not rows:
        return 0.0, []
    figure = ax.figure
    renderer = _renderer(figure)
    edge, between, apart = 0.07, 0.04, 0.20

    def wide(artists):
        """The width of the widest of these, in inches."""
        return max((artist.get_window_extent(renderer).width / figure.dpi
                    for artist in artists), default=0.0)

    def place(artists, x):
        for artist in artists:
            artist.set_x(x)

    drawn = []
    cursor = edge
    original_centre = None
    if dimensions:
        bracketed, stated = [], []
        for symbol, dimension in dimensions.items():
            y = (origins or {}).get(symbol)
            if y is None:
                continue
            # Bracketed, because the column states a dimension rather than a quantity:
            # `t  time` reads as t being a time, and `[t] = time` says what is meant.
            bracketed.append(ax.text(
                0.0, y, '$\\left[{}\\right]$'.format(expr_to_tex(symbol)),
                fontsize=font, ha='right', va='center'))
            stated.append(ax.text(
                0.0, y,
                '$=\\;{}$'.format(display_tex(dimension_to_tex(dimension))),
                # The same size as the left side: `[t] = [t]` is one bracket twice.
                fontsize=font, ha='left', va='center'))
        start = cursor
        cursor += wide(bracketed)
        place(bracketed, cursor)
        cursor += between
        place(stated, cursor)
        cursor += wide(stated)
        original_centre = (start + cursor) / 2
        cursor += apart
        drawn.extend(bracketed + stated)

    symbols, denoted, arrows, meanings = [], [], [], []
    for index, row in enumerate(rows):
        symbol = row['symbol']
        y = heights.get(symbol)
        if y is None:
            y = 1.0 - index / max(len(rows) - 1, 1)
        unit = row['role'] == 'unit'
        colour = UNIT_COLOUR if unit else VALUE_COLOUR

        # What the symbol denotes, in the caller's symbols, for a unit and a value
        # alike.  Never `y1 = e / x1`: true after substituting what x0 denotes, but
        # written on a figure it says that a reduced symbol depends on another reduced
        # symbol, and the reduced system's invariants and auxiliaries are independent
        # coordinates.  The equations the next box holds are in the y's alone.
        #
        # A value is set as the quantity against the rest of its monomial, which is the
        # reciprocal of the unit: `e * k_0_1/k_2_3` rather than `e*k_0_1/k_2_3`, so the
        # unit is a factor a reader can see.  Both halves are computed -- `of` from the
        # table, the other by dividing -- and a reduction desr has never met factors the
        # same way or not at all.
        quantity = row['of'] if not unit else None
        per_unit = (1 / row['unit_monomial']
                    if row['unit_monomial'] is not None else 1)
        if quantity is None or per_unit == 1:
            reads_as = expr_to_tex(row['denotes'])
        else:
            # Juxtaposed, with a thin space.  A dot earns its place between dimension
            # *names*, which are words, and gets in the way between two symbols.
            reads_as = '{}\\,{}'.format(expr_to_tex(quantity),
                                        expr_to_tex(per_unit))

        # Where the monomial goes under the normalisation.  An arrow means the
        # normalisation reached this row, so the rows without one are the quadrature
        # and the values measured in it -- nothing normalises a unit the scaling moves.
        #
        # In a column of its own, just clear of the longest monomial, so the arrows
        # line up and the rows the normalisation reaches can be read down the table.
        after = row['adopted']
        symbols.append(ax.text(0.0, y, '$' + expr_to_tex(symbol) + '$', fontsize=font,
                               ha='right', va='center', color=colour))
        denoted.append(ax.text(0.0, y, '$=\\;' + display_tex(reads_as) + '$',
                               fontsize=font - 0.5, ha='left', va='center'))
        if after != row['denotes']:
            arrows.append(ax.text(
                0.0, y, '$\\mapsto\\;' + display_tex(expr_to_tex(after)) + '$',
                fontsize=font - 0.5, ha='left', va='center'))

        if unit:
            # "unit of concentration", not "*the* unit of concentration": the model has
            # three of them, and the article says there is one.
            dimension = dimension_to_tex(row['dimension'])
            # A unit's own row states its dimension, in the brackets and the form the
            # left half uses: `[x_0] = [t]`.  That is what connects the two halves --
            # the dimensions on the left are the ones these symbols are the units of.
            bracketed_unit = '\\left[{}\\right]'.format(expr_to_tex(symbol))
            meaning = '${} = {}$'.format(bracketed_unit, display_tex(dimension))
            if after == row['denotes']:
                # The one kind of unit nobody can adopt: the scaling moves it, so it is
                # a quadrature rather than a number to set, and nothing sends it to 1.
                meaning += ', a quadrature'
        elif row['measured_in'] is not None:
            # In words, which is what keeps it out of the equations.  y1 is the value of
            # e in the unit x1, and saying so is a remark about two symbols rather than
            # a relation the reduced system holds.
            # Only the unit: the column is headed "unit", and the quantity measured is
            # the one the arrow beside it has just named.
            meaning = '${}$'.format(display_tex(expr_to_tex(row['measured_in'])))
        else:
            meaning = None

        if meaning:
            meanings.append((y, meaning, colour))

    # Each column a fixed gap after the widest entry of the one before.
    start = cursor
    cursor += wide(symbols)
    place(symbols, cursor)
    cursor += 0.03
    place(denoted, cursor)
    cursor += wide(denoted)
    if arrows:
        cursor += 0.07
        place(arrows, cursor)
        cursor += wide(arrows)
    reduced_centre = (start + cursor) / 2
    drawn.extend(symbols + denoted + arrows)

    # One column, flush left, clear of the longest expression.
    if meanings:
        cursor += 0.12
        written = [ax.text(cursor, y, meaning, fontsize=font - 1, ha='center',
                           va='center', color=colour)
                   for y, meaning, colour in meanings]
        # Centred under the heading, in a column as wide as its widest entry.
        unit_centre = cursor + wide(written) / 2
        place(written, unit_centre)
        cursor += wide(written)
        drawn.extend(written)
    else:
        unit_centre = None
    width = cursor + edge

    # The title is part of the table, and a narrow table is narrower than its title.
    # Then the columns move right by half the difference, so they stay under its middle.
    probe = ax.text(0.0, 0.0, title, fontsize=font + 1.5)
    short = wide([probe]) + 0.30 - width
    probe.remove()
    if short > 0:
        for artist in drawn:
            artist.set_x(artist.get_position()[0] + short / 2)
        reduced_centre += short / 2
        if original_centre is not None:
            original_centre += short / 2
        if unit_centre is not None:
            unit_centre += short / 2
        width += short

    # From here on a horizontal unit is an inch.
    ax.set_xlim(0.0, width)
    position = ax.get_position()
    centre = (position.x0 + position.x1) / 2
    share = width / figure.get_size_inches()[0]
    ax.set_position([centre - share / 2, position.y0, share, position.height])
    renderer = _renderer(figure)

    # Which half is which, in the same words as the boxes below, so the two halves are
    # the two systems the first graph runs between and not two unnamed lists.
    # Above the tallest row rather than at a height picked in advance: a display-size
    # fraction is twice the height of a bare symbol, so where the top row ends is
    # something to measure.
    rows_extent = Bbox.union([artist.get_window_extent(renderer)
                              for artist in drawn])
    rows_extent = rows_extent.transformed(ax.transData.inverted())
    if original_centre is not None:
        drawn.append(ax.text(original_centre, rows_extent.y1 + 0.07, 'original',
                             fontsize=font - 0.5, ha='center', va='bottom',
                             color=QUIET))
        drawn.append(ax.text(reduced_centre, rows_extent.y1 + 0.07, 'reduced',
                             fontsize=font - 0.5, ha='center', va='bottom',
                             color=QUIET))
    if unit_centre is not None:
        # A little lower than "reduced", so it reads as a heading within that half:
        # these are the reduced system's units, not the original's.
        drawn.append(ax.text(unit_centre, rows_extent.y1 + 0.025, 'unit',
                             fontsize=font - 1, ha='center', va='bottom',
                             color=QUIET))

    # Over the middle of the table.
    headed = Bbox.union([artist.get_window_extent(renderer) for artist in drawn])
    headed = headed.transformed(ax.transData.inverted())
    # Size alone, no weight: Computer Modern ships no bold face, so asking for one gets
    # the roman back with a warning.  The flowchart's panel headings are set the same way.
    drawn.append(ax.text(width / 2, headed.y1 + 0.06, title,
                         fontsize=font + 1.5, ha='center', va='bottom'))
    return width, drawn


def box_unit_table(ax, drawn):
    """
    The outline of the table :func:`draw_unit_table` drew, once its axes are in place.

    Separate because the outline is measured from the text as it finally sits, and the
    caller decides where that is after it knows how wide the table came out.
    """
    if drawn:
        _tight_box(ax, drawn, _renderer(ax.figure))
    return ax


def fit_graph(ax, separation, hang):
    """
    Make the axes of a drawn symbol graph exactly as wide as the graph needs.

    The two columns of vertices are put `separation` inches apart, and the axes reach as
    far to either side as the widest label hanging off that column, `hang` inches from
    it.  The labels are text of a fixed size, so their widths are measured, and what is
    left to choose is only how long the edges are.

    Returns:
        float: the width of the axes, in inches.
    """
    figure = ax.figure
    renderer = _renderer(figure)
    # How far each column's labels reach from where they are anchored, which is `hang`
    # outside the column.  A label of several pieces is anchored once and placed in
    # points from there, so each piece is measured from the anchor, not from itself.
    reach = {True: 0.0, False: 0.0}
    for label in ax.texts:
        anchor = getattr(label, 'xy', None) or label.get_position()
        on_the_left = anchor[0] < 0.5
        anchored = ax.transData.transform(anchor)[0]
        extent = label.get_window_extent(renderer)
        beyond = (anchored - extent.x0) if on_the_left else (extent.x1 - anchored)
        if label.get_bbox_patch() is not None:
            # The outline around an auxiliary's label, which stands this far outside
            # the text.  Its own extent is not known until it has been drawn.
            beyond += 0.3 * label.get_fontsize() * figure.dpi / 72.0
        reach[on_the_left] = max(reach[on_the_left], beyond / figure.dpi)
    edge = 0.03
    left = reach[True] + hang + edge
    right = reach[False] + hang + edge
    ax.set_xlim(-left / separation, 1.0 + right / separation)
    width = left + separation + right
    position = ax.get_position()
    centre = (position.x0 + position.x1) / 2
    share = width / figure.get_size_inches()[0]
    ax.set_position([centre - share / 2, position.y0, share, position.height])
    return width


def _renderer(figure):
    """The renderer, for asking how wide a piece of text came out."""
    try:
        return figure.canvas.get_renderer()
    except AttributeError:
        figure.canvas.draw()
        return figure.canvas.get_renderer()


def _tight_box(ax, drawn, renderer, pad=0.045):
    """
    A rounded box around the text that was drawn, `pad` outside its extent.
    """
    extent = Bbox.union([artist.get_window_extent(renderer) for artist in drawn])
    extent = extent.transformed(ax.transAxes.inverted())
    ax.add_patch(FancyBboxPatch(
        (extent.x0 - pad, extent.y0 - pad),
        extent.width + 2 * pad, extent.height + 2 * pad,
        boxstyle='round,pad=0.004,rounding_size=0.025', transform=ax.transAxes,
        facecolor='#fbfaf8', edgecolor=INERT, linewidth=0.7, zorder=0,
        clip_on=False))


def _system_box(figure, rect, title, system, font, note=None):
    ax = figure.add_axes(rect)
    ax.axis('off')
    # clip_on, because the box style's own padding puts the outline just outside the
    # axes, where matplotlib trims it and leaves the corners looking nibbled.
    ax.add_patch(FancyBboxPatch((0.012, 0.012), 0.976, 0.976,
                                boxstyle='round,pad=0.012,rounding_size=0.03',
                                transform=ax.transAxes, facecolor='#f6f5f2',
                                edgecolor=INERT, linewidth=0.7, zorder=0,
                                clip_on=False))
    ax.text(0.5, 0.95, title, fontsize=font + 1, ha='center',
            va='top', transform=ax.transAxes)
    lines = summary_lines(system)
    step = 0.66 / max(len(lines), 1)
    for index, line in enumerate(lines):
        ax.text(0.5, 0.74 - index * step, line, fontsize=font, ha='center', va='top',
                transform=ax.transAxes, color=PLAIN)
    if note:
        ax.text(0.5, 0.045, note, fontsize=font - 1, ha='center', va='bottom',
                transform=ax.transAxes, style='italic', color=QUIET, wrap=True)
    return ax


def step_labels(normalised, identified):
    """
    The words on the four arrows, worked out from what the two middle steps do.

    The outer two name the graph above them.  The inner two are the chain's own steps,
    and say what this reduction does at each: a reduction that normalises nothing is
    told apart from one that does.

    Args:
        normalised (set): The auxiliaries this reduction sets to 1.
        identified (dict): ``{invariant: the clock it becomes}``.

    Returns:
        tuple: four labels, one line each.
    """
    # Named for the system the step arrives at.  A label naming units would have to go
    # on and explain that the system it produces has none.
    units = 'normalise' if normalised else 'nothing to normalise'
    if identified:
        # The images, deduplicated: the step merges every symbol it touches into the one
        # clock, so naming that clock once is what the step does.  Several symbols reach
        # it whenever the reduced clock is renamed -- the system's own `t` goes there
        # alongside the invariant that is now the same quantity -- and listing one name
        # per source would print it twice.
        clock = 'clock is ' + ', '.join(
            sorted({'$' + expr_to_tex(image) + '$' for image in identified.values()}))
    else:
        clock = 'clock unchanged'
    return ('forward', units, clock, 'reverse')


#: How far an arrow between two boxes bends, as `arc3` measures it: the middle of the
#: arc stands half this fraction of the arrow's length above its ends.
STEP_BEND = 0.55


def _step_arrow(figure, start, end, y, label, font, colour):
    """
    One arrow from the top of a system box over to the top of the next, named above it
    with what the step does.

    Over the boxes rather than between them, so the boxes can stand close together and
    the label has the width of two half boxes to be written in.
    """
    figure.add_artist(FancyArrowPatch(
        (start, y), (end, y), transform=figure.transFigure,
        connectionstyle='arc3,rad={}'.format(-STEP_BEND),
        arrowstyle='-|>', mutation_scale=7, linewidth=0.9, color=colour,
        shrinkA=1.5, shrinkB=1.5))
    width, height = figure.get_size_inches()
    rise = 0.5 * STEP_BEND * (end - start) * width
    figure.text((start + end) / 2, y + (rise + 0.04) / height, label,
                fontsize=font + 0.5, ha='center', va='bottom', color=colour,
                linespacing=1.3)


#: The fixed parts of the height, in inches: the strip the title sits in, the gap under
#: the graphs, the system boxes, and the margin below them.  Everything left over is the
#: band the graphs and the table share, and that band is what grows with the system.
TITLE_IN, GAP_IN, BOXES_IN, FOOT_IN = 0.66, 0.62, 1.55, 0.50

#: The widths that are chosen rather than measured, in inches: between the two columns
#: of vertices in a graph, between two panels, at either edge of the figure, the least a
#: system box can be, and between two boxes.  The panels are as wide as a system's
#: monomials are long, which is measured, so the figure's own width is what these and
#: the measurements add up to.
SEPARATION_IN, BETWEEN_IN, MARGIN_IN, BOX_IN, BOX_GAP_IN = 0.65, 0.10, 0.05, 0.80, 0.14

#: How far from its vertex a graph's label is written, in inches.
HANG_IN = 0.11

#: How much height one row of symbols asks for.  Three rows is the floor, below which
#: the graphs crowd their own labels.
ROW_IN, MIN_ROWS = 0.42, 3


def symbol_map_figure(system, title=None, size=None, font=6.5, steps=None,
                      notes=None):
    """
    Draw the whole round trip for `system`.

    Everything here is read off the reduction, so any system draws without being told
    anything about itself.  `steps` and `notes` are where a caller says the part desr
    cannot: what this particular model's units mean, or why a particular stage is worth
    a second look.  The figures in this file use them, and they are the difference
    between a figure for a paper and the one a first call gives you.

    Args:
        system (ODESystem): The system to reduce.
        title (str, optional): Printed over the figure.  Left out there is none, which
            is what a journal wants: the caption says what the figure is, and a title
            says it a second time in a place the typesetter cannot reach.
        size (tuple, optional): Figure size in inches.  Left out, the height follows the
            number of symbols, which is the number of rows in the graphs and in the
            table alike, and the width is what the three panels measure, packed side by
            side, or what the five boxes need if that is more.  Given, the panels are
            spread over that width, and overlap if it is less than they need; a smaller
            `font` makes them need less.
        font (float): Base point size; everything else is set relative to it.
        steps (sequence of str, optional): The words on the four arrows, in order.
            Worked out from the reduction by :func:`step_labels` when left out.  A
            newline breaks the label over two lines, which the gaps between boxes are
            narrow enough to want.
        notes (dict, optional): ``{stage name: text}``, a line of the caller's own along
            the bottom of that system's box.  The stage names are ``original``,
            ``reduced``, ``normalised`` and ``nondimensionalised``; the chain ends on
            ``original`` a second time, which takes the same note.

    Returns:
        matplotlib.figure.Figure
    """
    # Every artist resolves its font when it is made, so the whole figure is drawn inside
    # the context rather than only created in it.  The drawing is a function of its own
    # so that the context wraps it without indenting a hundred lines under a `with`.
    with matplotlib.rc_context(FONT_RC):
        return _draw_symbol_map(system, title, size, font, steps, notes)


def _draw_symbol_map(system, title, size, font, steps, notes):
    """
    The body of :func:`symbol_map_figure`, called with the figure font in force.

    Args: as :func:`symbol_map_figure`, all of them positional, since nothing but that
        function calls this.

    :rtype: matplotlib.figure.Figure
    """
    notes = dict(notes or {})
    translator = ODETranslator(system)
    reduced = translator.translate_general(normalise_parameters=False,
                                           rescale_clock=False)
    # The step between: units chosen, the clock not yet rescaled.  Drawn so that the
    # two steps can be told apart.
    normalised = translator.translate_general(rescale_clock=False)
    final = translator.translate_general()
    forward = symbol_graph(translator)
    reverse = symbol_graph(translator, direction='reverse')

    # One row per symbol, in the graphs and in the table alike, so the band they share
    # is what the height is made of.  Everything else is a fixed number of inches, and
    # the fractions matplotlib wants are worked out from those.
    rows = max(len(system.symbols), MIN_ROWS)
    # Without a title the top strip still has to hold the panels' own names.
    title_in = TITLE_IN if title else TITLE_IN - 0.30
    fixed = title_in + GAP_IN + BOXES_IN + FOOT_IN
    given_width = None if size is None else size[0]
    height = fixed + rows * ROW_IN if size is None else size[1]
    band = max(height - fixed, ROW_IN)

    scale = font / 6.5
    figure = plt.figure(figsize=(12.0, height))
    bottom = (FOOT_IN + BOXES_IN + GAP_IN) / height
    top = bottom + band / height

    # The heights the forward graph gives each column of symbols, so the table lines up
    # with it row for row and either half can be followed back to a vertex.
    places = node_positions(forward)
    heights = {symbol: place[1] for (layer, symbol), place in places.items()
               if layer == REDUCED}
    origins = {symbol: place[1] for (layer, symbol), place in places.items()
               if layer == ORIGINAL}

    table = translator.unit_table()
    normalised_symbols = {symbol for symbol, _, _, _
                          in normalise_pairs(translator, system, reduced, normalised)}
    identified = {symbol: image for symbol, image, _, _
                  in rescale_pairs(translator, system, reduced, normalised, final)}
    if steps is None:
        steps = step_labels(normalised_symbols, identified)

    # Three panels: the forward graph, the table, the reverse graph.  Each is drawn in
    # axes put anywhere, measured, and then moved: how wide a panel is depends on the
    # system's monomials, so the figure's width is an outcome and not a setting.
    names = ('forward', 'dimensions, normalisations, and units', 'reverse')
    axes = [figure.add_axes([0.35, bottom, 0.30, top - bottom]) for _ in names]

    # The forward graph says what each new symbol *is* -- that is the whole content of a
    # change of coordinates, and without it the reduced column is a row of names nobody
    # can read.  The reverse graph says it the other way.
    #
    # Not the normalisation, though: that seam is the change of coordinates and the
    # normalisation happens at the next one.  Drawn here it would put a |-> 1 on the
    # caller's own k_0_1, which reads as a claim about a rate constant.  The table says
    # it, in the coordinates it acts in.
    #
    # One key for the figure, under everything, drawn by the forward graph because both
    # graphs mark the same things the same way and a key repeated is a key a reader has
    # to check is the same.
    widths = [0.0, 0.0, 0.0]
    for index, graph in ((0, forward), (2, reverse)):
        draw_symbol_graph(graph, ax=axes[index],
                          legend=index == 0, label_edges=False,
                          label_expressions=True, show_normalisation=False,
                          font_size=font, node_size=90,
                          label_padding=HANG_IN / SEPARATION_IN,
                          legend_kwargs={'loc': 'lower center', 'ncol': 4,
                                         'bbox_to_anchor': (0.5, 0.004),
                                         'bbox_transform': figure.transFigure,
                                         'columnspacing': 1.6})
        axes[index].set_title(names[index], fontsize=font + 2)
        widths[index] = fit_graph(axes[index], SEPARATION_IN * scale, HANG_IN * scale)
    # Its name goes inside its box, where it reads as the table's heading rather than as
    # a second subtitle under the figure's own.
    widths[1], table_text = draw_unit_table(
        axes[1], table, font + 0.5, heights, translator.normalisations(), names[1],
        dimensions=dimensions_shown(system, translator.scaling_matrix),
        origins=origins)

    # The five boxes need a width of their own: each wide enough for its longest line.
    # The arrows go over the top, from the middle of one box to the middle of the next,
    # so a label has that far to be written in and the boxes are wide enough for it.
    renderer = _renderer(figure)

    def written_width(words, fontsize):
        probe = figure.text(0, 0, words, fontsize=fontsize)
        measured = probe.get_window_extent(renderer).width / figure.dpi
        probe.remove()
        return measured

    stages = [('original', system), ('reduced', reduced),
              ('normalised', normalised), ('nondimensionalised', final),
              ('original', system)]
    box_in = max([BOX_IN * scale]
                 + [written_width(line, font) + 0.17
                    for _, each in stages for line in summary_lines(each)]
                 + [written_width(name, font + 1) + 0.17 for name, _ in stages])
    gap_in = BOX_GAP_IN * scale
    longest = max(written_width(line, font + 0.5)
                  for label in steps for line in label.split('\n'))
    box_in = max(box_in, longest + 0.08 - gap_in)
    chain_in = 5 * box_in + 4 * gap_in + 2 * MARGIN_IN

    packed = sum(widths) + 2 * BETWEEN_IN + 2 * MARGIN_IN
    if given_width is None:
        width = max(packed, chain_in)
    else:
        width = given_width
    # What the width has beyond the panels' own goes between them and at the edges, in
    # equal parts, so the panels stay evenly placed on a wider page.
    spare = max(width - packed, 0.0) / 4
    figure.set_size_inches(width, height)
    cursor = MARGIN_IN + spare
    for ax, panel in zip(axes, widths):
        ax.set_position([cursor / width, bottom, panel / width, top - bottom])
        cursor += panel + BETWEEN_IN + spare
    box_unit_table(axes[1], table_text)

    # The boxes stand close together in the middle of the page, whatever its width.
    box_width, gap = box_in / width, gap_in / width
    margin = (1.0 - 5 * box_width - 4 * gap) / 2
    box_left = [margin + index * (box_width + gap) for index in range(5)]

    # The boxes sit on top of the foot, which the key hangs in.
    box_bottom, box_height = FOOT_IN / height, BOXES_IN / height
    for index, (name, each) in enumerate(stages):
        _system_box(figure, [box_left[index], box_bottom, box_width, box_height],
                    name, each, font, notes.get(name))

    # The arrows, over the tops of the boxes they join.  The outer two are the graphs
    # above them and the inner two the chain's own steps, which is why they are coloured
    # as the table colours its two halves -- and greyed where the step does nothing, so
    # the colour means the same thing as the words.
    arrow_y = box_bottom + box_height + 0.03 / height
    colours = [QUIET,
               QUIET if normalised_symbols else INERT,
               QUIET if identified else INERT,
               QUIET]
    for index, (label, colour) in enumerate(zip(steps, colours)):
        _step_arrow(figure, box_left[index] + 0.72 * box_width,
                    box_left[index + 1] + 0.28 * box_width,
                    arrow_y, label, font, colour)

    if title:
        figure.suptitle(title, fontsize=font + 5, y=1 - 0.24 / height)
    return figure


