"""A reduction drawn as a chain of systems, with each map written out above it.

    original -> reduced -> normalised -> nondimensionalised -> original
      forward   adopt units   identify t        reverse

Five systems in a row, an arrow between each pair naming what that step does, and above
them the maps themselves.  The systems sit at the seams, because the column a map ends
on is the column the next one begins from, and that shared column *is* the system they
share.

The two outer maps are changes of coordinates, drawn as the exponent graphs that realise
them.  The two middle maps are one table, read downwards: adopting the units is what
makes the clock identifiable, so the two steps are one table and not two panels.

What each stage is
-------------------

**original** -- what you wrote down, with its dimensions.

**reduced** -- Hubert & Labahn's reduction, every column kept.  Its invariants are
scale-invariant and so dimensionless already; what it also carries is the auxiliaries,
which is where the dimensions went.  For a general-scheme system this stage reports
parameters its own equations never use, which is the thing that made them worth
normalising.

**normalised** -- the auxiliaries that are monomials in the parameters alone are
constant, so each is set to 1.  That is an adoption of units rather than an assumption:
:math:`x_0 = 1/k_{23}` becomes the unit of time, after which the model's own
characteristic time measures 1 by definition.

**nondimensionalised** -- an invariant that is the clock measured in a unit now equal to
1 is identified with the clock.  The second step, and a consequence of the first.

This system is nondimensionalised whichever scheme ran, because invariants are exactly
the dimensionless quantities and it is built from them.  What varies is whether anything
dimensional stands beside it: an auxiliary that survives is one the scaling moves, so it
carries dimension and has to be integrated.  That is a quadrature, and the box counts
them.  A reduction with none of them is the parameter scheme.

**original again** -- the round trip closes, given back the units it chose.  It will not
guess at them: a normalised auxiliary is 1 in the units the reduction chose, and taking
that for its value in the caller's units would hand back a rescaled system wearing the
original's symbols.

Run with::

    python examples/reduction_chain.py
"""

import matplotlib.pyplot as plt
import sympy
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.transforms import Bbox

from desr import ODESystem, ODETranslator
from desr.chemical_reaction_network import ChemicalReactionNetwork
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


def unused_parameters(system):
    """Parameters that occur in no right-hand side -- what issue #37 was about."""
    used = set()
    for variable, rate in system.derivative_dict.items():
        if rate is not None:
            used |= rate.free_symbols
    return [p for p in system.parameters if p not in used]



def auxiliary_split(system):
    """
    The carried auxiliaries, split by what they cost you.

    An auxiliary built from the original parameters alone cannot vary, so it is a *unit*
    still to be chosen -- set it to 1 once and it is gone.  One the scaling moves carries
    dimension and has to be integrated, which is Hubert & Labahn's *quadrature*.  Counting
    them together would say a reduction owes four integrals when it owes one.

    Returns:
        tuple: ``(units to choose, quadratures)``.
    """
    monomials, origin = system.named_symbols, system.original_system
    if not system.auxiliary_names or monomials is None or origin is None:
        return 0, 0
    moving = set(origin.dependents) | {origin.independent_variable}
    integrals = sum(1 for name in system.auxiliary_names
                    if monomials.get(name) is not None
                    and monomials[name].free_symbols & moving)
    return len(system.auxiliary_names) - integrals, integrals

#: What marks a line as one of the symbol count's parts rather than a fact of its own.
BULLET = '·  '


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
    parameters, idle = len(system.parameters), len(unused_parameters(system))
    lines = [plural(len(system.symbols), 'symbol'),
             BULLET + '1 independent variable',
             BULLET + plural(dependents, 'dependent variable')]
    if not parameters:
        lines.append(BULLET + 'no parameters')
    elif idle:
        # On its own line: side by side the two numbers run past the box, and the second
        # is the one worth stopping on.
        lines.append(BULLET + plural(parameters, 'parameter'))
        lines.append('{} of them in no equation'.format(idle))
    else:
        lines.append(BULLET + plural(parameters, 'parameter'))

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
    the arrow is a |-> and never an =.
    """
    import sympy

    monomials = reduced.named_symbols
    survivors = set(normalised.symbols)
    return [(symbol, sympy.Integer(1), monomials.get(symbol), 'normalise')
            for symbol in reduced.symbols if symbol not in survivors]


def identify_pairs(translator, system, reduced, normalised, final):
    """
    The second of the two middle maps: an invariant that *is* the clock is identified
    with it.

    **Nothing is ever asserted about a rate constant.**  It is tempting to say the
    invariant :math:`k_{23} t` becomes :math:`t` "because :math:`k_{23} = 1`", and that
    is false -- :math:`k_{23}` is whatever the model says it is.  What is true is this:

    .. math::  y_c = t / x_0

    where :math:`x_0` is an auxiliary with the dimensions of time, the system's
    **characteristic time**.  That relation holds from the start; it is what the
    invariant *is*, not something normalisation does to it.  The invariant is the clock
    measured in units of the characteristic time, and a quantity measured in units of
    itself is 1, so choosing :math:`x_0` as the unit makes the two numbers agree.

    Checked across five systems, ``t`` divided by the clock invariant is always a
    monomial in the auxiliaries alone -- never one that touches an invariant -- which is
    why the identification can be read off the units and never needs the parameters.

    Returns:
        list: ``(symbol, image, what it equals, 'identify')``, where "what it equals" is
        written in the reduced symbols, because that is the coordinate system the
        normalisation acts in.
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
        pairs.append((symbol, clock, reads_as, 'identify'))
    return pairs



#: The whole figure's colour, in one place, so that a colour means one thing wherever it
#: appears.  Colour names a **kind of symbol** and nothing else: **rust** the
#: auxiliaries, which are the units, and **blue** the invariants, which are the values.
#:
#: Not the steps.  An arrow is grey, because the two middle ones do not divide the way
#: the colours do -- adopting units touches every rust symbol, while identifying the
#: clock touches one blue one and makes it the independent variable, which is the one
#: thing blue does not mean.
#:
#: Black is for what is neither kind -- an expression, a dimension, a count.  Grey
#: labels the figure rather than saying anything about the system, and the palest grey
#: is a step that does nothing.
UNIT_COLOUR, VALUE_COLOUR = '#8a4500', '#1f5c86'
PLAIN, QUIET, INERT = 'black', '#555555', '#aaaaaa'


def dimensions_shown(system, scaling_matrix):
    """
    The dimensions to put on a figure: the physical ones where any were declared.

    The same rule :meth:`~desr.ode_translation.ODETranslator.unit_table` follows, so the
    two halves of the panel are never reading from different tables.
    """
    return (system.physical_dimensions_of(scaling_matrix)
            or system.dimensions_of(scaling_matrix))


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

    Args:
        rows (list): :meth:`~desr.ode_translation.ODETranslator.unit_table`.
        normalisations (dict):
            :meth:`~desr.ode_translation.ODETranslator.normalisations`.
        heights, origins (dict): the heights the graph alongside gives the reduced and
            the original symbols, so either can be followed across the figure.
        dimensions (dict, optional): ``{symbol: dimension}`` for the original system,
            from :func:`dimensions_shown`.  Left out, only the right half is drawn.
    """
    ax.set_xlim(-1.30, 1.30)
    ax.set_ylim(-0.18, 1.42)
    ax.axis('off')
    if not rows:
        return ax
    renderer = _renderer(ax.figure)

    drawn = []
    if dimensions:
        for symbol, dimension in dimensions.items():
            y = (origins or {}).get(symbol)
            if y is None:
                continue
            # Bracketed, because the column states a dimension rather than a quantity:
            # `t  time` reads as t being a time, and `[t] = time` says what is meant.
            written = ax.text(
                -0.84, y, '$\\left[{}\\right]$'.format(expr_to_tex(symbol)),
                fontsize=font, ha='right', va='center')
            drawn.append(written)
            drawn.append(ax.text(
                _right_edge(ax, written, renderer) + 0.05, y,
                '$=\\;{}$'.format(display_tex(dimension_to_tex(dimension)))
                if dimension is not None else '=  not declared',
                fontsize=font - 1, ha='left', va='center'))

    meanings, widest = [], 0.03
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
        per_unit = (sympy.simplify(row['denotes'] / quantity)
                    if quantity is not None else 1)
        if quantity is None or per_unit == 1:
            reads_as = expr_to_tex(row['denotes'])
        else:
            # Juxtaposed, with a thin space.  A dot earns its place between dimension
            # *names*, which are words, and gets in the way between two symbols.
            reads_as = '{}\\,{}'.format(expr_to_tex(quantity),
                                        expr_to_tex(per_unit))

        # Where the monomial goes once the units are adopted, for every row by the one
        # rule: substitute, and say so where anything changed.  `x0 |-> 1` and
        # `y1 |-> e` are then the same statement, which is what they are -- the second
        # is not an identification and the arrow does not claim one, any more than the
        # first claims that 1/k_2_3 equals 1.
        #
        # Part of the same string, not a column of its own: set apart it reads as a
        # later step rather than as what that line of the table says.
        after = sympy.simplify(row['denotes'].subs(normalisations, simultaneous=True))
        if after != row['denotes']:
            reads_as += '\\;\\;\\mapsto\\;' + expr_to_tex(after)

        drawn.append(ax.text(0.0, y, '$' + expr_to_tex(symbol) + '$', fontsize=font,
                             ha='right', va='center', color=colour))
        written = ax.text(0.03, y, '$=\\;' + display_tex(reads_as) + '$',
                          fontsize=font - 0.5, ha='left', va='center')
        drawn.append(written)
        widest = max(widest, _right_edge(ax, written, renderer))

        if unit:
            # "unit of concentration", not "*the* unit of concentration": the model has
            # three of them, and the article says there is one.
            dimension = (dimension_to_tex(row['dimension'])
                         if row['dimension'] is not None else None)
            meaning = ('unit of ${}$'.format(display_tex(dimension)) if dimension
                       else 'unit of a dimension not declared')
            if after == row['denotes']:
                # The one kind of unit nobody can adopt: the scaling moves it, so it is
                # a quadrature rather than a number to set, and nothing sends it to 1.
                meaning += ', a quadrature'
        elif row['measured_in'] is not None:
            # In words, which is what keeps it out of the equations.  y1 is the value of
            # e in the unit x1, and saying so is a remark about two symbols rather than
            # a relation the reduced system holds.
            meaning = '${}$ in units of ${}$'.format(
                expr_to_tex(row['of']), expr_to_tex(row['measured_in']))
        else:
            meaning = None

        if meaning:
            meanings.append((y, meaning, colour))

    # One column, flush left, clear of the longest expression.  The expressions are
    # monomials of no fixed width, so where that is lies in the measuring; set each
    # against its own row instead and the column comes out ragged.
    for y, meaning, colour in meanings:
        drawn.append(ax.text(widest + 0.10, y, meaning, fontsize=font - 1,
                             ha='left', va='center', color=colour))

    # Which half is which, in the same words as the boxes below, so the two halves are
    # the two systems the first graph runs between and not two unnamed lists.
    # Above the tallest row rather than at a height picked in advance: a display-size
    # fraction is twice the height of a bare symbol, so where the top row ends is
    # something to measure.
    rows_extent = Bbox.union([artist.get_window_extent(renderer)
                              for artist in drawn])
    rows_extent = rows_extent.transformed(ax.transData.inverted())
    if dimensions:
        drawn.append(ax.text(-0.62, rows_extent.y1 + 0.07, 'original',
                             fontsize=font - 0.5, ha='center', va='bottom',
                             color=QUIET))
        drawn.append(ax.text(0.34, rows_extent.y1 + 0.07, 'reduced',
                             fontsize=font - 0.5, ha='center', va='bottom',
                             color=QUIET))

    # Over the middle of the rows, which is where the heading of a table goes and is not
    # a column position any of them share.
    headed = Bbox.union([artist.get_window_extent(renderer) for artist in drawn])
    headed = headed.transformed(ax.transData.inverted())
    drawn.append(ax.text((headed.x0 + headed.x1) / 2, headed.y1 + 0.06, title,
                         fontsize=font + 1.5, weight='bold', ha='center',
                         va='bottom'))

    _tight_box(ax, drawn, renderer)
    return ax


def _renderer(figure):
    """The renderer, for asking how wide a piece of text came out."""
    try:
        return figure.canvas.get_renderer()
    except AttributeError:
        figure.canvas.draw()
        return figure.canvas.get_renderer()


def _right_edge(ax, artist, renderer):
    """Where a drawn text ends, in the axes' own coordinates."""
    box = artist.get_window_extent(renderer)
    return ax.transData.inverted().transform((box.x1, box.y1))[0]


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
    ax.text(0.5, 0.95, title, fontsize=font + 1, weight='bold', ha='center',
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
        identified (dict): ``{invariant: the symbol it is identified with}``.

    Returns:
        tuple: four labels, one line each.
    """
    # "adopt", not "choose".  Nobody picks these: the reduction coins them out of the
    # equations and then measures in them, and the only choice anyone made was to write
    # the model down.
    # No count: the table beside it lists them, and the arrow is for the step.
    units = 'adopt units' if normalised else 'no units to adopt'
    if identified:
        clock = 'identify ' + ', '.join(
            sorted('$' + expr_to_tex(image) + '$' for image in identified.values()))
    else:
        clock = 'no clock to identify'
    return ('forward', units, clock, 'reverse')


def _step_arrow(figure, start, end, y, label, font, colour):
    """
    One arrow between two system boxes, named above it with what the step does.
    """
    figure.add_artist(FancyArrowPatch(
        (start, y), (end, y), transform=figure.transFigure,
        arrowstyle='-|>', mutation_scale=7, linewidth=0.9, color=colour,
        shrinkA=0, shrinkB=0))
    figure.text((start + end) / 2, y + 0.020, label, fontsize=font + 0.5,
                ha='center', va='bottom', color=colour, linespacing=1.3)


#: The fixed parts of the height, in inches: the strip the title sits in, the gap under
#: the graphs, the system boxes, and the margin below them.  Everything left over is the
#: band the graphs and the table share, and that band is what grows with the system.
TITLE_IN, GAP_IN, BOXES_IN, FOOT_IN = 0.66, 0.22, 1.75, 0.50

#: The width.  Three panels of expressions side by side is what sets it, and they are
#: as wide as a system's monomials are long, so there is no number to compute it from --
#: pass `size` for a journal column.
WIDTH_IN = 10.8

#: How much height one row of symbols asks for.  Three rows is the floor, below which
#: the graphs crowd their own labels.
ROW_IN, MIN_ROWS = 0.42, 3


def chain_figure(system, title=None, size=None, font=6.5, steps=None, notes=None):
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
            table alike, and the width is :data:`WIDTH_IN`.  A journal column wants
            something nearer 5.5 wide and a smaller `font`.
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
    notes = dict(notes or {})
    translator = ODETranslator.from_ode_system(system)
    reduced = translator.translate_general(system, normalise_parameters=False,
                                           identify_indep_var=False)
    # The step between: units chosen, but the clock not yet identified.  Skipping it
    # would make the identification look like a second kind of substitution rather than
    # what it is -- a consequence of the first.
    normalised = translator.translate_general(system, identify_indep_var=False)
    final = translator.translate_general(system)
    forward = symbol_graph(translator, system)
    reverse = symbol_graph(translator, system, direction='reverse')

    # One row per symbol, in the graphs and in the table alike, so the band they share
    # is what the height is made of.  Everything else is a fixed number of inches, and
    # the fractions matplotlib wants are worked out from those.
    rows = max(len(system.symbols), MIN_ROWS)
    # Without a title the top strip still has to hold the panels' own names.
    title_in = TITLE_IN if title else TITLE_IN - 0.30
    fixed = title_in + GAP_IN + BOXES_IN + FOOT_IN
    if size is None:
        size = (WIDTH_IN, fixed + rows * ROW_IN)
    height = size[1]
    band = max(height - fixed, ROW_IN)

    figure = plt.figure(figsize=size)
    bottom = (FOOT_IN + BOXES_IN + GAP_IN) / height
    top = bottom + band / height

    # Five boxes distributed evenly across the full width, so the gaps between systems
    # are equal and the chain has a steady beat.  Everything else is placed from these:
    # each map is centred between the two boxes it maps between, and each arrow runs
    # across the gap between them.
    margin, box_width = 0.012, 0.116
    gap = (1.0 - 2 * margin - 5 * box_width) / 4
    box_left = [margin + index * (box_width + gap) for index in range(5)]
    box_centre = [place + box_width / 2 for place in box_left]

    # The heights the forward graph gives each column of symbols, so the table lines up
    # with it row for row and either half can be followed back to a vertex.
    places = node_positions(forward)
    heights = {symbol: place[1] for (layer, symbol), place in places.items()
               if layer == REDUCED}
    origins = {symbol: place[1] for (layer, symbol), place in places.items()
               if layer == ORIGINAL}

    table = translator.unit_table(system)
    normalised_symbols = {symbol for symbol, _, _, _
                          in normalise_pairs(translator, system, reduced, normalised)}
    identified = {symbol: image for symbol, image, _, _
                  in identify_pairs(translator, system, reduced, normalised, final)}

    # The forward graph spans the first seam and the reverse graph the last.  The table
    # spans both middle seams, because choosing the units is what makes the clock
    # identifiable: the two steps are one table read downwards, not two panels.
    panels = [('forward', 0, 1), ('dimensions, units, and normalization', 1, 3),
              ('reverse', 3, 4)]
    for name, start, end in panels:
        centre = (box_centre[start] + box_centre[end]) / 2
        reach = box_centre[end] - box_centre[start]
        # The table holds two tables side by side and wants the room; the graphs give it
        # up, since their own width is mostly the labels hanging off their columns.
        graph_width = reach * (1.04 if start == 1 else 0.88)
        ax = figure.add_axes([centre - graph_width / 2, bottom,
                              graph_width, top - bottom])
        if start == 1:
            # Its name goes inside its box, where it reads as the table's heading rather
            # than as a second subtitle under the figure's own.
            draw_unit_table(ax, table, font + 0.5, heights,
                            translator.normalisations(system), name,
                            dimensions=dimensions_shown(system,
                                                        translator.scaling_matrix),
                            origins=origins)
        else:
            # The forward graph says what each new symbol *is* -- that is the whole
            # content of a change of coordinates, and without it the reduced column is a
            # row of names nobody can read.  The reverse graph says it the other way.
            #
            # Not the normalisation, though: this seam is the change of coordinates, and
            # the units are chosen at the next one.  Drawn here it would put a |-> 1 on
            # the caller's own k_0_1, which reads as a claim about a rate constant and is
            # not what happens.  The table says it, once, in the coordinates it acts in.
            #
            # One key for the figure, under everything, drawn by the forward graph
            # because both graphs mark the same things the same way and a key repeated
            # is a key a reader has to check is the same.
            draw_symbol_graph(forward if name == 'forward' else reverse, ax=ax,
                              legend=name == 'forward', label_edges=False,
                              label_expressions=True, show_normalisation=False,
                              font_size=font, node_size=90,
                              legend_kwargs={'loc': 'lower center', 'ncol': 4,
                                             'bbox_to_anchor': (0.5, 0.004),
                                             'bbox_transform': figure.transFigure,
                                             'columnspacing': 1.6})
            ax.set_title(name, fontsize=font + 2)

    stages = [('original', system), ('reduced', reduced),
              ('normalised', normalised), ('nondimensionalised', final),
              ('original', system)]
    # The boxes sit on top of the foot, which the key hangs in.
    box_bottom, box_height = FOOT_IN / height, BOXES_IN / height
    for index, (name, each) in enumerate(stages):
        _system_box(figure, [box_left[index], box_bottom, box_width, box_height],
                    name, each, font, notes.get(name))

    # The arrows, across the gaps, at the height of the boxes they join.  The outer two
    # are the graphs above them and the inner two the chain's own steps, which is why
    # they are coloured as the table colours its two halves -- and greyed where the step
    # does nothing, so the colour means the same thing as the words.
    arrow_y = box_bottom + box_height * 0.52
    if steps is None:
        steps = step_labels(normalised_symbols, identified)
    colours = [QUIET,
               QUIET if normalised_symbols else INERT,
               QUIET if identified else INERT,
               QUIET]
    for index, (label, colour) in enumerate(zip(steps, colours)):
        _step_arrow(figure, box_left[index] + box_width, box_left[index + 1],
                    arrow_y, label, font, colour)

    if title:
        figure.suptitle(title, fontsize=font + 5, y=1 - 0.24 / height)
    return figure


def declared(system, dimensions):
    """
    A system with its physical dimensions named.

    Naming them changes no matrix and no reduction; it changes what the table can say.
    Undeclared, the enzyme chain's three species sit in three structural dimensions and
    the table writes them ``[e]``, ``[f]``, ``[x]``.  Declared, all three are
    concentration, and the table shows one dimension carrying three units -- which is
    the fact about the model, and the one a reader asks about.
    """
    system.set_dimensions(dimensions)
    return system


SYSTEMS = [
    ('Michaelis-Menten, two variable', declared(ODESystem.from_equations(
        ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
         'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']),
        {'t': 'time', 'e_0': 'concentration'})),
    ('predator-prey', declared(ODESystem.from_equations(
        ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)']),
        {'t': 'time', 'n': 'prey', 'p': 'predators'})),
    ('enzyme chain, a reaction network', declared(
        ChemicalReactionNetwork.from_diagram(
            'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system(),
        {'t': 'time', 'e': 'concentration', 'f': 'concentration',
         'x': 'concentration'})),
    # No parameters at all, so no units to choose: nothing normalises, and nothing is
    # declared either, so the table falls back to the structural dimensions.  The figure
    # has to stay honest about both.
    ('Hubert & Labahn 6.6', ODESystem.from_equations(
        ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
         'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])),
]

_label, _system = SYSTEMS[2]
figure = chain_figure(_system, _label)


if __name__ == '__main__':
    for label, each in SYSTEMS:
        translator = ODETranslator.from_ode_system(each)
        final = translator.translate_general(each)
        counts = translator.units_and_dimensions(each)
        units = '{} unit{}'.format(counts['units'],
                                   '' if counts['units'] == 1 else 's')
        measuring = ('measuring {} dimensions'.format(counts['dimensions'])
                     if counts['dimensions'] is not None
                     else 'dimensions undeclared')
        print('{:<34} {} -> {} symbols, {}, {}'.format(
            label, len(each.symbols), len(final.symbols), units, measuring))
        chain_figure(each, label)
    plt.show()
