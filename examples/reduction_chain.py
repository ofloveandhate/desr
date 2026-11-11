"""A reduction drawn as a chain of systems, with each map written out as a graph.

    original  --forward-->  reduced  --normalise-->  nondimensionalised  --reverse-->  original

Four systems and the three maps between them.  The maps run across the top, drawn as the
exponent graphs that realise them; the systems sit below at the seams, because the column
a graph ends on is the column the next one begins from, and that shared column *is* the
system they share.

No arrows are drawn between the systems.  The graph above a seam is already the map from
the box on its left to the box on its right, and drawing an arrow as well would state it
twice.

What each stage is
-------------------

**original** -- what you wrote down, with its dimensions.

**reduced** -- Hubert & Labahn's reduction, every column kept.  Its invariants are
scale-invariant and so dimensionless already; what it also carries is the auxiliaries,
which is where the dimensions went.  For a general-scheme system this stage reports
parameters its own equations never use, which is the thing that made them worth
normalising.

**nondimensionalised** -- the auxiliaries that are monomials in the parameters alone are
constant, so each is set to 1, which is a choice of units rather than an assumption.  An
invariant that has thereby become the independent variable is identified with it.

This system is nondimensionalised whichever scheme ran, because invariants are exactly
the dimensionless quantities and it is built from them.  What varies is whether anything
dimensional stands beside it: an auxiliary that survives is one the scaling moves, so it
carries dimension and has to be integrated.  That is a quadrature, and the box counts
them.  A reduction with none of them is the parameter scheme.

**original again** -- the round trip closes, given back the units it spent.  It will not
guess at them: a normalised auxiliary is 1 in the units the reduction chose, and taking
that for its value in the caller's units would hand back a rescaled system wearing the
original's symbols.

Run with::

    python examples/reduction_chain.py
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from desr import ODESystem, ODETranslator
from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.symbol_graph import (symbol_graph, draw_symbol_graph, node_positions,
                               REDUCED)


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

def summary_lines(system):
    """
    What a box says about a system.

    **No symbol names.**  A reduction's symbols are generated -- `x0`, `y3` -- and those
    names are desr's rather than the reader's; they will not be the names in a paper, and
    the graph beside the box is already where a symbol is named.  What a box is for is the
    shape of the system, which is what changes along the chain.
    """
    shape = (r'$\dot{\mathbf{y}} = F(\mathbf{y})$' if is_autonomous(system)
             else r'$\dot{\mathbf{y}} = F(\mathbf{y},\ '
                  + str(system.independent_variable) + ')$')
    lines = [shape, '{} symbols, {} dependent'.format(len(system.symbols),
                                                      len(system.dependents))]

    count, idle = len(system.parameters), len(unused_parameters(system))
    if not count:
        lines.append('no parameters')
    elif idle:
        lines.append('{} parameters, {} in no equation'.format(count, idle))
    else:
        lines.append('{} parameters'.format(count))

    units, integrals = auxiliary_split(system)
    if units:
        lines.append('{} unit{} to choose'.format(units, '' if units == 1 else 's'))
    if integrals:
        lines.append('{} quadrature{}'.format(integrals,
                                              '' if integrals == 1 else 's'))
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


def identify_pairs(translator, system, normalised, final):
    """
    The second: an invariant that has *become* the independent variable is identified
    with it.

    This is a consequence of the first and not a separate rule.  The monomial over the
    arrow is a rescaled clock -- :math:`k_{23} t` -- and the rescaling has just been set
    to 1, so the invariant is now :math:`t` itself, which the system already carries.  Its
    equation reads :math:`dt/dt = 1`, which the system already has.
    """
    clock = final.independent_variable
    monomials = normalised.named_symbols
    survivors = set(final.symbols)
    return [(symbol, clock, monomials.get(symbol), 'identify')
            for symbol in normalised.symbols if symbol not in survivors]


NORMALISE_COLOUR, IDENTIFY_COLOUR = '#8a4500', '#1f5c86'


def draw_substitution(ax, pairs, font, heights):
    """
    Two columns, one arrow each, with the monomial written over every arrow that moves.

    **Two steps, and the second follows from the first.**  Normalising sends a parameter
    auxiliary to 1, which is the choice of units.  Identifying sends an invariant to the
    independent variable -- and it can, only because the monomial over that arrow is a
    rescaled clock whose rescaling has just become 1.  Without the monomials the second
    step looks like a separate rule rather than a consequence, so they are what makes the
    panel say anything.
    """
    import sympy

    # Only the rows that move.  A symbol the normalisation does not touch would draw as
    # an arrow from itself to itself, which is a line of the panel spent saying nothing.
    pairs = [entry for entry in pairs if entry[3] is not None]
    if not pairs:
        ax.text(0.5, 0.5, 'nothing to do', fontsize=font, style='italic',
                ha='center', va='center', color='#777777', transform=ax.transAxes)
        ax.axis('off')
        return ax

    count = len(pairs)
    for index, (symbol, image, monomial, step) in enumerate(pairs):
        # At the height the symbol has in the graph to the left, so a reader can carry it
        # across the panels.  Only a symbol the graph does not show falls back to being
        # spread down the panel.
        y = heights.get(symbol)
        if y is None:
            y = 1.0 - index / max(count - 1, 1)
        colour = {'normalise': NORMALISE_COLOUR,
                  'identify': IDENTIFY_COLOUR}.get(step, '#aaaaaa')
        ax.text(-0.08, y, '$' + sympy.latex(symbol) + '$', fontsize=font,
                ha='right', va='center', color='black')
        ax.annotate('', xy=(0.92, y), xytext=(0.0, y),
                    arrowprops=dict(arrowstyle='-|>', color=colour,
                                    linewidth=1.4 if step else 0.7,
                                    shrinkA=0, shrinkB=0))
        if monomial is not None:
            ax.text(0.46, y + 0.030, '$' + sympy.latex(monomial) + '$',
                    fontsize=font - 1.5, ha='center', va='bottom', color=colour)
        ax.text(1.0, y, '$' + sympy.latex(image) + '$', fontsize=font,
                ha='left', va='center', color=colour if step else 'black',
                weight='bold' if step else 'normal')

    words = {'normalise': 'a choice of units',
             'identify': 'already the clock, once those 1s are made'}
    step = pairs[0][3]
    # In the panel's own coordinates: in axes coordinates this sat underneath the axes
    # and collided with the system boxes.
    ax.text(0.46, -0.14, words[step], fontsize=font - 2, ha='center', va='top',
            color=NORMALISE_COLOUR if step == 'normalise' else IDENTIFY_COLOUR)

    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.42, 1.12)
    ax.axis('off')
    return ax

def _system_box(figure, rect, title, system, font):
    ax = figure.add_axes(rect)
    ax.axis('off')
    ax.add_patch(FancyBboxPatch((0.01, 0.01), 0.98, 0.98,
                                boxstyle='round,pad=0.012,rounding_size=0.03',
                                transform=ax.transAxes, facecolor='#f6f5f2',
                                edgecolor='#aaaaaa', linewidth=0.7, zorder=0))
    ax.text(0.5, 0.95, title, fontsize=font + 1, weight='bold', ha='center',
            va='top', transform=ax.transAxes)
    lines = summary_lines(system)
    step = 0.72 / max(len(lines), 1)
    for index, line in enumerate(lines):
        ax.text(0.5, 0.72 - index * step, line, fontsize=font, ha='center', va='top',
                transform=ax.transAxes,
                color='#8a4500' if line.startswith(('parameters', 'no param',
                                                    'left over')) else 'black')
    return ax


def chain_figure(system, title, size=(13.3, 7.5), font=6.5):
    """
    Draw the whole round trip for `system`.

    Args:
        system (ODESystem): The system to reduce.
        title (str): Printed over the figure.
        size (tuple): Figure size in inches.  The default is a 16:9 slide; a journal
            page wants something nearer ``(5.5, 7)`` and a smaller `font`.
        font (float): Base point size; everything else is set relative to it.

    Returns:
        matplotlib.figure.Figure
    """
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

    figure = plt.figure(figsize=size)
    bottom, top = 0.355, 0.925

    # Five boxes distributed evenly across the full width, so the gaps between systems
    # are equal and the chain has a steady beat.  Everything else is placed from these:
    # each map is centred between the two boxes it maps between.
    margin, box_width = 0.012, 0.172
    gap = (1.0 - 2 * margin - 5 * box_width) / 4
    box_left = [margin + index * (box_width + gap) for index in range(5)]
    box_centre = [place + box_width / 2 for place in box_left]

    # The heights the forward graph gives the reduced symbols, so the two
    # substitution panels line up with it row for row.
    heights = {symbol: place[1] for (layer, symbol), place
               in node_positions(forward).items() if layer == REDUCED}

    maps = [('forward', None), ('normalise', normalise_pairs(
                translator, system, reduced, normalised)),
            ('identify', identify_pairs(translator, system, normalised, final)),
            ('reverse', None)]

    # Centred between the two boxes it maps between, so each reads as the arrow from one
    # to the next.  The two exponent graphs carry expressions and want the width; a
    # substitution is two short columns and does not.
    for index, (name, pairs) in enumerate(maps):
        centre = (box_centre[index] + box_centre[index + 1]) / 2
        reach = box_centre[index + 1] - box_centre[index]
        graph_width = reach * (0.60 if pairs is not None else 0.96)
        ax = figure.add_axes([centre - graph_width / 2, bottom,
                              graph_width, top - bottom])
        if pairs is not None:
            draw_substitution(ax, pairs, font + 1, heights)
        else:
            # The forward graph says what each new symbol *is* -- that is the whole
            # content of a change of coordinates, and without it the reduced column is a
            # row of names nobody can read.  The reverse graph says it the other way.
            draw_symbol_graph(forward if name == 'forward' else reverse, ax=ax,
                              legend=False, label_edges=False,
                              label_expressions=True, font_size=font, node_size=90)
        ax.set_title(name, fontsize=font + 2)

    stages = [('original', system), ('reduced', reduced),
              ('normalised', normalised), ('nondimensionalised', final),
              ('original', system)]
    for index, (name, each) in enumerate(stages):
        _system_box(figure, [box_left[index], 0.035, box_width, 0.265],
                    name, each, font)

    figure.suptitle(title, fontsize=font + 5, y=0.975)
    return figure


SYSTEMS = [
    ('Michaelis-Menten, two variable', ODESystem.from_equations(
        ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
         'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])),
    ('predator-prey', ODESystem.from_equations(
        ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)'])),
    ('enzyme chain, a reaction network', ChemicalReactionNetwork.from_diagram(
        'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()),
    # No parameters at all, so no units to spend: nothing normalises, the middle map is
    # the identity, and the ring does not draw.  The figure has to stay honest about that.
    ('Hubert & Labahn 6.6', ODESystem.from_equations(
        ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
         'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])),
]

_label, _system = SYSTEMS[2]
figure = chain_figure(_system, _label)


if __name__ == '__main__':
    for label, each in SYSTEMS:
        final = ODETranslator.from_ode_system(each).translate_general(each)
        print('{:<34} {} -> {} symbols, normalises {}, left over {}'.format(
            label, len(each.symbols), len(final.symbols),
            dict(final.normalisations) or '{}',
            tuple(map(str, final.auxiliary_names)) or '()'))
        chain_figure(each, label)
    plt.show()
