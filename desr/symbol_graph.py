'''
The exponent graph: the Hermite multiplier read as a bipartite graph.

The multiplier already *is* a weighted biadjacency matrix, so nothing here is computed that
the reduction did not already know.  What this module adds is a name for the object and a
way to look at it.

One layer of vertices is the original symbols, the other is the reduced ones -- the
:math:`r` auxiliary columns followed by the :math:`n - r` invariants -- and an edge carries
the exponent with which its source appears in its target.  Reading it:

- the **sign** of an edge is numerator versus denominator.  A negative edge into an
  invariant means that invariant divides by that symbol;
- :meth:`~desr.ode_translation.ODETranslator.symbols_that_must_not_be_zero` is an edge
  filter: a symbol qualifies when it has any edge into an auxiliary, which the reduction
  sets to 1, or a negative edge into an invariant.  A *positive* edge into an invariant and
  nothing into an auxiliary is the harmless case, where the invariant is simply zero, and is
  how a zero initial condition survives a reduction;
- the **degree** of an original symbol is how many reduced quantities it takes part in;
- an edge into the auxiliary block says which parameter was absorbed, and into which
  auxiliary.

Both directions are available, and they are different graphs.  ``'forward'`` is the
multiplier, original to reduced; ``'reverse'`` is its inverse, reduced to original.  They
compose to the identity, which in graph terms is the matrix product as weighted
path-counting: along a two-hop path the two exponents multiply, parallel paths sum, and the
result is 1 from a symbol back to itself and 0 to every other.

This module needs :mod:`networkx`, and :func:`draw_symbol_graph` also needs
:mod:`matplotlib`.  Like :mod:`desr.numerics` it is deliberately not imported by
``desr/__init__.py``, so :mod:`desr` itself keeps neither dependency.

Data and drawing are separate, as :meth:`invariant_monomials` and :meth:`invariant_names`
are: :func:`symbol_graph` returns a graph to analyse, and :func:`draw_symbol_graph` draws
one.  A graph can be relabelled with :func:`networkx.relabel_nodes`, written out as DOT
with :mod:`networkx.drawing.nx_pydot`, or measured with
:mod:`networkx.algorithms.bipartite`; a drawing can only be looked at.
'''

import networkx
import sympy

from desr.ode_translation import ReductionScheme
from desr.tex_tools import display_tex, expr_to_tex

#: The two vertex layers.  A node is a ``(layer, symbol)`` pair, because the two layers can
#: contain the *same* sympy symbol -- a naming scheme may reuse the original letters, and
#: the general and dependent-variable schemes keep the system's own clock -- and bare
#: symbols as nodes would silently merge those vertices.
ORIGINAL = 'original'
REDUCED = 'reduced'

#: The directions :func:`symbol_graph` accepts.
FORWARD = 'forward'
REVERSE = 'reverse'
DIRECTIONS = (FORWARD, REVERSE)


def _resolved_direction(direction):
    '''
    Check `direction` against :data:`DIRECTIONS`, naming the alternatives if it is wrong.
    '''
    if direction not in DIRECTIONS:
        raise ValueError(
            '{!r} is not a direction: use {}.  The forward graph is the multiplier, '
            'original to reduced; the reverse graph is its inverse.'.format(
                direction, ' or '.join(repr(name) for name in DIRECTIONS)))
    return direction


def _original_kinds(system):
    '''
    What each original symbol is: the independent variable, a dependent, or a parameter.
    '''
    dependents = set(system.dependents)
    kinds = {}
    for symbol in system.symbols:
        if symbol == system.independent_variable:
            kinds[symbol] = 'independent'
        elif symbol in dependents:
            kinds[symbol] = 'dependent'
        else:
            kinds[symbol] = 'parameter'
    return kinds


def _auxiliary_kind(monomial, system):
    '''
    Whether an auxiliary sits still or has to be integrated.

    Hubert & Labahn's auxiliary system is :math:`dx/dt = x \\star F(y^{W_\\partial})
    \\cdot V_i`, and they call it a quadrature: given the reduced system's solution, each
    auxiliary follows by one integration.  Its integrand is
    :math:`(z^{-1} \\star dz/dt) \\cdot V_i` -- the logarithmic derivatives of the original
    symbols, combined by a column of :math:`V_i` -- and a **parameter contributes zero to
    it**, having no derivative.

    So an auxiliary that is a monomial in the parameters alone has
    :math:`d(\\log x)/dt = 0`: it is constant, and recovering the original system needs a
    *constant of integration* rather than an integral.  That is the whole of why the
    parameter scheme can normalise its auxiliaries to 1 and drop them: its block form puts
    every column of :math:`V_i` on parameter rows, so every one of them is such a monomial.
    The paper states the conclusion -- constants there, quadratures in general -- without
    walking this step.

    Read per auxiliary rather than per scheme, which is finer than the paper needs to be
    and is what a reduction like the enzyme chain's calls for: three of its auxiliaries are
    monomials in rate constants and sit still, while the fourth is a species and moves.

    Returns:
        str: ``'parameter'`` where the quadrature is trivial, ``'dependent'`` where it is
        not.  The independent variable counts as something that moves, since a monomial
        containing it is not constant either.
    '''
    moving = set(system.variables) | {system.independent_variable}
    return 'dependent' if monomial.free_symbols & moving else 'parameter'


def _reduced_kinds(translator, system, monomials):
    '''
    What each column is, in column order.

    The invariants take the **role** that
    :meth:`~desr.ode_translation.ODETranslator.symbol_kinds` gives them in the reduced
    system -- not the structure of their monomial, because the two differ and the role is
    what a reader of the figure is after.  Under the general and dependent-variable schemes
    the reduction keeps the system's own independent variable, so an invariant that is a
    rescaled one -- the enzyme chain's :math:`k_{23} t`, Hubert & Labahn 6.4's :math:`t` --
    is a variable with an equation rather than the thing being integrated against.  Calling
    it independent would say the reduced system has two independent variables.

    The auxiliaries cannot come from there: :meth:`symbol_kinds` reports every one of them
    as ``'auxiliary'``, which is a fact about the scaling rather than about what the symbol
    does, and under the parameter scheme it does not report them at all.  Being an
    auxiliary is orthogonal to being a parameter or a dependent -- the reduced system says
    so itself, listing the enzyme chain's ``x0``, ``x1`` and ``x2`` among its parameters --
    so the two are carried separately here and :func:`_auxiliary_kind` decides the second.
    '''
    roles = translator.column_kinds(system)
    return [_auxiliary_kind(monomials[column], system) if column < translator.r
            else role
            for column, role in enumerate(roles)]


def _propagated_values(monomials, names, values):
    '''
    The values the known original ones force on the reduced symbols.

    One substitution is enough, and is **complete** rather than best-effort: no reduced
    symbol is a monomial in the other reduced symbols, because the columns of the
    multiplier are multiplicatively independent -- the multiplier is invertible -- so there
    is no redundancy on the reduced side to iterate over.

    A monomial need not have all of its symbols known to have a known value.  A zero in a
    numerator settles it whatever the rest are, which is the case that matters: `c_0 / s_0`
    is zero for every `s_0` once `c_0` is.  A zero in a *denominator* settles nothing, and
    is excluded by its value being non-finite rather than by inspecting the exponents --
    though a reduction that divides by a known zero is refused before it gets here.
    '''
    if not values:
        return {}

    propagated = {}
    for name, monomial in zip(names, monomials):
        substituted = monomial.subs(values)
        if substituted.is_number and substituted.is_finite:
            propagated[name] = substituted
    return propagated


def symbol_graph(translator, system, direction=FORWARD):
    '''
    The exponent graph of `translator` acting on `system`.

    Args:
        translator (desr.ode_translation.ODETranslator): The reduction to read.  Its
            multiplier is read as it currently stands, so a
            :meth:`~desr.ode_translation.ODETranslator.multiplier_add_columns` changes the
            graph -- which is the point, since the graph is what a column operation is for
            looking at.
        system (desr.ode_system.ODESystem): The system being reduced.  Required: the
            reduced symbols' names depend on which scheme runs, and their kinds and values
            on the system itself.
        direction (str): ``'forward'`` for the multiplier, original to reduced, or
            ``'reverse'`` for its inverse.

    **There is no scheme argument, and nothing is remembered.**  Which scheme applies is a
    feature of the multiplier:
    :meth:`~desr.ode_translation.ODETranslator._is_translate_parameter_compatible` is a
    condition on :math:`V_i`, :math:`V_n` and :math:`V^{-1}`, so a column operation can
    change it.  Michaelis-Menten starts on the parameter scheme and multiplying its clock
    invariant by a dependent one moves it to the general scheme.  A scheme recorded when
    something was last translated would therefore be stale as soon as a column moved, and a
    scheme passed in by the caller would be a claim this cannot check.  Asked afresh, the
    graph always describes the reduction this multiplier would give now, which is the only
    reading that cannot go out of date.

    Returns:
        networkx.DiGraph: Nodes are ``(layer, symbol)`` pairs over both layers, so the
        graph holds all :math:`2n` vertices whichever direction is asked for, and only the
        edges differ.

    Node attributes:
        ``layer``
            :data:`ORIGINAL` or :data:`REDUCED`.
        ``bipartite``
            ``0`` or ``1``, for :mod:`networkx.algorithms.bipartite`.
        ``symbol``
            the :class:`sympy.Symbol` itself, the node key's other half.
        ``index``
            its position within its layer: the variable order on the original side, the
            column order on the reduced side.  What a deterministic layout is drawn from.
        ``kind``
            ``'independent'``, ``'dependent'`` or ``'parameter'`` on the original side; the
            the same three on the reduced side.  **Not** a place where ``'auxiliary'``
            appears: being an auxiliary is orthogonal to this, and has its own flag.
        ``auxiliary``
            whether the scaling moves this symbol -- whether its exponent vector is a
            column of :math:`V_i` rather than of :math:`V_n`.  That is a fact about the
            group action and says nothing about whether the symbol evolves, which is what
            ``kind`` is for.  The two together are what a reader needs: a monomial in the
            parameters alone is constant whether or not the scaling moves it, so an
            auxiliary with ``kind`` ``'parameter'`` yields a constant of integration where
            one with ``kind`` ``'dependent'`` yields a quadrature.
        ``normalised``
            whether the reduction sets this one to 1 and discards it.  True only for the
            auxiliary columns, and only under the parameter scheme -- the general and
            dependent-variable schemes keep their auxiliaries as variables.
        ``must_be_supplied``
            whether this is one of the :math:`r` values the **reverse** direction has to be
            told.  A discarded auxiliary is an unknown coming back, and recovering the
            original values needs :math:`r` of them from somewhere -- which is what
            :class:`desr.numerics.InsufficientKnownValues` is raised about.

            A normalised auxiliary never carries a ``value``, in either direction.  The
            reduction reaches 1 by choosing a group element, which rescales the original
            symbols, so that 1 holds of the rescaled coordinates and not of the caller's:
            written beside the monomial it would read as a constraint on symbols the
            caller chose freely.
        ``must_not_be_zero``
            whether this symbol is one the translation in **this** direction divides by.
        ``monomial``
            what this symbol equals in terms of the **other** layer: a reduced symbol as a
            monomial in the original ones, an original symbol as a monomial in the reduced
            ones.  Present on both layers, and symmetric, so whichever direction is drawn
            the arrows point at nodes that say what they stand for.
        ``value``
            present only when known: given for the original symbols, derived for the
            reduced ones.

    Edge attributes:
        ``exponent``
            the exponent with which the source appears in the target.  Not ``weight``,
            which :mod:`networkx` reads as a cost by default -- a negative exponent is not
            one, and a shortest path over these numbers would mean nothing.

    Michaelis-Menten, the two-variable form with both initial conditions.  Twenty vertices
    and twenty-two edges, every exponent :math:`\\pm 1`:

    >>> from desr import ODESystem, ODETranslator
    >>> from desr.symbol_graph import symbol_graph
    >>> equations = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...              'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
    >>> system = ODESystem.from_equations(equations)
    >>> translator = ODETranslator.from_ode_system(system)
    >>> graph = symbol_graph(translator, system)
    >>> graph.number_of_nodes(), graph.number_of_edges()
    (14, 16)

    The reduction divides by what has an edge into the auxiliary block or a negative edge
    into an invariant, which is what the node attribute says:

    >>> sorted(str(symbol) for (layer, symbol), data in graph.nodes(data=True)
    ...        if layer == 'original' and data['must_not_be_zero'])
    ['e_0', 'k_1']
    >>> translator.symbols_that_must_not_be_zero()
    (e_0, k_1)

    An edge is one exponent, and the invariant it points into is the monomial it belongs
    to:

    >>> graph[('original', sympy.Symbol('k_2'))][('reduced', sympy.Symbol('kappa0'))]
    {'exponent': 1}
    >>> graph.nodes[('reduced', sympy.Symbol('kappa0'))]['monomial']
    k_2/(e_0*k_1)

    The reverse graph has the same vertices and different edges, and names a different set
    of symbols -- here both auxiliaries, which the parameter scheme sets to 1:

    >>> back = symbol_graph(translator, system, direction='reverse')
    >>> back.number_of_nodes()
    14
    >>> sorted(str(symbol) for (layer, symbol), data in back.nodes(data=True)
    ...        if layer == 'reduced' and data['must_not_be_zero'])
    ['x0', 'x1']

    Composing the two is the identity, so a two-hop round trip from a symbol to itself has
    exponent products summing to 1, and to 0 for any other pair:

    >>> k_1 = ('original', sympy.Symbol('k_1'))
    >>> def round_trip(start, finish):
    ...     return sum(graph[start][middle]['exponent'] * back[middle][finish]['exponent']
    ...                for middle in graph[start] if finish in back[middle])
    >>> round_trip(k_1, k_1)
    1
    >>> round_trip(k_1, ('original', sympy.Symbol('k_2')))
    0
    '''
    direction = _resolved_direction(direction)

    originals = list(translator.symbols_domain)
    names = list(translator.column_names(system))
    monomials = list(translator.column_monomials())

    original_kinds = _original_kinds(system)
    reduced_kinds = _reduced_kinds(translator, system, monomials)

    # What each original symbol equals in terms of the reduced ones, the mirror of
    # `column_monomials`.  Taken from `reverse_translate_parameter` -- which is
    # `scale_action(variables, inv_herm_mult)` and nothing else, so it applies whatever
    # scheme runs -- rather than indexing the inverse again here.  One place decides that
    # orientation, and it is the one that was always right.
    recoveries = list(translator.reverse_translate_parameter(names))

    original_values = dict(system.parameter_values)
    reduced_values = _propagated_values(monomials, names, original_values)

    # Which auxiliary columns the reduction normalises, **one at a time**.  It is not an
    # all-or-nothing property of the scheme: an auxiliary is normalised exactly when it is
    # a monomial in the original parameters alone, because only then is it constant, and a
    # reduction can have some such auxiliaries and some that move.  The parameter scheme is
    # the case where every one of them is, which is why asking per column covers it too
    # without a special case.
    normalised_monomials = set(translator.parameter_auxiliaries(system))

    # Which of the caller's own symbols that spends.  Asked of `normalisations`, which is
    # what the reduction records as having done, so the figure cannot disagree with it
    # about where the 1s land.
    #
    # The reduction adopts those symbols as its units -- `s_0 |-> 1` is "measure
    # concentration in units of s_0", which is why `u = s/s_0` has no `s_0` left in it.
    sent_to_one = set(translator.normalisations(system))

    # Note that this, `column_names` and `reduced_symbols_that_must_not_be_zero` are all
    # asked afresh, and none of them is told a scheme.  See the class docstring's note on
    # why there is no `scheme` argument.

    divides_by = {
        ORIGINAL: set(translator.symbols_that_must_not_be_zero(system)),
        REDUCED: set(translator.reduced_symbols_that_must_not_be_zero(system)),
    }

    graph = networkx.DiGraph()
    graph.graph['direction'] = direction

    for index, symbol in enumerate(originals):
        attributes = {
            'layer': ORIGINAL,
            'bipartite': 0,
            'symbol': symbol,
            'index': index,
            'kind': original_kinds[symbol],
            'monomial': recoveries[index],
            'auxiliary': False,
            'normalised': symbol in sent_to_one,
            'must_be_supplied': False,
            'must_not_be_zero': (direction == FORWARD
                                 and symbol in divides_by[ORIGINAL]),
        }
        if symbol in original_values:
            attributes['value'] = original_values[symbol]
        graph.add_node((ORIGINAL, symbol), **attributes)

    for index, symbol in enumerate(names):
        auxiliary = index < translator.r
        normalised = auxiliary and monomials[index] in normalised_monomials
        attributes = {
            'layer': REDUCED,
            'bipartite': 1,
            'symbol': symbol,
            'index': index,
            'kind': reduced_kinds[index],
            'monomial': monomials[index],
            'auxiliary': auxiliary,
            'normalised': normalised,
            'must_not_be_zero': (direction == REVERSE
                                 and symbol in divides_by[REDUCED]),
        }
        attributes['must_be_supplied'] = normalised and direction == REVERSE

        # A normalised auxiliary carries no `value`, in either direction.  Writing 1 on it
        # would be reading one coordinate system in another: the reduction reaches that 1
        # by choosing a group element, which rescales the original symbols, so the
        # statement holds of the rescaled coordinates and not of the caller's.  Labelled
        # `x_0 = 1/(k_1 s_0) = 1` it reads as a constraint on `k_1` and `s_0`, which are
        # the caller's to choose and are under no such constraint.  What is true of it --
        # that the reduction is free to normalise it away -- is said by `normalised`, and
        # shown by the hatch rather than by an equation.
        if not normalised and symbol in reduced_values:
            attributes['value'] = reduced_values[symbol]
        graph.add_node((REDUCED, symbol), **attributes)

    # The multiplier is indexed [original, reduced] and its inverse [reduced, original], so
    # each puts its own direction's source index first and neither needs transposing.
    if direction == FORWARD:
        exponents = translator.herm_mult
        sources, targets = (ORIGINAL, originals), (REDUCED, names)
    else:
        exponents = translator.inv_herm_mult
        sources, targets = (REDUCED, names), (ORIGINAL, originals)

    source_layer, source_symbols = sources
    target_layer, target_symbols = targets
    for row, source in enumerate(source_symbols):
        for column, target in enumerate(target_symbols):
            exponent = exponents[row, column]
            if exponent != 0:
                graph.add_edge((source_layer, source), (target_layer, target),
                               exponent=int(exponent))

    return graph


def _positions(graph, width=1.0, height=1.0):
    '''
    A fixed two-column layout, read off the ``index`` attributes.

    **The source layer is on the left, whichever direction this is**, so the arrows always
    run left to right.  Draw the forward graph beside the reverse one and the reduced
    column ends the first panel and begins the second, in the same vertical order, so a
    round trip reads as one continuous path across the pair: original, reduced, original.

    Deliberately not a force-directed layout.  Two of these get drawn to be compared --
    the pair above, or before a column operation against after -- and a layout that moves
    a vertex between the panels destroys exactly that.  Here a symbol's place depends only
    on the variable and column orders.
    '''
    by_layer = {ORIGINAL: [], REDUCED: []}
    for node, data in graph.nodes(data=True):
        by_layer[data['layer']].append((data['index'], node))

    if graph.graph.get('direction', FORWARD) == FORWARD:
        columns = ((ORIGINAL, 0.0), (REDUCED, width))
    else:
        columns = ((REDUCED, 0.0), (ORIGINAL, width))

    positions = {}
    for layer, column in columns:
        nodes = sorted(by_layer[layer])
        count = len(nodes)
        for place, (_, node) in enumerate(nodes):
            offset = 0.5 if count == 1 else place / (count - 1)
            positions[node] = (column, height * (1.0 - offset))
    return positions


def _node_label(data, with_monomial=False, latex=True, show_normalisation=True):
    '''
    The symbol; what it equals, where the arrows point at it; its value when one is known,
    so a zero is on the figure rather than in a caption; and an arrow where the reduction
    sets it to 1.

    The chain before the arrow is joined with ``=`` because every part of it genuinely is
    equal: :math:`\\kappa_4 = c_0/s_0 = 0` says what the invariant is and what it comes to,
    which is a zero initial condition carried through a reduction in one line.

    **The arrow is not an equals sign**, and the distinction is the point.  A value is a
    fact about the system the caller wrote down.  A :math:`\\mapsto 1` is something the
    reduction *does*: it chooses the group element that normalises the auxiliaries, and
    that choice rescales the original symbols.  So :math:`s_0 \\mapsto 1` says the reduced
    system is written in coordinates where :math:`s_0` is 1 -- that the reduction measures
    concentration in units of :math:`s_0` -- and not that :math:`s_0` was ever 1.  Written
    with ``=`` beside a monomial in the caller's own parameters it would read as a
    constraint on parameters that are theirs to choose.

    `show_normalisation` is off coming back, where the arrow would be describing the wrong
    direction: the reduction normalises on the way out, and on the way back those same
    auxiliaries are the values a recovery has to be told.

    Set in mathtext by default, which :mod:`matplotlib` renders itself -- no LaTeX
    installation -- so ``kappa0`` prints as :math:`\\kappa_0` and a quotient as a built-up
    fraction rather than with a slash.
    '''
    parts = [data['symbol']]
    if with_monomial:
        parts.append(data['monomial'])
    if 'value' in data:
        parts.append(data['value'])

    normalised = data['normalised'] and show_normalisation
    if latex:
        # expr_to_tex rather than sympy.latex: desr writes a negative subscript `k_m1`,
        # which sympy sets as a literal m.  display_tex because this is a figure, where
        # an inline fraction comes out half height.
        chain = ' = '.join(display_tex(expr_to_tex(part)) for part in parts)
        return '${}{}$'.format(chain, r' \mapsto 1' if normalised else '')
    chain = ' = '.join(str(part) for part in parts)
    return chain + (' |-> 1' if normalised else '')


#: Where along an edge its exponent is written.  Every edge of a two-column layout has its
#: midpoint at the same `x`, so a single label position stacks all of them into one
#: illegible column; cycling through a few spreads them out.  Deterministic, so the same
#: edge is labelled in the same place in both panels of a comparison.
_LABEL_POSITIONS = (0.30, 0.45, 0.60, 0.75)

#: Node fills, named once so the key and the figure cannot drift apart.  The auxiliary
#: grey is dark enough to read as a different thing beside the other two rather than as a
#: paler version of them, and it is carried by the label as well as the marker: what the
#: grey marks is a whole quantity, not a dot.
_ORIGINAL_FILL = '#dce6f1'
_INVARIANT_FILL = '#e8e0ef'
_AUXILIARY_FILL = '#dcdcdc'
_BORDER = '#333333'

#: The auxiliaries are hatched as well as filled, for the same reason the negative edges
#: are dashed: a pattern survives greyscale printing and colour blindness, where a fill
#: alone can come down to "slightly different grey".  The hatch takes its colour from the
#: edge, so it stays legible however pale the fill is.
_AUXILIARY_HATCH = '///'

#: Marker shapes, by what the symbol **is**.  Shape is free to carry this because the
#: layer is already said by which column a vertex sits in -- and said reliably, since the
#: source layer is always on the left.  So shape and fill carry the two independent facts:
#: a square is a parameter and a circle a dependent variable whichever side it is on.
#:
#: Being an auxiliary is **not** one of these.  An auxiliary is a parameter or a dependent
#: like any other symbol, so it keeps its shape and the hatch says the rest: a monomial in
#: the parameters alone sits still whether or not the scaling moves it.
_KIND_SHAPES = (('parameter', 's'),
                ('dependent', 'o'),
                ('independent', 'D'))

#: How to say each kind in a key.
_KIND_LABELS = {'parameter': 'parameter',
                'dependent': 'dependent variable',
                'independent': 'independent variable'}

#: Behind the label it is the fill alone, no hatch.  The hatch has to be dense to read as
#: a pattern at marker size, and at that density it competes with the glyphs of a label
#: like :math:`x_0 = 1/(k_1 s_0) = 1`.  The marker carries the pattern, the label carries
#: the colour, and the two together still mark the whole quantity rather than a dot.
_AUXILIARY_LABEL_HATCH = None

#: How faint a node drawn as removed is.  Pale enough to read as gone, dark enough that a
#: reader can still see *what* went -- the count of what the normalisation spends is part
#: of what the figure is for.
_GONE_ALPHA = 0.22


def _swatch_side(legend, fontsize):
    '''
    The size in pixels of a marker of `fontsize` points, which is what a swatch is drawn
    at so that it matches the plain markers beside it.
    '''
    return fontsize * legend.figure.dpi / 72.0


def _square_swatch(legend, orig_handle, xdescent, ydescent, width, height, fontsize):
    '''A square swatch, matching the marker a parameter is drawn with.'''
    import matplotlib.patches

    side = _swatch_side(legend, fontsize)
    return matplotlib.patches.Rectangle(
        (width / 2.0 - xdescent - side / 2.0, height / 2.0 - ydescent - side / 2.0),
        side, side,
        facecolor=orig_handle.get_facecolor(), edgecolor=orig_handle.get_edgecolor(),
        linewidth=orig_handle.get_linewidth(), hatch=orig_handle.get_hatch())


def _circle_swatch(legend, orig_handle, xdescent, ydescent, width, height, fontsize):
    '''A round swatch, matching the marker a dependent variable is drawn with.'''
    import matplotlib.patches

    return matplotlib.patches.Circle(
        (width / 2.0 - xdescent, height / 2.0 - ydescent),
        radius=0.5 * _swatch_side(legend, fontsize),
        facecolor=orig_handle.get_facecolor(), edgecolor=orig_handle.get_edgecolor(),
        linewidth=orig_handle.get_linewidth(), hatch=orig_handle.get_hatch())


#: A key entry for an auxiliary has to be hatched, and a marker cannot be hatched, so it is
#: a patch -- and :mod:`matplotlib` draws every patch handle as a wide rectangle, which
#: would say nothing about which shape the entry is describing.  These draw the shape the
#: figure actually uses.
_KIND_SWATCHES = {'parameter': _square_swatch, 'dependent': _circle_swatch}


def draw_symbol_graph(graph, ax=None, label_edges=True, label_expressions=True,
                      latex=True, legend=True, width=1.0, height=1.0, font_size=10,
                      node_size=320, faded=None, show_normalisation=None,
                      legend_kwargs=None):
    '''
    Draw an exponent graph.

    Positive exponents are solid and negative ones dashed, which survives greyscale
    printing and colour blindness; colour reinforces that rather than carrying it.  A
    symbol the translation divides by is drawn with a heavy border, and a known value is
    written into the label, so the figure needs no key to the zero propagation.

    **Shape says what a symbol is** -- a square is a parameter, a circle a dependent
    variable, a diamond the independent one -- and fill says which layer it belongs to.
    Shape is free to carry the kind because the layer is already given by which column a
    vertex sits in, and given reliably, since the source layer is always on the left.

    **Being an auxiliary is a hatch, not a shape**, because it is the other axis: an
    auxiliary is a parameter or a dependent like any other symbol, and giving it a shape of
    its own would put two independent questions on one channel and leave no way to say
    which kind it is.

    Args:
        graph (networkx.DiGraph): From :func:`symbol_graph`.
        ax (matplotlib.axes.Axes, optional): Where to draw.  A new figure is made when this
            is omitted.  Pass one per panel to put two graphs side by side -- the layout is
            deterministic, so a symbol sits in the same place in both.
        label_edges (bool): Write the exponent on each edge.
        label_expressions (bool): On the layer the arrows point at, write what each symbol
            equals in terms of the other layer.  So the forward graph says what each
            invariant *is*, and the reverse one says how each original symbol comes back.
        latex (bool): Set the labels in mathtext, so ``kappa0`` reads as
            :math:`\\kappa_0` and a quotient is a built-up fraction.  :mod:`matplotlib`
            renders this itself and no LaTeX installation is needed.
        legend (bool): A key to solid against dashed, and to the heavy border, which is
            only described when something carries it.
        width, height (float): The box the layout is drawn in.
        font_size (int): For the node and edge labels.
        node_size (int): Passed to :mod:`networkx`.
        faded (set, optional): Nodes to draw as removed -- pale, and every edge touching
            one pale with them.  They keep their places, which is the point: drawing the
            same graph twice with a set faded the second time shows what deleting it
            costs, and nothing else moves, so the two panels can be read against each
            other.
        show_normalisation (bool, optional): Whether a symbol the reduction sets to 1 says
            so with :math:`\\mapsto 1`.  Defaults to the forward direction alone, where the
            normalisation is performed; coming back, those same auxiliaries are the values
            a recovery has to be told, and an arrow there describes the wrong direction.
            Set it explicitly to draw the normalisation against a *reverse* graph, which
            is where the substitution can be watched being made.
        legend_kwargs (dict, optional): Passed to :meth:`~matplotlib.axes.Axes.legend`,
            over the defaults.  The key hangs under the axes, which is where the room is
            for a graph drawn on its own; a caller building a figure from several axes
            has the room somewhere else and places it with ``bbox_to_anchor`` and
            ``bbox_transform``, usually the figure's.

    Returns:
        matplotlib.axes.Axes: The axes drawn on.

    >>> import matplotlib
    >>> matplotlib.use('Agg')
    >>> from desr import ODESystem, ODETranslator
    >>> from desr.symbol_graph import symbol_graph, draw_symbol_graph
    >>> system = ODESystem.from_equations(['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...                                    'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
    >>> translator = ODETranslator.from_ode_system(system)
    >>> axes = draw_symbol_graph(symbol_graph(translator, system))
    >>> axes.get_xlim() is not None
    True
    '''
    import matplotlib.legend_handler
    import matplotlib.lines
    import matplotlib.patches
    import matplotlib.pyplot

    if ax is None:
        _, ax = matplotlib.pyplot.subplots(figsize=(7.0, 7.0))

    positions = _positions(graph, width=width, height=height)

    # Shape tracks the layer, so it still says which side a vertex belongs to when the
    # reverse graph puts the reduced layer on the left.  Grey is the auxiliary block,
    # which is the part of the reduced side that is not an invariant.
    faded = set() if faded is None else set(faded)
    normalisation = (graph.graph.get('direction', FORWARD) == FORWARD
                     if show_normalisation is None else show_normalisation)

    for layer, fill in ((ORIGINAL, _ORIGINAL_FILL), (REDUCED, _INVARIANT_FILL)):
        for kind, shape in _KIND_SHAPES:
            for auxiliary in (False, True):
                for gone in (False, True):
                    nodes = [node for node, data in graph.nodes(data=True)
                             if data['layer'] == layer and data['kind'] == kind
                             and data['auxiliary'] == auxiliary
                             and (node in faded) == gone]
                    if not nodes:
                        continue
                    borders = [2.5 if graph.nodes[node]['must_not_be_zero'] else 0.8
                               for node in nodes]
                    drawn = networkx.draw_networkx_nodes(
                        graph, positions, nodelist=nodes, ax=ax, node_shape=shape,
                        node_color=_AUXILIARY_FILL if auxiliary else fill,
                        node_size=node_size, edgecolors=_BORDER, linewidths=borders,
                        alpha=_GONE_ALPHA if gone else None)
                    if auxiliary:
                        drawn.set_hatch(_AUXILIARY_HATCH)

    for style, colour, keep in (('solid', '#1f3b63', lambda e: e > 0),
                                ('dashed', '#8c2f39', lambda e: e < 0)):
        for gone in (False, True):
            edges = [(tail, head) for tail, head, exponent
                     in graph.edges(data='exponent')
                     if keep(exponent)
                     and bool({tail, head} & faded) == gone]
            if not edges:
                continue
            networkx.draw_networkx_edges(graph, positions, edgelist=edges, ax=ax,
                                         style='dotted' if gone else style,
                                         edge_color=colour, node_size=node_size,
                                         arrowsize=12,
                                         alpha=_GONE_ALPHA if gone else None)

    # Labels go *outside* the columns, not on the markers: a marker big enough to hold
    # `kappa0` or `c_0 = 0` would be a box of text rather than a vertex, and the middle of
    # the figure is needed for the edges.
    padding = 0.11 * width
    target_layer = REDUCED if graph.graph.get('direction', FORWARD) == FORWARD else ORIGINAL
    for node, data in graph.nodes(data=True):
        x, y = positions[node]
        on_the_left = x < 0.5 * width
        # An auxiliary's label is set in the same grey as its marker, so the whole
        # quantity reads as the one thing the reduction is about to absorb.
        box = ({'boxstyle': 'round,pad=0.3', 'facecolor': _AUXILIARY_FILL,
                'edgecolor': _BORDER, 'linewidth': 0.8,
                'hatch': _AUXILIARY_LABEL_HATCH}
               if data['auxiliary'] else None)
        ax.text(x + (-padding if on_the_left else padding), y,
                _node_label(data,
                            with_monomial=label_expressions
                            and data['layer'] == target_layer,
                            latex=latex, show_normalisation=normalisation),
                fontsize=font_size, ha='right' if on_the_left else 'left', va='center',
                bbox=box, alpha=_GONE_ALPHA if node in faded else None)

    if label_edges:
        grouped = {}
        for tail, head, exponent in graph.edges(data='exponent'):
            # An edge drawn as removed keeps its line, faintly, so the reader can see what
            # went; its exponent is dropped, because a number is read before it is seen to
            # be pale and the panel fills up with arithmetic that no longer applies.
            if {tail, head} & faded:
                continue
            spread = (graph.nodes[tail]['index'] + graph.nodes[head]['index'])
            where = _LABEL_POSITIONS[spread % len(_LABEL_POSITIONS)]
            grouped.setdefault(where, {})[(tail, head)] = (
                '${}$'.format(exponent) if latex else str(exponent))
        for where, labels in sorted(grouped.items()):
            networkx.draw_networkx_edge_labels(
                graph, positions, ax=ax, label_pos=where,
                font_size=max(font_size - 2, 5), edge_labels=labels,
                bbox={'boxstyle': 'round,pad=0.15', 'facecolor': 'white',
                      'edgecolor': 'none', 'alpha': 0.85})

    if legend:
        handler_map = {tuple: matplotlib.legend_handler.HandlerTuple(
            ndivide=None, pad=1.1)}

        # With the text.  A marker set in points while the label shrinks with the figure
        # is a key whose shapes grow as the words get smaller, and a diamond -- drawn to
        # the same width as a circle and so taller -- shows it first.
        key_font = max(font_size - 1, 5)

        def marker(shape, face, border, scale=1.0):
            return matplotlib.lines.Line2D(
                [], [], linestyle='none', marker=shape,
                markersize=scale * key_font, markerfacecolor=face,
                markeredgecolor=_BORDER, markeredgewidth=border)

        handles = [matplotlib.lines.Line2D([], [], color='#1f3b63', linestyle='solid'),
                   matplotlib.lines.Line2D([], [], color='#8c2f39', linestyle='dashed')]
        labels = ['positive exponent', 'negative exponent']

        # Each of the rest only when something actually carries it, so the key never
        # describes an absent thing.  A graph can have nothing marked -- the reverse
        # direction of a reduction that divides by no reduced symbol, for one.
        # What a symbol is, in the shapes actually on the figure.  A neutral fill, since
        # the fill says which layer and these entries are about kind.
        present = {data['kind'] for _, data in graph.nodes(data=True)}
        for kind, shape in _KIND_SHAPES:
            if kind in _KIND_LABELS and kind in present:
                handles.append(marker(shape, 'white', 0.8))
                labels.append(_KIND_LABELS[kind])

        # Here rather than at the end, because it is another mark on those same shapes
        # and because it leaves the two auxiliary entries adjacent -- and a legend set
        # in columns fills them top to bottom, so adjacent is what puts the pair in a
        # column of its own.
        if any(data['must_not_be_zero'] for _, data in graph.nodes(data=True)):
            # Both shapes under the one label.  Either layer can carry the mark depending
            # on the direction drawn, so showing a single shape would read as a claim
            # about which layer it applies to.
            #
            # Smaller than the rest: two markers have to share one handle's width, and
            # at full size they touch.
            handles.append((marker('o', 'white', 2.0, scale=0.78),
                            marker('s', 'white', 2.0, scale=0.78)))
            labels.append('must not be zero')

        auxiliaries = [data for _, data in graph.nodes(data=True) if data['auxiliary']]
        if auxiliaries:
            def hatched(kind, text):
                '''One key entry, drawn in the shape the figure gives that kind.'''
                handle = matplotlib.patches.Patch(
                    facecolor=_AUXILIARY_FILL, edgecolor=_BORDER, linewidth=0.8,
                    hatch=_AUXILIARY_HATCH)
                handler_map[handle] = matplotlib.legend_handler.HandlerPatch(
                    patch_func=_KIND_SWATCHES[kind])
                handles.append(handle)
                labels.append('auxiliary, {}'.format(text))

            # Named by what each one costs to recover, which is Hubert & Labahn's own
            # division and holds under every scheme: a monomial in the parameters alone is
            # constant and supplies a constant of integration, while one that moves has to
            # be integrated.  Only the cases present are named.
            #
            # What then *becomes* of them -- normalised away, or supplied coming back --
            # is a property of the scheme and the direction rather than of the symbol, and
            # belongs in prose.  A key is for reading the marks.
            for kind, word in (('parameter', 'a constant of integration'),
                               ('dependent', 'a quadrature')):
                if any(data['kind'] == kind for data in auxiliaries):
                    hatched(kind, word)

        # One column.  Two splits the three kinds across the break -- and where they land
        # depends on how many entries there are, which varies with the figure -- so the
        # group that most wants reading together is the group that gets broken up.
        #
        # Under the axes, which is where the room is for a graph drawn on its own.  A
        # caller composing a figure out of several axes has the room somewhere else, and
        # says so with `legend_kwargs`.
        settings = dict(loc='upper center', bbox_to_anchor=(0.5, 0.0), frameon=False,
                        ncol=1, fontsize=key_font, handleheight=1.0,
                        handletextpad=0.7)
        settings.update(legend_kwargs or {})
        ax.legend(handles=handles, labels=labels, handler_map=handler_map, **settings)

    ax.set_axis_off()
    # Room for the labels outside the columns.  Text does not take part in autoscaling, so
    # a long symbol can still reach past this -- save with ``bbox_inches='tight'`` and it
    # will not be cut off.
    ax.set_xlim(-0.45 * width, 1.45 * width)
    ax.set_ylim(-0.12 * height, 1.12 * height)
    return ax


def node_positions(graph, width=1.0, height=1.0):
    """
    Where :func:`draw_symbol_graph` puts each vertex, as ``{node: (x, y)}``.

    For drawing something else beside a panel and keeping it in register.  A figure that
    puts a second picture next to a graph -- a substitution, a table of values -- wants
    its rows at the heights the symbols already have, so the eye can carry a symbol from
    one panel to the next; working those heights out again from the layout rule would be
    a second copy of it, free to drift.

    Args:
        graph (networkx.DiGraph): From :func:`symbol_graph`.
        width, height (float): As passed to :func:`draw_symbol_graph`.

    Returns:
        dict: ``{(layer, symbol): (x, y)}``.
    """
    return _positions(graph, width, height)
