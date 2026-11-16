"""The reduction as a flowchart, with the matrices that carry out each step.

The chain figure in ``reduction_chain.py`` shows five systems and names the map between
each pair.  This one is the other half: the same reduction, but with the matrices
themselves on the page and a step named on every arrow between them.

    system -> dimensionless -> K -> U K = H -> A -> V = [V_i | V_n] -> W -> reduced

**dimensionless.**  Every equation is first multiplied through by :math:`t/x`, which is
the step that makes it dimensionless: :math:`dx/dt` scales the way :math:`x/t` does, so
:math:`(t/x)\\,dx/dt` scales not at all.  A scaling leaves the system alone exactly when
it leaves every term of these expressions alone, and that -- not the original equation --
is what the exponents are taken from.

**K**, the exponent matrix, is one column per term of those expressions, one row per
symbol.  Everything after it is integer linear algebra, and the equations are not
consulted again until the last panel.

**U K = H**, the row Hermite normal form of K.  The zero rows at the bottom of H say
which combinations of K's rows vanish, and those are the scalings; the combinations
themselves are the matching rows of **U**.  H's whole job is to point at them, which is
why both matrices are on the page and only part of one is used.

**A**, the scaling matrix, is those rows put in row Hermite normal form, so the answer
does not depend on how the elimination ran.  Its :math:`r` rows are the scaling
symmetries of the system -- the panel writes each one out as the substitution it is --
and :math:`A K = 0`.

**V**, the Hermite multiplier of A, satisfies :math:`A V = [I_r | 0]` and so splits in
two: :math:`A V_i = I_r`, each of whose columns carries exactly one dimension, and
:math:`A V_n = 0`, whose columns are dimensionless.  Read a column of V down the old
symbols and it is a new symbol -- an auxiliary from :math:`V_i`, an invariant from
:math:`V_n`.

**W**, its inverse, is read down its columns too, and gives each old symbol back as a
monomial in the new ones.  Its top block is the scaling matrix again, entry for entry:
:math:`A V = [I_r | 0]` forces :math:`W_u = A`.  The normalisation is in this panel and
needs no solving -- set the auxiliary columns to 1 and read off what is left.

The last panel is the reduced system, where the equations come back.

Run with::

    python examples/reduction_flowchart.py
"""

import re

import matplotlib
import matplotlib.pyplot as plt
import sympy
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, PathPatch
from matplotlib.path import Path

from desr import ODESystem, ODETranslator
from desr.chemical_reaction_network import ChemicalReactionNetwork
# expr_to_tex rather than sympy.latex: desr writes a negative subscript `k_m1`, and
# sympy sets that as a literal m.
from desr.tex_tools import dimension_to_tex, display_tex, expr_to_tex


#: Colour means the same thing here as in the chain figure: **rust** the auxiliaries,
#: which are the units, and **blue** the invariants, which are the values.  A block of a
#: matrix takes the colour of the symbols it makes.
UNIT_COLOUR, VALUE_COLOUR = '#8a4500', '#1f5c86'
PLAIN, QUIET = 'black', '#555555'

#: The wash behind a highlighted row, and the panel fills.
PICKED = '#f0e2d2'
PANEL_FILL, PANEL_EDGE = '#fbfaf8', '#d9d6d0'

#: Text metrics, as multiples of the font size.  Layout happens before anything is
#: drawn, so these estimate what the renderer will do; both run a little wide, which
#: leaves the panels loose rather than clipped.
CHAR_W, LINE_H = 0.58, 1.55

#: Computer Modern, which matplotlib ships, so the figure sets in the same face as the
#: paper it belongs in.  The sans default also makes a :math:`\\kappa` hard to tell
#: from a :math:`k`, which this figure has both of.  ``axes.unicode_minus`` is off
#: because cmr10 has no U+2212 and matplotlib says so, loudly, once per label.
FONT_RC = {
    'font.family': 'serif',
    'font.serif': ['cmr10', 'DejaVu Serif'],
    'mathtext.fontset': 'cm',
    'mathtext.rm': 'serif',
    'axes.unicode_minus': False,
    # This figure has no tick labels at all, but matplotlib warns once per axes
    # unless the formatter is told it may set numbers as mathtext.
    'axes.formatter.use_mathtext': True,
}

#: Inches.  The width of the text block on a page -- pass `width` for a journal's own.
PAGE_IN = 6.5

#: Inches.  How wide a note or a title is allowed to run before it wraps.  Keeping the
#: words narrow is what lets two panels share a row, which halves the height.
TEXT_IN = 2.3

#: Inches.  Padding inside a panel, the gap between panels, and the bracket on a matrix.
PAD_IN, GAP_IN, BRACKET_IN = 0.09, 0.22, 0.055


def _plain(text):
    """
    A string as the characters it will actually set, for estimating its width.

    Mathtext markup is not drawn: ``$V_i$`` sets two glyphs, not six, and a command like
    ``\\lambda`` sets one.  Estimating from the raw string makes every panel that holds
    mathematics far too wide.
    """
    text = re.sub(r'\\[a-zA-Z]+', 'x', text)
    for mark in ('$', '{', '}', '\\,', '\\;', '^', '_'):
        text = text.replace(mark, '')
    return text


def _width(text, font):
    """How wide a string will come out, in inches."""
    return len(_plain(str(text))) * CHAR_W * font / 72.0


def _line(font):
    """One line of text, in inches."""
    return LINE_H * font / 72.0


def _wrap(text, font, limit=TEXT_IN):
    """`text` broken into lines that each fit `limit` inches."""
    lines, current = [], ''
    for word in text.split():
        candidate = (current + ' ' + word).strip()
        if current and _width(candidate, font) > limit:
            lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(current)
    return lines


def _math(tex):
    """A LaTeX fragment as mathtext, set the way a display rather than a line is."""
    return '$' + display_tex(tex) + '$'


def _tex(expression):
    """A sympy expression as mathtext."""
    return _math(expr_to_tex(expression))


class Matrix(object):
    """
    One integer matrix, with its brackets, its labels and any highlighted rows.

    Row and column labels are the symbols the axis is indexed by, and are what makes a
    matrix readable: K's rows are the system's symbols, V's columns are the reduction's.
    A divider draws the partition the name carries, so that :math:`[\\,V_i \\mid V_n\\,]`
    is a line on the page and not only a line in the caption.
    """

    def __init__(self, matrix, font, row_labels=None, col_labels=None,
                 highlight=(), col_colours=None, row_colours=None, label=None,
                 col_divider=None, row_divider=None, brace=None):
        self.matrix = matrix
        self.font = font
        self.row_labels = list(row_labels) if row_labels else None
        self.col_labels = list(col_labels) if col_labels else None
        self.highlight = set(highlight)
        self.col_colours = col_colours or {}
        self.row_colours = row_colours or {}
        self.label = label
        self.col_divider = col_divider
        self.row_divider = row_divider
        #: Words for a brace gathering the highlighted rows, drawn to their left.
        self.brace = brace

        self.entries = [[str(entry) for entry in row] for row in matrix.tolist()]
        self.rows, self.columns = matrix.shape

        widest = max((_width(entry, font)
                      for row in self.entries for entry in row), default=0.0)
        if self.col_labels:
            widest = max(widest, max(_width(label, font)
                                     for label in self.col_labels))
        self.cell_w = widest + 0.055
        self.cell_h = _line(font)
        self.label_w = (max(_width(label, font) for label in self.row_labels)
                        + 0.05) if self.row_labels else 0.0
        self.label_h = self.cell_h if self.col_labels else 0.0
        self.label_left = (_width(label, font) + 0.07) if label else 0.0
        self.brace_w = ((_width(brace, font) + 0.05 + BRACE_IN * 1.6 + 0.04)
                        if brace else 0.0)

    @property
    def size(self):
        width = (self.brace_w + self.label_left + self.label_w + 2 * BRACKET_IN
                 + self.columns * self.cell_w)
        return width, self.label_h + self.rows * self.cell_h

    def draw(self, ax, x, y):
        """Draw with the top-left of the whole block at (x, y), inches, y downward."""
        left = x + self.brace_w + self.label_left + self.label_w + BRACKET_IN
        top = y + self.label_h
        height = self.rows * self.cell_h
        right = left + self.columns * self.cell_w

        if self.label:
            ax.text(x + self.brace_w, top + height / 2.0, self.label,
                    fontsize=self.font + 0.5, ha='left', va='center', color=PLAIN)

        if self.brace and self.highlight:
            first, last = min(self.highlight), max(self.highlight)
            span = (top + first * self.cell_h, top + (last + 1) * self.cell_h)
            inner = left - BRACKET_IN - 0.03
            _brace(ax, inner, span[0], span[1])
            # Against the point of the brace, not against the panel's edge.
            ax.text(inner - BRACE_IN * 1.6 - 0.03, (span[0] + span[1]) / 2.0,
                    self.brace, fontsize=self.font - 0.5, ha='right', va='center',
                    color=QUIET)

        for index in sorted(self.highlight):
            ax.add_patch(plt.Rectangle(
                (left - BRACKET_IN * 0.6, top + index * self.cell_h),
                self.columns * self.cell_w + 1.2 * BRACKET_IN, self.cell_h,
                facecolor=PICKED, edgecolor='none', zorder=0))

        if self.col_labels:
            for column, label in enumerate(self.col_labels):
                ax.text(left + (column + 0.5) * self.cell_w, y + self.label_h * 0.5,
                        label, fontsize=self.font - 0.5, ha='center', va='center',
                        color=self.col_colours.get(column, QUIET))

        for row in range(self.rows):
            if self.row_labels:
                ax.text(x + self.label_left + self.label_w - 0.05,
                        top + (row + 0.5) * self.cell_h, self.row_labels[row],
                        fontsize=self.font - 0.5, ha='right', va='center',
                        color=self.row_colours.get(row, QUIET))
            for column in range(self.columns):
                ax.text(left + (column + 0.5) * self.cell_w,
                        top + (row + 0.5) * self.cell_h,
                        self.entries[row][column], fontsize=self.font,
                        ha='center', va='center', color=PLAIN)

        if self.col_divider is not None:
            at = left + self.col_divider * self.cell_w
            ax.plot([at, at], [top + 0.012, top + height - 0.012],
                    color=QUIET, linewidth=0.6)
        if self.row_divider is not None:
            at = top + self.row_divider * self.cell_h
            ax.plot([left + 0.012, right - 0.012], [at, at],
                    color=QUIET, linewidth=0.6)

        for edge, inward in ((left, 1), (right, -1)):
            ax.plot([edge, edge], [top, top + height], color=PLAIN, linewidth=0.7)
            for end in (top, top + height):
                ax.plot([edge, edge + inward * BRACKET_IN * 0.7], [end, end],
                        color=PLAIN, linewidth=0.7)


#: Inches.  How tall a column of lines may stand before the list is set in more
#: columns.  A height rather than a count, because a line holding a built fraction is
#: twice the height of a line of symbols and it is the inches that crowd the page.
COLUMN_IN = 1.0


class Lines(object):
    """
    A list of mathtext lines, left aligned, in as many columns as it takes.

    A long list -- every symbol's dimension, every column of V -- is set in columns
    rather than run down the page, because the page has width to spare and the panel
    beside it does not want to be stretched to match.
    """

    def __init__(self, lines, font, colours=None, column_in=COLUMN_IN):
        self.lines = list(lines)
        self.font = font
        self.colours = colours or {}
        tall = sum(self._height_of(drawn) for drawn in self.lines)
        self.columns = max(1, int(-(-tall // column_in))) if self.lines else 1
        self.per_column = -(-len(self.lines) // self.columns) if self.lines else 0

    def _column_of(self, index):
        return index // self.per_column, index % self.per_column

    def _height_of(self, drawn):
        """
        How tall one line stands.

        A built fraction is set a whole line above the baseline and another below it,
        so a line holding one needs about twice the room of a line of symbols, and a
        fraction inside a fraction needs more again.  Giving every line the same height
        is what made them overlap.
        """
        depth = drawn.count('\\dfrac') + drawn.count('\\frac')
        if not depth:
            return _line(self.font)
        return _line(self.font) * (1.95 if depth < 3 else 2.5)

    def _geometry(self):
        """``(column widths, column heights, per-line heights)``."""
        widths = [0.0] * self.columns
        heights = [0.0] * self.columns
        line_heights = []
        for index, drawn in enumerate(self.lines):
            column, _ = self._column_of(index)
            widths[column] = max(widths[column], _width(drawn, self.font))
            height = self._height_of(drawn)
            heights[column] += height
            line_heights.append(height)
        return widths, heights, line_heights

    @property
    def size(self):
        if not self.lines:
            return 0.0, 0.0
        widths, heights, _ = self._geometry()
        return sum(widths) + 0.18 * (self.columns - 1), max(heights)

    def draw(self, ax, x, y):
        widths, _, line_heights = self._geometry()
        starts, cursor = [], x
        for width in widths:
            starts.append(cursor)
            cursor += width + 0.18

        tops = [y] * self.columns
        for index, drawn in enumerate(self.lines):
            column, _ = self._column_of(index)
            height = line_heights[index]
            ax.text(starts[column], tops[column] + height / 2.0, drawn,
                    fontsize=self.font, ha='left', va='center',
                    color=self.colours.get(index, PLAIN))
            tops[column] += height


class Stack(object):
    """
    Pieces one above the next, centred, with a word between them.

    What it is for: a step that turns one matrix into another, where putting the two
    side by side would run off the page.
    """

    def __init__(self, pieces, font, joins=None, gap=0.10):
        self.pieces = list(pieces)
        self.font = font
        self.joins = list(joins) if joins else []
        self.gap = gap

    def _join_height(self, index):
        if index >= len(self.joins) or not self.joins[index]:
            return 0.0
        return _line(self.font) * 1.4

    @property
    def size(self):
        width = max(piece.size[0] for piece in self.pieces)
        # A join's words hang to the right of its arrow, which stands at the centre,
        # so the stack has to be wide enough for half of itself plus the words.
        width = max([width] + [width / 2.0 + 0.06 + _width(join, self.font)
                               for join in self.joins if join])
        height = sum(piece.size[1] for piece in self.pieces)
        height += self.gap * (len(self.pieces) - 1)
        height += sum(self._join_height(index)
                      for index in range(len(self.pieces) - 1))
        return width, height

    def draw(self, ax, x, y):
        width = self.size[0]
        cursor = y
        for index, piece in enumerate(self.pieces):
            piece_w, piece_h = piece.size
            piece.draw(ax, x + (width - piece_w) / 2.0, cursor)
            cursor += piece_h
            join_h = self._join_height(index)
            if join_h:
                # A step inside a panel is a step, so it gets an arrow of its own
                # rather than a word sitting in the gap.
                top = cursor + self.gap / 2.0
                ax.add_patch(FancyArrowPatch(
                    (x + width / 2.0, top), (x + width / 2.0, top + join_h),
                    arrowstyle='-|>', mutation_scale=7, linewidth=0.9, color=QUIET,
                    shrinkA=0, shrinkB=0, zorder=1))
                ax.text(x + width / 2.0 + 0.06, top + join_h / 2.0,
                        self.joins[index], fontsize=self.font - 0.5, ha='left',
                        va='center', color=QUIET)
            if index < len(self.pieces) - 1:
                cursor += self.gap + join_h


class Panel(object):
    """
    One station on the loop: a title, some pieces side by side, and a note underneath.

    The pieces are laid out in a row and centred on each other, so a matrix and the list
    of what its columns mean stand at the same height.
    """

    def __init__(self, title, pieces, font, note=None, piece_gap=0.22, number=None):
        #: Its place in the flow, so that a box can be referred to by a number.
        self.number = number
        self.title = title
        self.pieces = list(pieces)
        self.font = font
        self.note = _wrap(note, font - 0.5) if note else []
        self.piece_gap = piece_gap
        #: Pieces go across by default, which is short and wide.  :func:`_place` turns
        #: this on for a panel whose row would otherwise run off the page.
        self.stacked = False

    @property
    def can_stack(self):
        return len(self.pieces) > 1 and not self.stacked

    def _body(self):
        """The pieces' extent, laid out the way this panel is currently set to."""
        if self.stacked:
            return (max(piece.size[0] for piece in self.pieces),
                    sum(piece.size[1] for piece in self.pieces)
                    + self.piece_gap * 0.5 * (len(self.pieces) - 1))
        return (sum(piece.size[0] for piece in self.pieces)
                + self.piece_gap * (len(self.pieces) - 1),
                max(piece.size[1] for piece in self.pieces))

    @property
    def size(self):
        width, height = self._body()
        width = max(width, _width(self.title, self.font + 0.5))
        if self.note:
            width = max(width, max(_width(line, self.font - 0.5)
                                   for line in self.note))
            height += len(self.note) * _line(self.font - 0.5) + 0.04
        height += _line(self.font + 0.5) + 0.05
        return width + 2 * PAD_IN, height + 2 * PAD_IN

    def draw(self, ax, x, y):
        width, height = self.size
        ax.add_patch(FancyBboxPatch(
            (x + PAD_IN * 0.5, y + PAD_IN * 0.5),
            width - PAD_IN, height - PAD_IN,
            boxstyle='round,pad={},rounding_size=0.05'.format(PAD_IN * 0.5),
            facecolor=PANEL_FILL, edgecolor=PANEL_EDGE, linewidth=0.7, zorder=0,
            clip_on=False))

        title_h = _line(self.font + 0.5)
        ax.text(x + width / 2.0, y + PAD_IN + title_h * 0.5, self.title,
                fontsize=self.font + 0.5, ha='center', va='center', color=PLAIN)

        # In the corner, small and grey, taking no room of its own: it is there to be
        # pointed at, not to be read along with the title.
        if self.number is not None:
            ax.text(x + PAD_IN * 0.75, y + PAD_IN * 0.75, str(self.number),
                    fontsize=self.font - 1, ha='left', va='top', color=QUIET)

        body_top = y + PAD_IN + title_h + 0.05
        body_w, body_h = self._body()
        if self.stacked:
            cursor = body_top
            for piece in self.pieces:
                piece_w, piece_h = piece.size
                piece.draw(ax, x + (width - piece_w) / 2.0, cursor)
                cursor += piece_h + self.piece_gap * 0.5
        else:
            cursor = x + (width - body_w) / 2.0
            for piece in self.pieces:
                piece_w, piece_h = piece.size
                piece.draw(ax, cursor, body_top + (body_h - piece_h) / 2.0)
                cursor += piece_w + self.piece_gap

        step = _line(self.font - 0.5)
        for index, line in enumerate(self.note):
            ax.text(x + width / 2.0, body_top + body_h + 0.04 + (index + 0.5) * step,
                    line, fontsize=self.font - 0.5, ha='center', va='center',
                    color=QUIET)


def _place(rows, width):
    """
    Position the rows, each running the opposite way to the one above it.

    The rows are given rather than packed, because which panels belong together is a
    statement about the reduction: the exponent matrix beside the Hermite form it goes
    into, and V beside the inverse that undoes it.

    Returns:
        tuple: ``[(panel, x, y)]`` in flow order, and the total height.
    """
    placed, y = [], 0.0
    for index, row in enumerate(rows):
        # Across is short and wide, stacked is tall and narrow, and the page has a
        # width but no height.  So every panel lays its pieces across, and only a
        # panel in a row that overflows is made to stack -- the widest first, since
        # that is the one buying the most room.
        while sum(panel.size[0] for panel in row) \
                + GAP_IN * (len(row) - 1) > width:
            candidates = [panel for panel in row if panel.can_stack]
            if not candidates:
                break
            max(candidates, key=lambda panel: panel.size[0]).stacked = True

        row_h = max(panel.size[1] for panel in row)
        row_w = sum(panel.size[0] for panel in row) + GAP_IN * (len(row) - 1)
        # Every row reads left to right, because the order of the matrices is the
        # order of the steps and a reader should not have to work out which way a
        # given row runs.
        cursor, positions = (width - row_w) / 2.0, []
        for panel in row:
            positions.append(cursor)
            cursor += panel.size[0] + GAP_IN
        for panel, start in zip(row, positions):
            placed.append((panel, start, y + (row_h - panel.size[1]) / 2.0))
        y += row_h + GAP_IN * 1.1
    return placed, y - GAP_IN * 1.1


#: Inches.  The radius a corner of a connector is rounded by, and how far a brace
#: reaches out from the rows it gathers.
CORNER_IN, BRACE_IN = 0.09, 0.08


def _brace(ax, x, top, bottom, colour=QUIET, width=BRACE_IN):
    """
    A curly brace gathering the rows between `top` and `bottom`, pointing left.

    Drawn rather than set as a glyph, because a `}` scaled by its font size to span
    four rows of a matrix comes out far too heavy.
    """
    middle = (top + bottom) / 2.0
    inner, point = x, x - width
    vertices = [(inner, top),
                (point, top), (point, (top + middle) / 2.0),
                (point, middle), (point - width * 0.6, middle),
                (point, middle), (point, (middle + bottom) / 2.0),
                (point, bottom), (inner, bottom)]
    codes = [Path.MOVETO,
             Path.CURVE3, Path.CURVE3,
             Path.CURVE3, Path.CURVE3,
             Path.CURVE3, Path.CURVE3,
             Path.CURVE3, Path.CURVE3]
    ax.add_patch(PathPatch(Path(vertices, codes), facecolor='none',
                           edgecolor=colour, linewidth=0.8, zorder=2))


def _rounded(points, radius=CORNER_IN):
    """
    A path through `points`, square except that its corners are rounded off.

    The corners are quadratic Béziers whose control point is the corner itself, which
    is what makes the turn look drawn rather than mitred.
    """
    vertices, codes = [points[0]], [Path.MOVETO]
    for index in range(1, len(points) - 1):
        before, corner, after = points[index - 1], points[index], points[index + 1]
        into = _step_towards(corner, before, radius)
        away = _step_towards(corner, after, radius)
        vertices.extend([into, corner, away])
        codes.extend([Path.LINETO, Path.CURVE3, Path.CURVE3])
    vertices.append(points[-1])
    codes.append(Path.LINETO)
    return Path(vertices, codes)


def _step_towards(start, end, distance):
    """`distance` inches from `start` in the direction of `end`, or half way."""
    span = ((end[0] - start[0]) ** 2 + (end[1] - start[1]) ** 2) ** 0.5
    if not span:
        return start
    step = min(distance, span / 2.0) / span
    return (start[0] + (end[0] - start[0]) * step,
            start[1] + (end[1] - start[1]) * step)


def _connect(ax, placed, edges, font):
    """
    An arrow from each panel to the next, named with the step it performs.

    Every arrow leaves one box's border square to it and enters the next's square to
    that, turning at rounded right angles in the space between.  A curve drawn
    straight from box to box cuts across whatever lies between them; a route that only
    ever travels along the gaps cannot.
    """
    for source, target in edges:
        panel, x, y = placed[source]
        width, height = panel.size
        next_panel, next_x, next_y = placed[target]
        next_w, next_h = next_panel.size

        # Two panels that overlap vertically are side by side in a row; anything else
        # is a drop to the row below.
        drops = next_y >= y + height or y >= next_y + next_h

        if drops:
            # Down out of this row, along the empty band between the rows, and down
            # into the next.
            start = (x + width / 2.0, y + height)
            end = (next_x + next_w / 2.0, next_y)
            band = (start[1] + end[1]) / 2.0
            points = [start, (start[0], band), (end[0], band), end]
        else:
            # Across, stepping through the middle of the gap if the two boxes are not
            # centred on the same line.
            forwards = next_x > x
            start = (x + width if forwards else x, y + height / 2.0)
            end = (next_x if forwards else next_x + next_w, next_y + next_h / 2.0)
            if abs(start[1] - end[1]) < 0.01:
                points = [start, end]
            else:
                middle = (start[0] + end[0]) / 2.0
                points = [start, (middle, start[1]), (middle, end[1]), end]

        ax.add_patch(FancyArrowPatch(
            path=_rounded(points), arrowstyle='-|>', mutation_scale=8,
            linewidth=1.0, color=QUIET, shrinkA=0, shrinkB=0, zorder=1))


def equation_lines(system):
    """Every equation of a system, as mathtext."""
    lines = []
    clock = system.indep_var
    for symbol, derivative in system.derivative_dict.items():
        if derivative is None or symbol == clock:
            continue
        # A slash rather than a built fraction: at this size a numerator sets at about
        # four points, and the derivative is the part nobody needs to read closely.
        lines.append(_math('d{}/d{} = {}'.format(
            expr_to_tex(symbol), expr_to_tex(clock), expr_to_tex(derivative))))
    if system.parameters:
        lines.append('parameters: $' + '$, $'.join(
            expr_to_tex(parameter) for parameter in system.parameters) + '$')
    return lines


def dimensionless_lines(system):
    """
    Each equation multiplied through by :math:`t/x`, expanded into its terms.

    This is what the exponents are read from, and the panel that shows it is the only
    place a viewer can learn why :math:`t` turns up in rows of K that the original
    equations never mention.
    """
    clock = system.indep_var
    return [_math('({}/{})\\,d{}/d{} = {}'.format(
        expr_to_tex(clock), expr_to_tex(symbol), expr_to_tex(symbol),
        expr_to_tex(clock), expr_to_tex(expression)))
        for symbol, expression in system.dimensionless_equations().items()]


def symmetry_lines(system, scaling_matrix):
    """
    Each row of the scaling matrix written out as the substitution it stands for.

    A row holds powers, and a weight of :math:`-1` means *divided by*; reading that off
    a matrix is a step worth not making the reader take.
    """
    # One line per row, and the symbols it moves taken together: a row is one action,
    # not a handful of separate substitutions that happen to share a factor.  The
    # symbols it fixes are left out, which is what makes the line short enough to
    # read, and the factor is whatever desr named it.
    lines = []
    for moved, factor in zip(system.scaling_symmetries(scaling_matrix),
                             scaling_factors(system, scaling_matrix)):
        # \left( and \right), because the images hold built fractions and a pair of
        # parentheses set at the running size sits inside them rather than round them.
        lines.append(_math('{}:\\;\\left({}\\right) \\mapsto \\left({}\\right)'.format(
            expr_to_tex(factor),
            ',\\,'.join(expr_to_tex(symbol) for symbol in moved),
            ',\\,'.join(expr_to_tex(image) for image in moved.values()))))
    return lines


def scaling_factors(system, scaling_matrix):
    """The factor each symmetry is written in, read back off the symmetry itself."""
    factors = []
    for moved in system.scaling_symmetries(scaling_matrix):
        appearing = set()
        for image in moved.values():
            appearing |= image.free_symbols
        factors.append((appearing - set(system.symbols)).pop())
    return factors


def _bracket(symbol):
    """The dimension of `symbol`, as the name desr writes for it."""
    return sympy.Symbol('[{}]'.format(symbol))


def dimension_lines(system, scaling_matrix):
    """
    The dimensions, one line per set of symbols that share one.

    Symbols with the same column of the scaling matrix have the same dimension, so
    they belong on one line: that equality is the fact, and a line each would bury it
    in repetition.  A group whose dimension is the name of a symbol in it says all it
    has to say by the equalities alone, so the right-hand side is dropped rather than
    written as a tautology.
    """
    names, _ = system.dimension_matrix(scaling_matrix)
    lines = []
    for dimension, symbols in system.symbols_by_dimension(scaling_matrix):
        left = ' = '.join(dimension_to_tex(_bracket(symbol)) for symbol in symbols)
        # A group whose dimension is one of the names is the group that name came
        # from, so the equalities have already said it.
        if dimension in names:
            if len(symbols) == 1:
                continue
            lines.append('$' + left + '$')
        else:
            lines.append(_math(left + ' = ' + dimension_to_tex(dimension)))
    return lines


def _monomial_tex(monomial):
    """
    A monomial as one line, with a slash rather than a built fraction.

    These are single ratios, not fractions of sums, so a slash says the whole of it.
    Setting them as display fractions costs a line of height and a good deal of width
    each, and there are one of them per symbol.
    """
    numerator, denominator = sympy.fraction(sympy.together(monomial))
    if denominator == 1:
        return expr_to_tex(numerator)
    below = expr_to_tex(denominator)
    if not (denominator.is_Symbol or denominator.is_Number):
        below = '({})'.format(below)
    return '{}/{}'.format(expr_to_tex(numerator), below)


def monomial_lines(names, monomials, images=None):
    """
    ``name = monomial`` for each column of V or W, and where the monomial goes.

    The arrow appears on a row the normalisation reached and nowhere else, which is
    the rule the chain figure and the symbol graph both follow.  It is not an equals:
    :math:`t = \\tau x_0 \\mapsto \\tau` claims nothing about :math:`x_0` being 1, only
    that this is where the monomial goes once the unit is adopted.
    """
    lines = []
    for name, monomial in zip(names, monomials):
        line = '{} = {}'.format(expr_to_tex(name), _monomial_tex(monomial))
        image = (images or {}).get(name)
        if image is not None and image != monomial:
            line += ' \\mapsto ' + _monomial_tex(image)
        lines.append('$' + line + '$')
    return lines


def flowchart_panels(system, font):
    """
    The panels, in the order the reduction performs them, and the words between them.

    Returns:
        tuple: ``(panels, labels)``, one label per arrow.
    """
    translator = ODETranslator.from_ode_system(system)
    old = [_tex(symbol) for symbol in system.symbols]
    names = list(translator.column_names(system))
    new = [_tex(name) for name in names]
    r = translator.r

    K = system.exponent_matrix()
    H, U = system.exponent_matrix_hnf()
    picked, A = system.scaling_symmetry_rows(), system.maximal_scaling_matrix()
    zero_rows = range(H.shape[0] - picked.shape[0], H.shape[0])
    V, W = translator.herm_mult, translator.inv_herm_mult

    kinds = {index: (UNIT_COLOUR if index < r else VALUE_COLOUR)
             for index in range(len(new))}
    reduced = translator.translate(system)
    # The rows of A are the scalings, and naming them after their factors is what ties
    # the matrix to the substitutions written out beside it.
    lamdas = ['$\\lambda$'] if r == 1 else \
        ['$\\lambda_{{{}}}$'.format(index + 1) for index in range(r)]
    # A column of D is a symbol's dimension, read as exponents of these names.
    names_of_dimensions, D = system.dimension_matrix(A)
    dimension_names = ['$' + dimension_to_tex(name) + '$'
                       for name in names_of_dimensions]

    # Which panels share a row is a statement about the reduction, so the rows are
    # written down rather than packed by width: the exponent matrix beside the Hermite
    # form that selects from it, and V beside the inverse that undoes it.
    panels = [
        [Panel('Original system', [Lines(equation_lines(system), font)], font,
               note='{} symbols'.format(len(old))),
         Panel('Make equations dimensionless',
               [Lines(dimensionless_lines(system), font)], font)],

        [Panel('Form exponent matrix $K$', [Matrix(K, font, row_labels=old)], font),
         Panel('Row Hermite normal form $U K = H_K$ finds {}'.format(
                   '1 scaling' if r == 1 else '{} scalings'.format(r)),
               # The brace gathers the selected rows and counts them, which is the
               # one place they are pointed out.
               [Matrix(U, font, highlight=zero_rows, label='$U =$',
                       brace='$r = {}$'.format(r)),
                Matrix(H, font, highlight=zero_rows, label='$H_K =$')], font)],

        [Panel('Scaling matrix $A$ gives scalings',
               [Matrix(A, font, col_labels=old, row_labels=lamdas, label='$A =$'),
                Lines(symmetry_lines(system, A), font)], font),
         Panel('Dimension matrix $D$ gives dimensions',
               [Matrix(D, font, col_labels=old, row_labels=dimension_names,
                       label='$D =$'),
                Lines(dimension_lines(system, A), font, column_in=0.5)], font)],

        [Panel('Normal Hermite multiplier $V$: $A V = H_A$ defines invariants',
               [Matrix(V, font, row_labels=old, col_labels=new,
                       col_colours=kinds, col_divider=r),
                Lines(monomial_lines(names, translator.forward_monomials()),
                      font, colours=kinds)], font),
         Panel('$W = V^{-1}$ gives substitutions',
               [Matrix(W, font, row_labels=new, col_labels=old, row_divider=r,
                       row_colours=kinds),
                Lines(monomial_lines(system.symbols,
                                     translator.reverse_monomials(system),
                                     translator.substitutions(system)),
                      font)], font)],

        [Panel('Reduced system', [Lines(equation_lines(reduced), font)], font,
               note='{} symbols'.format(len(reduced.symbols)))],
    ]
    # The spine is the reduction, and the titles say what each step is, so the arrows
    # carry no words of their own.  The one arrow that leaves the spine goes to D,
    # because the dimensions are something A is read for and not something the
    # reduction does; the column Hermite normal form is taken of A, so that arrow
    # leaves A and not D.
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (4, 6), (6, 7), (7, 8)]

    # Numbered in flow order, so a box can be pointed at by number.  The numbers run
    # down the spine and the branch takes the next one, which is the order the panels
    # are written in above.
    for number, panel in enumerate(
            [panel for row in panels for panel in row], start=1):
        panel.number = number
    return panels, edges


def flowchart_figure(system, title=None, width=PAGE_IN, font=6.0):
    """
    The reduction as a loop of matrices, one page-width wide and as tall as it needs.

    Args:
        system (ODESystem): The system to reduce.
        title (str, optional): A heading.  A journal caption says this instead.
        width (float, optional): Inches.  Defaults to a page's text width.
        font (float, optional): Point size of the body text.

    Returns:
        matplotlib.figure.Figure
    """
    # Every artist resolves its font when it is made, so the whole figure is built
    # inside the context rather than only created in it.
    with matplotlib.rc_context(FONT_RC):
        rows, edges = flowchart_panels(system, font)
        placed, height = _place(rows, width - 2 * PAD_IN)

        title_h = _line(font + 3) + 0.12 if title else 0.0
        total = height + title_h + 2 * PAD_IN

        figure = plt.figure(figsize=(width, total))
        ax = figure.add_axes([0.0, 0.0, 1.0, 1.0])
        ax.set_xlim(0.0, width)
        ax.set_ylim(total, 0.0)
        ax.axis('off')

        if title:
            ax.text(width / 2.0, PAD_IN + title_h * 0.5, title, fontsize=font + 3,
                    weight='bold', ha='center', va='center', color=PLAIN)

        offset = PAD_IN + title_h
        shifted = [(panel, x + PAD_IN, y + offset) for panel, x, y in placed]
        for panel, x, y in shifted:
            panel.draw(ax, x, y)
        _connect(ax, shifted, edges, font)
    return figure


SYSTEMS = [
    ('Michaelis-Menten, two variable', ODESystem.from_equations(
        ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
         'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])),
    ('Hubert & Labahn 6.6', ODESystem.from_equations(
        ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
         'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])),
    ('enzyme chain, a reaction network', ChemicalReactionNetwork.from_diagram(
        'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()),
]

_label, _system = SYSTEMS[0]
figure = flowchart_figure(_system, _label)


if __name__ == '__main__':
    plt.show()
