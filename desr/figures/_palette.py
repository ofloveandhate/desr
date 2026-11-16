"""
The colours and the face the figures are drawn in, in one place.

Colour names a **kind of symbol** and nothing else: rust the auxiliaries, which are
the units, and blue the invariants, which are the values.  Not the steps -- a step is
grey, because the steps do not divide the way the symbols do.  Black is for what is
neither kind: an expression, a dimension, a count.
"""

UNIT_COLOUR, VALUE_COLOUR = '#8a4500', '#1f5c86'

#: Black for a statement, grey for a label about the figure rather than about the
#: system, and the palest grey for a step that does nothing.
PLAIN, QUIET, INERT = 'black', '#555555', '#aaaaaa'

#: The wash behind a highlighted row, and the panel fills.
PICKED = '#f0e2d2'
PANEL_FILL, PANEL_EDGE = '#fbfaf8', '#d9d6d0'

#: Computer Modern, which matplotlib ships, so a figure sets in the same face as the
#: paper it belongs in and a kappa is not a k.  ``axes.unicode_minus`` is off because
#: cmr10 has no U+2212 and matplotlib says so, loudly, once per label;
#: ``axes.formatter.use_mathtext`` silences a second warning about the same font.
FONT_RC = {
    'font.family': 'serif',
    'font.serif': ['cmr10', 'DejaVu Serif'],
    'mathtext.fontset': 'cm',
    'mathtext.rm': 'serif',
    'axes.unicode_minus': False,
    'axes.formatter.use_mathtext': True,
}
