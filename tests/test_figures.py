"""
The figures, drawn for systems they have not been tuned on.

What is worth testing here is not how anything looks.  It is that a figure drawn for a
model desr has never seen is either right or absent: that it fits the page it claims to
fit, that it says so rather than overflowing when it cannot, and that what it draws is
what the library returns.
"""

from unittest import TestCase

import matplotlib
matplotlib.use('Agg')

from desr.chemical_reaction_network import ChemicalReactionNetwork    # noqa: E402
from desr.figures import flowchart_figure, symbol_map_figure          # noqa: E402
from desr.figures.flowchart import (FONT_FLOOR, PAD_IN, _fitted,      # noqa: E402
                                    _overflow, _place, flowchart_panels)
from desr.ode_system import ODESystem                                 # noqa: E402
from desr.ode_translation import ODETranslator, ReductionScheme       # noqa: E402


def decay_chain(length):
    """A linear chain of `length` species, each decaying into the next."""
    equations = ['da0/dt = -r0*a0']
    for index in range(1, length):
        equations.append('da{}/dt = r{}*a{} - r{}*a{}'.format(
            index, index - 1, index - 1, index, index))
    return ODESystem.from_equations(equations)


class TestTheFiguresDrawUnfamiliarSystems(TestCase):

    def systems(self):
        return (('Michaelis-Menten', ODESystem.from_equations(
                    ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                     'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])),
                ('enzyme chain', ChemicalReactionNetwork.from_diagram(
                    'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()),
                ('HL 6.6', ODESystem.from_equations(
                    ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
                     'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])),
                ('HL 6.4, which keeps its clock', ODESystem.from_equations(
                    ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)'])),
                ('a decay chain', decay_chain(4)))

    @staticmethod
    def with_initial():
        system = ODESystem.from_equations(
            ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
             'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
        system.update_initial_conditions({'s': 's_0'})
        return system

    def test_the_flowchart_draws_and_keeps_to_its_width(self):
        for label, system in self.systems():
            with self.subTest(system=label):
                figure = flowchart_figure(system)
                self.assertEqual(6.5, round(figure.get_size_inches()[0], 2))

    def test_the_symbol_map_draws(self):
        for label, system in self.systems():
            with self.subTest(system=label):
                self.assertIsNotNone(symbol_map_figure(system))

    def test_the_symbol_map_is_as_wide_as_its_panels_or_its_chain(self):
        # The three panels are measured and packed side by side, so nothing drawn is
        # off the page and none of the three overlaps the next.
        for label, system in self.systems():
            with self.subTest(system=label):
                figure = symbol_map_figure(system)
                figure.canvas.draw()
                renderer = figure.canvas.get_renderer()
                page = figure.bbox
                panels = [ax for ax in figure.axes if ax.texts][:3]
                extents = []
                for ax in panels:
                    box = ax.get_tightbbox(renderer)
                    self.assertGreaterEqual(box.x0, page.x0 - 1)
                    self.assertLessEqual(box.x1, page.x1 + 1)
                    extents.append((ax.get_position().x0, ax.get_position().x1))
                extents.sort()
                for (_, ends), (starts, _) in zip(extents, extents[1:]):
                    self.assertLessEqual(ends, starts + 1e-9)
                matplotlib.pyplot.close(figure)

    def test_the_arrows_of_the_table_share_a_column(self):
        figure = symbol_map_figure(self.with_initial())
        starts = {round(text.get_position()[0], 6)
                  for ax in figure.axes for text in ax.texts
                  if text.get_text().startswith('$\\mapsto')}
        matplotlib.pyplot.close(figure)
        self.assertEqual(1, len(starts))

    def test_the_table_states_a_dimension_at_one_size(self):
        # `[t] = [t]` is one bracket twice, so both sides are set alike.
        figure = symbol_map_figure(self.with_initial())
        # The left half: a bracketed symbol alone, and the `= dimension` beside it.
        sizes = {text.get_fontsize() for ax in figure.axes for text in ax.texts
                 if (text.get_text().startswith('$\\left[') and '=' not in text.get_text())
                 or text.get_text().startswith(('$=\\;\\left[',
                                                '$=\\;\\dfrac{1}{\\left['))}
        matplotlib.pyplot.close(figure)
        self.assertEqual(1, len(sizes))

    def test_the_unit_column_gives_a_unit_its_dimension_and_a_value_its_unit(self):
        # Michaelis-Menten: x0 = 1/(e_0 k_1) is the unit of time and x1 = e_0 the unit
        # of concentration.  Their rows state their dimensions; every other row names
        # the unit it is measured in, and nothing in the column is a sentence.
        system = self.with_initial()
        translator = ODETranslator(system)
        figure = symbol_map_figure(system)
        written = [text.get_text() for ax in figure.axes for text in ax.texts]
        matplotlib.pyplot.close(figure)
        self.assertIn('unit', written)
        self.assertIn('dimensions, normalisations, and units', written)
        units = [row for row in translator.unit_table() if row['role'] == 'unit']
        self.assertEqual(2, len(units))
        for row in units:
            stated = [words for words in written
                      if words.startswith('$\\left[{}\\right] = '.format(
                          str(row['symbol'])[0] + '_{' + str(row['symbol'])[1:] + '}'))]
            self.assertEqual(1, len(stated), row['symbol'])
        self.assertFalse([words for words in written if 'in units of' in words])
        self.assertFalse([words for words in written if words.startswith('unit of')])

    def test_the_boxes_stand_together_under_arrows_that_arc_over_them(self):
        from matplotlib.patches import FancyArrowPatch
        figure = symbol_map_figure(self.with_initial())
        figure.canvas.draw()
        width = figure.get_size_inches()[0]
        boxes = sorted((ax.get_position() for ax in figure.axes if not ax.texts
                        or ax.texts[0].get_text() in ('original', 'reduced', 'normalised',
                                                      'nondimensionalised')),
                       key=lambda box: box.x0)[-5:]
        self.assertEqual(5, len(boxes))
        gaps = [(right.x0 - left.x1) * width for left, right in zip(boxes, boxes[1:])]
        for gap in gaps:
            self.assertLess(gap, 0.2)
        # Centred on the page.
        self.assertAlmostEqual(boxes[0].x0, 1.0 - boxes[-1].x1, places=6)
        # Four arrows, each starting and ending above the boxes, bent.
        arrows = [artist for artist in figure.artists
                  if isinstance(artist, FancyArrowPatch)]
        self.assertEqual(4, len(arrows))
        top = max(box.y1 for box in boxes)
        for arrow in arrows:
            (start_x, start_y), (end_x, end_y) = arrow._posA_posB
            self.assertGreater(start_y, top)
            self.assertEqual(start_y, end_y)
            self.assertLess(start_x, end_x)
            self.assertIn('arc3', str(arrow.get_connectionstyle().__class__.__name__).lower())
        matplotlib.pyplot.close(figure)

    def test_the_clock_step_says_what_the_clock_is(self):
        import sympy
        from desr.figures.symbol_map import step_labels
        tau = sympy.Symbol('tau')
        self.assertEqual(('forward', 'normalise', 'clock is $\\tau$', 'reverse'),
                         step_labels({sympy.Symbol('x0')}, {sympy.Symbol('y3'): tau}))
        self.assertEqual(('forward', 'nothing to normalise', 'clock unchanged', 'reverse'),
                         step_labels(set(), {}))

    def test_a_box_counts_its_symbols_in_short_lines(self):
        from desr.figures.symbol_map import summary_lines
        lines = summary_lines(self.with_initial())
        self.assertEqual(['8 symbols', '1 independent', '2 dependents', '5 parameters',
                          'autonomous'], lines)
        # The three kinds add up to the count.
        self.assertEqual(8, 1 + 2 + 5)
        self.assertEqual(8, len(self.with_initial().symbols))

    def test_a_given_size_is_the_size(self):
        figure = symbol_map_figure(self.with_initial(), size=(12.0, 7.0))
        self.assertEqual([12.0, 7.0], [round(side, 2) for side in figure.get_size_inches()])
        matplotlib.pyplot.close(figure)

    def test_no_panel_runs_off_the_page(self):
        # The silent failure this replaces: a panel wider than the page used to be
        # clamped to the margin and allowed to overflow, so the figure looked finished
        # and was not.
        for label, system in self.systems():
            with self.subTest(system=label):
                _, _, placed, _, _ = _fitted(system, 6.5 - 2 * PAD_IN, 6.0)
                self.assertEqual(0.0, _overflow(placed, 6.5 - 2 * PAD_IN))

    def test_a_wide_system_is_set_smaller_rather_than_overflowing(self):
        system = decay_chain(8)
        _, _, _, _, font = _fitted(system, 6.5 - 2 * PAD_IN, 6.0)
        self.assertLess(font, 6.0)
        self.assertGreaterEqual(font, FONT_FLOOR)
        self.assertEqual(6.5, round(flowchart_figure(system).get_size_inches()[0], 2))

    def test_a_system_too_wide_to_set_is_refused(self):
        # Past a certain size a page-wide figure of every matrix entry does not exist,
        # and a smaller font is not the answer.  The message says what it needed.
        with self.assertRaises(ValueError) as raised:
            flowchart_figure(decay_chain(10))
        message = str(raised.exception)
        self.assertIn('does not fit', message)
        self.assertIn('block structure', message)

    def test_a_wider_page_takes_what_the_page_width_refused(self):
        # The knob the message names actually works.
        system = decay_chain(10)
        figure = flowchart_figure(system, width=11.0)
        self.assertEqual(11.0, round(figure.get_size_inches()[0], 2))

    def test_a_reduction_that_keeps_the_clock_is_drawn_too(self):
        # Hubert & Labahn 6.4 reduces under the dependent-variable scheme: the scaling
        # leaves t alone, so t is already invariant and the reduction keeps it rather
        # than coining a symbol for it.  One column of the multiplier is the clock's,
        # and labelling the columns has to put it back or everything after it is
        # paired with the wrong symbol.
        system = ODESystem.from_equations(
            ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)'])
        translator = ODETranslator(system)
        self.assertEqual(ReductionScheme.DEP_VAR, translator.scheme())
        self.assertEqual(6.5, round(flowchart_figure(system).get_size_inches()[0], 2))
        self.assertIsNotNone(symbol_map_figure(system))

    def test_both_figures_draw_a_system_with_an_initial_condition(self):
        system = self.with_initial()
        self.assertIsNotNone(symbol_map_figure(system))
        self.assertIsNotNone(flowchart_figure(system))

    def test_the_boxes_are_numbered_in_flow_order(self):
        rows, _ = flowchart_panels(decay_chain(3), 6.0)
        numbers = [panel.number for row in rows for panel in row]
        self.assertEqual(list(range(1, len(numbers) + 1)), numbers)

    def test_every_row_fits_after_placement(self):
        rows, _ = flowchart_panels(decay_chain(4), 6.0)
        placed, height = _place(rows, 6.5 - 2 * PAD_IN)
        self.assertGreater(height, 0.0)
        self.assertEqual(0.0, _overflow(placed, 6.5 - 2 * PAD_IN))
