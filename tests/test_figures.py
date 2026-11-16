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
        translator = ODETranslator.from_ode_system(system)
        self.assertEqual(ReductionScheme.DEP_VAR, translator.scheme(system))
        self.assertEqual(6.5, round(flowchart_figure(system).get_size_inches()[0], 2))
        self.assertIsNotNone(symbol_map_figure(system))

    def test_an_initial_condition_is_refused_by_the_symbol_map(self):
        # The general scheme does not take initial conditions, and the symbol map
        # needs it.  The flowchart draws such a system, since it never translates.
        system = self.with_initial()
        with self.assertRaises(NotImplementedError):
            symbol_map_figure(system)
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
