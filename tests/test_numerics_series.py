'''
A time series taken through a reduction and back, by integrating the original system.

The other numeric tests build their series by hand.  These integrate the original
equations, push the trajectory forward and bring it back, which is the workflow the
module exists for -- and the one in which three defects sat unnoticed, because no test
ran a series through a reduction that carries a moving auxiliary, or started one at zero.
'''
import contextlib
import signal
from unittest import TestCase, skipUnless

import numpy as np
import sympy
from scipy.integrate import solve_ivp

from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.numerics import InsufficientKnownValues, NumericTranslator
from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator

S = sympy.Symbol


@contextlib.contextmanager
def time_limit(seconds):
    '''
    Fail rather than hang.  A rate that is not finite makes `solve_ivp` shrink its step
    for ever, so the regression these tests guard against does not fail -- it never
    returns -- and without a limit it would take the whole suite with it.
    '''
    def expired(signum, frame):
        raise AssertionError('did not finish within {} s'.format(seconds))
    previous = signal.signal(signal.SIGALRM, expired)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def integrate(system, parameters, start, times):
    '''
    A system integrated through the given samples of its own independent variable:
    ``{symbol: array, or scalar for a parameter}``.  Nothing of desr's numerics is used,
    so what comes back is an oracle and not an echo.
    '''
    clock = system.independent_variable
    dependents = list(system.dependents)
    rates = sympy.lambdify(
        [clock, dependents],
        [system.derivative_dict[v].subs(parameters) for v in dependents], modules='numpy')
    times = np.asarray(times, dtype=float)
    solved = solve_ivp(lambda t, y: np.broadcast_to(np.asarray(rates(t, y), dtype=float),
                                                    np.shape(y)),
                       (times[0], times[-1]), [start[v] for v in dependents],
                       t_eval=times, rtol=1e-11, atol=1e-13)
    values = {clock: solved.t}
    values.update({v: solved.y[row] for row, v in enumerate(dependents)})
    values.update(parameters)
    return values


def trajectory(system, parameters, start, span, points=40):
    return integrate(system, parameters, start, np.linspace(span[0], span[1], points))


def first_point(series):
    return {symbol: (float(value[0]) if np.ndim(value) else value)
            for symbol, value in series.items()}


def along(expression, series):
    ''' A monomial in a system's symbols, evaluated along a trajectory of that system. '''
    symbols = sorted(expression.free_symbols, key=str)
    value = sympy.lambdify(symbols, expression, modules='numpy')(
        *[series[symbol] for symbol in symbols])
    return np.asarray(value, dtype=float)


class SeriesCase(TestCase):
    def assert_recovers(self, original, recovered, rtol):
        self.assertEqual(set(original), set(recovered))
        for symbol, value in original.items():
            np.testing.assert_allclose(
                np.asarray(recovered[symbol], dtype=float), np.asarray(value, dtype=float),
                rtol=rtol, atol=0, err_msg=str(symbol))


class TestTheReducedSystemIsTheDynamicsOfItsMonomials(SeriesCase):
    '''
    The square that has to commute, and the one statement here that does not lean on
    desr's numerics to check desr's numerics.

    Integrate the original system.  Each reduced symbol denotes a monomial in the original
    ones, so evaluating those monomials along the trajectory gives a curve in the reduced
    symbols.  Integrate the reduced system separately, from that curve's first point.  The
    two curves must be the same curve.

    A wrong multiplier, a wrong equation, a clock that denotes something other than what
    its record says, or a pairing of names with the wrong monomials each break it.  A
    forward-then-reverse round trip breaks for none of them, as long as both directions
    are wrong the same way.
    '''

    CASES = {
        'enzyme chain, general, three normalised and the clock identified': (
            lambda: ChemicalReactionNetwork.from_diagram(
                'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system(),
            {S('k_0_1'): 1.3, S('k_2_3'): 0.7, S('k_4_5'): 2.1},
            {S('e'): 0.5, S('f'): 0.8, S('x'): 1.0, S('y'): 0.4, S('z'): 0.2},
            (0.5, 3.0)),
        'Michaelis-Menten, parameter scheme': (
            lambda: ODESystem.from_equations(
                ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                 'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']),
            {S('e_0'): 0.6, S('k_1'): 1.7, S('k_2'): 0.9, S('k_m1'): 0.4},
            {S('c'): 0.05, S('s'): 1.2}, (0.5, 4.0)),
        'Hubert & Labahn 6.4, dependent-variable scheme': (
            lambda: ODESystem.from_equations(
                ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)']),
            {}, {S('z1'): 0.3, S('z2'): 0.5}, (1.0, 1.6)),
        'Hubert & Labahn 6.6, general scheme, nothing identified': (
            lambda: ODESystem.from_equations(
                ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
                 'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)']),
            {}, {S('z1'): 0.9, S('z2'): 0.6}, (1.0, 2.0)),
    }

    def reductions(self, system):
        translator = ODETranslator.from_ode_system(system)
        yield 'translate()', translator.translate(system)
        yield ("Hubert & Labahn's form",
               translator.translate_general(system, normalise_parameters=False,
                                            identify_indep_var=False))

    @staticmethod
    def disagreement(original, reduced, denotes):
        '''
        How far the two curves are apart, as the worst relative gap over every variable.

        `denotes` is what each reduced symbol is taken to stand for.  It is an argument so
        that a test can hand in a wrong one and see that this notices.
        '''
        image = {symbol: along(monomial, original) for symbol, monomial in denotes.items()}
        constant = {p: float(np.ravel(image[p])[0]) for p in reduced.parameters}
        worst = 0.0
        for p in reduced.parameters:
            spread = np.max(np.abs(image[p] - constant[p])) / max(abs(constant[p]), 1e-300)
            worst = max(worst, float(spread))
        begin = {v: float(np.ravel(image[v])[0]) for v in reduced.dependents}
        solved = integrate(reduced, constant, begin, image[reduced.independent_variable])
        for variable in reduced.dependents:
            expected = np.broadcast_to(image[variable], solved[variable].shape)
            gap = np.max(np.abs(solved[variable] - expected)
                         / np.maximum(np.abs(expected), 1e-12))
            worst = max(worst, float(gap))
        return worst

    def test_the_square_commutes(self):
        for label, (build, parameters, start, span) in self.CASES.items():
            system = build()
            original = trajectory(system, parameters, start, span)
            for how, reduced in self.reductions(system):
                with self.subTest(system=label, reduction=how):
                    self.assertLess(
                        self.disagreement(original, reduced, reduced.named_symbols), 1e-6)

    def test_the_square_notices_a_symbol_paired_with_the_wrong_monomial(self):
        # Otherwise the test above could pass for a reason that has nothing to do with
        # the reduction.  Swap what two dependents are said to denote -- the pairing bug
        # that shipped in three releases -- and the curves must come apart.
        build, parameters, start, span = self.CASES[
            'enzyme chain, general, three normalised and the clock identified']
        system = build()
        original = trajectory(system, parameters, start, span)
        reduced = ODETranslator.from_ode_system(system).translate(system)

        swapped = dict(reduced.named_symbols)
        first, second = S('y3'), S('y4')
        swapped[first], swapped[second] = swapped[second], swapped[first]
        self.assertGreater(self.disagreement(original, reduced, swapped), 1e-2)

        # And a clock said to be the original one, when it is a rescaled time.
        mistimed = dict(reduced.named_symbols)
        mistimed[reduced.independent_variable] = system.independent_variable
        self.assertGreater(self.disagreement(original, reduced, mistimed), 1e-2)

    def test_forward_is_the_monomials(self):
        # `forward` rebuilds the multiplier from the record; this is the same statement
        # made directly from what each symbol is said to denote.
        for label, (build, parameters, start, span) in self.CASES.items():
            system = build()
            original = trajectory(system, parameters, start, span)
            for how, reduced in self.reductions(system):
                with self.subTest(system=label, reduction=how):
                    forward = NumericTranslator(reduced).forward(original)
                    self.assertEqual(set(reduced.symbols), set(forward))
                    for symbol, monomial in reduced.named_symbols.items():
                        np.testing.assert_allclose(
                            np.broadcast_to(forward[symbol],
                                            np.shape(along(monomial, original))),
                            along(monomial, original), rtol=1e-12, err_msg=str(symbol))


class TestASeriesThroughAReductionThatCarriesAnAuxiliary(SeriesCase):
    '''
    The enzyme chain's default reduction normalises three constant auxiliaries and keeps
    the one that moves, `x3`.  Its series is in the solution, so recovering the original
    system needs the three rate constants and nothing else.
    '''

    PARAMETERS = {S('k_0_1'): 1.3, S('k_2_3'): 0.7, S('k_4_5'): 2.1}
    START = {S('e'): 0.5, S('f'): 0.8, S('x'): 1.0, S('y'): 0.4, S('z'): 0.2}

    def setUp(self):
        self.system = ChemicalReactionNetwork.from_diagram(
            'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
        self.translator = ODETranslator.from_ode_system(self.system)
        self.original = trajectory(self.system, self.PARAMETERS, self.START, (0.0, 3.0))
        self.reduced = self.translator.translate(self.system)
        self.numeric = NumericTranslator(self.reduced)
        self.forward = self.numeric.forward(self.original)

    def test_the_premise(self):
        # Three normalised, one carried: neither "carries its auxiliaries" nor "carries
        # none", which is the case the series path did not have.
        self.assertEqual(3, len(self.reduced.normalisations))
        self.assertEqual((S('x3'),), tuple(self.reduced.auxiliary_names))
        self.assertFalse(self.numeric.carries_auxiliaries)
        self.assertEqual(1, np.ndim(self.forward[S('x3')]))

    def test_the_series_comes_back_from_the_constants_alone(self):
        recovered = self.numeric.reverse(self.forward, known_values=self.PARAMETERS)
        self.assert_recovers(self.original, recovered, rtol=1e-12)

    def test_the_recovered_auxiliaries_are_what_they_denote(self):
        # Each auxiliary is a monomial in the original symbols, so along the trajectory
        # its value is that monomial's -- known without asking the code under test.  The
        # three constant ones are 1/k_2_3, k_2_3/k_0_1 and k_2_3/k_4_5; the one that moves
        # is x itself, and it is exact because it is read off the solution rather than
        # integrated a second time.
        recovered = self.numeric.recover_auxiliaries(
            self.forward, known_values=self.PARAMETERS)
        monomials = list(self.translator.auxiliary_monomials())
        self.assertEqual(4, len(recovered))
        for value, monomial in zip(recovered, monomials):
            np.testing.assert_allclose(
                np.broadcast_to(value, np.shape(along(monomial, self.original))),
                along(monomial, self.original), rtol=1e-12, err_msg=str(monomial))
        self.assertEqual(S('x'), monomials[3])
        np.testing.assert_array_equal(recovered[3], self.original[S('x')])

    def test_a_constant_comes_back_as_a_scalar(self):
        recovered = self.numeric.reverse(self.forward, known_values=self.PARAMETERS)
        for constant in self.PARAMETERS:
            self.assertEqual(0, np.ndim(recovered[constant]), str(constant))

    def test_it_asks_for_three_values_not_four(self):
        # The carried auxiliary is one equation the reduced system states outright.
        with self.assertRaises(InsufficientKnownValues) as caught:
            self.numeric.reverse(self.forward)
        self.assertIn('3 are needed', str(caught.exception))
        self.assertIn('k_0_1, k_2_3, k_4_5', str(caught.exception))

    def test_the_point_and_the_series_agree_about_what_is_needed(self):
        point = first_point(self.forward)
        with self.assertRaises(InsufficientKnownValues) as at_a_point:
            self.numeric.reverse(point)
        with self.assertRaises(InsufficientKnownValues) as along_a_series:
            self.numeric.reverse(self.forward)
        self.assertIn('3 are needed', str(at_a_point.exception))
        self.assertIn('3 are needed', str(along_a_series.exception))


class TestWhatASeriesIsToldWhenNoConstantCanHelp(SeriesCase):
    '''
    In a series only constants may be known, and sometimes what is still free rescales
    only things that vary.  The message used to end "Supply the value of one of: ." --
    an empty list -- and now says where a value of anything may be given instead.
    '''

    def test_the_enzyme_chain_with_its_auxiliaries_dropped(self):
        system = ChemicalReactionNetwork.from_diagram(
            'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
        parameters = TestASeriesThroughAReductionThatCarriesAnAuxiliary.PARAMETERS
        original = trajectory(
            system, parameters,
            TestASeriesThroughAReductionThatCarriesAnAuxiliary.START, (0.0, 1.0))
        numeric = NumericTranslator(ODETranslator.from_ode_system(
            system).translate_general(system, include_aux_vars=False))
        with self.assertRaises(InsufficientKnownValues) as caught:
            numeric.reverse(numeric.forward(original), known_values=parameters)
        message = str(caught.exception)
        self.assertNotIn('one of: .', message)
        self.assertIn('auxiliaries_from()', message)
        self.assertIn('auxiliaries=', message)

    def test_a_system_with_no_constants_at_all(self):
        system = ODESystem.from_equations(
            ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
             'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])
        original = trajectory(system, {}, {S('z1'): 0.9, S('z2'): 0.6}, (1.0, 2.0))
        numeric = NumericTranslator(ODETranslator.from_ode_system(
            system).translate_general(system, include_aux_vars=False))
        with self.assertRaises(InsufficientKnownValues) as caught:
            numeric.reverse(numeric.forward(original))
        message = str(caught.exception)
        self.assertNotIn('one of: .', message)
        self.assertIn('auxiliaries_from()', message)


@skipUnless(hasattr(signal, 'SIGALRM'), 'needs SIGALRM to bound a run that may not return')
class TestASeriesThatStartsAtZero(SeriesCase):
    '''
    Hubert & Labahn write a system as `t dz/dt = z * F`, so a reduction kept in their form
    has a `1/t` in every rate.  For the enzyme chain that is only how the form writes it:
    the invariant `k_2_3*t` is a clock, and the reduction takes it as one unless it is told
    not to, which leaves no `1/t` at all.  Told not to, the invariant sits beside `t` and
    the quotient is a constant.  For Hubert & Labahn 6.6 the `1/t` is in the equations and
    the system is singular at zero.

    Either way the solver must not be handed a rate that is not finite, because it does
    not fail on one -- it never returns.
    '''

    def test_a_removable_singularity_is_removed(self):
        system = ChemicalReactionNetwork.from_diagram(
            'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
        case = TestASeriesThroughAReductionThatCarriesAnAuxiliary
        original = trajectory(system, case.PARAMETERS, case.START, (0.0, 3.0))
        numeric = NumericTranslator(ODETranslator.from_ode_system(
            system).translate_general(system, include_aux_vars=False,
                                      identify_indep_var=False))

        # The rate really does divide by the clock, and the series really does start there.
        clock = numeric.reduced.independent_variable
        self.assertEqual(system.independent_variable, clock)
        self.assertTrue(any(sympy.denom(sympy.together(rate)).has(clock)
                            for rate in numeric.auxiliary_growth_rates()))
        self.assertEqual(0.0, original[system.independent_variable][0])

        known = dict(case.PARAMETERS)
        known[S('x')] = case.START[S('x')]
        forward = numeric.forward(original)
        auxiliaries = numeric.auxiliaries_from(first_point(forward), known)
        with time_limit(60):
            recovered = numeric.reverse(forward, auxiliaries=auxiliaries)
        self.assert_recovers(original, recovered, rtol=1e-6)

    def test_with_the_change_of_clock_there_is_nothing_to_remove(self):
        # The default: k_2_3*t is the clock, no rate divides by it, and a series that
        # starts at zero comes back.
        system = ChemicalReactionNetwork.from_diagram(
            'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
        case = TestASeriesThroughAReductionThatCarriesAnAuxiliary
        original = trajectory(system, case.PARAMETERS, case.START, (0.0, 3.0))
        numeric = NumericTranslator(ODETranslator.from_ode_system(
            system).translate_general(system, include_aux_vars=False))

        clock = numeric.reduced.independent_variable
        self.assertEqual(S('k_2_3') * system.independent_variable,
                         numeric.reduced.named_symbols[clock])
        self.assertFalse(any(sympy.denom(sympy.together(rate)).has(clock)
                             for rate in numeric.auxiliary_growth_rates()))
        self.assertEqual(0.0, original[system.independent_variable][0])

        known = dict(case.PARAMETERS)
        known[S('x')] = case.START[S('x')]
        forward = numeric.forward(original)
        auxiliaries = numeric.auxiliaries_from(first_point(forward), known)
        with time_limit(60):
            recovered = numeric.reverse(forward, auxiliaries=auxiliaries)
        self.assert_recovers(original, recovered, rtol=1e-6)

    def test_a_real_singularity_is_refused(self):
        system = ODESystem.from_equations(
            ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
             'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])
        numeric = NumericTranslator(ODETranslator.from_ode_system(
            system).translate_general(system, include_aux_vars=False))
        through_zero = {S('t'): np.linspace(0.0, 1.0, 5),
                        S('y0'): np.ones(5), S('y1'): np.ones(5)}
        with time_limit(60):
            with self.assertRaises(ValueError) as caught:
                numeric.recover_auxiliaries(through_zero, auxiliaries=[1.0])
        message = str(caught.exception)
        self.assertIn('divides by', message)
        self.assertIn('t = 0.0', message)

    def test_the_same_system_away_from_zero_is_untouched(self):
        system = ODESystem.from_equations(
            ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
             'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])
        original = trajectory(system, {}, {S('z1'): 0.9, S('z2'): 0.6}, (1.0, 2.0))
        numeric = NumericTranslator(ODETranslator.from_ode_system(system).translate(system))
        with time_limit(60):
            recovered = numeric.reverse(numeric.forward(original))
        self.assert_recovers(original, recovered, rtol=1e-12)


class TestASeriesUnderTheDependentVariableScheme(SeriesCase):
    ''' Hubert & Labahn 6.4.  No numeric test took this scheme at all. '''

    def test_it_round_trips_and_the_clock_is_untouched(self):
        system = ODESystem.from_equations(
            ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)'])
        original = trajectory(system, {}, {S('z1'): 0.3, S('z2'): 0.5}, (1.0, 1.6))
        reduced = ODETranslator.from_ode_system(system).translate(system)
        numeric = NumericTranslator(reduced)
        forward = numeric.forward(original)

        # The reduced clock is the original one, sample for sample.
        np.testing.assert_array_equal(forward[reduced.independent_variable],
                                      original[system.independent_variable])
        self.assert_recovers(original, numeric.reverse(forward), rtol=1e-12)
