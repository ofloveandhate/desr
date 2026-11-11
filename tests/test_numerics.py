from unittest import TestCase

import sympy

from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator

# numpy is optional for using desr, but not for testing it: a test suite that quietly skips
# a module when its dependency is missing reports success it has not earned.
import numpy
from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.numerics import NumericTranslator, InsufficientKnownValues, ConflictingKnownValues


class TestNumericTranslator(TestCase):
    '''
    Carrying numbers across the Michaelis-Menten reduction, and the ways of getting it wrong.

    No integrator is involved: the "series" here is any set of arrays run through `forward`,
    which is enough to exercise every rule about what may be a known value.
    '''

    def setUp(self):
        system_tex = '''\\frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\\\
                     \\frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c'''
        self.system = ODESystem.from_tex(system_tex)
        self.system.update_initial_conditions({'s': 's_0'})
        # Parameters first, then the independent variable, then the dependents -- Hubert &
        # Labahn's (c, t, z); see ODESystem.default_order_variables.  Within the parameters
        # the order decides which are normalised to 1: the leading ones are, so s_0 and e_0
        # go first and k_2, k_m1 survive as the reduced constants.
        self.system.reorder_variables(['s_0', 'e_0', 'k_1', 'k_2', 'k_m1', 't', 's', 'c'])
        self.translator = ODETranslator.from_ode_system(
            self.system, naming_scheme={'parameters': 'c', 'independent': 'tau', 'dependents': ['u', 'v']})
        self.numeric = NumericTranslator(self.translator.translate())

        # Fetched by name rather than unpacked positionally, so that a change of variable
        # order is not a change to every test in this class.
        original = {str(v): v for v in self.system.symbols}
        self.t, self.s, self.c = original['t'], original['s'], original['c']
        self.k_m1, self.k_2 = original['k_m1'], original['k_2']
        self.k_1, self.e_0, self.s_0 = original['k_1'], original['e_0'], original['s_0']

        reduced = {str(v): v for v in self.numeric.reduced.symbols}
        self.tau, self.u, self.v = reduced['tau'], reduced['u'], reduced['v']
        self.c0, self.c1, self.c2 = reduced['c0'], reduced['c1'], reduced['c2']

        self.parameters = {self.k_1: 1.5, self.k_m1: 0.9, self.k_2: 0.4,
                           self.e_0: 0.3, self.s_0: 2.0}
        # A "trajectory": the values need not solve anything to test the bookkeeping.
        self.times = numpy.array([0.0, 1.0, 2.0, 3.0])
        self.s_values = numpy.array([2.0, 1.7, 1.5, 1.4])
        self.c_values = numpy.array([0.0, 0.2, 0.19, 0.18])
        self.series = self.numeric.forward({self.t: self.times, self.s: self.s_values,
                                            self.c: self.c_values, **self.parameters})
        self.start = self.numeric.forward({self.t: 0.0, self.s: 2.0, self.c: 0.0,
                                           **self.parameters})
        self.later = self.numeric.forward({self.t: 2.0, self.s: 1.5, self.c: 0.19,
                                           **self.parameters})

    def assertRecovers(self, recovered, **expected):
        for name, value in expected.items():
            symbol = getattr(self, name)
            numpy.testing.assert_allclose(recovered[symbol], value, rtol=1e-12, atol=0)

    # ---------------------------------------------------------------- round trips
    def test_point_round_trip(self):
        recovered = self.numeric.reverse(self.later,
                                         known_values={self.k_1: 1.5, self.s_0: 2.0})
        self.assertRecovers(recovered, t=2.0, s=1.5, c=0.19, k_m1=0.9, k_2=0.4, e_0=0.3)
        for symbol in (self.k_m1, self.k_2, self.e_0):
            self.assertEqual(numpy.ndim(recovered[symbol]), 0)

    def test_series_round_trip(self):
        recovered = self.numeric.reverse(self.series,
                                         known_values={self.k_1: 1.5, self.s_0: 2.0})
        self.assertRecovers(recovered, t=self.times, s=self.s_values, c=self.c_values,
                            k_m1=0.9, k_2=0.4, e_0=0.3)

    def test_series_recovers_constants_as_scalars(self):
        recovered = self.numeric.reverse(self.series,
                                         known_values={self.k_1: 1.5, self.s_0: 2.0})
        for symbol in (self.k_m1, self.k_2, self.e_0, self.k_1, self.s_0):
            self.assertEqual(numpy.ndim(recovered[symbol]), 0, symbol)
        self.assertEqual(recovered[self.s].shape, self.times.shape)

    def test_auxiliaries_from_a_point_translate_a_series(self):
        auxiliaries = self.numeric.auxiliaries_from(self.start,
                                                    known_values={self.k_1: 1.5, self.s: 2.0})
        recovered = self.numeric.reverse(self.series, auxiliaries=auxiliaries)
        self.assertRecovers(recovered, s=self.s_values, c=self.c_values, s_0=2.0, k_m1=0.9)

    # ------------------------------------------------- symbols that do not belong
    def test_foreign_symbol_as_known_value(self):
        zeta = sympy.Symbol('zeta')
        with self.assertRaises(InsufficientKnownValues) as caught:
            self.numeric.reverse(self.later, known_values={self.k_1: 1.5, zeta: 1.0})
        self.assertIn('zeta', str(caught.exception))

    def test_foreign_string_as_known_value(self):
        with self.assertRaises(InsufficientKnownValues) as caught:
            self.numeric.reverse(self.later, known_values={'k_1': 1.5, 'no_such_thing': 1.0})
        self.assertIn('no_such_thing', str(caught.exception))

    def test_string_that_sympifies_to_something_else(self):
        # 'S' sympifies to sympy's singleton registry, not a symbol, and must not slip through.
        with self.assertRaises(InsufficientKnownValues):
            self.numeric.reverse(self.later, known_values={'k_1': 1.5, 'S': 2.0})

    def test_foreign_symbol_in_forward_values(self):
        with self.assertRaises(ValueError) as caught:
            self.numeric.forward({self.t: 0.0, sympy.Symbol('zeta'): 1.0, **self.parameters})
        self.assertIn('zeta', str(caught.exception))

    def test_foreign_symbol_in_reverse_values(self):
        values = dict(self.later)
        values[sympy.Symbol('zeta')] = 1.0
        with self.assertRaises(ValueError) as caught:
            self.numeric.reverse(values, known_values={self.k_1: 1.5, self.s_0: 2.0})
        self.assertIn('zeta', str(caught.exception))

    def test_original_variable_offered_as_reduced_value(self):
        values = dict(self.later)
        values[self.s] = 1.5
        with self.assertRaises(ValueError) as caught:
            self.numeric.reverse(values, known_values={self.k_1: 1.5, self.s_0: 2.0})
        self.assertIn('the reduced system', str(caught.exception))

    # ------------------------------------------ point versus series: what may be known
    def test_series_rejects_dependent_variable(self):
        with self.assertRaises(ValueError) as caught:
            self.numeric.reverse(self.series, known_values={self.k_1: 1.5, self.s: 2.0})
        message = str(caught.exception)
        self.assertIn('s is a function of time', message)
        self.assertIn('s_0', message)

    def test_series_rejects_independent_variable(self):
        with self.assertRaises(ValueError) as caught:
            self.numeric.reverse(self.series, known_values={self.k_1: 1.5, self.t: 2.0})
        self.assertIn('t is a function of time', str(caught.exception))

    def test_point_accepts_dependent_variable(self):
        recovered = self.numeric.reverse(self.start, known_values={self.k_1: 1.5, self.s: 2.0})
        self.assertRecovers(recovered, s_0=2.0, k_m1=0.9, k_2=0.4, e_0=0.3)

    def test_point_dependent_variable_away_from_the_start(self):
        # At t = 2, s = 1.5 while s_0 = 2.0: knowing s here still fixes s_0 correctly.
        recovered = self.numeric.reverse(self.later, known_values={self.k_1: 1.5, self.s: 1.5})
        self.assertRecovers(recovered, s_0=2.0, k_m1=0.9)

    def test_dependent_variable_and_its_constant_constrain_the_same_scale(self):
        with self.assertRaises(InsufficientKnownValues) as caught:
            self.numeric.reverse(self.later, known_values={self.s_0: 2.0, self.s: 1.5})
        self.assertIn('1 more value', str(caught.exception))

    def test_dependent_variable_disagreeing_with_its_constant(self):
        # s = 1.5 at t = 2 implies s_0 = 2.0; s_0 = 3.0 contradicts it.
        with self.assertRaises(ConflictingKnownValues):
            self.numeric.reverse(self.later, known_values={self.s_0: 3.0, self.s: 1.5})

    def test_series_suggestions_are_constants_only(self):
        with self.assertRaises(InsufficientKnownValues) as caught:
            self.numeric.reverse(self.series, known_values={self.k_1: 1.5})
        suggestion = str(caught.exception).splitlines()[-1]
        self.assertIn('s_0', suggestion)
        for wrong in (' s,', ' c,', ' t,', ' s.', ' c.', ' t.'):
            self.assertNotIn(wrong, suggestion)

    def test_a_one_element_array_is_a_series(self):
        values = {x: numpy.array([float(self.start[x])]) for x in (self.tau, self.u, self.v)}
        values.update({x: self.start[x] for x in (self.c0, self.c1, self.c2)})
        with self.assertRaises(ValueError):
            self.numeric.reverse(values, known_values={self.k_1: 1.5, self.s: 2.0})

    # ---------------------------------------------------- how many, and which ones
    def test_none_at_all(self):
        with self.assertRaises(InsufficientKnownValues) as caught:
            self.numeric.reverse(self.later)
        self.assertIn('2 are needed', str(caught.exception))

    def test_too_few_says_what_is_free(self):
        with self.assertRaises(InsufficientKnownValues) as caught:
            self.numeric.reverse(self.later, known_values={self.k_1: 1.5})
        message = str(caught.exception)
        self.assertIn('1 more value is needed', message)
        self.assertIn('leaves k_1 unchanged', message)

    def test_degenerate_pair_is_insufficient_not_conflicting(self):
        # k_m1 and k_2 only fix k_1*s_0 between them.
        with self.assertRaises(InsufficientKnownValues) as caught:
            self.numeric.reverse(self.later, known_values={self.k_m1: 0.9, self.k_2: 0.4})
        suggestion = str(caught.exception).splitlines()[-1]
        self.assertIn('k_1', suggestion)
        self.assertNotIn('k_m1', suggestion)
        self.assertNotIn('k_2', suggestion)

    def test_invariant_ratio_pair_is_insufficient(self):
        # e_0 / s_0 is itself an invariant, so the pair fixes only one scale.
        with self.assertRaises(InsufficientKnownValues):
            self.numeric.reverse(self.later, known_values={self.e_0: 0.3, self.s_0: 2.0})

    def test_too_many_but_consistent(self):
        recovered = self.numeric.reverse(
            self.later, known_values={self.k_1: 1.5, self.s_0: 2.0, self.k_2: 0.4,
                                      self.e_0: 0.3})
        self.assertRecovers(recovered, k_m1=0.9)

    def test_too_many_and_conflicting(self):
        with self.assertRaises(ConflictingKnownValues) as caught:
            self.numeric.reverse(self.later,
                                 known_values={self.k_1: 1.5, self.s_0: 2.0, self.k_2: 99.0})
        message = str(caught.exception)
        self.assertIn('k_2 was given as 99', message)
        self.assertIn('imply 0.4', message)

    # ------------------------------------------------------------------ bad values
    def test_zero_known_value(self):
        with self.assertRaises(ValueError) as caught:
            self.numeric.reverse(self.start, known_values={self.k_1: 1.5, self.c: 0.0})
        self.assertIn('zero', str(caught.exception))

    def test_zero_reduced_value_behind_a_known_value(self):
        # c(0) = 0 so v = 0 at the start; a nonzero c there is impossible to reconcile.
        with self.assertRaises(ValueError) as caught:
            self.numeric.reverse(self.start, known_values={self.k_1: 1.5, self.c: 0.1})
        self.assertIn('zero', str(caught.exception))

    def test_array_as_known_value(self):
        with self.assertRaises(ValueError) as caught:
            self.numeric.reverse(self.later,
                                 known_values={self.k_1: numpy.array([1.5]), self.s_0: 2.0})
        self.assertIn('single number', str(caught.exception))

    def test_known_values_and_auxiliaries_together(self):
        with self.assertRaises(ValueError):
            self.numeric.reverse(self.later, known_values={self.k_1: 1.5, self.s_0: 2.0},
                                 auxiliaries=[1.0, 2.0])

    def test_wrong_number_of_auxiliaries(self):
        with self.assertRaises(ValueError):
            self.numeric.reverse(self.later, auxiliaries=[1.0])

    def test_known_values_that_fix_the_scale_only_through_a_square(self):
        # dz/dt = a z - b z**3 is odd in z, and b scales as 1/lambda**2: knowing a and b
        # fixes the size of z0 but not its sign, and that is refused rather than guessed.
        system = ODESystem.from_equations(['dz/dt = a*z - b*z**3'])
        system.update_initial_conditions({'z': 'z0'})
        # update_initial_conditions appends z0 at the end, which leaves the system out of
        # canonical order now that the parameters come first, so put it back.
        system.reorder_variables(['a', 'b', 'z0', 't', 'z'])
        numeric = NumericTranslator(ODETranslator.from_ode_system(system).translate())
        a, b, z0, t, z = system.symbols
        point = numeric.forward({t: 0.0, z: 1.5, a: 2.0, b: 0.5, z0: 1.5})
        with self.assertRaises(NotImplementedError) as caught:
            numeric.reverse(point, known_values={a: 2.0, b: 0.5})
        self.assertIn('roots of degree 2', str(caught.exception))
        # A variable the scaling acts on with power 1 fixes it outright.
        recovered = numeric.reverse(point, known_values={a: 2.0, z0: 1.5})
        numpy.testing.assert_allclose(float(recovered[b]), 0.5, rtol=1e-12)

    def test_known_values_refused_when_reduced_system_carries_auxiliaries(self):
        equations = ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)']
        system = ODESystem.from_equations(equations)
        system.reorder_variables(['t', 'z1', 'z2'])
        translator = ODETranslator.from_ode_system(system)
        numeric = NumericTranslator(translator.translate())
        self.assertTrue(numeric.carries_auxiliaries)
        t, z1, z2 = system.symbols
        point = numeric.forward({t: 1.0, z1: 0.4, z2: 0.3})
        with self.assertRaises(ValueError) as caught:
            numeric.reverse(point, known_values={z1: 0.4})
        self.assertIn('already determined', str(caught.exception))
        recovered = numeric.reverse(point)
        numpy.testing.assert_allclose(float(recovered[z1]), 0.4, rtol=1e-12)


class TestTheNumericTranslatorNeedsOnlyTheReduction(TestCase):
    '''
    It used to take a system, a translator and a reduced system: three things that had to
    agree, with nothing making them agree.  A reduced system now carries the system it came
    from and the whole change of variables, so one argument says all of it.
    '''

    EQUATIONS = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)']

    def reduction(self):
        system = ODESystem.from_equations(self.EQUATIONS)
        return ODETranslator.from_ode_system(system).translate()

    def test_the_reduced_system_makes_one(self):
        reduced = self.reduction()
        numeric = reduced.numeric_translator()
        self.assertIs(reduced, numeric.reduced)
        self.assertEqual(tuple(reduced.original_system.symbols),
                         tuple(numeric.system.symbols))
        self.assertEqual(reduced.scheme, numeric.scheme)

    def test_it_is_the_same_as_the_constructor(self):
        reduced = self.reduction()
        values = {sympy.Symbol(name): value for name, value in
                  (('K', 5.0), ('d', 0.7), ('h', 1.3), ('k', 0.9), ('r', 1.1),
                   ('s', 0.4), ('t', 2.0), ('n', 1.0), ('p', 0.5))}
        self.assertEqual(NumericTranslator(reduced).forward(values),
                         reduced.numeric_translator().forward(values))

    def test_the_matrices_come_from_the_record(self):
        '''
        The exponent vectors of the recorded monomials are the multiplier, so what is
        rebuilt here is the reduction that was performed, not a fresh guess at it.
        '''
        system = ODESystem.from_equations(self.EQUATIONS)
        translator = ODETranslator.from_ode_system(system)
        numeric = translator.translate().numeric_translator()
        self.assertEqual(translator.herm_mult_i, numeric._herm_mult_i)
        self.assertEqual(translator.herm_mult_n, numeric._herm_mult_n)
        self.assertEqual(translator.inv_herm_mult, numeric._inv_herm_mult)

    def test_a_translator_adjusted_afterwards_cannot_reach_it(self):
        '''
        The hazard this shape removes.  A column operation changes which invariants a
        translator produces; the reduced system it made earlier does not change with it.
        There is no longer anywhere to hand the moved-on translator in.
        '''
        system = ODESystem.from_equations(self.EQUATIONS)
        translator = ODETranslator.from_ode_system(system)
        reduced = translator.translate()
        before = reduced.numeric_translator()._herm_mult_n

        translator.multiplier_add_columns(translator.r + 0, translator.r + 3, 1)
        self.assertNotEqual(translator.herm_mult_n, before)
        # the reduction in hand is unmoved, and still describes itself
        self.assertEqual(before, reduced.numeric_translator()._herm_mult_n)

    def test_a_system_with_no_record_is_refused(self):
        '''
        A reduced system assembled by hand has no record of what its symbols mean, so
        nothing can work out what to translate.  Refused, rather than guessed at by
        matching column counts -- which is how `forward` once reported one column's value
        against another column's symbol.
        '''
        by_hand = ODESystem.from_equations(['dy0/dt = y0'], is_reduced=True)
        self.assertIsNone(by_hand.named_symbols)
        with self.assertRaises(ValueError) as caught:
            by_hand.numeric_translator()
        self.assertIn('no record', str(caught.exception))

    def test_a_partially_normalised_reduction_round_trips(self):
        '''
        The case the old count-matching got silently wrong: six kept symbols over nine
        columns.  Given back the parameters it normalised, it recovers the original
        exactly.
        '''
        enzyme = ChemicalReactionNetwork.from_diagram(
            'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
        reduced = ODETranslator.from_ode_system(enzyme).translate()
        numeric = reduced.numeric_translator()

        values = {sympy.Symbol(name): value for name, value in
                  (('k_0_1', 2.0), ('k_2_3', 3.0), ('k_4_5', 5.0), ('t', 1.5),
                   ('e', 0.7), ('f', 1.1), ('x', 0.9), ('y', 0.4), ('z', 0.2))}
        forward = numeric.forward(values)
        self.assertAlmostEqual(0.9, float(forward[sympy.Symbol('x3')]))

        recovered = numeric.reverse(forward, known_values={
            sympy.Symbol('k_0_1'): 2.0, sympy.Symbol('k_2_3'): 3.0,
            sympy.Symbol('k_4_5'): 5.0})
        for symbol, value in values.items():
            with self.subTest(symbol=str(symbol)):
                self.assertAlmostEqual(value, float(recovered[symbol]))
