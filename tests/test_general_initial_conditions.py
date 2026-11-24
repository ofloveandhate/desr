'''
Initial conditions under the general scheme.

A reduced symbol denotes a monomial in the original symbols.  Its initial value is that
monomial at the start: the independent variable at 0, each dependent variable at its
initial value.  The tests state that, and the cases where no such value exists.
'''
from unittest import TestCase

import sympy

from desr import ODESystem, ODETranslator


def michaelis_menten(initial_conditions):
    system = ODESystem.from_equations([
        'ds/dt = -k_1*e_0*s + (k_1*s + k_m1)*c',
        'dc/dt = k_1*e_0*s - (k_1*s + k_m1 + k_2)*c'])
    system.update_initial_conditions(initial_conditions)
    return system


def wrong_initial_values(system, reduced):
    '''
    The reduced initial conditions that are not what the definition gives.

    For a reduced symbol `y` that denotes the monomial `m`, the definition gives
    `m` with `t = 0` and each dependent at its initial value.  The reduced system states
    `y(0)` in its own parameters, so those are written back in the original symbols
    before the two are compared.  A dependent with no initial condition is left as it is,
    so it cannot be mistaken for agreement.
    '''
    denotes = reduced.named_symbols
    at_the_start = dict(system.initial_conditions)
    at_the_start[system.independent_variable] = 0
    wrong = []
    for symbol, stated in reduced.initial_conditions.items():
        defined = denotes[symbol].subs(at_the_start, simultaneous=True)
        in_original_symbols = sympy.sympify(stated).subs(denotes, simultaneous=True)
        if sympy.simplify(defined - in_original_symbols) != 0:
            wrong.append((symbol, stated, defined))
    return wrong


class TestTheInitialValueIsTheMonomialAtTheStart(TestCase):

    OPTIONS = ({},
               {'rescale_clock': False},
               {'normalise_parameters': False},
               {'normalise_parameters': False, 'rescale_clock': False},
               {'include_aux_vars': False})

    def systems(self):
        yield 'michaelis-menten, both', michaelis_menten({'s': 's_0', 'c': 'c_0'})
        yield 'michaelis-menten, one zero', michaelis_menten({'s': 's_0', 'c': 0})
        yield 'michaelis-menten, one', michaelis_menten({'s': 's_0'})
        decay = ODESystem.from_equations(['dz/dt = -k*z'])
        decay.update_initial_conditions({'z': 'z_0'})
        yield 'decay', decay
        # The scaling acts on t, and the initial value is the only parameter there is.
        square = ODESystem.from_equations(['dx/dt = x**2'])
        square.update_initial_conditions({'x': 'x_0'})
        yield 'dx/dt = x**2', square
        # Hubert & Labahn's example 6.6 shape: the scaling does not act on t.
        pair = ODESystem.from_equations(
            ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)'])
        pair.update_initial_conditions({'z1': 'a', 'z2': 'b'})
        yield 'two dependents scaled against each other', pair

    def test_every_stated_initial_value_is_the_one_the_definition_gives(self):
        for label, system in self.systems():
            translator = ODETranslator.from_ode_system(system)
            for options in self.OPTIONS:
                with self.subTest(system=label, **options):
                    reduced = translator.translate_general(system, **options)
                    self.assertEqual([], wrong_initial_values(system, reduced))

    def test_every_dependent_with_an_initial_condition_reaches_the_reduced_system(self):
        # Each initial value given appears in what some reduced initial condition was
        # worked out from: none is dropped on the way.
        for label, system in self.systems():
            translator = ODETranslator.from_ode_system(system)
            with self.subTest(system=label):
                reduced = translator.translate_general(system)
                denotes = reduced.named_symbols
                reached = set()
                for symbol in reduced.initial_conditions:
                    reached |= denotes[symbol].free_symbols
                self.assertLessEqual(set(system.initial_conditions), reached)

    def test_the_oracle_can_fail(self):
        system = michaelis_menten({'s': 's_0', 'c': 'c_0'})
        reduced = ODETranslator.from_ode_system(system).translate_general(system)
        symbol = next(iter(reduced.initial_conditions))
        reduced.update_initial_conditions({symbol: 2 * reduced.initial_conditions[symbol]})
        self.assertEqual(1, len(wrong_initial_values(system, reduced)))

    def test_a_clock_invariant_kept_beside_t_starts_at_zero(self):
        # With the clock not rescaled, e_0 k_1 t stays as an invariant with an equation.
        # It is t times parameters, so at t = 0 it is 0.
        system = michaelis_menten({'s': 's_0'})
        reduced = ODETranslator.from_ode_system(system).translate_general(
            system, rescale_clock=False)
        t = system.independent_variable
        clocks = [symbol for symbol in reduced.dependents
                  if reduced.named_symbols[symbol].has(t)]
        self.assertEqual(1, len(clocks))
        self.assertEqual(0, reduced.initial_conditions[clocks[0]])

    def test_a_zero_initial_value_stays_a_value(self):
        # c(0) = 0 is the value of c_0.  The reduced initial condition is a reduced
        # parameter, and that parameter's value is 0.
        system = michaelis_menten({'s': 's_0', 'c': 0})
        reduced = ODETranslator.from_ode_system(system).translate_general(system)
        c = sympy.Symbol('c')
        symbol = next(s for s in reduced.dependents if reduced.named_symbols[s].has(c))
        stated = reduced.initial_conditions[symbol]
        self.assertEqual(0, reduced.parameter_values[stated])

    def test_no_initial_conditions_gives_none(self):
        system = ODESystem.from_equations(['dz/dt = -k*z'])
        reduced = ODETranslator.from_ode_system(system).translate_general(system)
        self.assertEqual({}, reduced.initial_conditions)


class TestItAgreesWithTheParameterScheme(TestCase):
    '''
    Where both schemes apply they make the same change of variables under different
    names, so they state the same initial values.
    '''

    def test_michaelis_menten(self):
        for initial_conditions in ({'s': 's_0', 'c': 'c_0'}, {'s': 's_0', 'c': 0},
                                   {'s': 's_0'}):
            with self.subTest(initial_conditions=initial_conditions):
                system = michaelis_menten(initial_conditions)
                translator = ODETranslator.from_ode_system(system)
                general = translator.translate_general(system)
                parameter = translator.translate_parameter(system)

                def in_original_symbols(reduced):
                    denotes = reduced.named_symbols
                    return {denotes[symbol]: sympy.sympify(value).subs(denotes)
                            for symbol, value in reduced.initial_conditions.items()}

                self.assertEqual(in_original_symbols(parameter),
                                 in_original_symbols(general))
                self.assertEqual(len(initial_conditions), len(general.initial_conditions))


class TestSeveralDependentsInOneInvariant(TestCase):

    def mixed(self, initial_conditions):
        # Multiply the invariant c/e_0 by the invariant s/e_0.
        system = michaelis_menten(initial_conditions)
        translator = ODETranslator.from_ode_system(system)
        c, s, e_0 = sympy.symbols('c s e_0')
        translator.multiplier_add_columns(c / e_0, s / e_0, power=1)
        self.assertIn(c * s / e_0 ** 2, list(translator.invariant_monomials()))
        return system, translator

    def test_with_all_of_them_given_it_has_a_value(self):
        system, translator = self.mixed({'s': 's_0', 'c': 'c_0'})
        reduced = translator.translate_general(system)
        self.assertEqual([], wrong_initial_values(system, reduced))
        self.assertEqual(2, len(reduced.initial_conditions))

    def test_with_only_some_given_it_is_refused(self):
        system, translator = self.mixed({'s': 's_0'})
        with self.assertRaises(ValueError) as caught:
            translator.translate_general(system)
        message = str(caught.exception)
        self.assertIn('s, which has an initial condition', message)
        self.assertIn('c, which does not', message)

    def test_a_zero_in_a_denominator_is_refused(self):
        # s/c is not defined where c is 0.
        system = michaelis_menten({'s': 's_0', 'c': 0})
        translator = ODETranslator.from_ode_system(system)
        c, s, e_0 = sympy.symbols('c s e_0')
        translator.multiplier_add_columns(s / e_0, c / e_0, power=-1)
        with self.assertRaises(ValueError) as caught:
            translator.translate_general(system)
        self.assertIn('divides by c, whose initial value is 0', str(caught.exception))

    def test_a_zero_in_a_numerator_is_a_zero(self):
        system = michaelis_menten({'s': 's_0', 'c': 0})
        translator = ODETranslator.from_ode_system(system)
        c, s, e_0 = sympy.symbols('c s e_0')
        translator.multiplier_add_columns(c / e_0, s / e_0, power=1)
        reduced = translator.translate_general(system)
        self.assertEqual([], wrong_initial_values(system, reduced))


class TestWhatDoesNotCarryOver(TestCase):

    def test_the_independent_variable_in_an_invariant_with_a_dependent(self):
        # Multiply c/e_0 by the invariant e_0 k_1 t.  The product is 0 at t = 0 whatever
        # c(0) is.
        system = michaelis_menten({'s': 's_0', 'c': 'c_0'})
        translator = ODETranslator.from_ode_system(system)
        c, e_0, k_1, t = sympy.symbols('c e_0 k_1 t')
        translator.multiplier_add_columns(c / e_0, e_0 * k_1 * t, power=1)
        self.assertIn(c * k_1 * t, list(translator.invariant_monomials()))
        for options in ({}, {'rescale_clock': False}):
            with self.subTest(**options):
                with self.assertRaises(ValueError) as caught:
                    translator.translate_general(system, **options)
                message = str(caught.exception)
                self.assertIn('initial condition for c', message)
                self.assertIn('involves t as well', message)
                self.assertIn('isolate_clock', message)

    def test_isolating_the_clock_is_the_fix_the_message_names(self):
        system = michaelis_menten({'s': 's_0', 'c': 'c_0'})
        translator = ODETranslator.from_ode_system(system)
        c, e_0, k_1, t = sympy.symbols('c e_0 k_1 t')
        translator.multiplier_add_columns(c / e_0, e_0 * k_1 * t, power=1)
        translator.isolate_clock(system)
        reduced = translator.translate_general(system)
        self.assertEqual([], wrong_initial_values(system, reduced))
        self.assertEqual(2, len(reduced.initial_conditions))

    def test_ordering_the_initial_value_first_is_the_other_fix(self):
        # dx/dt = x**2 with t ahead of x_0: the invariants are t x and t x_0.
        system = ODESystem.from_equations(['dx/dt = x**2'])
        system.update_initial_conditions({'x': 'x_0'})
        system.reorder_variables(['t', 'x', 'x_0'])
        with self.assertRaises(ValueError) as caught:
            ODETranslator.from_ode_system(system).translate_general(system)
        self.assertIn('put x_0 ahead of t', str(caught.exception))
        system.reorder_variables(['x_0', 't', 'x'])
        reduced = ODETranslator.from_ode_system(system).translate_general(system)
        self.assertEqual([], wrong_initial_values(system, reduced))
        self.assertEqual(1, len(reduced.initial_conditions))

    def test_the_same_invariant_is_no_obstacle_without_that_initial_condition(self):
        system = michaelis_menten({'s': 's_0'})
        translator = ODETranslator.from_ode_system(system)
        c, e_0, k_1, t = sympy.symbols('c e_0 k_1 t')
        translator.multiplier_add_columns(c / e_0, e_0 * k_1 * t, power=1)
        reduced = translator.translate_general(system)
        self.assertEqual([], wrong_initial_values(system, reduced))
        self.assertEqual(1, len(reduced.initial_conditions))

    def test_a_dependent_used_as_a_unit(self):
        # dz1/dt = z1 (1 + z1 z2), dz2/dt = z2 (1/t - z1 z2): the one symmetry scales z1
        # against z2 and there is no parameter, so one of the two is the unit.  With the
        # initial values' symbols put after the dependents, that stays so when they are
        # added.
        system = ODESystem.from_equations(
            ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)'])
        bare = ODETranslator.from_ode_system(system)
        self.assertTrue(any(monomial.free_symbols & set(system.dependents)
                            for monomial in bare.auxiliary_monomials()))
        system.update_initial_conditions({'z1': 'a', 'z2': 'b'})
        system.reorder_variables(['t', 'z1', 'z2', 'a', 'b'])
        translator = ODETranslator.from_ode_system(system)
        self.assertTrue(any(monomial.free_symbols & set(system.dependents)
                            for monomial in translator.auxiliary_monomials()))
        with self.assertRaises(ValueError) as caught:
            translator.translate_general(system)
        self.assertIn('uses z1 as a unit', str(caught.exception))
        self.assertIn('Put a ahead of it', str(caught.exception))
        # The fix the message names.
        system.reorder_variables(['a', 'b', 't', 'z1', 'z2'])
        reduced = ODETranslator.from_ode_system(system).translate_general(system)
        self.assertEqual([], wrong_initial_values(system, reduced))
        self.assertEqual(2, len(reduced.initial_conditions))


class TestTheDependentVariableScheme(TestCase):
    '''
    The scaling does not act on the independent variable.  The reduced symbols are
    monomials in the other symbols, so they start by the same rule.
    '''

    EQUATIONS = ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)']

    def system(self, initial_conditions, order=None):
        system = ODESystem.from_equations(self.EQUATIONS)
        system.update_initial_conditions(initial_conditions)
        if order:
            system.reorder_variables(order)
        return system

    def test_every_stated_initial_value_is_the_one_the_definition_gives(self):
        for initial_conditions in ({'z1': 'a', 'z2': 'b'}, {'z1': 'a'}, {'z2': 'b'}):
            for options in ({}, {'normalise_parameters': False}):
                with self.subTest(initial_conditions=initial_conditions, **options):
                    system = self.system(initial_conditions)
                    reduced = ODETranslator.from_ode_system(system).translate_dep_var(
                        system, **options)
                    self.assertEqual([], wrong_initial_values(system, reduced))
                    reached = set()
                    for symbol in reduced.initial_conditions:
                        reached |= reduced.named_symbols[symbol].free_symbols
                    self.assertLessEqual(set(system.initial_conditions), reached)

    def test_it_states_what_the_general_scheme_states(self):
        system = self.system({'z1': 'a', 'z2': 'b'})
        translator = ODETranslator.from_ode_system(system)

        def in_original_symbols(reduced):
            denotes = reduced.named_symbols
            return {denotes[symbol]: sympy.sympify(value).subs(denotes)
                    for symbol, value in reduced.initial_conditions.items()}

        dep_var = translator.translate_dep_var(system)
        self.assertEqual(2, len(dep_var.initial_conditions))
        self.assertEqual(in_original_symbols(translator.translate_general(system)),
                         in_original_symbols(dep_var))

    def test_an_initial_value_takes_the_dependents_place_as_the_unit(self):
        # Without initial conditions the one symmetry scales z1 against z2 and there is
        # no parameter, so this is the dependent-variable scheme.  An initial value is a
        # parameter that scales as its dependent does, so with one given the default
        # order makes it the unit and `translate` takes the parameter scheme.
        from desr import ReductionScheme
        bare = ODESystem.from_equations(self.EQUATIONS)
        self.assertEqual(ReductionScheme.DEP_VAR,
                         ODETranslator.from_ode_system(bare).scheme())
        system = self.system({'z1': 'a', 'z2': 'b'})
        translator = ODETranslator.from_ode_system(system)
        self.assertEqual(ReductionScheme.PARAMETER, translator.scheme())
        reduced = translator.translate(system)
        self.assertEqual([], wrong_initial_values(system, reduced))
        self.assertEqual(2, len(reduced.initial_conditions))

    def test_a_dependent_used_as_a_unit_is_refused_and_the_fix_works(self):
        system = self.system({'z1': 'a', 'z2': 'b'}, order=['t', 'z1', 'z2', 'a', 'b'])
        with self.assertRaises(ValueError) as caught:
            ODETranslator.from_ode_system(system).translate_dep_var(system)
        self.assertIn('uses z1 as a unit', str(caught.exception))
        system.reorder_variables(['a', 'b', 't', 'z1', 'z2'])
        reduced = ODETranslator.from_ode_system(system).translate_dep_var(system)
        self.assertEqual([], wrong_initial_values(system, reduced))
