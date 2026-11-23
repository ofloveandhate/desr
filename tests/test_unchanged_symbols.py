'''
A symbol the reduction leaves unchanged keeps its name.

An invariant is a monomial in the original symbols.  When that monomial is a single
symbol, the scaling does not act on that symbol -- its column of the scaling matrix is
zero -- and the reduction carries it through as it stands.  A generated name would report
a change that did not happen.

An auxiliary is different, and is always given a new name: being an auxiliary is a role
the symbol did not have in the original system.

The statements here are about what a name denotes, read from `named_symbols`, and about
the scaling matrix.  Whether the reduced equations are true is asked of the chain rule.
'''
from unittest import TestCase

import sympy

from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator, ReductionScheme

from oracles import false_equations, forget

S = sympy.Symbol

SYSTEMS = {
    'no symmetry': lambda: ODESystem.from_equations(['dz/dt = z + 1 + t']),
    'a dimensionless parameter': lambda: ODESystem.from_equations(['dp/dt = p*(q - 1)']),
    'dependent-variable scheme with a parameter': lambda: ODESystem.from_equations(
        ['dz1/dt = a*z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)']),
    'Michaelis-Menten': lambda: ODESystem.from_equations(
        ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
         'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']),
    'predator-prey': lambda: ODESystem.from_equations(
        ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)']),
    'Hubert & Labahn 6.4': lambda: ODESystem.from_equations(
        ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)']),
    'Hubert & Labahn 6.6': lambda: ODESystem.from_equations(
        ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
         'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)']),
    'enzyme chain': lambda: ChemicalReactionNetwork.from_diagram(
        'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system(),
    'decay': lambda: ODESystem.from_equations(['dn/dt = -n/tc']),
}

RENAMING = {'rename_unchanged_symbols': True}


class TestTheRule(TestCase):

    def reductions(self):
        for label, build in SYSTEMS.items():
            system = build()
            translator = ODETranslator.from_ode_system(system)
            yield label, system, translator, translator.translate(system)

    def test_a_reduced_symbol_shares_a_name_only_with_the_symbol_it_denotes(self):
        # In both directions.  A name the two systems share means one quantity, and a
        # non-auxiliary reduced symbol that denotes an original symbol has its name.
        for label, system, translator, reduced in self.reductions():
            with self.subTest(system=label):
                original = set(system.symbols)
                auxiliaries = set(reduced.auxiliary_names)
                for symbol, monomial in reduced.named_symbols.items():
                    if symbol in original:
                        self.assertEqual(symbol, monomial)
                    if monomial in original and symbol not in auxiliaries:
                        self.assertEqual(monomial, symbol)

    def test_a_symbol_keeps_its_name_exactly_when_the_scaling_leaves_it_alone(self):
        # The premise, read off the scaling matrix rather than off the names: a kept
        # name belongs to a symbol whose column of A is zero.
        #
        # Of every symbol but the independent variable.  The general scheme retains the
        # original clock to integrate against, so its `t` is the original `t` even where
        # the scaling acts on it -- Hubert & Labahn 6.6, weight 3.  That is the clock's
        # own rule, tested with the clock.
        for label, system, translator, reduced in self.reductions():
            with self.subTest(system=label):
                A = translator.scaling_matrix
                clock = system.independent_variable
                unscaled = {symbol for column, symbol in enumerate(system.symbols)
                            if A[:, column].is_zero_matrix}
                kept = set(system.symbols) & set(reduced.symbols)
                self.assertEqual(unscaled - {clock}, kept - {clock})

    def test_an_auxiliary_is_always_given_a_new_name(self):
        for label, system, translator, reduced in self.reductions():
            with self.subTest(system=label):
                self.assertEqual(set(), set(reduced.auxiliary_names) & set(system.symbols))

    def test_the_equations_are_true_under_either_naming(self):
        for label, build in SYSTEMS.items():
            for naming in (None, RENAMING):
                with self.subTest(system=label, naming=naming):
                    system = build()
                    reduced = ODETranslator.from_ode_system(
                        system, naming_scheme=naming).translate(system)
                    self.assertEqual([], false_equations(system, reduced))

    def test_the_two_namings_differ_only_in_names(self):
        # Same monomials in the same order; the switch changes spelling and nothing else.
        for label, build in SYSTEMS.items():
            with self.subTest(system=label):
                system = build()
                kept = ODETranslator.from_ode_system(system).translate(system)
                renamed = ODETranslator.from_ode_system(
                    system, naming_scheme=RENAMING).translate(system)
                self.assertEqual(list(kept.named_symbols.values()),
                                 list(renamed.named_symbols.values()))


class TestWhatKeepsItsName(TestCase):

    def test_every_symbol_of_a_system_with_no_symmetry(self):
        system = SYSTEMS['no symmetry']()
        reduced = ODETranslator.from_ode_system(system).translate(system)
        self.assertEqual({symbol: symbol for symbol in system.symbols},
                         dict(reduced.named_symbols))
        self.assertEqual(system, forget(reduced))

    def test_a_dimensionless_parameter(self):
        system = SYSTEMS['a dimensionless parameter']()
        translator = ODETranslator.from_ode_system(system)
        q = S('q')
        # Its column of the scaling matrix is zero, which is what dimensionless means.
        column = list(system.symbols).index(q)
        self.assertTrue(translator.scaling_matrix[:, column].is_zero_matrix)

        reduced = translator.translate(system)
        self.assertEqual(q, reduced.named_symbols[q])
        self.assertIn(q, reduced.parameters)

    def test_the_parameter_schemes_clock_when_it_is_not_rescaled(self):
        # With no symmetry there is nothing to rescale time by, and the clock is `t`.
        system = SYSTEMS['no symmetry']()
        translator = ODETranslator.from_ode_system(system)
        self.assertEqual(ReductionScheme.PARAMETER, translator.scheme())
        self.assertEqual(system.independent_variable,
                         translator.translate(system).independent_variable)

    def test_the_parameter_schemes_clock_when_it_is_rescaled(self):
        # Michaelis-Menten's clock is e_0*k_1*t, a different quantity, and is `tau`.
        system = SYSTEMS['Michaelis-Menten']()
        reduced = ODETranslator.from_ode_system(system).translate(system)
        self.assertEqual(S('tau'), reduced.independent_variable)
        self.assertNotEqual(system.independent_variable,
                            reduced.named_symbols[reduced.independent_variable])

    def test_a_generated_name_keeps_its_index(self):
        # `a` keeps its name and the next invariant is still `y1`: an index is a column.
        system = SYSTEMS['dependent-variable scheme with a parameter']()
        reduced = ODETranslator.from_ode_system(system).translate(system)
        self.assertEqual((S('t'), S('x0'), S('a'), S('y1')), tuple(reduced.symbols))


class TestWhatDoesNot(TestCase):

    def test_an_auxiliary_that_is_one_of_the_original_symbols(self):
        # The enzyme chain's surviving auxiliary denotes `x` and is not called `x`.
        system = SYSTEMS['enzyme chain']()
        reduced = ODETranslator.from_ode_system(system).translate(system)
        auxiliary, = reduced.auxiliary_names
        self.assertEqual(S('x'), reduced.named_symbols[auxiliary])
        self.assertNotEqual(S('x'), auxiliary)

    def test_an_invariant_that_is_the_clock_while_the_clock_is_also_carried(self):
        # Hubert & Labahn 6.4 through the general scheme, with identification declined,
        # carries the clock and an invariant whose monomial is `t`.  One name cannot be
        # two variables, so that invariant keeps its generated name.
        system = SYSTEMS['Hubert & Labahn 6.4']()
        translator = ODETranslator.from_ode_system(system)
        reduced = translator.translate_general(system, rescale_clock=False)
        t = system.independent_variable
        denoting_t = [symbol for symbol, monomial in reduced.named_symbols.items()
                      if monomial == t]
        self.assertEqual(2, len(denoting_t))
        self.assertEqual(len(reduced.symbols), len(set(reduced.symbols)))
        self.assertEqual([], false_equations(system, reduced))

    def test_names_the_caller_wrote_out(self):
        system = SYSTEMS['a dimensionless parameter']()
        reduced = ODETranslator.from_ode_system(
            system, naming_scheme={'dependents': ['Q']}).translate(system)
        self.assertEqual(S('q'), reduced.named_symbols[S('Q')])
        self.assertNotIn(S('q'), reduced.symbols)


class TestTheSwitch(TestCase):

    def test_renaming_gives_every_invariant_a_generated_name(self):
        for label in ('no symmetry', 'a dimensionless parameter',
                      'dependent-variable scheme with a parameter'):
            with self.subTest(system=label):
                system = SYSTEMS[label]()
                reduced = ODETranslator.from_ode_system(
                    system, naming_scheme=RENAMING).translate(system)
                shared = (set(system.symbols) & set(reduced.symbols)
                          - {system.independent_variable})
                self.assertEqual(set(), shared)

    def test_the_general_schemes_clock_follows_its_denotation_either_way(self):
        # The switch is about invariants the scheme names.  The general and
        # dependent-variable schemes' clock keeps `t` whenever it denotes `t`,
        # which is a statement about the clock and was not made optional.
        system = SYSTEMS['Hubert & Labahn 6.6']()
        reduced = ODETranslator.from_ode_system(
            system, naming_scheme=RENAMING).translate(system)
        self.assertEqual(system.independent_variable, reduced.independent_variable)

    def test_it_must_be_true_or_false(self):
        system = SYSTEMS['no symmetry']()
        with self.assertRaises(ValueError) as caught:
            ODETranslator.from_ode_system(
                system, naming_scheme={'rename_unchanged_symbols': 'yes'})
        self.assertIn('True or False', str(caught.exception))

    def test_it_is_not_mistaken_for_a_name(self):
        # It sits in the naming scheme beside the names, and the check for names the
        # caller wrote out must pass over it.
        system = SYSTEMS['Michaelis-Menten']()
        for value in (True, False):
            with self.subTest(rename_unchanged_symbols=value):
                reduced = ODETranslator.from_ode_system(
                    system,
                    naming_scheme={'rename_unchanged_symbols': value}).translate(system)
                self.assertEqual([], false_equations(system, reduced))


class TestAKeptNameMeetingAGeneratedOne(TestCase):
    '''
    A model whose own dimensionless parameter is called `y1`, reduced by a scheme that
    generates a `y1`.  The reduced system would have two symbols with one name.
    '''

    EQUATIONS = ['dz1/dt = y1*z1*(z1**5*z2 - 2)/(3*t)',
                 'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)']

    def test_it_is_refused(self):
        system = ODESystem.from_equations(self.EQUATIONS)
        with self.assertRaises(ValueError) as caught:
            ODETranslator.from_ode_system(system).translate(system)
        message = str(caught.exception)
        self.assertIn('y1 twice', message)
        self.assertIn('patterns avoid y1', message)

    def test_a_pattern_that_avoids_the_name_is_the_way_out(self):
        system = ODESystem.from_equations(self.EQUATIONS)
        reduced = ODETranslator.from_ode_system(
            system, naming_scheme={'dependents': 'w'}).translate(system)
        self.assertEqual(len(reduced.symbols), len(set(reduced.symbols)))
        self.assertEqual(S('y1'), reduced.named_symbols[S('y1')])
        self.assertEqual([], false_equations(system, reduced))

    def test_renaming_it_is_not(self):
        # Renamed, the parameter no longer holds the name `y1` -- but the scheme still
        # generates a `y1` for another quantity, and the original still has one.  That
        # is the other refusal: one spelling across two systems for two quantities.
        system = ODESystem.from_equations(self.EQUATIONS)
        with self.assertRaises(ValueError) as caught:
            ODETranslator.from_ode_system(system, naming_scheme=RENAMING).translate(system)
        self.assertIn('two different quantities', str(caught.exception))


class TestReducingTwiceWithNoNamingScheme(TestCase):
    '''
    With unchanged symbols keeping their names, a reduced system can be reduced again
    without being told what to call anything.
    '''

    def test_the_second_reduction_keeps_every_invariant(self):
        for label in ('Hubert & Labahn 6.4', 'Hubert & Labahn 6.6', 'enzyme chain',
                      'decay', 'Michaelis-Menten', 'predator-prey'):
            with self.subTest(system=label):
                system = SYSTEMS[label]()
                first = ODETranslator.from_ode_system(system).translate(system)
                once = forget(first)
                second = ODETranslator.from_ode_system(once).translate(once)
                self.assertEqual([], false_equations(once, second))
                self.assertEqual(len(once.symbols), len(second.symbols))
                carried = set(first.auxiliary_names)
                for symbol in once.symbols:
                    if symbol not in carried:
                        self.assertIn(symbol, second.symbols)
                        self.assertEqual(symbol, second.named_symbols[symbol])

    def test_the_third_is_the_second(self):
        for label in ('Hubert & Labahn 6.4', 'Hubert & Labahn 6.6', 'enzyme chain',
                      'decay', 'Michaelis-Menten', 'predator-prey'):
            with self.subTest(system=label):
                system = SYSTEMS[label]()
                once = forget(ODETranslator.from_ode_system(system).translate(system))
                twice = forget(ODETranslator.from_ode_system(once).translate(once))
                thrice = forget(ODETranslator.from_ode_system(twice).translate(twice))
                self.assertEqual(twice, thrice)
