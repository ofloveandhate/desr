'''
The three reduction schemes agree where more than one applies.

The general scheme always applies.  Its two further steps -- normalising the constant
auxiliaries and identifying a clock invariant with the clock -- are what the other
schemes do by construction, so on a system one of those fits, the general scheme with
both steps gives the same reduced system up to the names of its symbols.  The second
step is a change of independent variable: an invariant that is the clock times
parameters is taken as the clock, whatever the first step did.  With both steps off it
is Hubert & Labahn's reduction, checked here against their printed examples.

"The same up to names" is decided through the monomials: two reduced symbols are the
same quantity when they denote the same monomial in the original symbols.
'''
from unittest import TestCase

import sympy

from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator, ReductionScheme

from oracles import false_equations


MICHAELIS_MENTEN = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                    'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
HL_6_4 = ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)']


def disagreements(first, second):
    '''
    Where two reductions of one system differ, compared quantity by quantity.

    Returns a list of descriptions; empty when they are the same system up to names.
    '''
    name_in_second = {monomial: name for name, monomial in second.named_symbols.items()}
    found = []
    if set(first.named_symbols.values()) != set(name_in_second):
        return ['they do not denote the same monomials: {} against {}'.format(
            sorted(map(str, first.named_symbols.values())),
            sorted(map(str, name_in_second)))]
    renaming = {name: name_in_second[monomial]
                for name, monomial in first.named_symbols.items()}
    second_rates = dict(zip(second.symbols, second.derivatives))
    for name, rate in zip(first.symbols, first.derivatives):
        other = second_rates[renaming[name]]
        if (rate is None) != (other is None):
            found.append('{} has an equation in one and not the other'.format(name))
        elif rate is not None and sympy.simplify(
                rate.subs(renaming, simultaneous=True) - other) != 0:
            found.append('the equation for {} differs'.format(name))
    return found


class TestTheGeneralSchemeWithBothStepsIsTheParameterScheme(TestCase):

    SYSTEMS = {
        'Michaelis-Menten': ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                             'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'],
        'predator-prey': ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )',
                          'dp/dt = s*p*(1 - h*p / n)'],
        'logistic growth': ['dn/dt = r*n*(1 - n/K)'],
    }

    def test_they_are_the_same_system_up_to_names(self):
        for label, equations in self.SYSTEMS.items():
            with self.subTest(system=label):
                system = ODESystem.from_equations(equations)
                translator = ODETranslator.from_ode_system(system)
                self.assertEqual(ReductionScheme.PARAMETER, translator.scheme())
                by_parameter = translator.translate_parameter(system)
                by_general = translator.translate_general(system)
                self.assertEqual([], disagreements(by_general, by_parameter))
                self.assertEqual([], false_equations(system, by_general))
                self.assertEqual([], false_equations(system, by_parameter))

    def test_with_a_step_declined_they_differ(self):
        # So the test above is not satisfied by any two reductions whatever.
        system = ODESystem.from_equations(self.SYSTEMS['Michaelis-Menten'])
        translator = ODETranslator.from_ode_system(system)
        by_parameter = translator.translate_parameter(system)
        unnormalised = translator.translate_general(system, normalise_parameters=False)
        self.assertNotEqual([], disagreements(unnormalised, by_parameter))


class TestTheGeneralSchemeOnADependentVariableSystem(TestCase):

    SYSTEMS = {
        'Hubert & Labahn 6.4': ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)'],
        'with a parameter': ['dz1/dt = a*z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)'],
    }

    def test_it_is_the_dependent_variable_scheme_up_to_names(self):
        for label, equations in self.SYSTEMS.items():
            with self.subTest(system=label):
                system = ODESystem.from_equations(equations)
                translator = ODETranslator.from_ode_system(system)
                self.assertEqual(ReductionScheme.DEP_VAR, translator.scheme())
                by_dep_var = translator.translate_dep_var(system)
                by_general = translator.translate_general(system)
                self.assertEqual([], false_equations(system, by_dep_var))
                self.assertEqual([], false_equations(system, by_general))
                self.assertEqual([], disagreements(by_general, by_dep_var))


class TestWithBothStepsOffItIsHubertAndLabahn(TestCase):
    ''' Their printed reduced systems, in their Examples 6.4 and 6.6. '''

    def test_example_6_6(self):
        # "t dy1/dt = y1 (y1 y2 - 1)  and  t dy2/dt = y2 (1 + y2)", for the invariants
        # t z1^3 and z1^2 z2 / t.
        system = ODESystem.from_equations(
            ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
             'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])
        reduced = ODETranslator.from_ode_system(system).translate_general(
            system, normalise_parameters=False, rescale_clock=False)
        t, z1, z2 = system.symbols
        name = {monomial: symbol for symbol, monomial in reduced.named_symbols.items()}
        y1, y2 = name[t * z1 ** 3], name[z1 ** 2 * z2 / t]
        rates = dict(zip(reduced.symbols, reduced.derivatives))
        self.assertEqual(0, sympy.simplify(t * rates[y1] - y1 * (y1 * y2 - 1)))
        self.assertEqual(0, sympy.simplify(t * rates[y2] - y2 * (1 + y2)))

    def test_example_6_4(self):
        # "dy/dt = y (1 + 1/t)" for the invariant y = z1 z2.
        system = ODESystem.from_equations(
            ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)'])
        reduced = ODETranslator.from_ode_system(system).translate_dep_var(
            system, normalise_parameters=False)
        t, z1, z2 = system.symbols
        name = {monomial: symbol for symbol, monomial in reduced.named_symbols.items()}
        y = name[z1 * z2]
        rates = dict(zip(reduced.symbols, reduced.derivatives))
        self.assertEqual(0, sympy.simplify(rates[y] - y * (1 + 1 / t)))


class TestATranslatorBuiltWithoutTheClocksColumn(TestCase):
    '''
    The dependent-variable scheme's scaling leaves the independent variable alone, and a
    scaling matrix written over the dependents only -- Hubert & Labahn's `A = [1, -1]`
    for their 6.4 -- says the same thing without the clock's zero column.

    `translate_dep_var` takes such a matrix to be the full one with that zero column.
    '''

    HL_6_4 = ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)']

    def test_hubert_and_labahns_own_matrix_for_6_4(self):
        # They give A = [1, -1] on (z1, z2), and the reduced equation dy/dt = y (1 + 1/t)
        # for the invariant y = z1 z2.
        system = ODESystem.from_equations(self.HL_6_4)
        reduced = ODETranslator(sympy.Matrix([[1, -1]])).translate_dep_var(
            system, normalise_parameters=False)
        self.assertEqual([], false_equations(system, reduced))
        t, z1, z2 = system.symbols
        name = {monomial: symbol for symbol, monomial in reduced.named_symbols.items()}
        y = name[z1 * z2]
        rates = dict(zip(reduced.symbols, reduced.derivatives))
        self.assertEqual(0, sympy.simplify(rates[y] - y * (1 + 1 / t)))

    def test_the_clock_need_not_be_the_first_symbol(self):
        # The zero column goes where the independent variable is, not at the front.
        # The auxiliary may differ from the full translator's -- z1 or 1/z1 both carry
        # the scaling -- so what is required is a true reduction with n - r invariants,
        # the clock among them.
        system = ODESystem.from_equations(self.HL_6_4)
        system.reorder_variables('z1 t z2')
        self.assertEqual(1, system.independent_variable_index)
        full = ODETranslator.from_ode_system(system)
        reduced = ODETranslator(sympy.Matrix([[1, -1]])).translate_dep_var(system)
        self.assertEqual([], false_equations(system, reduced))
        monomials = list(reduced.named_symbols.values())
        invariant = [m for m in monomials if full.is_invariant_expr(m)]
        self.assertEqual(len(system.symbols) - full.r, len(invariant))
        self.assertIn(system.independent_variable, invariant)

    def test_a_multiplier_given_with_the_short_matrix_is_the_one_used(self):
        # V = [[0, 1], [-1, 1]] is unimodular with A V = [1, 0] for A = [1, -1]: auxiliary
        # 1/z2, invariant z1 z2.  The normal multiplier takes z1 as the auxiliary.
        system = ODESystem.from_equations(self.HL_6_4)
        t, z1, z2 = system.symbols
        short = sympy.Matrix([[1, -1]])
        by_default = ODETranslator(short).translate_dep_var(system)
        self.assertIn(z1, by_default.named_symbols.values())

        chosen = ODETranslator(short, hermite_multiplier=sympy.Matrix([[0, 1], [-1, 1]]))
        reduced = chosen.translate_dep_var(system)
        self.assertEqual([], false_equations(system, reduced))
        self.assertIn(1 / z2, reduced.named_symbols.values())
        self.assertNotIn(z1, reduced.named_symbols.values())

    def test_a_matrix_of_any_other_width_is_refused(self):
        system = ODESystem.from_equations(self.HL_6_4)
        for width in (1, 4):
            with self.subTest(width=width):
                with self.assertRaises(ValueError) as caught:
                    ODETranslator(sympy.ones(1, width)).translate_dep_var(system)
                self.assertIn('{} columns and the system has 3 symbols'.format(width),
                              str(caught.exception))

    def test_it_reduces_to_the_same_system(self):
        system = ODESystem.from_equations(
            ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)'])
        full = ODETranslator.from_ode_system(system)
        clock = system.independent_variable_index
        without_clock = full.scaling_matrix.copy()
        without_clock.col_del(clock)

        short = ODETranslator(without_clock)
        self.assertEqual(len(system.symbols) - 1, short.scaling_matrix.cols)
        reduced = short.translate_dep_var(system)
        self.assertEqual([], false_equations(system, reduced))
        self.assertEqual([], disagreements(reduced, full.translate_dep_var(system)))


class TestAShortScalingMatrixWithItsSymbolsNamed(TestCase):

    def test_the_independent_variable_is_added_to_the_symbols_too(self):
        system = ODESystem.from_equations(HL_6_4)
        t, z1, z2 = system.symbols
        short = ODETranslator(sympy.Matrix([[1, -1]]), symbols_domain=[z1, z2])
        full = short._with_independent_variable(system)
        self.assertEqual((t, z1, z2), tuple(full.symbols_domain))
        self.assertEqual([], false_equations(system, short.translate_dep_var(system)))


class TestTheOrderOfTheSymbolsDoesNotChangeTheTruth(TestCase):
    '''
    The parameter scheme needs the parameters first.  In another order the system is the
    same system, `translate` takes another scheme, and the reduction is still true.
    '''

    def test_michaelis_menten_with_a_dependent_first(self):
        system = ODESystem.from_equations(MICHAELIS_MENTEN)
        self.assertEqual(ReductionScheme.PARAMETER,
                         ODETranslator.from_ode_system(system).scheme())
        symbols = list(system.symbols)
        system.reorder_variables([symbols[-1]] + symbols[:-1])
        self.assertIn(system.symbols[0], system.dependents)
        translator = ODETranslator.from_ode_system(system)
        self.assertNotEqual(ReductionScheme.PARAMETER, translator.scheme())
        self.assertEqual([], false_equations(system, translator.translate(system)))
