'''
The three reduction schemes agree where more than one applies.

The general scheme always applies.  Its two further steps -- normalising the constant
auxiliaries and identifying a clock invariant with the clock -- are what the other
schemes do by construction, so on a system one of those fits, the general scheme with
both steps must give the same reduced system up to the names of its symbols.  With both
steps off it is Hubert & Labahn's reduction, checked here against their printed examples.

"The same up to names" is decided through the monomials: two reduced symbols are the
same quantity when they denote the same monomial in the original symbols.
'''
from unittest import TestCase, expectedFailure

import sympy

from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator, ReductionScheme

from oracles import false_equations


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
            system, normalise_parameters=False, identify_indep_var=False)
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

    `translate_dep_var` has a branch for exactly this, and it does not work: it raises
    an IndexError from `_dep_var_translator`.  Nothing had ever run it.  Marked as an
    expected failure, in the convention of `test_edge_cases.py`, so that it turns into
    a failure when the branch is fixed or removed.
    '''

    @expectedFailure
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
