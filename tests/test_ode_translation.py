from unittest import TestCase

import sympy

from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator


class TestODETranslatorTraps(TestCase):
    ''' Things that used to fail silently, or somewhere other than where the mistake was. '''

    def test_non_unimodular_multiplier_is_refused_at_construction(self):
        scaling = sympy.Matrix([[1, 0, 1, 1, -1], [0, 1, 0, -1, 1]])
        doubled = sympy.Matrix([[0, 0, 2, 0, 0], [0, 0, 0, 1, 0], [1, 1, -2, -1, 0],
                                [0, 0, 0, 0, 1], [0, 1, 0, -1, 1]])
        with self.assertRaises(ValueError) as caught:
            ODETranslator(scaling, hermite_multiplier=doubled)
        self.assertIn('not unimodular', str(caught.exception))
        self.assertIn('roots of degree 2', str(caught.exception))

    def test_a_conserved_variable_does_not_depend_on_its_name(self):
        # A catalyst gets a derivative of identically zero straight from the stoichiometry.
        # That expression has no monomials, so it used to contribute a 0x0 block to the
        # exponent matrix; sympy's hstack absorbs such a block on the left but raises on
        # the right, so the same network reduced or crashed according to where the catalyst
        # fell in the variable ordering -- that is, according to what it was called.
        scalings = {}
        for catalyst in ('e', 'z'):
            network = ChemicalReactionNetwork.from_diagram(
                'x + {c} -> y + {c}'.format(c=catalyst))
            system = network.to_ode_system()
            self.assertEqual(system.derivative_dict[sympy.Symbol(catalyst)],
                             sympy.sympify(0))
            scalings[catalyst] = (list(system.variables), system.maximal_scaling_matrix())

        # The two differ only by the name of the catalyst, so permuting the columns of one
        # into the other's variable order must give the same scaling action.
        first_variables, first = scalings['e']
        other_variables, other = scalings['z']
        renamed = [sympy.Symbol('z') if v == sympy.Symbol('e') else v
                   for v in first_variables]
        order = [other_variables.index(v) for v in renamed]
        permuted = other[:, order]
        self.assertEqual(first.rref()[0], permuted.rref()[0])

    def test_invariant_constant_that_cancels_is_still_a_variable(self):
        # Hubert & Labahn 6.6, with initial-condition constants adjoined so that each z_i0
        # scales like its z_i.  kappa0 = t0*z10**3 is the reduced initial time: it appears
        # in no right-hand side, but it is one of the n - r coordinates of the reduced
        # system, and translate_parameter_substitutions names it.  Building the reduced
        # system from the symbols surviving into the equations used to lose it.
        eq1 = 'dz1/dt = z1/t * (-2/3 + 1/3 *z1**5 *z2 )'
        eq2 = 'dz2/dt = z2/t * ( 10/3 - 2/3 *z1**5 *z2 + z1**2*z2/t )'
        system = ODESystem.from_equations([eq1, eq2])
        t0, z10, z20 = sympy.symbols('t0 z10 z20')
        extended = ODESystem(list(system.variables) + [t0, z10, z20],
                             list(system.derivatives) + [None, None, None],
                             indep_var=system.indep_var)
        # Parameters, then the independent variable, then the dependents.  The scaling
        # matrix below is unchanged by the move because its pattern repeats with period 3.
        extended.reorder_variables(['t0', 'z10', 'z20', 't', 'z1', 'z2'])
        translator = ODETranslator(sympy.Matrix([[3, -1, 5, 3, -1, 5]]),
                                     variables_domain=extended.variables)
        reduced = translator.translate(extended)

        self.assertEqual(len(reduced.variables),
                         len(extended.variables) - translator.r)
        # Every symbol the substitutions mention must exist in the reduced system.
        mentioned = set()
        for value in translator.translate_parameter_substitutions(extended).values():
            mentioned |= value.free_symbols
        self.assertFalse(mentioned - set(reduced.variables))

    def test_rational_section_returns_its_result(self):
        equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)']
        system = ODESystem.from_equations(equations)
        translator = ODETranslator.from_ode_system(system)
        section = translator.rational_section()
        self.assertIsNotNone(section)
        self.assertEqual(section.shape, (1, translator.r))

    def test_reverse_translate_dispatches_by_count(self):
        translator = ODETranslator(scaling_matrix=sympy.Matrix([[0, 1, -1]]),
                                     hermite_multiplier=sympy.Matrix([[0, 1, 0],
                                                                      [1, 0, 1],
                                                                      [0, 0, 1]]))
        t, x0, y0 = sympy.symbols('t x0 y0')
        self.assertEqual(translator.reverse_translate((x0, y0), 0), (x0, y0 / x0))
        self.assertEqual(translator.reverse_translate((x0, t, y0), 0), (t, x0, y0 / x0))
        # The general branch checks its inputs are genuine solutions, so use the exact one.
        general = ODETranslator(scaling_matrix=sympy.Matrix([[3, -1, 5]]),
                                 hermite_multiplier=sympy.Matrix([[1, 1, -1],
                                                                  [2, 3, 2],
                                                                  [0, 0, 1]]))
        exact = (t,
                 sympy.sympify('c3/(t**(1/3)*(ln(t-c1)-ln(t)+c2)**(2/3))'),
                 sympy.sympify('c1/(t*(ln(t-c1)-ln(t)+c2))'),
                 sympy.sympify('t/(c1 - t)'))
        self.assertEqual(general.reverse_translate(exact),
                         general.reverse_translate_general(exact))


class TestNamingScheme(TestCase):
    '''
    The naming scheme is a dict, and every reduction scheme honours it.

    Until 2.0 it was a positional triple, and the general scheme accepted it and threw it
    away -- a scheme passed to `translate_general` changed nothing at all, not even the
    independent variable's name.  Both halves are pinned here.
    '''

    def setUp(self):
        # Hubert & Labahn 6.6: the general scheme, with one auxiliary and two invariants.
        self.general_system = ODESystem.from_equations(
            ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
             'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])
        self.general_system.reorder_variables(['t', 'z1', 'z2'])
        # Michaelis-Menten: the parameter scheme.
        self.parameter_system = ODESystem.from_equations(
            ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
             'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
        self.parameter_system.reorder_variables(['e_0', 'k_1', 'k_2', 'k_m1', 't', 's', 'c'])

    def names(self, system):
        return tuple(map(str, system.variables))

    def test_a_tuple_is_refused_and_says_what_to_write(self):
        with self.assertRaises(ValueError) as caught:
            ODETranslator.from_ode_system(self.parameter_system,
                                          naming_scheme=('tau', 'nu', 'kappa'))
        message = str(caught.exception)
        self.assertIn('dict', message)
        for key in ('parameters', 'independent', 'dependents'):
            self.assertIn(key, message)

    def test_an_unknown_key_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            ODETranslator.from_ode_system(self.parameter_system,
                                          naming_scheme={'indep': 'tau'})
        self.assertIn('indep', str(caught.exception))

    def test_a_list_where_a_pattern_is_needed_is_refused(self):
        for key in ('parameters', 'auxiliaries'):
            with self.assertRaises(ValueError):
                ODETranslator.from_ode_system(self.parameter_system,
                                              naming_scheme={key: ['a', 'b']})

    def test_parameter_scheme_defaults_are_unchanged(self):
        translator = ODETranslator.from_ode_system(self.parameter_system)
        reduced = translator.translate(self.parameter_system)
        self.assertEqual(('kappa0', 'kappa1', 'tau', 'nu0', 'nu1'), self.names(reduced))

    def test_general_scheme_defaults_follow_hubert_labahn(self):
        translator = ODETranslator.from_ode_system(self.general_system)
        reduced = translator.translate_general(self.general_system)
        self.assertEqual(('t', 'x0', 'y0', 'y1'), self.names(reduced))

    def test_general_scheme_honours_the_scheme(self):
        ''' The regression: this used to be discarded in full. '''
        translator = ODETranslator.from_ode_system(
            self.general_system,
            naming_scheme={'independent': 'T', 'dependents': ['alpha', 'beta'],
                           'auxiliaries': 'aux'})
        reduced = translator.translate_general(self.general_system)
        self.assertEqual(('T', 'aux0', 'alpha', 'beta'), self.names(reduced))

    def test_general_scheme_renames_only_the_auxiliaries_if_that_is_all_that_is_set(self):
        translator = ODETranslator.from_ode_system(self.general_system,
                                                   naming_scheme={'auxiliaries': 'w'})
        reduced = translator.translate_general(self.general_system)
        self.assertEqual(('t', 'w0', 'y0', 'y1'), self.names(reduced))

    def test_parameter_scheme_takes_a_partial_scheme(self):
        translator = ODETranslator.from_ode_system(self.parameter_system,
                                                   naming_scheme={'parameters': 'c'})
        reduced = translator.translate(self.parameter_system)
        self.assertEqual(('c0', 'c1', 'tau', 'nu0', 'nu1'), self.names(reduced))

    def test_one_scheme_serves_both_reductions(self):
        '''
        Each reduction falls back to its own defaults for what the scheme leaves unset, so
        the general one keeps y for its invariants while both take T and aux.
        '''
        shared = {'parameters': 'p', 'independent': 'T', 'auxiliaries': 'aux'}
        reduced = ODETranslator.from_ode_system(
            self.parameter_system, naming_scheme=shared).translate(self.parameter_system)
        self.assertEqual(('p0', 'p1', 'T', 'nu0', 'nu1'), self.names(reduced))
        reduced = ODETranslator.from_ode_system(
            self.general_system, naming_scheme=shared).translate_general(self.general_system)
        self.assertEqual(('T', 'aux0', 'y0', 'y1'), self.names(reduced))

    def test_a_wrong_length_dependents_list_is_refused_with_both_counts(self):
        translator = ODETranslator.from_ode_system(
            self.general_system, naming_scheme={'dependents': ['only_one']})
        with self.assertRaises(ValueError) as caught:
            translator.translate_general(self.general_system)
        message = str(caught.exception)
        self.assertIn('1 names', message)
        self.assertIn('2 of', message)

    def test_the_scheme_as_set_is_readable(self):
        translator = ODETranslator.from_ode_system(self.general_system,
                                                   naming_scheme={'auxiliaries': 'w'})
        self.assertEqual({'auxiliaries': 'w'}, translator.naming_scheme)

    def test_setting_no_scheme_leaves_it_empty(self):
        translator = ODETranslator.from_ode_system(self.general_system)
        self.assertEqual({}, translator.naming_scheme)
