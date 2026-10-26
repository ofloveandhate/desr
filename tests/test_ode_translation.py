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
        extended.reorder_variables(['t', 'z1', 'z2', 't0', 'z10', 'z20'])
        translation = ODETranslator(sympy.Matrix([[3, -1, 5, 3, -1, 5]]),
                                     variables_domain=extended.variables)
        reduced = translation.translate(extended)

        self.assertEqual(len(reduced.variables),
                         len(extended.variables) - translation.r)
        # Every symbol the substitutions mention must exist in the reduced system.
        mentioned = set()
        for value in translation.translate_parameter_substitutions(extended).values():
            mentioned |= value.free_symbols
        self.assertFalse(mentioned - set(reduced.variables))

    def test_rational_section_returns_its_result(self):
        equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)']
        system = ODESystem.from_equations(equations)
        translation = ODETranslator.from_ode_system(system)
        section = translation.rational_section()
        self.assertIsNotNone(section)
        self.assertEqual(section.shape, (1, translation.r))

    def test_reverse_translate_dispatches_by_count(self):
        translation = ODETranslator(scaling_matrix=sympy.Matrix([[0, 1, -1]]),
                                     hermite_multiplier=sympy.Matrix([[0, 1, 0],
                                                                      [1, 0, 1],
                                                                      [0, 0, 1]]))
        t, x0, y0 = sympy.symbols('t x0 y0')
        self.assertEqual(translation.reverse_translate((x0, y0), 0), (x0, y0 / x0))
        self.assertEqual(translation.reverse_translate((x0, t, y0), 0), (t, x0, y0 / x0))
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
