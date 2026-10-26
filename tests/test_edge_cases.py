'''
Adversarial systems.

Nothing here is a realistic model.  These are the shapes that fall between the examples --
systems with symbols connected to nothing, systems with no parameters at all, names that
sympy reads as something other than a symbol -- and the invariants that should hold no
matter which shape a system has.

Tests marked `expectedFailure` record bugs that are known and not yet fixed.  They turn
into failures the moment the bug is, which is the point.
'''
import unittest
from unittest import TestCase

import sympy

from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslation


def expected_variable_count(system, translation, scheme=None, include_aux_vars=True):
    '''
    How many variables the reduced system should have.

    This is scheme-dependent, and the library does not expose it yet; when it does, this
    helper should be deleted in favour of asking.

      - parameter scheme: the n - r coordinates, and nothing else
      - dependent-variable scheme: the independent variable, r auxiliaries and the
        invariants of the remaining n - 1 variables, so n in total
      - general scheme: as above but the independent variable is outside the count, n + 1,
        or n - r + 1 when the auxiliaries are dropped

    `scheme` says which one, for a caller who picked it rather than letting `translate`
    choose.  Inferring it is only right when `translate` did the choosing -- calling
    `translate_general` on a system that `translate` would have given to another scheme is
    perfectly legal and gives a different count.
    '''
    n = len(system.variables)
    if scheme is None:
        if translation._is_translate_parameter_compatible(system):
            scheme = 'parameter'
        elif translation.scaling_matrix[:, system.indep_var_index].is_zero_matrix:
            scheme = 'dep_var'
        else:
            scheme = 'general'
    if scheme == 'parameter':
        return n - translation.r
    if scheme == 'dep_var':
        return n
    return n + 1 if include_aux_vars else n - translation.r + 1


class TestDisconnectedSymbols(TestCase):
    ''' Symbols that nothing in the dynamics refers to. '''

    @staticmethod
    def with_inert(names):
        '''`dz/dt = a*z - b*z**2` with extra parameters that appear in no equation.'''
        base = ODESystem.from_equations(['dz/dt = a*z - b*z**2'])
        extra = [sympy.Symbol(name) for name in names]
        return ODESystem(list(base.variables) + extra,
                         list(base.derivatives) + [None] * len(extra),
                         indep_var=base.indep_var)

    def test_one_inert_parameter(self):
        system = self.with_inert(['rho'])
        self.assertIn(sympy.Symbol('rho'), system.parameters)
        scaling = system.maximal_scaling_matrix()
        self.assertEqual(scaling.cols, len(system.variables))

    def test_many_inert_parameters(self):
        names = ['rho_{}'.format(i) for i in range(12)]
        system = self.with_inert(names)
        scaling = system.maximal_scaling_matrix()
        self.assertEqual(scaling.cols, len(system.variables))
        # Each one is free to scale on its own, so each adds a dimension to the action.
        plain = ODESystem.from_equations(['dz/dt = a*z - b*z**2'])
        self.assertEqual(scaling.rows,
                         plain.maximal_scaling_matrix().rows + len(names))

    def test_an_inert_parameter_does_not_change_the_rest_of_the_reduction(self):
        plain = ODESystem.from_equations(['dz/dt = a*z - b*z**2'])
        padded = self.with_inert(['rho'])
        plain_reduced = ODETranslation.from_ode_system(plain).translate(plain)
        padded_reduced = ODETranslation.from_ode_system(padded).translate(padded)
        # The extra symbol is free, so it is normalised away and the reduction is unchanged.
        self.assertEqual(len(padded_reduced.variables), len(plain_reduced.variables))


class TestDegenerateSystems(TestCase):
    ''' Systems at the edges of what counts as a system. '''

    def test_a_system_with_no_invariants(self):
        # `dz/dt = z` reduces completely: r = 1, the only invariant is t, and z is the
        # auxiliary.  A total reduction is a perfectly good answer; it used to raise
        #   ShapeError: Matrix shapes must agree (0, 0) != (1, 0)
        # from multiplying empty matrices.
        system = ODESystem.from_equations(['dz/dt = z'])
        translation = ODETranslation.from_ode_system(system)
        self.assertEqual(translation.r, 1)
        reduced = translation.translate(system)
        self.assertEqual(len(reduced.variables),
                         expected_variable_count(system, translation))
        # Nothing is left but the independent variable and the one auxiliary.
        self.assertEqual([str(v) for v in reduced.variables], ['t', 'x0'])

    def test_a_system_with_no_invariants_under_the_general_scheme(self):
        system = ODESystem.from_equations(['dz/dt = z'])
        translation = ODETranslation.from_ode_system(system)
        reduced = translation.translate_general(system)
        self.assertEqual(len(reduced.variables),
                         expected_variable_count(system, translation, scheme='general'))
        dropped = translation.translate_general(system, include_aux_vars=False)
        self.assertEqual(len(dropped.variables),
                         expected_variable_count(system, translation, scheme='general',
                                                 include_aux_vars=False))

    def test_single_dependent_variable_with_one_parameter(self):
        system = ODESystem.from_equations(['dz/dt = a*z'])
        translation = ODETranslation.from_ode_system(system)
        reduced = translation.translate(system)
        self.assertEqual(len(reduced.variables),
                         expected_variable_count(system, translation))

    def test_a_conserved_variable_keeps_its_equation(self):
        base = ODESystem.from_equations(['dz/dt = a*z'])
        conserved = ODESystem(list(base.variables) + [sympy.Symbol('w')],
                              list(base.derivatives) + [sympy.sympify(0)],
                              indep_var=base.indep_var)
        w = sympy.Symbol('w')
        # It has an equation, so it is not a parameter ...
        self.assertNotIn(w, conserved.parameters)
        # ... and that equation says it does not move.
        self.assertEqual(conserved.derivative_dict[w], sympy.sympify(0))
        # Whatever else happens, computing the scaling must not raise.
        conserved.maximal_scaling_matrix()


class TestInvariantsOfReduction(TestCase):
    '''
    Properties that should hold for every system, checked across a spread of them rather
    than asserted of one printed answer.
    '''

    SYSTEMS = {
        'verhulst': ['dn/dt = r*n*(1-n/k)'],
        'pred_prey': ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )',
                      'dp/dt = s*p*(1 - h*p / n)'],
        'one_parameter': ['dz/dt = a*z'],
        'cubic': ['dz/dt = a*z - b*z**3'],
        'two_coupled': ['dx/dt = a*x - b*x*y', 'dy/dt = c*x*y - d*y'],
    }

    def test_reduced_system_has_the_expected_number_of_variables(self):
        for name, equations in self.SYSTEMS.items():
            with self.subTest(system=name):
                system = ODESystem.from_equations(equations)
                translation = ODETranslation.from_ode_system(system)
                reduced = translation.translate(system)
                self.assertEqual(len(reduced.variables),
                                 expected_variable_count(system, translation))

    def test_every_invariant_is_actually_invariant(self):
        for name, equations in self.SYSTEMS.items():
            with self.subTest(system=name):
                system = ODESystem.from_equations(equations)
                translation = ODETranslation.from_ode_system(system)
                # A*V_n = 0 is what makes the columns of V_n invariants.
                product = translation.scaling_matrix * translation.herm_mult_n
                self.assertTrue(product.is_zero_matrix, name)

    def test_substitutions_only_mention_symbols_of_the_reduced_system(self):
        for name, equations in self.SYSTEMS.items():
            with self.subTest(system=name):
                system = ODESystem.from_equations(equations)
                translation = ODETranslation.from_ode_system(system)
                if not translation._is_translate_parameter_compatible(system):
                    continue
                reduced = translation.translate(system)
                mentioned = set()
                for value in translation.translate_parameter_substitutions(system).values():
                    mentioned |= value.free_symbols
                self.assertFalse(mentioned - set(reduced.variables), name)

    def test_renaming_a_variable_does_not_change_the_size_of_the_reduction(self):
        system = ODESystem.from_equations(['dn/dt = r*n*(1-n/k)'])
        renamed = ODESystem.from_equations(['dpopulation/dt = growth*population*'
                                            '(1-population/capacity)'])
        self.assertEqual(ODETranslation.from_ode_system(system).r,
                         ODETranslation.from_ode_system(renamed).r)
        self.assertEqual(len(system.variables), len(renamed.variables))

    def test_reordering_variables_does_not_change_the_number_of_symmetries(self):
        equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)']
        first = ODESystem.from_equations(equations)
        original_r = ODETranslation.from_ode_system(first).r
        for order in ('tnpKdhkrs', 'rhKskdtnp', 'srkhdKpnt'):
            with self.subTest(order=order):
                system = ODESystem.from_equations(equations)
                system.reorder_variables(order)
                self.assertEqual(ODETranslation.from_ode_system(system).r, original_r)


class TestAdversarialNames(TestCase):
    '''
    sympy reads several names as something other than a symbol.  desr parses with
    `_clash1`, which covers the single-letter cases `E I N O Q S`, but not sympy's
    `_clash2` -- `beta gamma pi zeta` -- so those are read as sympy's own objects.
    '''

    def test_single_letter_clashes_are_handled(self):
        for name in ('E', 'I', 'S', 'N', 'O', 'Q'):
            with self.subTest(name=name):
                system = ODESystem.from_equations(['dx/dt = {}*x'.format(name)])
                self.assertIn(sympy.Symbol(name), system.variables)

    def test_a_parameter_called_lamda(self):
        system = ODESystem.from_equations(['dz/dt = lamda*z'])
        self.assertIn(sympy.Symbol('lamda'), system.variables)

    def test_a_plain_sir_model(self):
        system = ODESystem.from_equations(['dS/dt = -beta*S*I',
                                           'dI/dt = beta*S*I - gamma*I',
                                           'dR/dt = gamma*I'])
        for name in ('S', 'I', 'R', 'beta', 'gamma'):
            self.assertIn(sympy.Symbol(name), system.variables)

    def test_a_parameter_called_pi_is_not_swallowed(self):
        # Left to sympy `pi` is the number, and the parameter would vanish without a word,
        # so the system reduced would not be the system written down.
        system = ODESystem.from_equations(['dx/dt = pi*x'])
        self.assertIn(sympy.Symbol('pi'), system.variables)

    def test_pi_can_still_be_asked_for_as_a_number(self):
        from desr.sympy_helper import sympify_expression
        self.assertIn(sympy.Symbol('pi'), sympify_expression('pi*x').free_symbols)
        self.assertNotIn(sympy.Symbol('pi'),
                         sympify_expression('pi*x', numeric=['pi']).free_symbols)

    def test_a_name_resolves_the_same_way_whichever_door_it_comes_in(self):
        # Both `from_equations` and desr.numerics go through sympy_helper, so 'S' is the
        # same symbol in each.  It used to be sympy's singleton registry in numerics, so
        # naming a variable by its own string did not find it.
        from desr.numerics import NumericTranslation
        system = ODESystem.from_equations(['dS/dt = -k*S'])
        system.reorder_variables(['t', 'S', 'k'])
        self.assertIn(sympy.Symbol('S'), system.variables)
        numeric = NumericTranslation(system, ODETranslation.from_ode_system(system))

        by_symbol = numeric.forward({sympy.Symbol('t'): 1.0, sympy.Symbol('S'): 2.0,
                                     sympy.Symbol('k'): 0.5})
        by_name = numeric.forward({'t': 1.0, 'S': 2.0, 'k': 0.5})
        self.assertEqual(by_symbol, by_name)
