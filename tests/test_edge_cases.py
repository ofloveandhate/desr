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


def parameters_of(system):
    '''The variables with no equation of their own.'''
    return system.constant_variables


def expected_variable_count(system, translation, include_aux_vars=True):
    '''
    How many variables the reduced system should have.

    This is scheme-dependent, and the library does not expose it yet; when it does, this
    helper should be deleted in favour of asking.

      - parameter scheme: the n - r coordinates, and nothing else
      - dependent-variable scheme: the independent variable, r auxiliaries and the
        invariants of the remaining n - 1 variables, so n in total
      - general scheme: as above but the independent variable is outside the count, n + 1,
        or n - r + 1 when the auxiliaries are dropped
    '''
    n = len(system.variables)
    if translation._is_translate_parameter_compatible(system):
        return n - translation.r
    if translation.scaling_matrix[:, system.indep_var_index].is_zero_matrix:
        return n
    return n + 1 if include_aux_vars else n - translation.r + 1


class TestDisconnectedSymbols(TestCase):
    ''' Symbols that nothing in the dynamics refers to. '''

    @staticmethod
    def with_inert(names):
        '''`dz/dt = a*z - b*z**2` with extra parameters that appear in no equation.'''
        base = ODESystem.from_equations(['dz/dt = a*z - b*z**2'])
        derivatives = [None if v in base.constant_variables else base.derivative_dict[v]
                       for v in base.variables]
        extra = [sympy.Symbol(name) for name in names]
        return ODESystem(list(base.variables) + extra,
                         derivatives + [None] * len(extra),
                         indep_var=base.indep_var)

    def test_one_inert_parameter(self):
        system = self.with_inert(['rho'])
        self.assertIn(sympy.Symbol('rho'), parameters_of(system))
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

    @unittest.expectedFailure
    def test_a_system_with_no_invariants(self):
        # `dz/dt = z` reduces completely: r = 1, the only invariant is t, and z is the
        # auxiliary.  A total reduction is a perfectly good answer, but translate_dep_var
        # multiplies empty matrices and raises
        #   ShapeError: Matrix shapes must agree (0, 0) != (1, 0)
        system = ODESystem.from_equations(['dz/dt = z'])
        translation = ODETranslation.from_ode_system(system)
        self.assertEqual(translation.r, 1)
        reduced = translation.translate(system)
        self.assertEqual(len(reduced.variables),
                         expected_variable_count(system, translation))

    def test_single_dependent_variable_with_one_parameter(self):
        system = ODESystem.from_equations(['dz/dt = a*z'])
        translation = ODETranslation.from_ode_system(system)
        reduced = translation.translate(system)
        self.assertEqual(len(reduced.variables),
                         expected_variable_count(system, translation))

    def test_a_conserved_variable_keeps_its_equation(self):
        base = ODESystem.from_equations(['dz/dt = a*z'])
        derivatives = [None if v in base.constant_variables else base.derivative_dict[v]
                       for v in base.variables]
        conserved = ODESystem(list(base.variables) + [sympy.Symbol('w')],
                              derivatives + [sympy.sympify(0)],
                              indep_var=base.indep_var)
        w = sympy.Symbol('w')
        # It has an equation, so it is not a parameter ...
        self.assertNotIn(w, parameters_of(conserved))
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

    @unittest.expectedFailure
    def test_a_plain_sir_model(self):
        # beta and gamma are sympy functions, so this raises
        #   TypeError: bad operand type for unary -: 'FunctionClass'
        system = ODESystem.from_equations(['dS/dt = -beta*S*I',
                                           'dI/dt = beta*S*I - gamma*I',
                                           'dR/dt = gamma*I'])
        for name in ('S', 'I', 'R', 'beta', 'gamma'):
            self.assertIn(sympy.Symbol(name), system.variables)

    @unittest.expectedFailure
    def test_a_parameter_called_pi_is_not_swallowed(self):
        # `pi` sympifies to the number, so the parameter vanishes without a word and the
        # system reduced is not the system written down.
        system = ODESystem.from_equations(['dx/dt = pi*x'])
        self.assertIn(sympy.Symbol('pi'), system.variables)

    @unittest.expectedFailure
    def test_a_name_resolves_the_same_way_whichever_door_it_comes_in(self):
        # `from_equations` reads 'S' as a symbol, via _clash1.  desr.numerics sympifies
        # user-supplied names without it, so the same string is sympy's singleton registry
        # there.  One convention, or the same name means two things.
        from desr.numerics import NumericTranslation
        system = ODESystem.from_equations(['dS/dt = -k*S'])
        system.reorder_variables(['t', 'S', 'k'])
        self.assertIn(sympy.Symbol('S'), system.variables)
        translation = ODETranslation.from_ode_system(system)
        numeric = NumericTranslation(system, translation)
        point = numeric.forward({sympy.Symbol('t'): 1.0, sympy.Symbol('S'): 2.0,
                                 sympy.Symbol('k'): 0.5})
        # Naming the same variable by its string should work here too.
        numeric.reverse(point, known_values={'S': 2.0})
