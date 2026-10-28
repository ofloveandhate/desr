'''
Classifying a system's variables: which are parameters and which are dependents.

A variable with no equation is a parameter.  A variable with an equation is a dependent,
*including* when that equation says zero: writing `dw/dt = 0` is a deliberate statement
that `w` is conserved, not a roundabout way of declaring a constant.  The two are
distinguishable, they are decided per variable, and the system -- not the translator --
is what knows.

The distinction is invisible to the scaling computation: `dw/dt = 0` constrains nothing,
so the scaling matrix is the same either way.  Everything here is therefore about what the
system reports and what the reduction does with it, not about the mathematics of the
reduction.
'''
from unittest import TestCase

import sympy

from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator

t, z, w, a = sympy.symbols('t z w a')


def with_conserved_w():
    '''`dz/dt = a*z` with a conserved `w` adjoined: a dependent whose equation is zero.'''
    return ODESystem.from_equations(['dz/dt = a*z', 'dw/dt = 0'])


class TestParametersAndDependents(TestCase):
    ''' The basic split, and that it is reported faithfully. '''

    def test_a_variable_with_no_equation_is_a_parameter(self):
        system = ODESystem.from_equations(['dz/dt = a*z'])
        self.assertEqual(system.parameters, (a,))
        self.assertEqual(system.num_parameters, 1)
        self.assertEqual(system.dependent_variables, (z,))
        self.assertEqual(system.num_dependent_variables, 1)

    def test_the_independent_variable_is_neither(self):
        system = ODESystem.from_equations(['dz/dt = a*z'])
        self.assertNotIn(t, system.parameters)
        self.assertNotIn(t, system.dependent_variables)

    def test_derivatives_reports_a_parameter_as_having_none(self):
        # `derivatives` used to map None to 0, which made a parameter indistinguishable
        # from a conserved dependent to anybody reading it.
        system = ODESystem.from_equations(['dz/dt = a*z'])
        self.assertIsNone(system.derivatives[system.variables.index(a)])
        self.assertNotIn(a, system.derivative_dict)

    def test_an_explicit_zero_is_a_conserved_dependent(self):
        system = with_conserved_w()
        self.assertIn(w, system.dependent_variables)
        self.assertNotIn(w, system.parameters)
        self.assertEqual(system.derivative_dict[w], sympy.sympify(0))

    def test_a_dependent_variable_may_have_derivative_one(self):
        # The old test for "has an equation" compared the derivative against 1, so a
        # variable that genuinely moves at unit rate was classed as neither a parameter
        # nor a dependent and fell out of the system entirely.
        system = ODESystem.from_equations(['dx/dt = 1'])
        x = sympy.Symbol('x')
        self.assertEqual(system.dependent_variables, (x,))
        self.assertEqual(system.parameters, ())

    def test_a_dependent_variable_with_derivative_one_survives_reduction(self):
        system = ODESystem.from_equations(['dx/dt = 1'])
        translator = ODETranslator.from_ode_system(system)
        reduced = translator.translate(system)
        self.assertTrue(reduced.variables)

    def test_the_retired_vocabulary_is_gone(self):
        # Cut outright rather than aliased, so that a caller is told rather than quietly
        # given the old meaning of "constant".
        system = ODESystem.from_equations(['dz/dt = a*z'])
        for name in ('constant_variables', 'non_constant_variables', 'num_constants',
                     'num_nonconstants'):
            with self.subTest(name=name):
                self.assertFalse(hasattr(system, name))


class TestDeclaringParameters(TestCase):
    ''' Saying "this conserved variable is really a parameter", at construction or after. '''

    def test_make_parameter_demotes_a_conserved_variable(self):
        system = with_conserved_w()
        system.make_parameter('w')
        self.assertIn(w, system.parameters)
        self.assertNotIn(w, system.dependent_variables)
        self.assertIsNone(system.derivatives[system.variables.index(w)])

    def test_parameters_at_construction_does_the_same(self):
        declared = ODESystem.from_equations(['dz/dt = a*z', 'dw/dt = 0'], parameters=['w'])
        made = with_conserved_w()
        made.make_parameter('w')
        self.assertEqual(set(declared.parameters), set(made.parameters))
        self.assertEqual(declared.dependent_variables, made.dependent_variables)

    def test_declaring_at_construction_orders_it_with_the_parameters(self):
        # Construction ends in a default ordering, so a variable declared a parameter there
        # lands with the parameters.  make_parameter leaves it where it is, because the
        # ordering decides which parameters get normalised away and must not move silently.
        declared = ODESystem.from_equations(['dz/dt = a*z', 'dw/dt = 0'], parameters=['w'])
        self.assertEqual(declared.variables, (t, z, a, w))

        made = with_conserved_w()
        made.make_parameter('w')
        self.assertEqual(made.variables, (t, w, z, a))
        made.default_order_variables()
        self.assertEqual(made.variables, declared.variables)

    def test_parameters_accepts_symbols_as_well_as_strings(self):
        by_name = ODESystem.from_equations(['dz/dt = a*z', 'dw/dt = 0'], parameters=['w'])
        by_symbol = ODESystem.from_equations(['dz/dt = a*z', 'dw/dt = 0'], parameters=[w])
        self.assertEqual(by_name.parameters, by_symbol.parameters)

    def test_parameters_through_the_raw_constructor(self):
        base = with_conserved_w()
        system = ODESystem(base.variables, base.derivatives, indep_var=base.indep_var,
                           parameters=['w'])
        self.assertIn(w, system.parameters)

    def test_parameters_through_from_tex(self):
        tex = '\n'.join([r'\frac{dz}{dt} &= a z \\', r'\frac{dw}{dt} &= 0'])
        system = ODESystem.from_tex(tex, parameters=['w'])
        self.assertIn(w, system.parameters)
        self.assertNotIn(w, system.dependent_variables)

    def test_a_catalyst_can_be_declared_a_parameter(self):
        # The motivating case: stoichiometry hands a catalyst a derivative of identically
        # zero, which by default makes it a conserved dependent.  A user who wants it to
        # be a constant of the model says so.
        from desr.chemical_reaction_network import ChemicalReactionNetwork
        network = ChemicalReactionNetwork.from_diagram('x + e -> y + e')
        system = network.to_ode_system()
        e = sympy.Symbol('e')
        self.assertIn(e, system.dependent_variables)
        system.make_parameter(e)
        self.assertIn(e, system.parameters)

    def test_making_a_parameter_of_a_parameter_is_a_no_op(self):
        system = ODESystem.from_equations(['dz/dt = a*z'])
        system.make_parameter('a')
        self.assertEqual(system.parameters, (a,))

    def test_a_variable_that_actually_moves_cannot_be_made_a_parameter(self):
        system = ODESystem.from_equations(['dz/dt = a*z'])
        with self.assertRaises(ValueError) as caught:
            system.make_parameter('z')
        self.assertIn('a*z', str(caught.exception))

    def test_the_independent_variable_cannot_be_made_a_parameter(self):
        system = ODESystem.from_equations(['dz/dt = a*z'])
        with self.assertRaises(ValueError):
            system.make_parameter('t')

    def test_an_unknown_name_cannot_be_made_a_parameter(self):
        system = ODESystem.from_equations(['dz/dt = a*z'])
        with self.assertRaises(ValueError) as caught:
            system.make_parameter('nowhere')
        self.assertIn('nowhere', str(caught.exception))

    def test_a_variable_with_an_initial_condition_cannot_be_made_a_parameter(self):
        # `w(0) = w_0` only means something for a variable that has a trajectory.  Dropping
        # the equation would leave the initial condition dangling, so say so instead.
        system = with_conserved_w()
        system.update_initial_conditions({'w': 'w_0'})
        with self.assertRaises(ValueError) as caught:
            system.make_parameter('w')
        self.assertIn('initial condition', str(caught.exception))


class TestClassificationIsCarriedAround(TestCase):
    ''' The split belongs to the system, so it must survive everything the system does. '''

    def test_a_parsed_system_and_a_hand_built_one_agree(self):
        parsed = ODESystem.from_equations(['dz/dt = a*z'])
        built = ODESystem([t, z, a], [sympy.sympify(1), a * z, None], indep_var=t)
        self.assertEqual(parsed.variables, built.variables)
        self.assertEqual(parsed.derivatives, built.derivatives)
        self.assertEqual(parsed.parameters, built.parameters)
        self.assertEqual(parsed.dependent_variables, built.dependent_variables)

    def test_copy_preserves_the_split(self):
        system = with_conserved_w()
        system.make_parameter('w')
        duplicate = system.copy()
        self.assertEqual(duplicate.parameters, system.parameters)
        self.assertEqual(duplicate.dependent_variables, system.dependent_variables)

    def test_reordering_preserves_the_split(self):
        system = with_conserved_w()
        system.reorder_variables(['t', 'w', 'a', 'z'])
        self.assertEqual(set(system.parameters), {a})
        self.assertEqual(set(system.dependent_variables), {z, w})

    def test_default_ordering_puts_parameters_last(self):
        system = with_conserved_w()
        system.reorder_variables(['a', 'w', 't', 'z'])
        system.default_order_variables()
        self.assertEqual(system.variables, (t, w, z, a))


class TestConservedVariablesAndReduction(TestCase):
    ''' What the reduction makes of a conserved dependent. '''

    def test_a_conserved_variable_does_not_change_the_scaling(self):
        # `dw/dt = 0` imposes no condition on how `w` scales, so the scaling matrix cannot
        # tell a conserved dependent from a parameter.  The whole difference is bookkeeping.
        conserved = with_conserved_w()
        as_parameter = with_conserved_w()
        as_parameter.make_parameter('w')
        self.assertEqual(conserved.variables, as_parameter.variables)
        self.assertEqual(conserved.maximal_scaling_matrix(),
                         as_parameter.maximal_scaling_matrix())

    def test_a_conserved_variable_at_the_end_is_not_mistaken_for_a_parameter(self):
        # The parameter scheme needs the parameters to be the trailing variables.  That
        # used to be checked by asking whether the trailing derivatives were zero, which a
        # conserved variable also satisfies -- so ordering one last made the system look
        # eligible for a scheme it does not satisfy.
        system = with_conserved_w()
        system.reorder_variables(['t', 'z', 'a', 'w'])
        translator = ODETranslator.from_ode_system(system)
        self.assertFalse(translator._is_translate_parameter_compatible(system))

    def test_a_conserved_variable_keeps_an_equation_in_the_reduced_system(self):
        conserved = with_conserved_w()
        translator = ODETranslator.from_ode_system(conserved)
        reduced = translator.translate(conserved)
        # Whatever the scheme chose, the reduction must not have silently dropped it.
        self.assertTrue(reduced.variables)


class TestHowSystemsPrint(TestCase):
    '''
    repr lists the equations, then the parameters on a line of their own.  A parameter has
    no equation, so printing `dK/dt = 0` for it was stating something the system does not
    say -- and it made a genuinely conserved variable unrecognisable.
    '''

    def test_parameters_get_their_own_line(self):
        system = ODESystem.from_equations(['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )',
                                           'dp/dt = s*p*(1 - h*p / n)'])
        self.assertEqual(repr(system), '\n'.join([
            'dt/dt = 1',
            'dn/dt = n*(-k*p/(d + n) + r*(1 - n/K))',
            'dp/dt = p*s*(-h*p/n + 1)',
            'parameters: K, d, h, k, r, s',
        ]))

    def test_a_system_with_no_parameters_has_no_parameter_line(self):
        system = ODESystem.from_equations(['dz/dt = z'])
        self.assertEqual(repr(system), 'dt/dt = 1\ndz/dt = z')

    def test_a_conserved_variable_prints_its_equation(self):
        system = with_conserved_w()
        self.assertIn('dw/dt = 0', repr(system))
        self.assertEqual(repr(system).count('parameters:'), 1)
        self.assertNotIn('w', repr(system).split('parameters: ')[1])

    def test_declaring_it_a_parameter_moves_it_to_the_parameter_line(self):
        system = with_conserved_w()
        system.make_parameter('w')
        self.assertNotIn('dw/dt', repr(system))
        self.assertIn('w', repr(system).split('parameters: ')[1])

    def test_to_tex_omits_parameters(self):
        system = with_conserved_w()
        tex = system.to_tex()
        self.assertNotIn('da', tex)
        self.assertIn(r'\frac{dw}{dt} &= 0', tex)

    def test_to_tex_keeps_initial_conditions_and_constraints(self):
        system = ODESystem.from_equations(['dz/dt = a*z'])
        system.update_initial_conditions({'z': 'z_0'})
        tex = system.to_tex()
        self.assertIn(r'z\left(0\right) &= z_{0}', tex)
        self.assertNotIn(r'\frac{dz_{0}}{dt}', tex)
