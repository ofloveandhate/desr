from unittest import TestCase

import sympy

from desr.matrix_normal_forms import hnf_row
from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator
from desr.sympy_helper import monomial_to_powers


class TestODESystemScaling(TestCase):
    ''' Test ode_system.py scaling methods '''

    def test_example_pred_prey_hub_lab(self):
        ''' Predator prey model from :cite:`Hubert2013c`
            dn/dt = n( r(1 - n/K) - kp/(n+d) )
            dp/dt = sp(1 - hp / n)
        '''

        equations = '''dn/dt = n*( r*(1 - n/K) - k*p/(n+d) );dp/dt = s*p*(1 - h*p / n)'''.split(';')

        system = ODESystem.from_equations(equations)

        # Swap around some columns so we agree on the order (quick and dirty)
        answer_var = 'rhKskdtnp'
        system.reorder_variables(answer_var)

        # Take the maximal scaling matrix
        max_scal = system.maximal_scaling_matrix()
        self.assertEqual(max_scal.shape, (3, 9))

        # # Now do a couple of row ops so we get exactly the same matrix. This amount to changing the scalar
        # # operations lambda1 -> lambda1^-1, swap lambda1<->lambda2 etc
        # max_scal = max_scal.extract([1, 0, 2], range(max_scal.shape[1]))
        # max_scal[1, :] -= max_scal[2, :]

        max_scal_ans = sympy.Matrix([[-1, 0, 0, -1, -1, 0, 1, 0, 0],
                                     [0, 1, 1, 0, 1, 1, 0, 1, 0],
                                     [0, -1, 0, 0, -1, 0, 0, 0, 1]])

        # Compare row HNFs
        self.assertEqual(max_scal, hnf_row(max_scal_ans)[0])


    def test_example_6_4_hub_lab(self):
        ''' Predator prey model from :cite:`Hubert2013c`, example 6.4
            dz1/dt = z1*(1+z1*z2)
            dz2/dt = z2*(1/t - z1*z2)
        '''
        equations = '''dz1/dt = z1*(1+z1*z2);dz2/dt = z2*(1/t - z1*z2)'''.split(';')

        system = ODESystem.from_equations(equations)
        var_order = ['t', 'z1', 'z2']
        system.reorder_variables(var_order)

        ## Check against the answer in the paper
        # Match the maximal scaling matrix.  The paper prints [[0, 1, -1]]; desr reports its
        # negation, because a canonical row Hermite normal form has a positive pivot and the
        # pivot is the *last* nonzero entry in this convention, not the first.  Same row
        # lattice, so the same scaling group.  See desr.matrix_normal_forms._CONVENTION.
        max_scal = system.maximal_scaling_matrix()
        self.assertEqual(max_scal, sympy.Matrix([[0, -1, 1]]))

        # The paper's own scaling matrix and Hermite multiplier, used together below.
        max_scal_paper = sympy.Matrix([[0, 1, -1]])
        hermite_multiplier_example = sympy.Matrix([[0, 1, 0],
                                                   [1, 0, 1],
                                                   [0, 0, 1]])

        translator = ODETranslator(max_scal_paper, hermite_multiplier=hermite_multiplier_example)

        translated = translator.translate_dep_var(system)
        answer = ODESystem.from_equations('dx0/dt = x0*(y0 + 1)\ndy0/dt = y0*(1 + 1/t)')
        self.assertEqual(translated, answer)
        self.assertEqual(translator.translate_dep_var(system), translator.translate(system))

        ## Check our answer using the general translation
        translated = translator.translate_general(system)
        answer = ODESystem.from_equations('dx0/dt = x0*(y0*y1 + y0)/t\ndy0/dt = y0/t\ndy1/dt = y1*(y0 + 1)/t')
        self.assertEqual(translated, answer)

        ## Test reverse translating
        t_var, c1_var, c2_var = sympy.var('t c1 c2')
        reduced_soln = (c2_var*sympy.exp(t_var+c1_var*(1-t_var)*sympy.exp(t_var)), c1_var * t_var * sympy.exp(t_var))
        orig_soln = translator.reverse_translate_dep_var(reduced_soln, system.indep_var_index)
        self.assertTupleEqual(orig_soln, (reduced_soln[0], reduced_soln[1] / reduced_soln[0]))

        ## Check our answer hasn't changed, using our own Hermite multiplier
        translator = ODETranslator.from_ode_system(system)
        translated = translator.translate_dep_var(system)
        answer = ODESystem.from_equations('dx0/dt = x0*(-y0 - 1)\ndy0/dt = y0*(1 + 1/t)')
        self.assertEqual(translated, answer)

    def test_example_6_6_hub_lab(self):
        ''' Example 6.6 from :cite:`Hubert2013c`, where we act on time
            dz1/dt = z1*(z1**5*z2 - 2)/(3*t)
            dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)
        '''
        equations = '''dz1/dt = z1*(z1**5*z2 - 2)/(3*t);dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'''.split(';')

        system = ODESystem.from_equations(equations)
        system.reorder_variables(['t', 'z1', 'z2'])

        ## Check against the answer in the paper
        # Match the maximal scaling matrix
        max_scal = system.maximal_scaling_matrix()
        max_scal_ans = sympy.Matrix([[3, -1, 5]])
        self.assertEqual(max_scal, max_scal_ans)

        # Give Hermite multiplier from the paper. Padd it with a row and column for t to work with current infrastructure
        hermite_multiplier_ans = sympy.Matrix([[1, 1, -1],
                                               [2, 3, 2],
                                               [0, 0, 1]])

        translator = ODETranslator(max_scal_ans, hermite_multiplier=hermite_multiplier_ans)

        translated = translator.translate(system)
        answer = ODESystem.from_equations('dx0/dt = x0*(2*y0*y1/3 - 1/3)/t\ndy0/dt = y0*(y0*y1 - 1)/t\ndy1/dt = y1*(y1 + 1)/t')
        self.assertEqual(translated, answer)

        ## Check reverse translation
        reduced_soln = (sympy.var('t'),
                        sympy.sympify('c3/(t**(1/3)*(ln(t-c1)-ln(t)+c2)**(2/3))'),  # x
                        sympy.sympify('c1/(t*(ln(t-c1)-ln(t)+c2))'),  # y1
                        sympy.sympify('t/(c1 - t)'))  # y2

        orig_soln = translator.reverse_translate_general(reduced_soln)

        # Solution from the paper
        orig_soln_paper = [reduced_soln[2] / reduced_soln[1],
                                          reduced_soln[1] ** 5 * reduced_soln[3] / reduced_soln[2] ** 4]
        orig_soln_paper = [soln.subs({sympy.var('t'): sympy.sympify('t / (c3**3 / c1**2)')}) for soln in orig_soln_paper]

        self.assertEqual(len(orig_soln), len(orig_soln_paper))
        for sol1, sol2 in zip(orig_soln, orig_soln_paper):
            self.assertEqual(sol1, sol2)

        ## desr's own Hermite multiplier is now the paper's.  The normal Hermite multiplier
        ## of Prop 2.3(ii) is unique, and desr reports it in the paper's convention, so
        ## there is nothing left for the default to disagree with.  This used to be a
        ## regression guard against a different answer of desr's own.
        translator = ODETranslator.from_ode_system(system)
        self.assertEqual(translator.herm_mult, hermite_multiplier_ans)
        self.assertEqual(translator.translate(system), answer)


    def test_verhulst_log_growth(self):
        ''' Verhult logistic growth model from :cite:`Hubert2013c`
            dn/dt = r*n*(1-n/k)
        '''
        equations = '''dn/dt = r*n*(1-n/k)'''

        system = ODESystem.from_equations(equations)
        system.reorder_variables('rktn')

        ## Check against the answer in the paper
        # Match the maximal scaling matrix
        # This now matches the paper directly.  It used to need one row negated first, which
        # was the convention difference showing up: a canonical row Hermite normal form has
        # a positive pivot, and the pivot is the last nonzero entry of the row rather than
        # the first.
        max_scal = system.maximal_scaling_matrix()
        self.assertEqual(max_scal, sympy.Matrix([[-1, 0, 1, 0], [0, 1, 0, 1]]))

        herm_mult_paper = sympy.Matrix([[-1, 0, 1, 0],
                                       [0, 1, 0, -1],
                                       [0, 0, 1, 0],
                                       [0, 0, 0, 1]])
        translator = ODETranslator(max_scal,
                                    symbols_domain=system.symbols,
                                    hermite_multiplier=herm_mult_paper)
        # The paper's answer carries auxiliary variables, so it is the general reduction and
        # is asked for by name.  `translate` would now choose the parameter scheme: the
        # variable order is canonical under this convention, which it was not before, so the
        # system is eligible for a scheme it used to be passed over for.
        translated = translator.translate_general(system)
        answer = ODESystem.from_equations('dx0/dt = 0\ndx1/dt = 0\ndy0/dt = y0/t\ndy1/dt = y1*(-y0*y1 + y0)/t')
        self.assertEqual(translated, answer)

        ## Check reverse translation
        reduced_solutions = (sympy.var('t'),  # indep var t
                     sympy.var('c1'),  # x0
                     sympy.var('c2'),  # x1
                     sympy.sympify('t'),  # y0
                     sympy.sympify('1/(1+c*exp(-t))'))  # y1

        general_soln = translator.reverse_translate_general(reduced_solutions, system_indep_var_index=2)
        paper_soln = tuple(map(sympy.sympify, ['1/c1', 'c2', 'c2/(c*exp(-t/c1) + 1)']))
        self.assertTupleEqual(general_soln, paper_soln)

        ## And desr's own Hermite multiplier is the paper's, with no reordering needed.
        ## This used to reorder to 'tnrk' first -- the independent variable at the front --
        ## which was the old convention's layout, and then compared against a cached value
        ## of desr's own rather than against the paper.  Neither is necessary now.
        translator = ODETranslator.from_ode_system(system)
        self.assertEqual(translator.herm_mult, herm_mult_paper)
        self.assertEqual(translator.translate_general(system), answer)
        general_soln = translator.reverse_translate_general(reduced_solutions,
                                                            system_indep_var_index=2)
        self.assertTupleEqual(general_soln, paper_soln)

    def test_example_pred_prey_choosing_invariants(self):
        ''' Predator prey model from :cite:`Hubert2013c`
            dn/dt = n( r(1 - n/K) - kp/(n+d) )
            dp/dt = sp(1 - hp / n)

            Rather than using the invariants suggested by the algorithm, we pick our own and extend it.
        '''
        equations = '''dn/dt = n*( r*(1 - n/K) - k*p/(n+d) );dp/dt = s*p*(1 - h*p / n)'''.split(';')
        system = ODESystem.from_equations(equations)

        # Take the maximal scaling matrix
        max_scal = ODETranslator.from_ode_system(system)

        # Fetched by name, not unpacked positionally: which column a variable occupies is a
        # property of the convention and not something a test should encode.
        named = {str(v): v for v in max_scal.symbols_domain}
        t, n, p = named['t'], named['n'], named['p']
        K, d, h, k, r, s = (named[x] for x in ('K', 'd', 'h', 'k', 'r', 's'))

        self.assertTupleEqual(tuple(max_scal.invariants()),
                              (d/K, h*r/k, h*s/k, k*t/h, n/K, h*p/K))

        def as_choice(monomials):
            '''The exponent matrix of the given monomials, one per column.

            extend_from_invariants takes exponents, which means knowing where each variable
            sits.  Naming the monomials and converting here keeps the test independent of
            that -- and is the shape the method itself ought to accept.
            '''
            variables = list(max_scal.symbols_domain)
            return sympy.Matrix([monomial_to_powers(m, variables) for m in monomials]).T

        # Choose a new invariant
        max_scal2 = max_scal.extend_from_invariants(as_choice([r*t, h*p/d]))
        self.assertTupleEqual(tuple(max_scal2.invariants()),
                              (r*t, h*p/d, h*s/k, k*t/h, n/K, h*p/K))

        # This should work even if we move the time about
        max_scal3 = max_scal.extend_from_invariants(as_choice([h*p/d]))
        self.assertTupleEqual(tuple(max_scal3.invariants()),
                              (h*p/d, h*r/k, h*s/k, k*t/h, n/K, h*p/K))


        reduced_system = max_scal.translate(system=system)
        reduced_system2 = max_scal2.translate(system=system)

        self.assertNotEqual(reduced_system, reduced_system2)


class TestInitialConditions(TestCase):
    ''' Test methods relating to initial conditions '''

    def test_reduced_michaelis_menten(self):
        '''
        Michaelis Menten reaction scheme as written in equation (1) in :cite:`Gunawardena2014`, originally from from :cite:`michaelis1913kinetik`

        '''

        system_tex = '''\frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + K c - k_2 c \\\\
                 \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - K c'''

        # First, walk through the standard analysis of the reduced MM equations without an initial condition.
        original_system = ODESystem.from_tex(system_tex)

        # Parameters, then the independent variable, then the dependents.  Within the
        # parameters the leading ones are normalised to 1, so this is the old order's
        # parameter list reversed -- which keeps the reduction it was chosen to produce.
        variables = ['e_0', 'k_1', 'k_2', 'K', 't', 's', 'c']
        original_system.reorder_variables(variables)
        self.assertEqual(variables, list(map(str, original_system.symbols)))
        self.assertEqual(original_system.exponent_matrix(),
                         sympy.Matrix([[ 0,  1,  0,  0,  0,  0,  1],
                                       [ 0,  1,  0,  1,  0,  1,  1],
                                       [ 1,  0,  0,  0,  0,  0,  0],
                                       [ 0,  0,  1,  0,  1,  0,  0],
                                       [ 1,  1,  1,  1,  1,  1,  1],
                                       [-1,  0, -1,  0,  0,  1,  1],
                                       [ 1,  0,  1,  1,  0,  0, -1]]))
        max_scal1 = ODETranslator.from_ode_system(original_system)
        self.assertEqual(max_scal1.scaling_matrix,
                         sympy.Matrix([[ 0, -1, -1, -1,  1,  0,  0],
                                       [ 1, -1,  0,  0,  0,  1,  1]]))

        # Add in the initial condition s_0
        original_system.update_initial_conditions({'s': 's_0'})
        self.assertSetEqual(set(original_system.initial_conditions.items()),
                            set(map(lambda t: (sympy.sympify(t[0]), sympy.sympify(t[1])),
                                    (('s', 's_0'),))))
        # s_0 is a parameter, so it lands at the end of the parameter block and the system
        # stays in canonical order -- parameters, independent variable, dependents.
        self.assertEqual(['e_0', 'k_1', 'k_2', 'K', 's_0', 't', 's', 'c'],
                         list(map(str, original_system.symbols)))
        self.assertEqual(original_system.indep_var_index, original_system.num_parameters)
        # Moving it to the front is now a choice, not a repair: leading parameters are the
        # ones the reduction normalises, and the rest of this test wants s_0 normalised.
        original_system.reorder_variables(['s_0'] + variables)

        self.assertEqual(original_system.exponent_matrix(),
                         sympy.Matrix([[ 0,  0,  0,  0,  0,  0,  0, -1],
                                       [ 0,  1,  0,  0,  0,  0,  1,  0],
                                       [ 0,  1,  0,  1,  0,  1,  1,  0],
                                       [ 1,  0,  0,  0,  0,  0,  0,  0],
                                       [ 0,  0,  1,  0,  1,  0,  0,  0],
                                       [ 1,  1,  1,  1,  1,  1,  1,  0],
                                       [-1,  0, -1,  0,  0,  1,  1,  1],
                                       [ 1,  0,  1,  1,  0,  0, -1,  0]]))

        max_scal1 = ODETranslator.from_ode_system(original_system)
        self.assertEqual(max_scal1.scaling_matrix,
                         sympy.Matrix([[ 0,  0, -1, -1, -1,  1,  0,  0],
                                       [ 1,  1, -1,  0,  0,  0,  1,  1]]))

        # The same reduction as before the convention changed, with kappa0 and kappa2
        # exchanged: the reduced constants are numbered in parameter order, and that order
        # is the reverse of the one this test used to give.
        reduced_system = max_scal1.translate(original_system)
        self.assertSetEqual(set(reduced_system.derivative_dict.items()),
                            set(map(lambda t: (sympy.sympify(t[0]), sympy.sympify(t[1])),
                                    (('nu0', '-kappa0*nu0 - kappa1*nu1 + kappa2*nu1 + nu0*nu1'),
                                     ('nu1', 'kappa0*nu0 - kappa2*nu1 - nu0*nu1'),
                                     ('tau', 1)))))
        self.assertSetEqual(set(reduced_system.initial_conditions.items()),
                            set(map(lambda t: (sympy.sympify(t[0]), sympy.sympify(t[1])),
                                    (('nu0', 1),))))


class TestCanonicalOrder(TestCase):
    '''
    Every way of adding a variable must leave the system in canonical order.

    The order is parameters, then the independent variable, then the dependents.  A new
    parameter appended to the *end* of the variable list lands after the dependents, which
    costs the parameter reduction scheme while leaving `repr` looking correct, since `repr`
    groups the variables by classification.  That is how this escaped notice until
    `examples/linear_compartment.py` tripped over it, so the invariant is asserted directly
    here rather than inferred from a reduction succeeding.
    '''

    def assertCanonicalOrder(self, system):
        self.assertEqual(system.indep_var_index, system.num_parameters,
                         'independent variable should sit at index num_parameters, after the '
                         'parameters and before the dependents.  Order is {}'.format(
                             system.symbols))
        self.assertEqual(tuple(system.parameters) + (system.indep_var,)
                         + tuple(system.dependents),
                         system.symbols)

    def test_from_equations_is_canonical(self):
        system = ODESystem.from_equations(['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)'])
        self.assertCanonicalOrder(system)

    def test_add_constraint_naming_a_new_parameter(self):
        ''' The `examples/linear_compartment.py` case: the constraint introduces `a_2`. '''
        system = ODESystem.from_equations(['dz_1/dt = -a_1*z_1 + a_12*z_2',
                                           'dz_2/dt = -b_22*z_2 + u_2'])
        self.assertCanonicalOrder(system)
        system.add_constraint('b_22', 'a_2 + a_12')
        self.assertIn(sympy.Symbol('a_2'), system.parameters)
        self.assertCanonicalOrder(system)

    def test_add_constraint_naming_several_new_parameters(self):
        system = ODESystem.from_equations(['dx/dt = c_0*x'])
        system.add_constraint('c_0', 'a*b')
        self.assertCanonicalOrder(system)
        for name in ('a', 'b'):
            self.assertIn(sympy.Symbol(name), system.parameters)

    def test_update_initial_conditions(self):
        system = ODESystem.from_equations(['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)'])
        system.update_initial_conditions({'x': 'x_0'})
        self.assertIn(sympy.Symbol('x_0'), system.parameters)
        self.assertCanonicalOrder(system)

    def test_subs_introducing_a_new_parameter(self):
        system = ODESystem.from_equations(['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)'])
        substituted = system.subs({'c_0': 'a*b'})
        self.assertIn(sympy.Symbol('a'), substituted.parameters)
        self.assertCanonicalOrder(substituted)

    def test_a_new_parameter_respects_the_chosen_within_block_order(self):
        ''' The dependents keep the order the caller asked for; c_2 goes at the block end. '''
        system = ODESystem.from_equations(['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)'])
        system.reorder_variables(['c_1', 'c_0', 't', 'y', 'x'])
        self.assertCanonicalOrder(system)
        system.add_constraint('c_2', 'c_0 + c_1')
        self.assertEqual(('c_1', 'c_0', 'c_2', 't', 'y', 'x'),
                         tuple(map(str, system.symbols)))
        self.assertCanonicalOrder(system)

    def test_a_non_canonical_order_is_not_silently_rearranged(self):
        '''
        `reorder_variables` takes any permutation, by design, so a caller may deliberately
        leave the system non-canonical.  Adding a parameter then puts it immediately before
        the independent variable and disturbs nothing else -- it does not quietly shuffle
        the system into canonical order behind the caller's back.
        '''
        system = ODESystem.from_equations(['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)'])
        system.reorder_variables(['c_0', 't', 'x', 'y', 'c_1'])  # an interleaved order
        system.add_constraint('c_2', 'c_0 + c_1')
        self.assertEqual(('c_0', 'c_2', 't', 'x', 'y', 'c_1'),
                         tuple(map(str, system.symbols)))
        self.assertIn(sympy.Symbol('c_2'), system.parameters)


class TestBlockReordering(TestCase):
    '''
    `reorder_parameters` and `reorder_dependents` let the caller say only what is
    meaningful: the order *within* a block. Which block comes first is forced by the
    Hermite convention and is not theirs to think about.

    `reorder_variables` stays as it was, taking a full permutation, so an interleaved order
    remains expressible -- the user was explicit that one function could not do both jobs.
    '''

    def setUp(self):
        self.system = ODESystem.from_equations(
            ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
             'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])

    def names(self, variables):
        return tuple(map(str, variables))

    def assertCanonical(self, system):
        self.assertEqual(system.indep_var_index, system.num_parameters)

    def test_a_partial_parameter_list_brings_those_to_the_front(self):
        self.system.reorder_parameters(['k_2', 'k_m1'])
        self.assertEqual(('k_2', 'k_m1', 'e_0', 'k_1'), self.names(self.system.parameters))
        self.assertCanonical(self.system)

    def test_a_single_parameter_need_not_be_in_a_list(self):
        self.system.reorder_parameters('k_2')
        self.assertEqual(('k_2', 'e_0', 'k_1', 'k_m1'), self.names(self.system.parameters))

    def test_a_full_parameter_list_is_taken_exactly(self):
        self.system.reorder_parameters(['k_m1', 'k_2', 'k_1', 'e_0'])
        self.assertEqual(('k_m1', 'k_2', 'k_1', 'e_0'), self.names(self.system.parameters))

    def test_an_empty_order_changes_nothing(self):
        before = self.system.symbols
        self.system.reorder_parameters([])
        self.assertEqual(before, self.system.symbols)

    def test_reordering_parameters_leaves_the_dependents_alone(self):
        before = self.system.dependents
        self.system.reorder_parameters('k_2')
        self.assertEqual(before, self.system.dependents)

    def test_reordering_dependents_leaves_the_parameters_alone(self):
        before = self.system.parameters
        self.system.reorder_dependents('s')
        self.assertEqual(('s', 'c'), self.names(self.system.dependents))
        self.assertEqual(before, self.system.parameters)
        self.assertCanonical(self.system)

    def test_the_derivatives_travel_with_their_variables(self):
        before = dict(self.system.derivative_dict)
        self.system.reorder_parameters('k_2')
        self.system.reorder_dependents('s')
        self.assertEqual(before, dict(self.system.derivative_dict))

    def test_a_dependent_named_as_a_parameter_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            self.system.reorder_parameters('s')
        self.assertIn('a dependent variable', str(caught.exception))

    def test_the_independent_variable_named_as_a_dependent_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            self.system.reorder_dependents('t')
        self.assertIn('the independent variable', str(caught.exception))

    def test_the_independent_variable_named_as_a_parameter_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            self.system.reorder_parameters('t')
        self.assertIn('the independent variable', str(caught.exception))

    def test_an_unknown_name_is_refused_and_says_what_the_block_holds(self):
        with self.assertRaises(ValueError) as caught:
            self.system.reorder_parameters('nonesuch')
        message = str(caught.exception)
        self.assertIn('not a variable of this system', message)
        self.assertIn('k_1', message)

    def test_a_repeated_name_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            self.system.reorder_parameters(['k_2', 'k_2'])
        self.assertIn('twice', str(caught.exception))

    def test_an_interleaved_order_is_canonicalised(self):
        '''
        `reorder_variables` permits a non-canonical layout on purpose. The block reorderers
        are the ergonomic layer, so they put the blocks back -- documented, and the reason
        they exist.
        '''
        self.system.reorder_variables(['k_1', 't', 's', 'c', 'e_0', 'k_2', 'k_m1'])
        self.assertNotEqual(self.system.indep_var_index, self.system.num_parameters)
        self.system.reorder_parameters('k_2')
        self.assertCanonical(self.system)
        self.assertEqual(('k_2', 'k_1', 'e_0', 'k_m1'), self.names(self.system.parameters))

    def normalised_parameters(self, system):
        ''' The parameters this system's reduction sends to 1. '''
        translator = ODETranslator.from_ode_system(system)
        substitutions = translator.translate_parameter_substitutions(system)
        return [str(v) for v in system.parameters if substitutions[v] == 1]

    def test_reordering_parameters_changes_which_are_normalised(self):
        ''' The whole point: moving a parameter to the front makes it a candidate. '''
        self.system.reorder_parameters(['e_0', 'k_1'])
        self.assertEqual(['e_0', 'k_1'], self.normalised_parameters(self.system))
        self.system.reorder_parameters('k_2')
        self.assertIn('k_2', self.normalised_parameters(self.system))

    def test_leading_position_is_a_preference_and_not_a_guarantee(self):
        '''
        Exactly `r` parameters are normalised, and the leading ones are only the natural
        candidates -- the Hermite structure decides which set actually *can* be, and skips
        one that cannot participate.  Here `k_2` and `k_m1` cannot both go, so asking for
        both gets `k_2` and then `e_0`, not `k_m1`.  `reorder_parameters`' docstring makes
        this claim, so it is pinned.
        '''
        self.system.reorder_parameters(['k_2', 'k_m1'])
        self.assertEqual(('k_2', 'k_m1', 'e_0', 'k_1'), self.names(self.system.parameters))
        normalised = self.normalised_parameters(self.system)
        self.assertEqual(2, len(normalised))
        self.assertIn('k_2', normalised)
        self.assertNotIn('k_m1', normalised)
        self.assertIn('e_0', normalised)


class TestSubsCallShapes(TestCase):
    ''' `subs` is named after sympy's, so it takes the argument shapes sympy's takes. '''

    def setUp(self):
        self.system = ODESystem.from_equations(['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)'])
        self.expected = self.system.subs({'c_0': 'a'})

    def test_two_arguments(self):
        self.assertEqual(self.expected, self.system.subs('c_0', 'a'))

    def test_iterable_of_pairs(self):
        self.assertEqual(self.expected, self.system.subs([('c_0', 'a')]))

    def test_tuple_of_pairs(self):
        self.assertEqual(self.expected, self.system.subs((('c_0', 'a'),)))

    def test_several_pairs_match_the_dict(self):
        expected = self.system.subs({'c_0': 'a', 'c_1': 'b'})
        self.assertEqual(expected, self.system.subs([('c_0', 'a'), ('c_1', 'b')]))

    def test_a_lone_symbol_is_refused(self):
        with self.assertRaises(ValueError):
            self.system.subs('c_0')

    def test_too_many_arguments_are_refused(self):
        with self.assertRaises(ValueError):
            self.system.subs('c_0', 'a', 'b')

    def test_a_group_that_is_not_a_pair_is_refused(self):
        with self.assertRaises(ValueError):
            self.system.subs([('c_0', 'a', 'b')])

    def test_lamda_is_the_spelling_that_parses(self):
        '''
        `lambda` is a reserved word, so sympy cannot parse it and uses `lamda`, which desr
        follows.  Worth pinning: `lam` also parses, but renders as literal `lam` in TeX.
        '''
        from desr.tex_tools import expr_to_tex
        substituted = self.system.subs({'c_0': 'lamda'})
        self.assertIn(sympy.Symbol('lamda'), substituted.parameters)
        self.assertIn(r'\lambda', expr_to_tex(substituted.derivative_dict[sympy.Symbol('x')]))
        with self.assertRaises(Exception):
            self.system.subs({'c_0': 'lambda'})
