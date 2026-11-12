from unittest import TestCase

import sympy

from desr.matrix_normal_forms import hnf_row
from desr.ode_system import ConflictingParameterValues, ODESystem
from desr.ode_translation import ODETranslator, ReductionScheme
from desr.sympy_helper import CONSISTENCY_TOLERANCE, monomial_to_powers


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

        ## Check our answer using the general translation.  The paper's answer, so
        ## both of the further steps are turned off: they are departures from it.
        translated = translator.translate_general(system, normalise_parameters=False,
                                                  identify_indep_var=False)
        answer = ODESystem.from_equations('dx0/dt = x0*(y0*y1 + y0)/t\ndy0/dt = y0/t\ndy1/dt = y1*(y0 + 1)/t')
        self.assertEqual(translated, answer)

        ## By default y0 is identified with t -- its monomial is t itself -- so the
        ## equation dy0/dt = y0/t becomes dt/dt = 1, which the system already has.
        default = translator.translate_general(system)
        self.assertEqual(('t', 'x0', 'y1'), tuple(map(str, default.symbols)))

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
        translated = translator.translate_general(system, normalise_parameters=False,
                                                  identify_indep_var=False)
        self.assertEqual(('t', 'x0', 'x1', 'y0', 'y1'),
                         tuple(map(str, translated.symbols)))
        # x0 and x1 are monomials in the original parameters alone, so they are parameters
        # of the reduced system -- not dependents whose equation says zero.  They are also
        # the auxiliaries, which is a separate fact and recorded separately.
        self.assertEqual(('x0', 'x1'), tuple(map(str, translated.parameters)))
        self.assertEqual(('x0', 'x1'), tuple(map(str, translated.auxiliary_names)))
        self.assertEqual((), translated.conserved_variables)
        self.assertEqual({'y0': 'y0/t', 'y1': 'y1*(-y0*y1 + y0)/t'},
                         {str(v): str(translated.derivative_dict[v])
                          for v in translated.dependents})

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
        # Same multiplier, so the same reduced system as the one checked above -- asked
        # for the same way, since `translated` above is the paper's reduction.
        self.assertEqual(translator.translate_general(system, normalise_parameters=False,
                                                      identify_indep_var=False),
                         translated)
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

        self.assertTupleEqual(tuple(max_scal.invariant_monomials()),
                              (d/K, h*r/k, h*s/k, k*t/h, n/K, h*p/K))

        # Choose a new invariant, by writing it down.  This used to need the exponent
        # matrix, which meant knowing where each variable sat in the order.
        max_scal2 = max_scal.extend_from_invariants([r*t, h*p/d])
        self.assertTupleEqual(tuple(max_scal2.invariant_monomials()),
                              (r*t, h*p/d, h*s/k, k*t/h, n/K, h*p/K))

        # This should work even if we move the time about
        max_scal3 = max_scal.extend_from_invariants(h*p/d)
        self.assertTupleEqual(tuple(max_scal3.invariant_monomials()),
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


class TestConstraintCallShapes(TestCase):
    '''
    A constraint is one equation on paper, so it can be given as one.
    '''

    EQUATIONS = ['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)']

    def system(self):
        return ODESystem.from_equations(self.EQUATIONS)

    def expected(self):
        c_0, c_1, c_2 = sympy.symbols('c_0 c_1 c_2')
        return sympy.Eq(c_2, c_0 + c_1)

    def test_two_sides_still_work(self):
        system = self.system()
        system.add_constraint('c_2', 'c_0 + c_1')
        self.assertEqual(system.constraints, [self.expected()])

    def test_one_string_says_the_same_thing(self):
        system = self.system()
        system.add_constraint('c_2 = c_0 + c_1')
        self.assertEqual(system.constraints, [self.expected()])

    def test_a_double_equals_reads_the_same_way(self):
        # Which is how a constraint prints, so retyping what was printed works.
        system = self.system()
        system.add_constraint('c_2 == c_0 + c_1')
        self.assertEqual(system.constraints, [self.expected()])

    def test_whitespace_around_the_equals_is_optional(self):
        system = self.system()
        system.add_constraint('c_2=c_0 + c_1')
        self.assertEqual(system.constraints, [self.expected()])

    def test_a_sympy_equality_is_taken_whole(self):
        system = self.system()
        system.add_constraint(self.expected())
        self.assertEqual(system.constraints, [self.expected()])

    def test_two_equals_signs_is_refused(self):
        system = self.system()
        with self.assertRaises(ValueError) as caught:
            system.add_constraint('c_2 = c_0 = c_1')
        self.assertIn('not one constraint', str(caught.exception))

    def test_one_side_alone_is_refused(self):
        system = self.system()
        with self.assertRaises(ValueError) as caught:
            system.add_constraint('c_2')
        self.assertIn('0 equals signs', str(caught.exception))

    def test_a_bare_expression_alone_is_refused(self):
        system = self.system()
        with self.assertRaises(ValueError) as caught:
            system.add_constraint(sympy.Symbol('c_2'))
        self.assertIn('one side of an equation', str(caught.exception))

    def test_the_whole_form_is_checked_like_the_two_sided_one(self):
        # Same path, so the dependent-variable refusal still applies.
        system = self.system()
        with self.assertRaises(ValueError) as caught:
            system.add_constraint('c_2 = c_0 + x')
        self.assertIn('dependent variables', str(caught.exception))


class TestParameterValues(TestCase):
    '''
    A parameter is a constant of the model, so its numeric value is true everywhere and is
    a fact the system can hold.  That is the whole difference from `desr.numerics`'
    `known_values`, which are values at one point on a solution and may name a dependent
    variable.
    '''

    EQUATIONS = ['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)']

    def system(self):
        return ODESystem.from_equations(self.EQUATIONS)

    def test_a_value_is_recorded_against_the_parameter(self):
        system = self.system()
        system.update_parameter_values({'c_0': 2})
        self.assertEqual(system.parameter_values, {sympy.Symbol('c_0'): sympy.Integer(2)})

    def test_the_accessor_hands_back_a_copy(self):
        system = self.system()
        system.update_parameter_values({'c_0': 2})
        system.parameter_values[sympy.Symbol('c_1')] = 9
        self.assertNotIn(sympy.Symbol('c_1'), system.parameter_values)

    def test_values_accumulate_across_calls(self):
        system = self.system()
        system.update_parameter_values({'c_0': 2})
        system.update_parameter_values({'c_1': 3})
        self.assertEqual(set(map(str, system.parameter_values)), {'c_0', 'c_1'})

    def test_setting_the_same_parameter_again_replaces_the_value(self):
        # An edit by the caller, not a conflict: there is only one of them, and the second
        # call is the later word.  A conflict is between two sources that both claim to be
        # right, which is what the constraint check and `known_values` are about.
        system = self.system()
        system.update_parameter_values({'c_0': 2})
        system.update_parameter_values({'c_0': 5})
        self.assertEqual(system.parameter_values[sympy.Symbol('c_0')], sympy.Integer(5))

    def test_a_dependent_variable_is_refused_and_told_where_to_go(self):
        system = self.system()
        with self.assertRaises(ValueError) as caught:
            system.update_parameter_values({'x': 2})
        message = str(caught.exception)
        self.assertIn('dependent variable', message)
        self.assertIn('known_values', message)

    def test_the_independent_variable_is_refused(self):
        system = self.system()
        with self.assertRaises(ValueError) as caught:
            system.update_parameter_values({'t': 0})
        self.assertIn('independent variable', str(caught.exception))

    def test_a_symbol_the_system_does_not_have_is_refused(self):
        system = self.system()
        with self.assertRaises(ValueError) as caught:
            system.update_parameter_values({'q': 1})
        self.assertIn('not a symbol of this system', str(caught.exception))

    def test_an_expression_is_refused_and_told_it_is_a_constraint(self):
        system = self.system()
        with self.assertRaises(ValueError) as caught:
            system.update_parameter_values({'c_1': 'c_0 + 1'})
        self.assertIn('add_constraint', str(caught.exception))

    def test_a_float_is_kept(self):
        system = self.system()
        system.update_parameter_values({'c_0': 0.25})
        self.assertAlmostEqual(float(system.parameter_values[sympy.Symbol('c_0')]), 0.25)

    def test_zero_is_an_ordinary_value_here(self):
        # Zero is refused as a *scale*; as a value it is a number like any other, and this
        # is where a zero initial condition's value will live.
        system = self.system()
        system.update_parameter_values({'c_0': 0})
        self.assertEqual(system.parameter_values[sympy.Symbol('c_0')], sympy.Integer(0))

    def test_a_rejected_key_leaves_no_partial_update(self):
        system = self.system()
        with self.assertRaises(ValueError):
            system.update_parameter_values({'c_0': 2, 'x': 3})
        self.assertEqual(system.parameter_values, {})


class TestParameterValuesAgainstConstraints(TestCase):
    '''
    A constraint says something true of the parameters, so values that fail it and a
    constraint the values fail are the same mistake from two sides.  Whichever arrives
    second is refused.
    '''

    EQUATIONS = ['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)']

    def system(self):
        return ODESystem.from_equations(self.EQUATIONS)

    def test_values_consistent_with_a_constraint_are_accepted(self):
        system = self.system()
        system.add_constraint('c_2 = c_0 + c_1')
        system.update_parameter_values({'c_0': 2, 'c_1': 3, 'c_2': 5})
        self.assertEqual(system.parameter_values[sympy.Symbol('c_2')], sympy.Integer(5))

    def test_values_inconsistent_with_a_constraint_are_refused(self):
        system = self.system()
        system.add_constraint('c_2 = c_0 + c_1')
        with self.assertRaises(ConflictingParameterValues) as caught:
            system.update_parameter_values({'c_0': 2, 'c_1': 3, 'c_2': 6})
        message = str(caught.exception)
        self.assertIn('c_2 == c_0 + c_1', message)
        self.assertIn('6 == 5', message)

    def test_a_constraint_inconsistent_with_standing_values_is_refused(self):
        system = self.system()
        system.update_parameter_values({'c_0': 2, 'c_1': 3})
        with self.assertRaises(ConflictingParameterValues):
            system.add_constraint('c_1 = c_0')

    def test_a_refused_constraint_is_not_added(self):
        system = self.system()
        system.update_parameter_values({'c_0': 2, 'c_1': 3})
        with self.assertRaises(ConflictingParameterValues):
            system.add_constraint('c_1 = c_0')
        self.assertEqual(system.constraints, [])

    def test_a_partly_valued_constraint_says_nothing_yet(self):
        # c_2 has no value, so the constraint can still be satisfied by choosing one.
        system = self.system()
        system.add_constraint('c_2 = c_0 + c_1')
        system.update_parameter_values({'c_0': 2, 'c_1': 3})
        self.assertNotIn(sympy.Symbol('c_2'), system.parameter_values)

    def test_floats_agree_within_the_shared_tolerance(self):
        system = self.system()
        system.add_constraint('c_2 = c_0 + c_1')
        system.update_parameter_values({'c_0': 1.0, 'c_1': 2.0,
                                        'c_2': 3.0 + CONSISTENCY_TOLERANCE / 10})

    def test_floats_further_apart_than_the_tolerance_are_refused(self):
        system = self.system()
        system.add_constraint('c_2 = c_0 + c_1')
        with self.assertRaises(ConflictingParameterValues):
            system.update_parameter_values({'c_0': 1.0, 'c_1': 2.0,
                                            'c_2': 3.0 + CONSISTENCY_TOLERANCE * 100})

    def test_the_tolerance_is_one_constant_for_both_halves_of_the_library(self):
        # ode_system cannot import numerics -- numpy would become mandatory for
        # `import desr` -- so the constant lives in sympy_helper and numerics re-exports
        # it.  A value the system accepted and the numerics then rejected would be a
        # contradiction the caller could do nothing about.
        from desr import numerics
        self.assertIs(numerics.CONSISTENCY_TOLERANCE, CONSISTENCY_TOLERANCE)


class TestCopyCarriesEverything(TestCase):
    '''
    `copy` is on the critical path of every reduction: a translator stores a copy of the
    system it was built for, and a reduced system stores a copy of the system it came
    from.  A field `copy` forgets is a field that silently empties, everywhere.

    So these compare the whole of `__dict__` rather than a list of fields anyone has to
    remember to extend.  A field added to `ODESystem` later and left out of `copy` fails
    here without being named.
    '''

    def fully_populated(self):
        '''An original system with every optional piece of state set.'''
        system = ODESystem.from_equations(['dx/dt = c_0*x*y',
                                           'dy/dt = c_1*(1-x)*(1-y)'])
        system.update_initial_conditions({'x': 'x_0'})
        system.add_constraint('c_0', 'c_1*2')
        system.update_parameter_values({'c_1': 3})
        system.require_same_unit({sympy.Symbol('x'): sympy.Symbol('y')})
        # Populated, not merely present: comparing two empty dicts passes while
        # checking nothing, which is how `copy` came to drop this one silently.
        system.set_dimensions({'x': 'concentration', 'y': 'concentration'})
        return system

    def a_reduction(self):
        '''A reduced system, which carries the provenance an original one does not.'''
        system = ODESystem.from_equations(
            ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)'])
        return ODETranslator.from_ode_system(system).translate_general(
            system, normalise_parameters=False)

    def assert_copy_is_complete(self, system):
        duplicate = system.copy()
        self.assertEqual(sorted(vars(system)), sorted(vars(duplicate)),
                         'copy produced an object with different state')
        for field in sorted(vars(system)):
            with self.subTest(field=field):
                self.assertEqual(vars(system)[field], vars(duplicate)[field])

    def test_an_original_system_copies_completely(self):
        self.assert_copy_is_complete(self.fully_populated())

    def test_a_reduced_system_copies_completely(self):
        self.assert_copy_is_complete(self.a_reduction())

    def test_the_provenance_a_reduction_carries_is_not_empty(self):
        '''
        Guards the test above: comparing two empty records passes while saying nothing, and
        prey-predator's general reduction keeps no auxiliaries once they are normalised,
        which is how a check like that goes quietly hollow.
        '''
        reduced = self.a_reduction()
        self.assertTrue(reduced.named_symbols)
        self.assertTrue(reduced.symbol_kinds)
        self.assertTrue(reduced.auxiliary_names)
        self.assertIsNotNone(reduced.original_system)
        self.assertIsNotNone(reduced.scheme)

    def test_a_copy_does_not_share_mutable_state(self):
        ''' A copy that shares a dict is not a copy; editing one would edit the other. '''
        system = self.fully_populated()
        duplicate = system.copy()
        duplicate.update_parameter_values({'c_0': 6})
        self.assertNotIn(sympy.Symbol('c_0'), system.parameter_values)

    def test_editing_the_original_does_not_reach_a_translator(self):
        '''
        A translator stores a copy for this reason: the scaling matrix and multiplier
        describe the system as it was, so a later edit is not one the translator could
        honour even if it saw it.
        '''
        system = ODESystem.from_equations(['dx/dt = c_0*x'])
        translator = ODETranslator.from_ode_system(system)
        system.update_parameter_values({'c_0': 5})
        self.assertEqual({}, translator.system.parameter_values)
        self.assertEqual(sympy.Integer(5), system.parameter_values[sympy.Symbol('c_0')])


class TestParameterValuesSurviveTheSystem(TestCase):
    ''' Values are part of the system, so they travel with it. '''

    EQUATIONS = ['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)']

    def system(self):
        return ODESystem.from_equations(self.EQUATIONS)

    def test_copy_carries_the_values(self):
        system = self.system()
        system.update_parameter_values({'c_0': 2})
        self.assertEqual(system.copy().parameter_values, system.parameter_values)

    def test_the_constructor_takes_them(self):
        system = self.system()
        rebuilt = ODESystem(system.symbols, system.derivatives,
                            parameter_values={'c_0': 2})
        self.assertEqual(rebuilt.parameter_values, {sympy.Symbol('c_0'): sympy.Integer(2)})

    def test_the_constructor_checks_them_against_its_constraints(self):
        system = self.system()
        system.add_constraint('c_1 = c_0')
        with self.assertRaises(ConflictingParameterValues):
            ODESystem(system.symbols, system.derivatives,
                      constraints=list(system.constraints),
                      parameter_values={'c_0': 2, 'c_1': 3})

    def test_repr_names_the_values(self):
        system = self.system()
        system.update_parameter_values({'c_0': 2})
        self.assertIn('parameter values: c_0 = 2', repr(system))

    def test_repr_is_unchanged_when_no_value_is_known(self):
        self.assertNotIn('parameter values', repr(self.system()))

    def test_a_renaming_substitution_carries_the_value_to_the_new_name(self):
        system = self.system()
        system.update_parameter_values({'c_0': 2})
        renamed = system.subs({'c_0': 'd_0'})
        self.assertEqual(renamed.parameter_values,
                         {sympy.Symbol('d_0'): sympy.Integer(2)})

    def test_substituting_a_parameter_to_its_known_value_is_fine(self):
        system = self.system()
        system.update_parameter_values({'c_0': 2})
        substituted = system.subs({'c_0': '2'})
        self.assertEqual(substituted.parameter_values, {})

    def test_substituting_a_parameter_to_a_different_number_is_refused(self):
        # Two statements about the same quantity.  Dropping either one quietly is how a
        # wrong number reaches a paper.
        system = self.system()
        system.update_parameter_values({'c_0': 2})
        with self.assertRaises(ConflictingParameterValues) as caught:
            system.subs({'c_0': '7'})
        self.assertIn('c_0', str(caught.exception))


class TestNumericInitialConditions(TestCase):
    '''
    A number is a value, and in a system that has not been reduced a value still has
    units: a `0` for a concentration is a concentration, not the dimensionless `0` of a
    reduced system.  So the number becomes the value of a symbol, and the symbol carries
    the units.  Which is what lets a zero initial condition be stated at all.
    '''

    EQUATIONS = ['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)']

    def system(self):
        return ODESystem.from_equations(self.EQUATIONS)

    def test_a_number_becomes_a_symbol_with_that_value(self):
        system = self.system()
        system.update_initial_conditions({'x': 3})
        self.assertEqual(system.initial_conditions,
                         {sympy.Symbol('x'): sympy.Symbol('x_0')})
        self.assertEqual(system.parameter_values,
                         {sympy.Symbol('x_0'): sympy.Integer(3)})

    def test_zero_is_allowed_and_is_the_whole_point(self):
        system = self.system()
        system.update_initial_conditions({'x': 0})
        self.assertEqual(system.parameter_values[sympy.Symbol('x_0')], sympy.Integer(0))

    def test_the_minted_symbol_is_a_parameter(self):
        system = self.system()
        system.update_initial_conditions({'x': 0})
        self.assertIn(sympy.Symbol('x_0'), system.parameters)

    def test_the_symbol_enters_the_scaling_as_a_units_statement(self):
        # `x/x_0` is true about units whatever x_0's value turns out to be, so the
        # symbolic side is identical to naming the symbol by hand.
        by_number, by_name = self.system(), self.system()
        by_number.update_initial_conditions({'x': 0})
        by_name.update_initial_conditions({'x': 'x_0'})
        self.assertEqual(by_number.maximal_scaling_matrix(),
                         by_name.maximal_scaling_matrix())

    def test_a_float_works_too(self):
        system = self.system()
        system.update_initial_conditions({'x': 0.5})
        self.assertAlmostEqual(
            float(system.parameter_values[sympy.Symbol('x_0')]), 0.5)

    def test_a_bare_symbol_still_carries_no_value(self):
        system = self.system()
        system.update_initial_conditions({'x': 'x_0'})
        self.assertEqual(system.parameter_values, {})

    def test_an_existing_parameter_of_that_name_is_reused(self):
        system = self.system()
        system.update_initial_conditions({'x': 'x_0'})
        system.update_initial_conditions({'x': 0})
        self.assertEqual(system.initial_conditions,
                         {sympy.Symbol('x'): sympy.Symbol('x_0')})
        self.assertEqual(system.parameter_values[sympy.Symbol('x_0')],
                         sympy.Integer(0))

    def test_a_name_already_another_variables_initial_condition_is_refused(self):
        system = self.system()
        system.update_initial_conditions({'y': 'x_0'})
        with self.assertRaises(ValueError) as caught:
            system.update_initial_conditions({'x': 0})
        self.assertIn('already the initial condition of y', str(caught.exception))

    def test_a_name_taken_by_a_variable_is_refused(self):
        system = ODESystem.from_equations(['dx/dt = c_0*x*x_0',
                                           'dx_0/dt = -c_0*x*x_0'])
        with self.assertRaises(ValueError) as caught:
            system.update_initial_conditions({'x': 0})
        self.assertIn('already a variable', str(caught.exception))

    def test_an_expression_is_still_refused(self):
        system = self.system()
        with self.assertRaises(ValueError) as caught:
            system.update_initial_conditions({'x': 'c_0 + c_1'})
        self.assertIn('c_0, c_1', str(caught.exception))

    def test_a_reduced_system_keeps_the_number_as_a_number(self):
        # Its symbols are dimensionless, so there is nothing for a symbol to carry.  This
        # is how `nu0(0) = 1` comes out of a reduction.
        system = ODESystem.from_equations(self.EQUATIONS, is_reduced=True)
        system.update_initial_conditions({'x': 1})
        self.assertEqual(system.initial_conditions,
                         {sympy.Symbol('x'): sympy.Integer(1)})
        self.assertEqual(system.parameter_values, {})


class TestZeroInitialConditionsThroughAReduction(TestCase):
    '''
    Michaelis-Menten with `C(0) = 0` and `P(0) = 0` stated rather than skipped, which is
    what the narrative has always had to gloss over.
    '''

    TEX = '\n'.join(['\\frac{dS}{dt} &= - k_1 E S + k_{-1} C \\\\',
                     '\\frac{dE}{dt} &= - k_1 E S + k_{-1} C + k_2 C \\\\',
                     '\\frac{dC}{dt} &= k_1 E S - k_{-1} C - k_2 C \\\\',
                     '\\frac{dP}{dt} &= k_2 C'])

    def system(self):
        system = ODESystem.from_tex(self.TEX)
        system.update_initial_conditions({'S': 'S_0', 'E': 'E_0', 'C': 0, 'P': 0})
        return system

    def test_the_zeros_can_be_stated(self):
        system = self.system()
        self.assertEqual(system.parameter_values,
                         {sympy.Symbol('C_0'): sympy.Integer(0),
                          sympy.Symbol('P_0'): sympy.Integer(0)})

    def test_the_default_order_does_not_scale_by_them(self):
        # New parameters land at the back of the parameter block, and the back is the keep
        # end, so a zero-valued initial condition is out of harm's way by default.
        system = self.system()
        translator = ODETranslator.from_ode_system(system)
        must_not_be_zero = translator.symbols_that_must_not_be_zero()
        self.assertNotIn(sympy.Symbol('C_0'), must_not_be_zero)
        self.assertNotIn(sympy.Symbol('P_0'), must_not_be_zero)

    def test_the_reduction_goes_through(self):
        system = self.system()
        reduced = ODETranslator.from_ode_system(system).translate(system)
        self.assertTrue(reduced.is_reduced)

    def test_the_zero_is_visible_in_the_reduced_system(self):
        # Two reduced parameters are zero, each the invariant holding one of C_0, P_0 in
        # its numerator.  Derived, not restated.
        system = self.system()
        reduced = ODETranslator.from_ode_system(system).translate(system)
        zero_valued = [parameter for parameter, value in reduced.parameter_values.items()
                       if value == 0]
        self.assertEqual(len(zero_valued), 2)

    def test_each_zero_reduced_parameter_is_an_initial_condition(self):
        system = self.system()
        reduced = ODETranslator.from_ode_system(system).translate(system)
        zero_valued = {parameter for parameter, value in reduced.parameter_values.items()
                       if value == 0}
        self.assertTrue(zero_valued <= set(reduced.initial_conditions.values()))

    def test_scaling_by_a_zero_is_refused_if_the_order_asks_for_it(self):
        system = self.system()
        system.reorder_parameters(['C_0', 'P_0', 'k_1', 'k_2', 'k_m1', 'S_0', 'E_0'])
        translator = ODETranslator.from_ode_system(system)
        with self.assertRaises(ValueError) as caught:
            translator.translate(system)
        self.assertIn('C_0', str(caught.exception))

    def test_naming_a_scale_that_already_exists_adds_no_symbol(self):
        # The other way to give a zero-start variable a scale: Segel & Slemrod measure the
        # complex in units of e_0, not of c(0).  No new symbol, so no extra reduced
        # parameter -- which is the cost of minting one.
        system = ODESystem.from_tex(self.TEX)
        system.update_initial_conditions({'S': 'S_0', 'E': 'E_0'})
        before = len(system.symbols)
        system.require_same_unit({'C': 'E_0'})
        self.assertEqual(len(system.symbols), before)


class TestInitialConditionCallShapes(TestCase):
    '''
    An initial condition is `C(0) = 0` on paper, so it can be written that way, and the
    bare `C` that only this method ever needed still works.  All of the shapes mean the
    same thing.
    '''

    EQUATIONS = ['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)']

    def system(self):
        return ODESystem.from_equations(self.EQUATIONS)

    def expected(self):
        return {sympy.Symbol('x'): sympy.Symbol('x_0')}

    def test_the_bare_variable_as_a_key(self):
        system = self.system()
        system.update_initial_conditions({'x': 0})
        self.assertEqual(system.initial_conditions, self.expected())

    def test_the_variable_at_zero_as_a_key(self):
        system = self.system()
        system.update_initial_conditions({'x(0)': 0})
        self.assertEqual(system.initial_conditions, self.expected())

    def test_one_equation_as_a_string(self):
        system = self.system()
        system.update_initial_conditions('x(0) = 0')
        self.assertEqual(system.initial_conditions, self.expected())

    def test_a_list_of_equations(self):
        system = self.system()
        system.update_initial_conditions(['x(0) = 0', 'y(0) = y_0'])
        self.assertEqual(set(map(str, system.initial_conditions)), {'x', 'y'})

    def test_a_set_of_equations(self):
        system = self.system()
        system.update_initial_conditions({'x(0) = 0'})
        self.assertEqual(system.initial_conditions, self.expected())

    def test_a_double_equals_reads_the_same_way(self):
        system = self.system()
        system.update_initial_conditions('x(0) == 0')
        self.assertEqual(system.initial_conditions, self.expected())

    def test_the_value_may_come_first(self):
        system = self.system()
        system.update_initial_conditions('x_0 = x(0)')
        self.assertEqual(system.initial_conditions, self.expected())

    def test_a_sympy_equality_works(self):
        system = self.system()
        system.update_initial_conditions(
            sympy.Eq(sympy.Function('x')(0), sympy.Integer(0)))
        self.assertEqual(system.initial_conditions, self.expected())

    def test_whitespace_inside_the_parenthesis_is_fine(self):
        system = self.system()
        system.update_initial_conditions({' x ( 0 ) ': 0})
        self.assertEqual(system.initial_conditions, self.expected())

    def test_all_the_shapes_agree(self):
        shapes = [{'x': 0}, {'x(0)': 0}, 'x(0) = 0', ['x(0) = 0'], {'x(0) = 0'}]
        results = []
        for shape in shapes:
            system = self.system()
            system.update_initial_conditions(shape)
            results.append((system.initial_conditions, system.parameter_values))
        for shape, result in zip(shapes[1:], results[1:]):
            with self.subTest(shape=shape):
                self.assertEqual(result, results[0])

    def test_an_equation_with_no_variable_at_zero_is_refused(self):
        system = self.system()
        with self.assertRaises(ValueError) as caught:
            system.update_initial_conditions('x = 0')
        self.assertIn('Neither side', str(caught.exception))

    def test_two_equals_signs_is_refused(self):
        system = self.system()
        with self.assertRaises(ValueError) as caught:
            system.update_initial_conditions('x(0) = 0 = 0')
        self.assertIn('not one initial condition', str(caught.exception))

    def test_one_side_alone_is_refused(self):
        system = self.system()
        with self.assertRaises(ValueError) as caught:
            system.update_initial_conditions(['x(0)'])
        self.assertIn('0 equals signs', str(caught.exception))


class TestTheTwoFacesOfAValue(TestCase):
    '''
    A reduced parameter is an invariant and an invariant is a monomial, so sometimes the
    expression is what you want to read and sometimes the number is.  Both are available,
    before and after a reduction.
    '''

    TEX = '\n'.join(['\\frac{dS}{dt} &= - k_1 E S + k_{-1} C \\\\',
                     '\\frac{dE}{dt} &= - k_1 E S + k_{-1} C + k_2 C \\\\',
                     '\\frac{dC}{dt} &= k_1 E S - k_{-1} C - k_2 C \\\\',
                     '\\frac{dP}{dt} &= k_2 C'])

    def system(self):
        system = ODESystem.from_tex(self.TEX)
        system.update_initial_conditions(['S(0) = S_0', 'E(0) = E_0',
                                          'C(0) = 0', 'P(0) = 0'])
        return system

    def test_the_expression_face_names_the_symbol(self):
        system = self.system()
        self.assertEqual(system.initial_conditions[sympy.Symbol('C')],
                         sympy.Symbol('C_0'))

    def test_the_value_face_gives_the_number(self):
        system = self.system()
        self.assertEqual(system.initial_condition_values[sympy.Symbol('C')],
                         sympy.Integer(0))

    def test_a_variable_whose_initial_value_is_unknown_is_absent(self):
        # Absent rather than present and meaningless.
        system = self.system()
        self.assertNotIn(sympy.Symbol('S'), system.initial_condition_values)
        self.assertIn(sympy.Symbol('S'), system.initial_conditions)

    def test_both_faces_survive_a_reduction(self):
        system = self.system()
        reduced = ODETranslator.from_ode_system(system).translate(system)
        # The expression face is whichever invariant C's scale became; the value face is
        # still the zero that was stated.
        valued = reduced.initial_condition_values
        self.assertEqual(set(valued.values()), {sympy.Integer(0)})
        self.assertEqual(len(valued), 2)
        for variable in valued:
            with self.subTest(variable=variable):
                self.assertIn(reduced.initial_conditions[variable],
                              reduced.parameters)

    def test_value_of_evaluates_any_expression(self):
        system = self.system()
        system.update_parameter_values({'k_1': 2, 'k_2': 3})
        self.assertEqual(system.value_of('k_1 * k_2'), sympy.Integer(6))

    def test_value_of_returns_none_when_a_symbol_has_no_value(self):
        system = self.system()
        self.assertIsNone(system.value_of('k_1 * k_2'))

    def test_value_of_returns_none_when_a_variable_is_left_undetermined(self):
        system = self.system()
        self.assertIsNone(system.value_of('C_0 + C'))

    def test_a_known_zero_can_settle_an_expression_naming_a_variable(self):
        # C_0 is zero, so `C_0 * C` is zero whatever C does.  The question is whether the
        # values determine the expression, not whether every symbol in it has one.
        system = self.system()
        self.assertEqual(system.value_of('C_0 * C'), sympy.Integer(0))

    def test_value_of_a_known_zero_is_zero_not_none(self):
        # `None` means "not settled"; a settled zero must not be confused with it.
        system = self.system()
        self.assertEqual(system.value_of('C_0'), sympy.Integer(0))
        self.assertIsNotNone(system.value_of('C_0'))


class TestScalingSymmetries(TestCase):
    '''
    A row of the scaling matrix is a one-parameter group of scalings the equations
    admit.  `scaling_symmetries` writes the row out with its factor in place, which is
    the form a reader can check against the equations.
    '''

    MICHAELIS_MENTEN = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                        'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
    HL64 = ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)']

    def systems(self):
        from desr.chemical_reaction_network import ChemicalReactionNetwork
        return (('Michaelis-Menten', ODESystem.from_equations(self.MICHAELIS_MENTEN)),
                ('prey-predator', ODESystem.from_equations(
                    ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )',
                     'dp/dt = s*p*(1 - h*p / n)'])),
                ('enzyme chain', ChemicalReactionNetwork.from_diagram(
                    'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()),
                ('HL 6.4', ODESystem.from_equations(self.HL64)))

    def test_every_symmetry_leaves_the_equations_invariant(self):
        # The claim the report makes, checked against the equations rather than against
        # the matrix it came from.  A derivative dv/dt scales by the weight of v less
        # the weight of t, so substituting the symmetry into a right-hand side gives
        # that power of the factor times what was there.
        factor = sympy.Symbol('lamda')
        for label, system in self.systems():
            matrix = system.maximal_scaling_matrix()
            symbols = list(system.symbols)
            clock = symbols.index(system.indep_var)
            for row, symmetry in enumerate(system.scaling_symmetries()):
                for variable, derivative in system.derivative_dict.items():
                    if variable == system.indep_var:
                        continue
                    with self.subTest(system=label, symmetry=row, variable=str(variable)):
                        weight = matrix[row, symbols.index(variable)] - matrix[row, clock]
                        scaled = derivative.subs(symmetry, simultaneous=True)
                        self.assertEqual(0, sympy.simplify(
                            scaled - factor ** weight * derivative))

    def test_a_symbol_the_symmetry_fixes_is_left_out(self):
        # Hubert & Labahn 6.4's scaling does not act on t, which is the condition the
        # dependent-variable scheme is chosen by.  A symbol multiplied by 1 is noise in
        # a report of what moves.
        system = ODESystem.from_equations(self.HL64)
        symmetry, = system.scaling_symmetries()
        self.assertNotIn(system.indep_var, symmetry)
        self.assertEqual({sympy.Symbol('z1'), sympy.Symbol('z2')}, set(symmetry))

    def test_the_images_read_as_a_change_of_unit(self):
        # t lengthens by lambda and the rate constants shorten by it, which is what
        # changing the unit of time does and what [t] names that direction after.
        system = ODESystem.from_equations(self.MICHAELIS_MENTEN)
        time, concentration = system.scaling_symmetries()
        factor = sympy.Symbol('lamda')
        self.assertEqual(factor * sympy.Symbol('t'), time[sympy.Symbol('t')])
        self.assertEqual(sympy.Symbol('k_2') / factor, time[sympy.Symbol('k_2')])
        self.assertEqual(factor * sympy.Symbol('s'),
                         concentration[sympy.Symbol('s')])

    def test_the_factor_can_be_renamed(self):
        system = ODESystem.from_equations(self.MICHAELIS_MENTEN)
        symmetry, _ = system.scaling_symmetries(factor='mu')
        self.assertEqual(sympy.Symbol('mu') * sympy.Symbol('t'),
                         symmetry[sympy.Symbol('t')])

    def test_the_rows_line_up_with_the_dimension_names(self):
        # The i-th symmetry is the change of unit of the i-th dimension, which is what
        # lets a reader pair the two reports.
        system = ODESystem.from_equations(self.MICHAELIS_MENTEN)
        dimensions = system.dimensions_of()
        for symmetry, (name, sign) in zip(
                system.scaling_symmetries(),
                system._dimension_names(system.maximal_scaling_matrix())):
            with self.subTest(dimension=str(name)):
                # The symbol the direction is named after is one the symmetry moves.
                named = sympy.Symbol(str(name)[1:-1])
                self.assertIn(named, symmetry)
                self.assertEqual(name ** sign, dimensions[named])


class TestStructuralDimensions(TestCase):
    '''
    A symbol's column of the scaling matrix says how it transforms under each scaling --
    its weight -- and weights written multiplicatively are dimensions.

    These are structural dimensions, from the equations.  Physical dimensions are
    contained among them, since a dimensionally consistent system is invariant under its
    own changes of unit, but the containment is usually strict and nothing here
    distinguishes the two.  Declaring physical dimensions changes none of this.
    '''

    MICHAELIS_MENTEN = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                        'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
    ENZYME_CHAIN = 'x + e -> y + e\ny -> z\nz + f -> x + f'
    HL66 = ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
            'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)']
    HL64 = ['dz1/dt = z1*(1+z1*z2)', 'dz2/dt = z2*(1/t - z1*z2)']

    def systems(self):
        from desr.chemical_reaction_network import ChemicalReactionNetwork
        return (('Michaelis-Menten', ODESystem.from_equations(self.MICHAELIS_MENTEN)),
                ('prey-predator', ODESystem.from_equations(
                    ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )',
                     'dp/dt = s*p*(1 - h*p / n)'])),
                ('enzyme chain', ChemicalReactionNetwork.from_diagram(
                    self.ENZYME_CHAIN).to_ode_system()),
                ('HL 6.6', ODESystem.from_equations(self.HL66)),
                ('HL 6.4', ODESystem.from_equations(self.HL64)))

    @staticmethod
    def dimension_of(monomial, dimensions, system):
        '''The dimension of a monomial, from the dimensions of its symbols.'''
        powers = monomial.as_powers_dict()
        return sympy.simplify(sympy.prod(
            [dimensions[s] ** powers.get(s, 0) for s in system.symbols]))

    def test_every_invariant_is_dimensionless(self):
        # A value is a number relative to a scale, so only at dimension 1 is it
        # independent of which scale was chosen.
        from desr.ode_translation import scale_action
        for label, system in self.systems():
            dimensions = system.dimensions_of()
            translator = ODETranslator.from_ode_system(system)
            for monomial in scale_action(system.symbols, translator.herm_mult_n):
                with self.subTest(system=label, invariant=str(monomial)):
                    self.assertEqual(
                        1, self.dimension_of(monomial, dimensions, system))

    def test_the_auxiliaries_carry_one_dimension_each(self):
        # A * V_i is the identity, so the auxiliaries are a basis of dimension space:
        # r directions and r auxiliaries, one apiece and no two alike.
        from desr.ode_translation import scale_action
        for label, system in self.systems():
            dimensions = system.dimensions_of()
            translator = ODETranslator.from_ode_system(system)
            carried = [self.dimension_of(monomial, dimensions, system)
                       for monomial in scale_action(system.symbols,
                                                    translator.herm_mult_i)]
            with self.subTest(system=label):
                self.assertEqual(translator.r, len(carried))
                self.assertEqual(len(carried), len(set(carried)))
                self.assertNotIn(1, carried)

    def test_symbols_the_equations_tie_together_share_a_dimension(self):
        # Nobody said c, e_0 and s are concentrations; the equations did.
        system = ODESystem.from_equations(self.MICHAELIS_MENTEN)
        dimensions = system.dimensions_of()
        shared = {dimensions[sympy.Symbol(name)] for name in ('c', 'e_0', 's')}
        self.assertEqual(1, len(shared))

    def test_a_rate_constant_comes_out_an_inverse_time(self):
        # Read off the equations rather than asserted: k_2 is first order and k_1 is
        # second, and the dimensions say so without anyone naming a unit.
        system = ODESystem.from_equations(self.MICHAELIS_MENTEN)
        dimensions = system.dimensions_of()
        time = dimensions[sympy.Symbol('t')]
        self.assertEqual(1 / time, dimensions[sympy.Symbol('k_2')])
        self.assertEqual(1 / (time * dimensions[sympy.Symbol('e_0')]),
                         dimensions[sympy.Symbol('k_1')])

    def test_a_direction_is_named_after_a_symbol_that_has_it(self):
        # Michaelis-Menten's first direction holds t at +1 and k_2 at -1.  Naming it
        # after the rate constant would be right and backwards to read, so +1 wins.
        system = ODESystem.from_equations(self.MICHAELIS_MENTEN)
        dimensions = system.dimensions_of()
        self.assertEqual(sympy.Symbol('[t]'), dimensions[sympy.Symbol('t')])

    def test_a_direction_with_nothing_at_plus_one_falls_back_to_minus_one(self):
        # Hubert & Labahn 6.6 has no symbol of weight +1 in its single direction, but
        # z1 sits there at -1.  Both generate the same lattice direction, so naming by
        # the other one is presentation and changes no matrix.
        system = ODESystem.from_equations(self.HL66)
        dimensions = system.dimensions_of()
        self.assertEqual(sympy.Symbol('[z1]'), dimensions[sympy.Symbol('z1')])
        self.assertEqual(sympy.Symbol('[z1]') ** -3, dimensions[sympy.Symbol('t')])

    def test_the_clock_can_be_dimensionless(self):
        # Hubert & Labahn 6.4's scaling does not act on t, which is the condition the
        # dependent-variable scheme is chosen by.  So t has dimension 1.
        system = ODESystem.from_equations(self.HL64)
        self.assertEqual(1, system.dimensions_of()[sympy.Symbol('t')])
        self.assertEqual(ODETranslator.from_ode_system(system).scheme(system),
                         ReductionScheme.DEP_VAR)

    def test_naming_physical_dimensions_would_change_none_of_this(self):
        # The enzyme chain has four structural dimensions.  A modeller would call three
        # of them concentration and one time, which says something about those
        # directions without changing them.
        from desr.chemical_reaction_network import ChemicalReactionNetwork
        system = ChemicalReactionNetwork.from_diagram(
            self.ENZYME_CHAIN).to_ode_system()
        dimensions = system.dimensions_of()
        species = [dimensions[sympy.Symbol(name)] for name in ('e', 'f', 'x')]
        # All concentrations physically, and three different structural dimensions.
        self.assertEqual(3, len(set(species)))

    def test_the_translator_reads_dimensions_in_its_own_scaling(self):
        system = ODESystem.from_equations(self.MICHAELIS_MENTEN)
        translator = ODETranslator.from_ode_system(system)
        self.assertEqual(system.dimensions_of(), translator.dimensions_of())

    def test_an_explicit_scaling_matrix_is_honoured(self):
        # A translator may hold a scaling other than the maximal one, and then the
        # dimensions are that one's.  Here the two agree, which is the point: the
        # argument is plumbed, not ignored.
        system = ODESystem.from_equations(self.MICHAELIS_MENTEN)
        matrix = system.maximal_scaling_matrix()
        self.assertEqual(system.dimensions_of(),
                         system.dimensions_of(scaling_matrix=matrix))


class TestPhysicalDimensions(TestCase):
    '''
    Physical dimensions are an overlay on the structural ones.  A dimensionally
    consistent system is invariant under its own changes of unit, so physical directions
    sit among the structural ones, usually strictly.  The gap is reduction that no choice
    of unit accounts for.
    '''

    MICHAELIS_MENTEN = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                        'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
    ENZYME_CHAIN = 'x + e -> y + e\ny -> z\nz + f -> x + f'

    def michaelis_menten(self):
        return ODESystem.from_equations(self.MICHAELIS_MENTEN)

    def enzyme(self):
        from desr.chemical_reaction_network import ChemicalReactionNetwork
        return ChemicalReactionNetwork.from_diagram(self.ENZYME_CHAIN).to_ode_system()

    def test_a_few_declarations_determine_the_rest(self):
        # Two declarations out of seven symbols.  The structural weights fix the
        # relations between symbols, so the others follow.
        system = self.michaelis_menten()
        system.set_dimensions({'t': 'time', 'e_0': 'concentration'})
        physical = system.physical_dimensions_of()

        time, concentration = sympy.Symbol('time'), sympy.Symbol('concentration')
        self.assertEqual(time, physical[sympy.Symbol('t')])
        self.assertEqual(concentration, physical[sympy.Symbol('s')])
        self.assertEqual(1 / time, physical[sympy.Symbol('k_2')])
        self.assertEqual(1 / (concentration * time), physical[sympy.Symbol('k_1')])

    def test_declaring_one_dimension_merges_structural_directions(self):
        # The enzyme chain gives e, f and x three different structural dimensions.
        # Declaring all three to be concentration identifies those directions, and
        # symbols that were never declared follow.
        system = self.enzyme()
        structural = system.dimensions_of()
        species = [structural[sympy.Symbol(name)] for name in ('e', 'f', 'x')]
        self.assertEqual(3, len(set(species)))

        system.set_dimensions({'t': 'time', 'e': 'concentration',
                               'f': 'concentration', 'x': 'concentration'})
        physical = system.physical_dimensions_of()
        self.assertEqual(
            1, len({physical[sympy.Symbol(name)] for name in ('e', 'f', 'x')}))
        self.assertEqual(physical[sympy.Symbol('e')], physical[sympy.Symbol('y')])
        self.assertEqual(physical[sympy.Symbol('e')], physical[sympy.Symbol('z')])

    def test_every_invariant_is_physically_dimensionless_too(self):
        # A change of unit is one of the scalings an invariant is invariant under, so an
        # invariant carries no dimension.  This is also the condition under which the
        # declarations are consistent.
        from desr.ode_translation import scale_action
        system = self.enzyme()
        system.set_dimensions({'t': 'time', 'e': 'concentration',
                               'f': 'concentration', 'x': 'concentration'})
        physical = system.physical_dimensions_of()
        translator = ODETranslator.from_ode_system(system)
        for monomial in scale_action(system.symbols, translator.herm_mult_n):
            with self.subTest(invariant=str(monomial)):
                powers = monomial.as_powers_dict()
                self.assertEqual(1, sympy.simplify(sympy.prod(
                    [physical[s] ** powers.get(s, 0) for s in system.symbols])))

    def test_a_declaration_the_equations_contradict_is_refused(self):
        # k_2 multiplies a concentration to give a rate of change of concentration, so
        # it is an inverse time.  Declaring it a concentration has no solution.
        system = self.michaelis_menten()
        system.set_dimensions({'t': 'time', 'k_2': 'concentration'})
        with self.assertRaises(ValueError) as caught:
            system.physical_dimensions_of()
        self.assertIn('dimensionally consistent', str(caught.exception))

    def test_a_sum_is_not_a_dimension(self):
        system = self.michaelis_menten()
        with self.assertRaises(ValueError) as caught:
            system.set_dimension('t', 'time + length')
        self.assertIn('not a dimension', str(caught.exception))

    def test_a_symbol_of_another_system_is_refused(self):
        system = self.michaelis_menten()
        with self.assertRaises(ValueError) as caught:
            system.set_dimension('nonsense', 'time')
        self.assertIn('not a symbol of this system', str(caught.exception))

    def test_a_symbol_the_declarations_do_not_reach_has_no_dimension(self):
        # Declaring t and x reaches two of the enzyme chain's four directions.  e is
        # left unsaid, and k_0_1 multiplies e, so neither follows.  Reporting 1 for
        # those would read as "dimensionless", which nobody said.
        system = self.enzyme()
        system.set_dimensions({'t': 'time', 'x': 'concentration'})
        physical = system.physical_dimensions_of()

        self.assertEqual(sympy.Symbol('concentration'), physical[sympy.Symbol('x')])
        self.assertEqual(sympy.Symbol('time'), physical[sympy.Symbol('t')])
        for name in ('e', 'f', 'k_0_1', 'k_4_5'):
            with self.subTest(symbol=name):
                self.assertIsNone(physical[sympy.Symbol(name)])

    def test_declaring_the_rest_turns_none_into_a_dimension(self):
        # The same system, said in full.  Nothing but a declaration can fill the gap.
        system = self.enzyme()
        system.set_dimensions({'t': 'time', 'x': 'concentration'})
        self.assertIsNone(system.physical_dimensions_of()[sympy.Symbol('e')])

        system.set_dimensions({'e': 'concentration', 'f': 'concentration'})
        physical = system.physical_dimensions_of()
        concentration, time = sympy.Symbol('concentration'), sympy.Symbol('time')
        self.assertEqual(concentration, physical[sympy.Symbol('e')])
        # k_0_1 is second order, k_2_3 first, read off the equations as always.
        self.assertEqual(1 / (concentration * time), physical[sympy.Symbol('k_0_1')])
        self.assertEqual(1 / time, physical[sympy.Symbol('k_2_3')])

    def test_nothing_declared_is_not_a_claim_that_there_are_none(self):
        system = self.michaelis_menten()
        self.assertEqual({}, system.physical_dimensions)
        self.assertEqual({}, system.physical_dimensions_of())

    def test_naming_changes_nothing(self):
        # Physical dimensions are an overlay.  Whatever is declared, the scaling matrix,
        # the multiplier, its inverse, the structural dimensions and the reduced system
        # are the same.
        plain, named = self.enzyme(), self.enzyme()
        named.set_dimensions({'t': 'time', 'e': 'concentration',
                              'f': 'concentration', 'x': 'concentration'})

        self.assertEqual(plain.maximal_scaling_matrix(),
                         named.maximal_scaling_matrix())
        self.assertEqual(plain.dimensions_of(), named.dimensions_of())

        before = ODETranslator.from_ode_system(plain)
        after = ODETranslator.from_ode_system(named)
        self.assertEqual(before.herm_mult, after.herm_mult)
        self.assertEqual(before.inv_herm_mult, after.inv_herm_mult)
        self.assertEqual(before.r, after.r)
        self.assertEqual(before.translate(), after.translate())
