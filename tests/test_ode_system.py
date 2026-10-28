from unittest import TestCase

import sympy

from desr.matrix_normal_forms import hnf_row
from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator


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
        # Match the maximal scaling matrix
        max_scal = system.maximal_scaling_matrix()
        self.assertTrue(max_scal == sympy.Matrix([[0, 1, -1]]))

        # Give Hermite multiplier from the paper. Padd it with a row and column for t to work with current infrastructure
        hermite_multiplier_example = sympy.Matrix([[0, 1, 0],
                                                   [1, 0, 1],
                                                   [0, 0, 1]])

        translator = ODETranslator(max_scal, hermite_multiplier=hermite_multiplier_example)

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
        answer = ODESystem.from_equations('dx0/dt = x0*(y0 - 1/t)\ndy0/dt = y0*(1 + 1/t)')
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

        ## Check our answer hasn't changed, using our own Hermite multiplier
        translator = ODETranslator.from_ode_system(system)
        translated = translator.translate(system)
        answer = ODESystem.from_equations('dx0/dt = x0*(2*y1/3 + 2/3 + y1/y0)/t\ndy0/dt = y0*(y1 - 1)/t\ndy1/dt = y1*(5*y1/3 - 10/3 + (2*y0*(-y1 + 5)/3 + y1)/y0)/t')
        self.assertEqual(translated, answer)


    def test_verhulst_log_growth(self):
        ''' Verhult logistic growth model from :cite:`Hubert2013c`
            dn/dt = r*n*(1-n/k)
        '''
        equations = '''dn/dt = r*n*(1-n/k)'''

        system = ODESystem.from_equations(equations)
        system.reorder_variables('rktn')

        ## Check against the answer in the paper
        # Match the maximal scaling matrix
        max_scal = system.maximal_scaling_matrix()
        # Compare to the paper by performing one row operation
        max_scal[0, :] *= -1
        self.assertEqual(max_scal, sympy.Matrix([[-1, 0, 1, 0], [0, 1, 0, 1]]))

        herm_mult_paper = sympy.Matrix([[-1, 0, 1, 0],
                                       [0, 1, 0, -1],
                                       [0, 0, 1, 0],
                                       [0, 0, 0, 1]])
        translator = ODETranslator(max_scal,
                                    variables_domain=system.variables,
                                    hermite_multiplier=herm_mult_paper)
        translated = translator.translate(system)
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

        ## Now do it with our own Hermite multiplier, except we need the variables in a different order
        system.reorder_variables('tnrk')
        translator = ODETranslator.from_ode_system(system)
        translated = translator.translate_general(system)
        self.assertEqual(translated, answer)
        general_soln = translator.reverse_translate_general(reduced_solutions, system_indep_var_index=0)
        saved_soln = tuple(map(sympy.sympify, ['c2/(c*exp(-t/c1) + 1)', '1/c1', 'c2']))  # Just a cached value
        self.assertTupleEqual(general_soln, saved_soln)

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

        t, n, p, K, d, h, k, r, s = max_scal.variables_domain

        self.assertTupleEqual(tuple(max_scal.invariants()), (s*t, n / d, k*p/(d*s), K/d, h*s/k, r/s))

        # Choose a new invariant
        invariant_choice = sympy.Matrix([[1, 0, 0, 0, 0, 0, 0, 1, 0],  # t * r
                                        [0, 0, 1, 0, -1, 1, 0, 0, 0],  # p * h /d
                                        ]).T  # Note the transpose! Each column expresses an invariant

        max_scal2 = max_scal.extend_from_invariants(invariant_choice=invariant_choice)
        self.assertTupleEqual(tuple(max_scal2.invariants()), (r*t, h*p/d, s*t, n/d, k*p/(d*s), K/d))

        # This should work even if we move the time about
        invariant_choice = sympy.Matrix([[0, 0, 1, 0, -1, 1, 0, 0, 0],  # p * h /d
                                        ]).T  # Note the transpose! Each column expresses an invariant
        max_scal3 = max_scal.extend_from_invariants(invariant_choice=invariant_choice)
        self.assertTupleEqual(tuple(max_scal3.invariants()), (h*p/d, s*t, n/d, k*p/(d*s), K/d, r/s))


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

        variables = ['t', 's', 'c', 'K', 'k_2', 'k_1', 'e_0']
        original_system.reorder_variables(variables)
        self.assertEqual(variables, list(map(str, original_system.variables)))
        self.assertEqual(original_system.exponent_matrix(), sympy.Matrix([[ 1, 1,  1, 1, 1, 1,  1],
                                                                       [-1, 0, -1, 0, 0, 1,  1],
                                                                       [ 1, 0,  1, 1, 0, 0, -1],
                                                                       [ 0, 0,  1, 0, 1, 0,  0],
                                                                       [ 1, 0,  0, 0, 0, 0,  0],
                                                                       [ 0, 1,  0, 1, 0, 1,  1],
                                                                       [ 0, 1,  0, 0, 0, 0,  1]]))
        max_scal1 = ODETranslator.from_ode_system(original_system)
        self.assertEqual(max_scal1.scaling_matrix, sympy.Matrix([[1, 0, 0, -1, -1, -1, 0],
                                                                 [0, 1, 1,  0,  0, -1, 1]]))

        # Add in the initial condition s_0
        original_system.update_initial_conditions({'s': 's_0'})
        self.assertSetEqual(set(original_system.initial_conditions.items()),
                            set(map(lambda t: (sympy.sympify(t[0]), sympy.sympify(t[1])),
                                    (('s', 's_0'),))))
        self.assertEqual(variables + ['s_0'], list(map(str, original_system.variables)))
        self.assertEqual(original_system.exponent_matrix(), sympy.Matrix([[ 1, 1,  1, 1, 1, 1,  1,  0],
                                                                       [-1, 0, -1, 0, 0, 1,  1,  1],
                                                                       [ 1, 0,  1, 1, 0, 0, -1,  0],
                                                                       [ 0, 0,  1, 0, 1, 0,  0,  0],
                                                                       [ 1, 0,  0, 0, 0, 0,  0,  0],
                                                                       [ 0, 1,  0, 1, 0, 1,  1,  0],
                                                                       [ 0, 1,  0, 0, 0, 0,  1,  0],
                                                                       [ 0, 0,  0, 0, 0, 0,  0, -1]]))

        max_scal1 = ODETranslator.from_ode_system(original_system)
        self.assertEqual(max_scal1.scaling_matrix, sympy.Matrix([[1, 0, 0, -1, -1, -1, 0, 0],
                                                                 [0, 1, 1,  0,  0, -1, 1, 1]]))

        reduced_system = max_scal1.translate(original_system)
        self.assertSetEqual(set(reduced_system.derivative_dict.items()),
                            set(map(lambda t: (sympy.sympify(t[0]), sympy.sympify(t[1])),
                                    (('nu0', 'kappa0*nu1 - kappa1*nu1 - kappa2*nu0 + nu0*nu1'),
                                     ('nu1', '-kappa0*nu1 + kappa2*nu0 - nu0*nu1'),
                                     ('tau', 1)))))
        self.assertSetEqual(set(reduced_system.initial_conditions.items()),
                            set(map(lambda t: (sympy.sympify(t[0]), sympy.sympify(t[1])),
                                    (('nu0', 1),))))
