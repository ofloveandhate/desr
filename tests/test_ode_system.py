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
                                    variables_domain=system.variables,
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
        named = {str(v): v for v in max_scal.variables_domain}
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
            variables = list(max_scal.variables_domain)
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
        self.assertEqual(variables, list(map(str, original_system.variables)))
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
        # s_0 is appended, which leaves the system out of canonical order: it is a parameter
        # and the parameters come first.  Putting it back is the caller's job for now; see
        # the note about update_initial_conditions in the task list.
        self.assertEqual(variables + ['s_0'], list(map(str, original_system.variables)))
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
