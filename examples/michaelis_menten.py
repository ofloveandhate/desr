"""Michaelis-Menten Equations

This example script walks through an analysis of the Michaelis Menten equations.

The original source for the equations is
L. Michaelis and M. Menten, Die kinetik der Invertinwirkung, Biochem Z., 49 (1913), pp. 333–369.

and for a modern source, see

J. Gunawardena, Time-scale separation–Michaelis and Menten’s old idea, still bearing
fruit., The FEBS journal, 281 (2014), pp. 473–88, https://doi.org/10.1111/febs.12532,
http://www.ncbi.nlm.nih.gov/pubmed/24103070http://www.pubmedcentral.nih.gov/
articlerender.fcgi?artid=PMC3991559.
"""

import sympy
from desr import ODESystem, ODETranslator, scale_action
from desr.matrix_normal_forms import smith_normal_form
from desr.tex_tools import expr_to_tex




def michaelis_menten_four_variables():
    """
    Example: Michaelis-Menten kinetics in four variables with initial conditions

    We start with ODEs.
    """

    print('\n\n----------------------\nmichaelis_menten_four_variables\n-----------------------\n\n')


    def system_four_variable_michaelis_menten():


        # make string to hold the system in latex form
        # note that, so make valid variable names,
        # minus signs  `-` in subscripts are changed to the character `m`
        system_tex = \
              '''\frac{dS}{dt} &= - k_1 E S + k_{-1} C \\\\
                 \frac{dE}{dt} &= - k_1 E S + k_{-1} C + k_2 C \\\\
                 \frac{dC}{dt} &=   k_1 E S - k_{-1} C - k_2 C \\\\
                 \frac{dP}{dt} &=   k_2 C'''

        system = ODESystem.from_tex(system_tex)

        return system


    # a 'lil helper function that moves values to the front.  ⚠️ assumes that `sorted` is a stable sort...
    put_at_front = lambda v, order: sorted(order, key=lambda n: str(n)!=v)



    def four_variable_with_ics():
        print('\n\nFour-variable with ICs\n-----------------------------------')

        system = system_four_variable_michaelis_menten()
        system.update_initial_conditions({'S': 'S_0', 'E':'E_0', 'C':'C_0', 'P':'P_0'})

        variable_order = list(system.symbols)



        # put s_0 and k_1 at the front so will likely get normalized out
        variable_order = put_at_front('k_1', variable_order)
        variable_order = put_at_front('k_2', variable_order)
        variable_order = put_at_front('k_m1', variable_order)
        variable_order = put_at_front('S_0', variable_order)
        print(f'{variable_order=}')

        # variable_order[-1], variable_order[-4] = variable_order[-4], variable_order[-1]
        system.reorder_variables(variable_order)
        print('system\n',system)


        translator = ODETranslator.from_ode_system(system)


        print('variable order:', translator.symbols_domain)

        # Print scaling matrices
        print('Scaling matrix:')
        print(translator.scaling_matrix.__repr__())

        # Print invariants
        print('Invariants: ', translator.invariant_monomials())

        # Print translated result
        print('Reduced system:')
        print(translator.translate(system))

        print('Substitutions:')
        print(translator.translate_parameter_substitutions(system=system))

    def four_variable_using_specific_row_ops():
        print('\n\nFour-variable without ICs:\nWhat happens when we make k_2 the last variable?\n-----------------------------------')

        # make the system
        system = system_four_variable_michaelis_menten()


        variable_order = list(system.symbols)
        put_at_front('k_2',variable_order)
        system.reorder_variables(variable_order)
        print('New variable order:',variable_order)


        system_from_ode = ODETranslator.from_ode_system(system)
        print('Invariants:', ', '.join(map(str, system_from_ode.invariant_monomials())))

        # Now do the reduction
        reduced_system = system_from_ode.translate(system)
        print('Reduced system:', reduced_system)


        print('Imposing a choice of invariants')
        print('--------------------------------')
        # Extend a choice of invariants   t  C  E  P  S k_1 k_2 k_{-1}
        invariant_choice = sympy.Matrix([[0, 1, 0, 0, 0, 1, -1, 0],
                                         [0, 0, 0, 1, 0, 1, 0, -1]]).T
        print('P =')
        print(invariant_choice.__repr__())
        print('Chosen invariants:')
        print(scale_action(system_from_ode.symbols_domain, invariant_choice))

        ## Method that does the extension automatically.
        max_scal2 = system_from_ode.extend_from_invariants(invariant_choice=invariant_choice)

        ## Stepping through the above function for the sake of the paper.
        ## Step 1: Check we have invariants
        choice_actions = system_from_ode.scaling_matrix * invariant_choice
        assert choice_actions.is_zero_matrix  # Else we have to stop.

        ## Step 2: Try and extend the choices by a basis of invariants
        ## Step 2a: Extend (W_d . invariant_choice) to a unimodular matrix
        WdP = system_from_ode.inv_herm_mult_d * invariant_choice
        normal_form, row_ops, col_ops = smith_normal_form(WdP)

        print('Smith normal form decomposition:')
        print('{}\n{}\n{}\n{}\n=\n{}'.format(*map(lambda x: x.__repr__(), (row_ops,
                                                                           system_from_ode.inv_herm_mult_d,
                                                                           invariant_choice,
                                                                           col_ops,
                                                                           smith_normal_form)))
             )

        print('U^{-1} = ',row_ops.inv())
        print(f'{WdP=}')


        col_ops = sympy.Matrix.hstack(WdP, row_ops.inv()[:, 2:])  # col_ops is the column operations we're going to apply to Vn
        print(f'{col_ops=}')


        print('New Vn\n',system_from_ode.herm_mult_n * col_ops)

        print('New invariants:\n',', '.join(map(str, scale_action(system_from_ode.symbols_domain, system_from_ode.herm_mult_n * col_ops))))

        hermite_mult = system_from_ode.herm_mult_n * col_ops
        # The permutation is (0 1 3 4 2 5)
        hermite_mult.col_swap(0, 1)
        hermite_mult.col_swap(0, 3)
        hermite_mult.col_swap(0, 4)
        hermite_mult.col_swap(0, 2)
        hermite_mult.col_swap(0, 5)

        print('Permuted Vn:\n',hermite_mult)


        # We need to add on Vi - the original will do.
        hermite_mult = sympy.Matrix.hstack(system_from_ode.herm_mult_i, hermite_mult)
        translator = ODETranslator(system_from_ode.scaling_matrix, hermite_multiplier=hermite_mult)

        print('Reduced system:')
        print(translator.translate(system))


    four_variable_with_ics()
    # four_variable_using_specific_row_ops()



def michaelis_menten_two_variable(verbose=True):
    """
    Example `sec:michaelis_menten_comparison_with_classical`
    """

    print('\n\n----------------------\nmichaelis_menten_two_variable\n-----------------------\n\n')

    # Enable pretty printing
    sympy.init_printing(pretty_print=True, use_latex=True)

    def system_tex_two_var():
        return '''\frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\\\
             \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c'''


    def with_initial_conditions_and_a_constant():
        '''
        this code reproduces the michaelis-menten example in the literature.
        '''
        system_tex = system_tex_two_var()


        system = ODESystem.from_tex(system_tex)
        system.add_constraint('Km', '(k_m1 + k_2) / k_1')




        print('\n\nAdding in the initial conditions')
        print('-------------------------------------')
        # The complex starts at nothing.  Saying so used to be impossible, and the
        # literature just glosses over it; `c(0) = 0` records the value against a symbol
        # `c_0`, which is what carries the units -- a zero concentration is still a
        # concentration, and the reduction stays entirely symbolic.
        system.update_initial_conditions(['s(0) = s0', 'c(0) = 0'])
        print('initial conditions:', system.initial_conditions)
        print('as values         :', system.initial_condition_values)

        # `c_0` goes towards the back, which is the end that is *kept*.  A zero must not be
        # scaled by, and the front of the order is the scale-by end.
        system.reorder_variables(['s0','e_0', 'k_1', 'k_2', 'k_m1', 'Km', 'c_0', 't', 'c', 's'])


        print('System:', system)
        print('Variable order: ', system.symbols)

        print('Exponent Matrix:', system.exponent_matrix().__repr__())
        # The dependents are in the order (c, s), so Murray's names go (v, u): his `u` is
        # the substrate and his `v` is the complex.
        #
        # `kappa` rather than `c` for the parameters.  They are substituted to Murray's
        # names below, so the letter is arbitrary -- but `c` is a dependent variable here,
        # whose initial value is the minted `c_0`, and printing `c_0` beside `c0, c1, ...`
        # asks a reader to tell apart two things that look the same.  They are different
        # symbols, so nothing breaks; it just cannot be read.
        translator = ODETranslator.from_ode_system(system, naming_scheme={'parameters': 'kappa', 'independent': 'tau', 'dependents': ['v','u']})
        print('symbols this reduction divides by, so none of them may be zero:',
              translator.symbols_that_must_not_be_zero())


        print('Scaling matrix:',translator.scaling_matrix.__repr__())

        print('Hermite multiplier', translator.herm_mult.__repr__())
        print('Invariants: ', translator.invariant_monomials())


        print('Reduced system:')
        print(translator.translate(system))#.to_tex()


        print('Changing invariants by column operations on the Hermite multiplier')

        # Two groups still differ from Murray's: his clock is k_1 e_0 t rather than
        # k_1 s0 t, and he measures the complex in units of the enzyme, c/e_0, rather than
        # of the substrate.  Both are a factor of epsilon = e_0/s0 away, so both are one
        # column operation.
        #
        # The column numbers come from `column_of`, so they are derived from the invariant
        # each one carries rather than counted off by hand.  They can still be printed and
        # quoted -- and they move when the variable order moves, which is exactly why they
        # should not be written as literals.
        clock = translator.column_of('k_1*s0*t')
        epsilon = translator.column_of('e_0/s0')
        complexed = translator.column_of('c/s0')
        print(f'{clock=}, {epsilon=}, {complexed=}')

        translator.multiplier_add_columns(clock, epsilon)              # tau -> k_1 e_0 t
        translator.multiplier_add_columns(complexed, epsilon, power=-1)  # v -> c/e_0

        print('Invariants: ',
            ['{}={}'.format(a,b) for a,b in zip(translator.invariant_names(system), translator.invariant_monomials())])
        print('Hermite multiplier', translator.herm_mult.__repr__())


        # print('\nMichaelis-Menten Reparametrisation 2')
        # print('Changing invariants by column operations on the Hermite multiplier')
        # # Divide time through by epsilon
        # translator.multiplier_add_columns(i=0, j=8, power=1)

        # print('Hermite multiplier', translator.herm_mult.__repr__())
        # print('Invariants: ', translator.invariant_monomials())

        print('W:', translator.inv_herm_mult.__repr__())
        # Print translated system

        print('substitutions:',translator.translate_parameter_substitutions(system))

        reduced_system = translator.translate(system)

        print('Reduced system, in arbitrary parameter names:')
        print(reduced_system)#.to_tex()

        print('the initial conditions, as values:',
              reduced_system.initial_condition_values)

        # Read off the monomials rather than guessed at, because the order of the kappa_i
        # follows the variable order and so is not something to remember:
        #
        #   kappa0 = e_0/s0        = epsilon
        #   kappa1 = k_2/(k_1 s0)  = lambda, which this script calls L
        #   kappa2 = k_m1/(k_1 s0) = K - L
        #   kappa3 = Km/s0         = K
        #   kappa4 = c_0/s0        = the initial complex, in substrate units
        #
        # Murray has no name for the last one because he sets it to zero, which is exactly
        # what we have recorded.  `gamma` until then.
        #
        # The check is the constraint.  `kappa3 == kappa1 + kappa2` has to become
        # `K == L + (K - L)`, a tautology, and so disappear.  If a constraint survives the
        # substitution, the mapping is wrong -- which is how this was caught: it used to
        # leave the false `epsilon == K` behind.
        to_sub = dict(zip(translator.invariant_names(system),
                          ['epsilon', 'L', 'K - L', 'K', 'gamma']))
        print('substituting Murray\'s names:', to_sub)
        reduced_system = reduced_system.subs(to_sub, subs_constraints=True)
        assert not reduced_system.constraints, (
            'a constraint survived the renaming, so the names do not match the monomials: '
            '{}'.format(reduced_system.constraints))

        print('after subbing into canonical names from Murray')
        print(reduced_system)#.to_tex()

        # `v(0) = gamma/epsilon`, and gamma is zero, so v starts at nothing however small
        # epsilon is.  Both faces of that are available: the expression says *what* the
        # initial value is, and the value says what it *is*.
        print('the initial conditions, as values:',
              reduced_system.initial_condition_values)

        # Murray writes three groups, not four, because he sets the initial complex to
        # zero.  Substituting the value the system already holds does that, and the group
        # goes with it -- derived from the stated initial condition rather than asserted.
        # Substituting anything else here is refused, since it would contradict `c(0) = 0`.
        print('\nSubstituting the initial complex, which is known to be zero')
        print('-----------------------------------------------------------')
        reduced_system = reduced_system.subs({'gamma': 0})
        print(reduced_system)
        assert set(map(str, reduced_system.parameters)) == {'K', 'L', 'epsilon'}, (
            'expected Murray (6.20) three groups, got {}'.format(
                reduced_system.parameters))

        return




    with_initial_conditions_and_a_constant()
    # another_way()

if __name__ == '__main__':
    # Enable pretty printing
    sympy.init_printing(pretty_print=True, use_latex=True)

    michaelis_menten_two_variable()

    # michaelis_menten_four_variables()
