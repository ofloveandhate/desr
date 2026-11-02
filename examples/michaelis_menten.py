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




        print('\n\nAdding in the initial condition for s')
        print('-------------------------------------')
        system.update_initial_conditions({'s': 's0'})

        system.reorder_variables(['s0','e_0', 'k_1', 'k_2', 'k_m1', 'Km',  't', 'c', 's'])


        print('System:', system)
        print('Variable order: ', system.symbols)

        print('Exponent Matrix:', system.exponent_matrix().__repr__())
        translator = ODETranslator.from_ode_system(system, naming_scheme={'parameters': 'c', 'independent': 'tau', 'dependents': ['u','v']})


        print('Scaling matrix:',translator.scaling_matrix.__repr__())

        print('Hermite multiplier', translator.herm_mult.__repr__())
        print('Invariants: ', translator.invariant_monomials())


        print('Reduced system:')
        print(translator.translate(system))#.to_tex()


        # print('\nMichaelis-Menten Reparametrisation 1')
        print('Changing invariants by column operations on the Hermite multiplier')
        # Python start at 0, but  the paper starts at 1, so these are down by 1 from the paper
        # translator.multiplier_add_columns(i=2, j=8, power=1)
        # translator.multiplier_add_columns(i=4, j=8, power=-1)

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

        to_sub =  {'c0': 'K', 'c1': 'K-L', 'c2':'L', 'c3':'epsilon'}
        reduced_system = reduced_system.subs(to_sub, subs_constraints=True)

        print('after subbing into canonical names from Murray')
        print(reduced_system)#.to_tex()

        return




    with_initial_conditions_and_a_constant()
    # another_way()

if __name__ == '__main__':
    # Enable pretty printing
    sympy.init_printing(pretty_print=True, use_latex=True)

    michaelis_menten_two_variable()

    # michaelis_menten_four_variables()
