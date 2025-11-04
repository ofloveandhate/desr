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
from desr.matrix_normal_forms import smf
from desr.ode_system import ODESystem
from desr.ode_system import maximal_scaling_matrix, rational_expr_to_power_matrix, hnf_col, hnf_row, normal_hnf_col
from desr.matrix_normal_forms import normal_hnf_row
from desr.ode_translation import ODETranslation, scale_action
from desr.tex_tools import expr_to_tex, matrix_to_tex




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

    
    # a 'lil helper function that moves values to the end.  assumes that `sorted` is a stable sort...
    put_at_end = lambda v, order: sorted(order, key=lambda n: str(n)==v)



    def four_variable_with_ics():
        print('\n\nFour-variable with ICs\n-----------------------------------')

        system = system_four_variable_michaelis_menten()
        system.update_initial_conditions({'S': 'S_0', 'E':'E_0', 'C':'C_0', 'P':'P_0'})

        variable_order = list(system.variables)

        

        # put s_0 and k_1 at the end so will likely get normalized out
        variable_order = put_at_end('k_1', variable_order)
        variable_order = put_at_end('k_2', variable_order)
        variable_order = put_at_end('k_m1', variable_order)
        variable_order = put_at_end('S_0', variable_order)
        print(f'{variable_order=}')

        # variable_order[-1], variable_order[-4] = variable_order[-4], variable_order[-1]
        system.reorder_variables(variable_order)
        print('system\n',system)


        translation = ODETranslation.from_ode_system(system)


        print('variable order:', translation.variables_domain)

        # Print scaling matrices
        print('Scaling matrix:')
        print(translation.scaling_matrix.__repr__())

        # Print invariants
        print('Invariants: ', translation.invariants())

        # Print translated result
        print('Reduced system:')
        print(translation.translate(system))

        print('Substitutions:')
        print(translation.translate_parameter_substitutions(system=system))

    def four_variable_using_specific_row_ops():
        print('\n\nFour-variable without ICs:\nWhat happens when we make k_2 the last variable?\n-----------------------------------')

        # make the system 
        system = system_four_variable_michaelis_menten()


        variable_order = list(system.variables)
        put_at_end('k_2',variable_order)
        system.reorder_variables(variable_order)
        print('New variable order:',variable_order)
        
        
        system_from_ode = ODETranslation.from_ode_system(system)
        print('Invariants:', ', '.join(map(str, system_from_ode.invariants())))
        
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
        print(scale_action(system_from_ode.variables_domain, invariant_choice))

        ## Method that does the extension automatically.
        max_scal2 = system_from_ode.extend_from_invariants(invariant_choice=invariant_choice)

        ## Stepping through the above function for the sake of the paper.
        ## Step 1: Check we have invariants
        choice_actions = system_from_ode.scaling_matrix * invariant_choice
        assert choice_actions.is_zero_matrix  # Else we have to stop.

        ## Step 2: Try and extend the choices by a basis of invariants
        ## Step 2a: Extend (W_d . invariant_choice) to a unimodular matrix
        WdP = system_from_ode.inv_herm_mult_d * invariant_choice
        smith_normal_form, row_ops, col_ops = smf(WdP)

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

        print('New invariants:\n',', '.join(map(str, scale_action(system_from_ode.variables_domain, system_from_ode.herm_mult_n * col_ops))))

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
        translation = ODETranslation(system_from_ode.scaling_matrix, hermite_multiplier=hermite_mult)

        print('Reduced system:')
        print(translation.translate(system))


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

    def with_intitial_conditions_and_a_constant_K_not_Km():
        system_tex = system_tex_two_var()
        system_tex_reduced_km1 = system_tex.replace('k_{-1}', '(K - k_2)')
        
        reduced_system_km1 = ODESystem.from_tex(system_tex_reduced_km1)
        
        reduced_system_km1.reorder_variables(['t', 's', 'c', 'K', 'k_2', 'k_1', 'e_0'])
        
        # Print variable order
        print('Variable order: ', reduced_system_km1.variables)
        
        print('Power Matrix:', reduced_system_km1.power_matrix().__repr__())
        
        # Print scaling matrices
        max_scal1 = ODETranslation.from_ode_system(reduced_system_km1)
        print('Scaling matrix:')
        print(max_scal1.scaling_matrix.__repr__())
        
        # Print invariants
        print('Invariants: ', max_scal1.invariants())
        # print(',\quad '.join(map(expr_to_tex, max_scal1.invariants())))
        
        # Print translated system
        print('Reduced system:')
        # print(max_scal1.translate(original_system))#.to_tex()
        
        print('Adding in the initial condition for s')
        print('-------------------------------------')
        reduced_system_km1.update_initial_conditions({'s': 's_0'})
        max_scal2 = ODETranslation.from_ode_system(reduced_system_km1)
        
        print('Invariants: ', max_scal2.invariants())
        print('Hermite multiplier', max_scal2.herm_mult.__repr__())


        print('\nMichaelis-Menten Reparametrisation 1')
        print('Changing invariants by column operations on the Hermite multiplier')
        max_scal2.multiplier_add_columns(2, -1, 1)
        max_scal2.multiplier_add_columns(4, -1, -1)
        
        print('Invariants: ', max_scal2.invariants())
        print('Hermite multiplier', max_scal2.herm_mult.__repr__())
        
        # Print translated system
        print('Reduced system:')
        print(max_scal2.translate(reduced_system_km1))#.to_tex()
        
        print('\nMichaelis-Menten Reparametrisation 2')
        print('Changing invariants by column operations on the Hermite multiplier')
        # Divide time through by epsilon
        max_scal2.multiplier_add_columns(2, -1, -1)
        
        
        print('Invariants: ', max_scal2.invariants())
        # print(max_scal2.herm_mult.__repr__())
        
        # Print translated system
        print('Reduced system:')
        print(max_scal2.translate(reduced_system_km1))#.to_tex()

        return


    def with_initial_conditions_and_a_constant():
        '''
        silviana wants this one to exactly match the example in the paper, in the section by the same name.
        '''
        system_tex = system_tex_two_var()
        # system_tex_reduced_km1 = system_tex.replace('k_{-1}', '(K - k_2)')
        
        system = ODESystem.from_tex(system_tex)
        system.add_constraint('Km', '(k_m1 + k_2) / k_1')

        system.reorder_variables(['t', 's', 'c', 'Km', 'k_m1', 'k_2', 'k_1', 'e_0'])
        
        # Print variable order

        
        # # Print scaling matrices
        # max_scal1 = ODETranslation.from_ode_system(system)
        # print('Scaling matrix:')
        # print(max_scal1.scaling_matrix.__repr__())
        
        # # Print invariants
        # print('Hermite multiplier', max_scal1.herm_mult.__repr__())
        # print('Invariants: ', max_scal1.invariants())
        # # print(',\quad '.join(map(expr_to_tex, max_scal1.invariants())))
        
        # # Print translated system
        # print('Reduced system:')
        # print(max_scal1.translate(system))#.to_tex()
        
        print('\n\nAdding in the initial condition for s')
        print('-------------------------------------')
        system.update_initial_conditions({'s': 's0'})


        print('Variable order: ', system.variables)
        
        print('Power Matrix:', system.power_matrix().__repr__())
        translation = ODETranslation.from_ode_system(system, renaming_scheme=('tau',['u','v'], 'c'))
        
        
        print('Scaling matrix:',translation.scaling_matrix.__repr__())

        print('Hermite multiplier', translation.herm_mult.__repr__())
        print('Invariants: ', translation.invariants())


        print('Reduced system:')
        print(translation.translate(system))#.to_tex()
        

        # print('\nMichaelis-Menten Reparametrisation 1')
        print('Changing invariants by column operations on the Hermite multiplier')
        # Python start at 0, but  the paper starts at 1, so these are down by 1 from the paper
        translation.multiplier_add_columns(i=2, j=8, alpha=1)
        translation.multiplier_add_columns(i=4, j=8, alpha=-1)
        
        print('Invariants: ', translation.invariants())
        print('Hermite multiplier', translation.herm_mult.__repr__())
        

        # print('\nMichaelis-Menten Reparametrisation 2')
        # print('Changing invariants by column operations on the Hermite multiplier')
        # # Divide time through by epsilon
        # translation.multiplier_add_columns(i=0, j=8, alpha=1)
        
        # print('Hermite multiplier', translation.herm_mult.__repr__())
        # print('Invariants: ', translation.invariants())
        
        print('W:', translation.inv_herm_mult.__repr__())
        # Print translated system

        print('substitutions:',translation.translate_parameter_substitutions(system))

        reduced_system = translation.translate(system)

        print('Reduced system, in arbitrary parameter names:')
        print(reduced_system)#.to_tex()

        to_sub =  {'c0': 'K', 'c1': 'K-L', 'c2':'L', 'c3':'epsilon'}
        reduced_system = reduced_system.diff_subs(to_sub, subs_constraints=True)

        print('after subbing into canonical names from Murray')
        print(reduced_system)#.to_tex()

        return


    def another_way():

        print('\n---------------\nSetting k_1 = (k_m1 + k_2) / Km, and seeing what happens\n-------\n')

        system_tex = system_tex_two_var()
        # system_tex_reduced_km1 = system_tex.replace('k_1', '(k_{-1} + k_2) / Km')
        
        reduced_system_km1 = ODESystem.from_tex(system_tex)
        reduced_system_km1.add_constraint('k_1', '(k_m1 + k_2) / Km')

        # print('\n\nAdding in the initial condition for s')
        # print('-------------------------------------')
        reduced_system_km1.update_initial_conditions({'s': 's_0'})

        reduced_system_km1.reorder_variables(['t', 's', 'c', 'Km', 'k_m1', 'k_2', 'k_1', 'e_0', 's_0'])
        
        # Print variable order
        
        print('Power Matrix:', reduced_system_km1.power_matrix().__repr__())
        print('constants: ',reduced_system_km1.constant_variables)

        max_scal2 = ODETranslation.from_ode_system(reduced_system_km1)
        
        print('Scaling matrix:',max_scal2.scaling_matrix.__repr__())
        print('Variable order: ', reduced_system_km1.variables)
        print('Hermite multiplier', max_scal2.herm_mult.__repr__())
        print('Invariants: ', max_scal2.invariants())
        print('W:', max_scal2.inv_herm_mult.__repr__())
        print('substitutions:',max_scal2.translate_parameter_substitutions(reduced_system_km1))
        print('Reduced system:', max_scal2.translate(reduced_system_km1))#.to_tex()
        

        # # print('\nMichaelis-Menten Reparametrisation 1')
        # print('Changing invariants by column operations on the Hermite multiplier')
        # # Python start at 0, but  the paper starts at 1, so these are down by 1 from the paper
        # max_scal2.multiplier_add_columns(i=2, j=8, alpha=1)
        # max_scal2.multiplier_add_columns(i=4, j=8, alpha=-1)
        
        # print('Hermite multiplier', max_scal2.herm_mult.__repr__())
        # print('Invariants: ', max_scal2.invariants())
        
        
        # # print('\nMichaelis-Menten Reparametrisation 2')
        # # print('Changing invariants by column operations on the Hermite multiplier')
        # # # Divide time through by epsilon
        # # max_scal2.multiplier_add_columns(i=0, j=8, alpha=1)
        
        # # print('Hermite multiplier', max_scal2.herm_mult.__repr__())
        # # print('Invariants: ', max_scal2.invariants())
        
        # print('W:', max_scal2.inv_herm_mult.__repr__())

        
        # print('Reduced system:', max_scal2.translate(reduced_system_km1))#.to_tex()

        return

    def multiple_timescales():
        # print('Michaelis-Menten Reparametrisation 2')
        print('What if epsilon = e_0 / s_0 is not small?')
        
        system_tex = system_tex_two_var()
        # Substitute K_m into the equations
        system_tex_reduced_l = system_tex.replace('k_{-1}', '(K - k_2)').replace('K', 'K_m k_1')
        # Now set L = K_m + s_0
        # system_tex_reduced_l = system_tex_reduced_km.replace('K_m', '(L - s_0)')

        print('In order ')
        print(system_tex_reduced_l)

        # return
        # print system_tex_reduced_km
        reduced_system_l = ODESystem.from_tex(system_tex_reduced_l)
        reduced_system_l.update_initial_conditions({'s': 's_0'})
        reduced_system_l.add_constraint('L', 's_0 + K_m')
        reduced_system_l.reorder_variables(['t', 's', 'c', 'k_2', 'k_1', 'e_0', 's_0', 'L', 'K_m'])

        translation = ODETranslation.from_ode_system(reduced_system_l)

        # print translation.invariants()
        # Scale t correctly to t/t_C = k_1 L t
        translation.multiplier_add_columns(2, -1, 1)
        # Scale s correctly to s / s_0
        translation.multiplier_add_columns(3, -2, -1)

        # Scale c correctly to c / (e_0 s_0 / L)
        translation.multiplier_add_columns(4, 6, -1)
        translation.multiplier_add_columns(4, 7, -1)
        translation.multiplier_add_columns(4, -1, 1)
        # Find kappa = k_{-1} / k_2 = (K_m k_1 / k_2) - 1
        translation.multiplier_negate_column(5)
        # Find epsilon = e_0 / L
        translation.multiplier_add_columns(6, -1, -1)
        # Find sigma = s_0 / K_m

        print('invariants:', translation.invariants())

        # silviana says: i dislike that these are hardcoded.  these should be computed.
        print('We now have:')
        print('c_0 = kappa + 1')
        print('c_1 = epsilon')
        print('c_2 = sigma')
        print('c_3 = L / K_m')
        # print reduced_system_l
        print(translation.translate(reduced_system_l))

        return


    with_initial_conditions_and_a_constant()
    # another_way()

if __name__ == '__main__':
    # Enable pretty printing
    sympy.init_printing(pretty_print=True, use_latex=True)

    michaelis_menten_two_variable()

    # michaelis_menten_four_variables()
