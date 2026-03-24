# implements example from 
#
# R. P. Weissburg, P. W. Berman, J. L. Cleland, D. Eastman, F. Farina, S. Frie, A. Lim, J. Mordenti,
# M. R. Peterson, K. Yim, and M. F. Powell, Characterization of the MN gp120 HIV-1 vaccine:
# Antigen binding to alum, Pharmaceutical Research, 12 (1995), pp. 1439–1446, 
# https://doi.org/10.1023/A:1016266916893.
#
# and
# 
# example 13.6 in 
# J. DiStefano, Dynamic Systems Biology Modeling and Simulation, Elsevier Science, 2015, 
# https://books.google.co.uk/books?id=kkkhmwEACAAJ


import sympy
from desr.matrix_normal_forms import smf
from desr.ode_system import ODESystem
from desr.ode_system import maximal_scaling_matrix, rational_expr_to_power_matrix, hnf_col, hnf_row, normal_hnf_col
from desr.matrix_normal_forms import normal_hnf_row
from desr.ode_translation import ODETranslation, scale_action
from desr.tex_tools import expr_to_tex, matrix_to_tex


sympy.init_printing(pretty_print=True, use_latex=True)



def run_symmetry_reduction(sys):

    print(f'original variables: {sys.variables}')
    print(f'original system: {sys}')

    translation = ODETranslation.from_ode_system(sys)

    print('\nScaling matrix:')
    print(translation.scaling_matrix.__repr__())

    print('Dimension of scaling action: ', translation.r)
    # Print invariants
    print('Invariants: ', translation.invariants())
    print('Substitutions: ', translation.translate_parameter_substitutions(system=sys))

    # Print translated system
    reduced_system = translation.translate(sys)

    print(f'\nreduced variables: {reduced_system.variables}')
    print('\nReduced system:')
    print(reduced_system)


    return reduced_system

print('\n\nA\n\n')
tex_a =  \
r'''
\frac{dz_1}{dt} &= -a_1 z_1 + a_12 z_2\\
\frac{dz_2}{dt} &= -(a_2 + a_12) z_2 + u_2
'''

system_a = ODESystem.from_tex(tex_a)
# system_a.add_constraint(lhs='w_1',rhs='z_1')
# desr cannot have constraints on non-constant variables, so commented out
translation_a = run_symmetry_reduction(system_a)




print('\n\nB\n\n')
tex_b =  \
r'''
\frac{dz_1}{dt} &= -a_1 z_1 + a_12 z_2\\
\frac{dz_2}{dt} &= -b_22 z_2 + u_2
'''

system_b = ODESystem.from_tex(tex_b)
# system_b.add_constraint(lhs='w_1',rhs='z_1')
# desr cannot have constraints on non-constant variables, so commented out



translation_b = run_symmetry_reduction(system_b)


print('\n\nC\n\n')
tex_c =  \
r'''
\frac{dz_1}{dt} &= -a_1 z_1 + a_12 z_2\\
\frac{dz_2}{dt} &= -b_22 z_2 + u_2
'''

system_c = ODESystem.from_tex(tex_c)
# system_c.add_constraint(lhs='w_1',rhs='z_1')
# desr cannot have constraints on non-constant variables, so commented out
system_c.add_constraint(lhs='b_22',rhs='a_2 + a_12')


translation_c = run_symmetry_reduction(system_c)
