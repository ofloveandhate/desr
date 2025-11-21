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
\frac{dx_1}{dt} &=-(a_1+a_21) x_1+u_1\\
\frac{dx_2}{dt} &= a_21 - a_2 x_2
'''


system_a = ODESystem.from_tex(tex_a)
# system_a.add_constraint(lhs='y_2',rhs='x_2')


translation_a = run_symmetry_reduction(system_a)


print('\n\nB\n\n')
tex_b =  \
r'''
\frac{dx_1}{dt} &=-b_11 x_1+u_1\\
\frac{dx_2}{dt} &= a_21 - a_2 x_2
'''

system_b = ODESystem.from_tex(tex_b)
# system_a.add_constraint(lhs='y_2',rhs='x_2')


translation_b = run_symmetry_reduction(system_b)


print('\n\nC\n\n')
tex_c =  \
r'''
\frac{dx_1}{dt} &=-b_11 x_1+u_1\\
\frac{dx_2}{dt} &= a_21 - a_2 x_2
'''

system_c = ODESystem.from_tex(tex_c)
system_c.add_constraint(lhs='b_11',rhs='a_1 + a_21')


translation_c = run_symmetry_reduction(system_c)
