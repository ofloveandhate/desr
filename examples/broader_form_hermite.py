import sympy
from desr.matrix_normal_forms import smf
from desr.ode_system import ODESystem
from desr.ode_system import maximal_scaling_matrix, rational_expr_to_exponent_matrix, hnf_col, hnf_row, normal_hnf_col
from desr.matrix_normal_forms import normal_hnf_row
from desr.ode_translation import ODETranslation, scale_action
from desr.tex_tools import expr_to_tex, matrix_to_tex


sympy.init_printing(pretty_print=True, use_latex=True)



def run_symmetry_reduction(sys):




    translation = ODETranslation.from_ode_system(sys)

    print('Variable order: ', translation.variables_domain)

    K = sys.exponent_matrix()
    print('exponent matrix K=',
    K)

    k_hnf_row_form, k_hnf_row_mult = hnf_row(K)

    print(k_hnf_row_form, k_hnf_row_mult)


    print('hermite multiplier:',
    translation.herm_mult
    )

    print('Scaling matrix:')
    print(translation.scaling_matrix.__repr__())

    print('scaling_matrix_hnf',translation._scaling_matrix_hnf)

    # Print invariants
    print('Invariants: ', translation.invariants())

    try:
        print('substitutions:',translation.translate_parameter_substitutions(system))
    except Exception as e:
        print(e)

    # Print translated system
    print('Reduced system:')
    print(translation.translate(sys))
    return translation






system = ODESystem.from_tex('\frac{dz}{dt} = \frac{ c z^3}{t}')
translation = run_symmetry_reduction(system)


