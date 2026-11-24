"""
A two-compartment model of a vaccine injection.

The equations are Example 13.6 of J. DiStefano, Dynamic Systems Biology Modeling and
Simulation, Elsevier (2015), a model for the study of R. P. Weissburg et al.,
Characterization of the MN gp120 HIV-1 vaccine: antigen binding to alum,
Pharmaceutical Research 12 (1995) 1439-1446, https://doi.org/10.1023/A:1016266916893.
That paper gives its model as a diagram, its Scheme 1, and the equations here have not
been compared with the book.

Three versions are reduced: as written (A), with the total loss rate a_2 + a_12
replaced by b_22 (B), and with b_22 and the relation that defines it (C).
"""

import sympy
from desr import ODESystem, ODETranslator


sympy.init_printing(pretty_print=True, use_latex=True)



def run_symmetry_reduction(sys):

    print(f'original variables: {sys.symbols}')
    print(f'original system: {sys}')

    translator = ODETranslator.from_ode_system(sys)

    print('\nScaling matrix:')
    print(translator.scaling_matrix.__repr__())

    print('Dimension of scaling action: ', translator.r)
    # Print invariants
    print('Invariants: ', translator.invariant_monomials())
    print('Substitutions: ', translator.translate_parameter_substitutions(system=sys))

    # Print translated system
    reduced_system = translator.translate(sys)

    print(f'\nreduced variables: {reduced_system.symbols}')
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
translator_a = run_symmetry_reduction(system_a)




print('\n\nB\n\n')
tex_b =  \
r'''
\frac{dz_1}{dt} &= -a_1 z_1 + a_12 z_2\\
\frac{dz_2}{dt} &= -b_22 z_2 + u_2
'''

system_b = ODESystem.from_tex(tex_b)
# system_b.add_constraint(lhs='w_1',rhs='z_1')
# desr cannot have constraints on non-constant variables, so commented out



translator_b = run_symmetry_reduction(system_b)


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


translator_c = run_symmetry_reduction(system_c)
