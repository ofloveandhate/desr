import sympy
from desr import ODESystem, ODETranslator


sympy.init_printing(pretty_print=True, use_latex=True)



def run_symmetry_reduction(sys):




    translator = ODETranslator.from_ode_system(sys)

    print('Symbol order: ', translator.symbols_domain)

    print('Scaling matrix:')
    print(translator.scaling_matrix.__repr__())

    # Print invariants
    print('Invariants: ', translator.invariants())

    # Print translated system
    print('Reduced system:')
    print(translator.translate(sys))
    return translator




print("\n\n\nsystem a")

system_a = ODESystem.from_tex('\frac{dx}{dt} &= (a + b) x + b')
translator_a = run_symmetry_reduction(system_a)


# the same system, but naively with c = a+b

print("\n\n\nsystem b")
system_b = ODESystem.from_tex('\frac{dx}{dt} &= (c) x + b')
translator_b = run_symmetry_reduction(system_b)


# the same system, i already did substitution c = a+b
# but i'm also letting the software know that b,c have to have the same unit
print("\n\n\nsystem c")
system_c = ODESystem.from_tex('\frac{dx}{dt} &= (c) x + b')
system_c.require_same_unit({'c':'b'})

translator_c = run_symmetry_reduction(system_c)


# let the software do the substitution,
print("\n\n\nsystem d")
system_d = ODESystem.from_tex('\frac{dx}{dt} &= (a + b) x + b')
system_d.add_constraint('c', 'a+b')

translator_d = run_symmetry_reduction(system_d)

