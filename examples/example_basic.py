import sympy
from desr.chemical_reaction_network import ChemicalSpecies, ChemicalReactionNetwork, Complex, Reaction
from desr.ode_system import ODESystem

species = sympy.var('x1 x2')
species = list(map(ChemicalSpecies, species))
x1, x2 = species

complex0 = Complex({x1: 1, x2: 1})
complex1 = Complex({x2: 2})
complex2 = Complex({x1: 1})
complex3 = Complex({x2: 1})
complexes = (complex0, complex1, complex2, complex3)

r1 = Reaction(complex0, complex1)
r2 = Reaction(complex3, complex2)
reactions = [r1, r2]

reaction_network = ChemicalReactionNetwork(species, complexes, reactions)
system = reaction_network.to_ode_system()
print(system)

answer = ODESystem.from_equations('dx1/dt = -k_0_1*x1*x2 + k_3_2*x2\ndx2/dt = k_0_1*x1*x2 - k_3_2*x2')
print(answer)
print(system == answer)
