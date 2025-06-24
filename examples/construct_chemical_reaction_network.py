# an example that shows how to construct 
# a system of ordinary diffirential equations 
# from a chemical reaction network

# do the necessary imports
import sympy
from desr.chemical_reaction_network import ChemicalSpecies, ChemicalReactionNetwork, Complex, Reaction
from desr.ode_system import ODESystem


# make variables, then convert to `ChemicalSpecies` from `desr`.
species = sympy.var('x1 x2')
species = list(map(ChemicalSpecies, species))
x1, x2 = species

# form the complexes for each chemical reaction
complex0 = Complex({x1: 1, x2: 1})
complex1 = Complex({x2: 2})
complex2 = Complex({x1: 1})
complex3 = Complex({x2: 1})

# lump together so can later pass to `ChemicalReactionNetwork` constructor
complexes = (complex0, complex1, complex2, complex3)

# make two reactions, left hand and right hand side of each
r1 = Reaction(complex0, complex1)
r2 = Reaction(complex3, complex2)

# lump together
reactions = [r1, r2]


# construct the network from 
#   the individual species, 
#   the species in their complexes, and 
#   the reactions between complexes

reaction_network = ChemicalReactionNetwork(species, complexes, reactions)

# make an ODE system from the reaction network
system = reaction_network.to_ode_system()
print(system)

# check that the result is the same as an 
# already-known representation for this reaction network.
answer = ODESystem.from_equations('dx1/dt = -k_0_1*x1*x2 + k_3_2*x2\ndx2/dt = k_0_1*x1*x2 - k_3_2*x2')
print(answer)

assert system == answer
