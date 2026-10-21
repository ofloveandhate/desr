from unittest import TestCase

import sympy

from desr.chemical_reaction_network import (ChemicalReactionNetwork, ChemicalSpecies,
                                            Complex, Reaction)
from desr.ode_system import ODESystem


class TestChemicalReactionNetwork(TestCase):
    ''' Test cases for the ChemicalReactionNetwork class '''

    def test_crn_harrington(self):
        ''' Example 2.8 from Harrington - Joining and decomposing '''
        species = list(map(ChemicalSpecies, sympy.var('x1 x2')))
        x1, x2 = species

        complex0 = Complex()
        complex1 = Complex({x1: 1})
        complex2 = Complex({x2: 1})
        complexes = (complex0, complex1, complex2)

        r01 = Reaction(complex0, complex1)
        r12 = Reaction(complex1, complex2)
        r21 = Reaction(complex2, complex1)
        r20 = Reaction(complex2, complex0)
        reactions = [r01, r12, r21, r20]

        reaction_network = ChemicalReactionNetwork(species, complexes, reactions)

        system = reaction_network.to_ode_system()

        answer = ODESystem.from_equations('dx1/dt = k_0_1 - k_1_2*x1 + k_2_1*x2\ndx2/dt = k_1_2*x1 - k_2_0*x2 - k_2_1*x2')

        self.assertEqual(system, answer)

    def test_crn_harrington2(self):
        ''' Example 1 from Harrington board notes - Joining and decomposing '''
        species = list(map(ChemicalSpecies, sympy.var('x1 x2')))
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

        answer = ODESystem.from_equations('dx1/dt = -k_0_1*x1*x2 + k_3_2*x2\ndx2/dt = k_0_1*x1*x2 - k_3_2*x2')

        self.assertEqual(system, answer)
