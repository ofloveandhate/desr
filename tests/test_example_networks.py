"""
The six signalling networks in `example_networks.py`, reduced.

They are the largest systems in the repository.  These are properties rather than
printed output: that the reduction is the size the translator predicts, that the
invariants really are invariant, that the round trip returns the symbols it started
from, and that the numbers tabulated in
`docs/source/examples/signalling_networks.rst` are the ones the library gives.
"""

from unittest import TestCase

import sympy

from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.ode_translation import ODETranslator, ReductionScheme

from example_networks import NETWORKS
from oracles import false_equations


def michaelis_menten(substrate, enzyme, intermediate, product):
    ''' `S + E <-> C -> P + E`, as the three one-way reactions `from_diagram` reads. '''
    return ['{} + {} -> {}'.format(substrate, enzyme, intermediate),
            '{} -> {} + {}'.format(intermediate, substrate, enzyme),
            '{} -> {} + {}'.format(intermediate, product, enzyme)]


#: The reactions of each network.  For the five from Feliu & Wiuf, arXiv:1101.5799v2,
#: these are the reaction diagrams printed with each motif: (a) on p. 21, (e) on p. 27,
#: (g) on p. 34, (h) on p. 37, (j) on p. 41.  The kinase is `Ev` here where the paper has `E`.  The
#: phosphorelay's are the reactions its equations describe.
REACTIONS = {
    'one-site modification':
        michaelis_menten('S0', 'Ev', 'X', 'S1') + michaelis_menten('S1', 'F', 'Y', 'S0'),
    'two-site modification':
        michaelis_menten('S0', 'Ev1', 'X1', 'S1') + michaelis_menten('S1', 'Ev2', 'X2', 'S2')
        + michaelis_menten('S1', 'F1', 'Y1', 'S0') + michaelis_menten('S2', 'F2', 'Y2', 'S1'),
    'two-site modification, shared enzymes':
        michaelis_menten('S0', 'Ev', 'X1', 'S1') + michaelis_menten('S1', 'Ev', 'X2', 'S2')
        + michaelis_menten('S1', 'F', 'Y1', 'S0') + michaelis_menten('S2', 'F', 'Y2', 'S1'),
    'two-substrate modification':
        michaelis_menten('S0', 'Ev', 'X1', 'S1') + michaelis_menten('P0', 'Ev', 'X2', 'P1')
        + michaelis_menten('S1', 'F1', 'Y1', 'S0') + michaelis_menten('P1', 'F2', 'Y2', 'P0'),
    'two-layer cascade':
        michaelis_menten('S0', 'Ev', 'X1', 'S1') + michaelis_menten('P0', 'S1', 'X2', 'P1')
        + michaelis_menten('S1', 'F1', 'Y1', 'S0') + michaelis_menten('P1', 'F2', 'Y2', 'P0'),
    'two-step phosphorelay':
        michaelis_menten('X1', 'S0', 'Y1', 'Xs1') + michaelis_menten('Xs1', 'F1', 'Y2', 'X1')
        + ['Xs1 + X2 -> X1 + Xs2'] + michaelis_menten('Xs2', 'F2', 'Y3', 'X2'),
}


def where_it_appears(system, constant):
    '''
    Every term a rate constant is in: the species whose equation it is, and the term
    with the constant divided out.  Two rate constants that appear in the same places
    are the same constant under two names.
    '''
    found = set()
    for species, rate in system.derivative_dict.items():
        if species == system.independent_variable:
            continue
        for term in sympy.Add.make_args(sympy.expand(rate)):
            if term.has(constant):
                found.add((species, sympy.simplify(term / constant)))
    return frozenset(found)


class TestTheNetworksAreTheirReactions(TestCase):
    '''
    The equations in `example_networks.py` were written out by hand.  The same networks
    can be given as reactions and left to the library to turn into mass-action
    equations.  The two must be one system, up to what the rate constants are called.

    This is what proofreads the transcription.  It is how a wrong symbol in the
    phosphorelay was found: one term had `km5` for `km6`, a constant no reaction has.
    '''

    def test_every_network_has_its_reactions_listed(self):
        self.assertEqual(set(NETWORKS), set(REACTIONS))

    def test_the_hand_written_equations_are_the_mass_action_ones(self):
        for label, build in NETWORKS.items():
            with self.subTest(network=label):
                written = build()
                derived = ChemicalReactionNetwork.from_diagram(
                    '\n'.join(REACTIONS[label])).to_ode_system()
                self.assertEqual(set(written.dependents), set(derived.dependents))
                self.assertEqual(len(written.parameters), len(derived.parameters))

                places = {}
                for constant in derived.parameters:
                    places.setdefault(where_it_appears(derived, constant), []).append(constant)
                renaming = {}
                for constant in written.parameters:
                    match = places.get(where_it_appears(written, constant), [])
                    self.assertEqual(1, len(match),
                                     '{} has no single counterpart'.format(constant))
                    renaming[match[0]] = constant
                self.assertEqual(len(derived.parameters), len(renaming))

                for species in written.dependents:
                    difference = (derived.derivative_dict[species].subs(
                        renaming, simultaneous=True) - written.derivative_dict[species])
                    self.assertEqual(0, sympy.expand(difference), str(species))

    def test_a_wrong_symbol_is_noticed(self):
        # The phosphorelay as it was: km5 where km6 belongs.  No reaction's rate
        # constant appears where km5 then does, or where km6 then does.
        written = NETWORKS['two-step phosphorelay']()
        Xs2, Y3, km5, km6 = sympy.symbols('Xs2 Y3 km5 km6')
        rates = dict(written.derivative_dict)
        rates.pop(written.independent_variable)
        rates[Xs2] = rates[Xs2] - km6 * Y3 + km5 * Y3
        wrong = type(written).from_dict(rates)
        derived = ChemicalReactionNetwork.from_diagram(
            '\n'.join(REACTIONS['two-step phosphorelay'])).to_ode_system()
        places = {where_it_appears(derived, constant) for constant in derived.parameters}
        self.assertNotIn(where_it_appears(wrong, km5), places)
        self.assertNotIn(where_it_appears(wrong, km6), places)
        self.assertEqual(len(derived.parameters) + 1, len(wrong.parameters))


class TestTheExampleReactionNetworks(TestCase):

    def systems(self):
        for label, build in NETWORKS.items():
            yield label, build()

    def test_every_one_of_them_reduces(self):
        for label, system in self.systems():
            with self.subTest(network=label):
                translator = ODETranslator.from_ode_system(system)
                reduced = translator.translate(system)
                self.assertEqual(translator.expected_variable_count(system),
                                 len(reduced.symbols))
                self.assertLess(len(reduced.symbols), len(system.symbols))

    def test_their_invariants_are_dimensionless(self):
        # What makes an invariant an invariant.  Asked of the system rather than of the
        # matrix the invariants came from, so the two have to agree.
        for label, system in self.systems():
            translator = ODETranslator.from_ode_system(system)
            for invariant in translator.invariant_monomials():
                with self.subTest(network=label, invariant=str(invariant)):
                    self.assertEqual(1, system.dimension_of(invariant))

    def test_the_scaling_leaves_the_equations_alone(self):
        # Checked against the equations, not against the matrix that was solved from
        # them: every term of every dimensionless expression has weight zero.
        for label, system in self.systems():
            scaling = system.maximal_scaling_matrix()
            symbols = list(system.symbols)
            for row in system.dimensionless_expressions():
                for term in sympy.Add.make_args(row['expression']):
                    with self.subTest(network=label, term=str(term)):
                        powers = term.as_powers_dict()
                        weights = scaling * sympy.Matrix(
                            [powers.get(symbol, 0) for symbol in symbols])
                        self.assertTrue(weights.is_zero_matrix)

    def test_the_round_trip_returns_the_original_symbols(self):
        # Substituting the forward monomials into the reverse ones has to give the
        # symbols back, which is the whole of what it means for V and its inverse to
        # be inverse on these systems.
        for label, system in self.systems():
            with self.subTest(network=label):
                translator = ODETranslator.from_ode_system(system)
                forward = dict(zip(translator.column_names(),
                                   translator.forward_monomials()))
                recovered = [
                    sympy.simplify(monomial.subs(forward, simultaneous=True))
                    for monomial in translator.reverse_monomials()]
                self.assertEqual(list(system.symbols), recovered)

    def test_every_reduced_equation_is_true(self):
        for label, system in self.systems():
            with self.subTest(network=label):
                self.assertEqual([], false_equations(
                    system, system.reduce_scaling_symmetries()))

    #: network: (dependent variables, rate constants, rate constants left)
    TABULATED = {
        'one-site modification': (6, 6, 4),
        'two-step phosphorelay': (10, 10, 8),
        'two-site modification': (11, 12, 10),
        'two-site modification, shared enzymes': (9, 12, 10),
        'two-substrate modification': (11, 12, 10),
        'two-layer cascade': (11, 12, 10),
    }

    def test_the_table_on_the_documentation_page(self):
        self.assertEqual(set(NETWORKS), set(self.TABULATED))
        for label, system in self.systems():
            with self.subTest(network=label):
                translator = ODETranslator.from_ode_system(system)
                reduced = translator.translate(system)
                self.assertEqual(self.TABULATED[label],
                                 (len(system.dependents), len(system.parameters),
                                  len(reduced.parameters)))
                self.assertEqual(2, translator.num_scaling_symmetries)
                self.assertEqual(ReductionScheme.PARAMETER, translator.scheme())

    def test_they_all_choose_the_same_two_units(self):
        # Time in 1/k2, and every concentration in k2/k1.
        t, k1, k2 = sympy.symbols('t k1 k2')
        for label, system in self.systems():
            with self.subTest(network=label):
                reduced = system.reduce_scaling_symmetries()
                denotes = reduced.named_symbols
                self.assertEqual(k2 * t, denotes[reduced.independent_variable])
                for symbol in reduced.dependents:
                    species = sympy.simplify(denotes[symbol] * k2 / k1)
                    self.assertIn(species, system.dependents)
