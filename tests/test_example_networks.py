"""
The five signalling networks in `example_networks.py`, reduced.

They are the largest systems in the repository.  These are properties rather than
printed output: that the reduction is the size the translator predicts, that the
invariants really are invariant, that the round trip returns the symbols it started
from, and that the numbers tabulated in
`docs/source/examples/signalling_networks.rst` are the ones the library gives.
"""

from unittest import TestCase

import sympy

from desr.ode_translation import ODETranslator, ReductionScheme

from example_networks import NETWORKS
from oracles import false_equations


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
