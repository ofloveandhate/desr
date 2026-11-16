"""
The reaction networks in ``examples/chemical_reaction_networks.py``, reduced.

They are the largest systems in the repository -- the two-step phosphorelay has ten
symbols -- and until now nothing checked their reductions.  The script prints them, and
CI checked only that it exited zero, so a change of convention could have altered every
one of these answers without anything saying so.

These are properties rather than printed output: that the reduction is the size the
translator predicts, that the invariants really are invariant, and that the round trip
returns the symbols it started from.  A convention that drifts breaks at least one of
them, and none of them has to be rewritten when a naming scheme changes.
"""

import os
import sys
from unittest import TestCase

import sympy

from desr.ode_translation import ODETranslator

# examples/ is not part of the installed package, so it is reached by path rather than
# by import.  If those scripts move, this moves with them.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

from examples import chemical_reaction_networks as networks    # noqa: E402


class TestTheExampleReactionNetworks(TestCase):

    BUILDERS = (
        ('two step phosphorelay', networks.example_two_step_phosphorelay),
        ('one site modification', networks.example_one_site_modification),
        ('two site modification', networks.example_two_site_modification),
        ('two substrate modification', networks.example_two_substrate_modification),
        ('two layer cascade', networks.example_two_layer_cascade),
    )

    def systems(self):
        import warnings
        for label, build in self.BUILDERS:
            with warnings.catch_warnings():
                # The phosphorelay warns that conservation forces km5 == km6, which is
                # a statement about the model and not something to hear in a test run.
                warnings.simplefilter('ignore', UserWarning)
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
                forward = dict(zip(translator.column_names(system),
                                   translator.forward_monomials()))
                recovered = [
                    sympy.simplify(monomial.subs(forward, simultaneous=True))
                    for monomial in translator.reverse_monomials(system)]
                self.assertEqual(list(system.symbols), recovered)
