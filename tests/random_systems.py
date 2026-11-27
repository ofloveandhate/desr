'''
Small systems drawn from a seeded generator, for sweeping statements over more systems
than anyone would write by hand.

Each right-hand side is a sum of Laurent monomials with integer coefficients, which is
the class of system the scaling reduction is defined on.  A seed fixes the system, so a
failure names a seed and can be reproduced.

Not a test module itself.
'''
import random

import sympy

from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.ode_system import ODESystem


#: Shapes of system to draw, by name.
SHAPES = {
    'usual': {},
    'small': {'max_dependents': 1, 'max_parameters': 2, 'max_terms': 2},
    'large': {'max_dependents': 4, 'max_parameters': 6, 'max_terms': 4, 'max_power': 3},
    'non-autonomous': {'clock_chance': 0.7},
    # Mass-action kinetics of a random reaction network, built by the library.  Every
    # term has a rate constant of its own, and these are the systems the parameter
    # scheme applies to with symmetries to remove; among the other shapes it applies
    # almost only where there is no symmetry at all.
    'reaction network': {'reaction_network': True},
}


def random_system(seed, max_dependents=3, max_parameters=4, max_terms=3, max_power=2,
                  clock_chance=0.25, reaction_network=False):
    '''
    The system for `seed`.

    Args:
        seed (int): Fixes everything drawn.
        max_dependents (int): Between one and this many dependent variables.
        max_parameters (int): Between none and this many parameters.
        max_terms (int): Between one and this many monomials in each right-hand side.
        max_power (int): Exponents are drawn from ``-max_power .. max_power``.
        clock_chance (float): The chance that a monomial involves the independent
            variable at all.  Most models are autonomous, and a system whose every term
            involves the clock rarely has a symmetry worth reducing.
        reaction_network (bool): Ignore the rest, and return the mass-action system
            of a random reaction network instead: see :func:`random_network_diagram`.

    Returns:
        ODESystem: with every parameter drawn appearing in at least one equation.
    '''
    draw = random.Random(seed)
    if reaction_network:
        return _reaction_network(draw)
    t = sympy.Symbol('t')
    # `u` and `p`, because the reduction's own default names are x, y, tau, nu, kappa.
    dependents = sympy.symbols('u1:{}'.format(draw.randint(1, max_dependents) + 1))
    parameters = sympy.symbols('p1:{}'.format(draw.randint(0, max_parameters) + 1))

    def monomial():
        result = sympy.Integer(draw.choice([-3, -2, -1, 1, 1, 1, 2, 3]))
        # Few symbols in each term, as in a model: a term in everything at once makes
        # a system with no symmetry at all.
        involved = draw.sample(dependents + parameters,
                               draw.randint(1, min(3, len(dependents + parameters))))
        for symbol in involved:
            result *= symbol ** draw.choice(
                [power for power in range(-max_power, max_power + 1) if power != 0])
        if draw.random() < clock_chance:
            result *= t ** draw.choice([-1, 1])
        return result

    while True:
        rates = {dependent: sum(monomial() for _ in range(draw.randint(1, max_terms)))
                 for dependent in dependents}
        used = set().union(*(rate.free_symbols for rate in rates.values()))
        # A right-hand side that cancelled to zero, or a parameter in no equation, is a
        # different system from the one drawn; draw again.
        if all(rate != 0 for rate in rates.values()) and set(parameters) <= used:
            break

    symbols = [t] + list(dependents) + list(parameters)
    derivatives = ([sympy.Integer(1)] + [rates[dependent] for dependent in dependents]
                   + [None] * len(parameters))
    return ODESystem._from_variables(symbols, derivatives, indep_var=t)



def random_network_diagram(seed, max_species=4, max_reactions=5):
    '''
    A reaction diagram for `seed`, in the form
    :meth:`~desr.chemical_reaction_network.ChemicalReactionNetwork.from_diagram` reads.

    Each reaction has at most two species on a side, each once or twice.  A side may be
    empty, which is an inflow or a degradation.
    '''
    draw = random.Random(seed)
    species = ['s{}'.format(index + 1) for index in range(draw.randint(2, max_species))]

    def side():
        chosen = draw.sample(species, draw.choice([0, 1, 1, 1, 2, 2]))
        return ' + '.join(name if draw.random() < 0.8 else '2*' + name
                          for name in sorted(chosen))

    reactions = []
    while len(reactions) < draw.randint(2, max_reactions):
        left, right = side(), side()
        if left != right and (left, right) not in reactions:
            reactions.append((left, right))
    return '\n'.join('{} -> {}'.format(left, right) for left, right in reactions)


def _reaction_network(draw):
    ''' The mass-action system of a random reaction network. '''
    return ChemicalReactionNetwork.from_diagram(
        random_network_diagram(draw.randrange(10 ** 9))).to_ode_system()
