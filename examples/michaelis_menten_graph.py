"""Michaelis-Menten: the reduction drawn as a graph.

Companion figure to ``docs/source/examples/mm_init_cond.rst``.

The Hermite multiplier is the weighted biadjacency matrix of a bipartite graph between
the original symbols and the reduced ones, so the reduction can be drawn rather than only
printed.  Two panels: going out, and coming back.  The source layer is on the left in
both, so the reduced column ends the first panel and begins the second and a round trip
reads as one path across the pair.

What the figure says about this system, which is a page of prose otherwise:

- the reduction divides by ``s_0``, ``e_0`` and ``k_1`` -- the heavy borders -- and does
  not divide by ``c_0``, which is why a zero initial condition survives it;
- ``c_0`` has a single edge, positive, into an invariant rather than an auxiliary, so
  that invariant is simply zero: :math:`\\kappa_4 = c_0/s_0 = 0`;
- coming back, :math:`c_0 = \\kappa_4 x_1`, so the zero returns the way it came;
- the two auxiliaries are hatched.  Going out they are normalised to 1, which is how
  ``s_0`` and ``k_1`` are absorbed; coming back they are the two values a reverse
  translation has to be told.

Run with::

    python examples/michaelis_menten_graph.py
"""

import matplotlib.pyplot as plt

from desr import ODESystem, ODETranslator

system_tex = (r'\frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\'
              '\n'
              r'\frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c')

system = ODESystem.from_tex(system_tex)
system.add_constraint('K_m', '(k_m1 + k_2) / k_1')
system.update_initial_conditions(['s(0) = s_0', 'c(0) = 0'])
system.reorder_variables(['s_0', 'e_0', 'k_1', 'k_2', 'k_m1', 'K_m', 'c_0', 't', 's', 'c'])

translator = ODETranslator.from_ode_system(
    system,
    naming_scheme={'parameters': 'kappa', 'independent': 'tau', 'dependents': ['u', 'v']})

# The same two column operations the narrative makes, so the figure is of the reduction
# the page arrives at rather than of the one desr reaches first.
translator.multiplier_add_columns(translator.column_of('k_1*s_0*t'),
                                  translator.column_of('e_0/s_0'))
translator.multiplier_add_columns(translator.column_of('c/s_0'),
                                  translator.column_of('e_0/s_0'), power=-1)

figure, (going_out, coming_back) = plt.subplots(1, 2, figsize=(16, 9))

translator.draw_symbol_graph(system, ax=going_out, legend=False, font_size=12)
going_out.set_title('forward', fontsize=14)

translator.draw_symbol_graph(system, direction='reverse', ax=coming_back, font_size=12)
coming_back.set_title('reverse', fontsize=14)

figure.tight_layout()


if __name__ == '__main__':
    graph = translator.symbol_graph(system)
    print('divides by, going out :', translator.symbols_that_must_not_be_zero(system))
    print('divides by, coming back:',
          translator.reduced_symbols_that_must_not_be_zero(system))
    print('values the graph carries:')
    for (layer, symbol), data in graph.nodes(data=True):
        if 'value' in data:
            print('   {:<8} {} = {}'.format(layer, symbol, data['value']))
    plt.show()
