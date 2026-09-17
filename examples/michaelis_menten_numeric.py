"""Michaelis-Menten: solve the reduced system, recover the original.

Companion figure to ``docs/source/examples/mm_numeric_round_trip.rst``.

The reduced system is solved numerically, its solution is translated back into the
original variables with :class:`desr.numerics.NumericTranslation`, and the result is
compared against solving the original system directly.  The two agree to the tolerance
of the integrator, which is the point: the translation itself is exact.

Run with::

    python examples/michaelis_menten_numeric.py
"""

import matplotlib.pyplot as plt
import numpy as np
import sympy
from scipy.integrate import solve_ivp

from desr.numerics import NumericTranslation
from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslation

# Categorical slots 1 and 2 of the reference palette; identity is carried by hue and
# repeated as a direct label, never by colour alone.
SUBSTRATE = '#2a78d6'
COMPLEX = '#eb6834'
INK = '#0b0b0b'
INK_SOFT = '#52514e'
GRID = '#dcdbd7'


# ---------------------------------------------------------------- the system
SYSTEM_TEX = r'''\frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\
                 \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c'''

system = ODESystem.from_tex(SYSTEM_TEX)
system.update_initial_conditions({'s': 's_0'})
system.reorder_variables(['t', 's', 'c', 'k_m1', 'k_2', 'k_1', 'e_0', 's_0'])

translation = ODETranslation.from_ode_system(system, naming_scheme=('tau', ['u', 'v'], 'c'))
reduced_system = translation.translate(system)
numeric = NumericTranslation(system, translation, reduced_system)

t, s, c, k_m1, k_2, k_1, e_0, s_0 = system.variables
tau, u, v, c0, c1, c2 = reduced_system.variables

parameters = {k_1: 1.5, k_m1: 0.9, k_2: 0.4, e_0: 0.3, s_0: 2.0}
initial_state = [2.0, 0.0]
final_time = 70.0


def as_function(a_system, constants):
    """Turn a system with its constants fixed into a right-hand side for solve_ivp."""
    variables = list(a_system.non_constant_variables)
    rhs = [a_system.derivative_dict[x].subs(constants) for x in variables]
    return sympy.lambdify([a_system.indep_var, variables], rhs, modules='numpy')


# --------------------------------------------------------- solve both systems
times = np.linspace(0.0, final_time, 400)
reference_soln = solve_ivp(as_function(system, parameters), (0.0, final_time), initial_state,
                      t_eval=times, rtol=1e-10, atol=1e-12)

start = numeric.forward({t: 0.0, s: initial_state[0], c: initial_state[1], **parameters})
reduced_times = numeric.forward({t: times, **parameters})[tau]
reduced_parameters = {c0: start[c0], c1: start[c1], c2: start[c2]}

reduced_solution = solve_ivp(as_function(reduced_system, reduced_parameters),
                     (reduced_times[0], reduced_times[-1]),
                     [float(start[u]), float(start[v])],
                     t_eval=reduced_times, rtol=1e-10, atol=1e-12)

recovered_soln = numeric.reverse({tau: reduced_solution.t, u: reduced_solution.y[0], v: reduced_solution.y[1],
                             **reduced_parameters},
                            known_values={k_1: parameters[k_1], s_0: parameters[s_0]})


# ------------------------------------------------------------------- the plot
def style(axes, title, xlabel, ylabel):
    axes.set_title(title, color=INK, fontsize=11, pad=10)
    axes.set_xlabel(xlabel, color=INK_SOFT, fontsize=9)
    axes.set_ylabel(ylabel, color=INK_SOFT, fontsize=9)
    axes.grid(True, color=GRID, linewidth=0.6)
    axes.set_axisbelow(True)
    axes.tick_params(colors=INK_SOFT, labelsize=8)
    for edge in ('top', 'right'):
        axes.spines[edge].set_visible(False)
    for edge in ('left', 'bottom'):
        axes.spines[edge].set_color(GRID)


def label_end(axes, x, y, text, colour):
    axes.annotate(text, xy=(x[-1], y[-1]), xytext=(4, 0), textcoords='offset points',
                  color=colour, fontsize=9, va='center')


figure, (left, middle, right) = plt.subplots(1, 3, figsize=(13.0, 4.2))

# The reduced system, in its own variables.
left.plot(reduced_solution.t, reduced_solution.y[0], color=SUBSTRATE, linewidth=2)
left.plot(reduced_solution.t, reduced_solution.y[1], color=COMPLEX, linewidth=2)
label_end(left, reduced_solution.t, reduced_solution.y[0], r'$u$', SUBSTRATE)
label_end(left, reduced_solution.t, reduced_solution.y[1], r'$v$', COMPLEX)
style(left, 'Reduced system\n3 parameters', r'$\tau$', 'invariant')

# The original system: solved directly, and recovered from the reduced solution.
marks = slice(None, None, 25)
middle.plot(times, reference_soln.y[0], color=SUBSTRATE, linewidth=2, label='solved directly')
middle.plot(times, reference_soln.y[1], color=COMPLEX, linewidth=2)
middle.plot(recovered_soln[t][marks], recovered_soln[s][marks], linestyle='none', marker='o',
            markersize=8, markerfacecolor='none', markeredgewidth=1.6,
            color=SUBSTRATE, label='solution translation')
middle.plot(recovered_soln[t][marks], recovered_soln[c][marks], linestyle='none', marker='o',
            markersize=8, markerfacecolor='none', markeredgewidth=1.6, color=COMPLEX)
label_end(middle, times, reference_soln.y[0], r'$s$', SUBSTRATE)
label_end(middle, times, reference_soln.y[1], r'$c$', COMPLEX)
style(middle, 'Original system\n5 parameters', r'$t$', 'concentration')
legend = middle.legend(frameon=False, fontsize=9, loc='center right',
                       handler_map={}, labelcolor=INK_SOFT)
for handle in legend.legend_handles:
    handle.set_color(INK_SOFT)
    handle.set_markeredgecolor(INK_SOFT)

# What the round trip costs: nothing the integrator did not already cost.  The
# difference sits two orders of magnitude below the tolerance we asked the solver for,
# so the translation contributes nothing measurable of its own.
difference_s = np.abs(recovered_soln[s] - reference_soln.y[0])
difference_c = np.abs(recovered_soln[c] - reference_soln.y[1])
right.axhline(1e-10, color=INK_SOFT, linewidth=1, linestyle=(0, (4, 3)))
right.annotate('tolerance asked of the solver', xy=(0.2, 1e-10), xytext=(0, 5),
               textcoords='offset points', color=INK_SOFT, fontsize=8)
right.semilogy(times, difference_s, color=SUBSTRATE, linewidth=1.2, label=r'$s$')
right.semilogy(times, difference_c, color=COMPLEX, linewidth=1.2, label=r'$c$')
right.set_ylim(1e-15, 3e-10)
style(right, 'Difference between solutions to nondimensionalized model',
      r'$t$', 'absolute difference')
error_legend = right.legend(frameon=False, fontsize=9, loc='lower right',
                            labelcolor=INK_SOFT, ncols=2)

figure.tight_layout()


if __name__ == '__main__':
    print(reduced_solution)
    print('invariants  :', translation.invariants())
    print('recovered_soln   :', {str(x): round(float(recovered_soln[x]), 10)
                            for x in (k_m1, k_2, e_0)})
    print('largest difference in s: {:.2e}'.format(
        np.abs(recovered_soln[s] - reference_soln.y[0]).max()))
    print('largest difference in c: {:.2e}'.format(
        np.abs(recovered_soln[c] - reference_soln.y[1]).max()))
    plt.show()
