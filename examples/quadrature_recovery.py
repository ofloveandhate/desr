"""Recovering a dropped auxiliary variable by quadrature.

Companion figure to ``docs/source/examples/quadrature_round_trip.rst``.

Hubert & Labahn example 6.6 is reduced with ``include_aux_vars=False``, so the auxiliary
variable is never solved for.  It is recovered afterwards by integrating its own equation
along the invariants, and the original system is reconstructed from the result.

Run with::

    python examples/quadrature_recovery.py
"""

import matplotlib.pyplot as plt
import numpy as np
import sympy
from scipy.integrate import solve_ivp

from desr.numerics import NumericTranslation
from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslation

# Slots 1, 2 and 3 of the reference palette.  Identity is carried by hue and repeated as a
# direct label, never by colour alone.
FIRST = '#2a78d6'
SECOND = '#eb6834'
AUXILIARY = '#1baf7a'
INK = '#0b0b0b'
INK_SOFT = '#52514e'
GRID = '#dcdbd7'

EQUATIONS = ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
             'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)']

system = ODESystem.from_equations(EQUATIONS)
system.reorder_variables(['t', 'z1', 'z2'])
translation = ODETranslation.from_ode_system(system)
# Move onto the basis Hubert & Labahn use on p. 504, so the figure matches the paper.
translation.multiplier_add_columns(2, 1, -1)
translation.multiplier_add_columns(0, 2, -1)
with_aux = translation.translate_general(system)
reduced = translation.translate_general(system, include_aux_vars=False)
numeric = NumericTranslation(system, translation, reduced)

t, z1, z2 = system.variables
x0 = with_aux.variables[1]
y0, y1 = numeric.invariant_variables

start_time, final_time = 1.0, 2.0
initial_state = [0.5, 0.2]

# The paper's exact solution of the reduced system, with the constants that match the
# starting point above.  Everything numerical below is held to this, not to another solve.
c1, c2, c3 = sympy.symbols('c1 c2 c3')
L = sympy.log(c1 - t) - sympy.log(t) + c2
exact = {x0: c3 / (t ** sympy.Rational(1, 3) * L ** sympy.Rational(2, 3)),
         y0: c1 / (t * L),
         y1: t / (c1 - t)}
constants = {c1: 21, c2: 168 - sympy.log(20), c3: 21 ** sympy.Rational(2, 3)}
z1_exact, z2_exact = translation.reverse_translate_general((t, exact[x0], exact[y0], exact[y1]))
exact_z1 = sympy.lambdify(t, z1_exact.subs(constants), modules='numpy')
exact_z2 = sympy.lambdify(t, z2_exact.subs(constants), modules='numpy')
exact_x0 = sympy.lambdify(t, exact[x0].subs(constants), modules='numpy')


def as_function(a_system):
    variables = list(a_system.non_constant_variables)
    rhs = [a_system.derivative_dict[x] for x in variables]
    return sympy.lambdify([a_system.indep_var, variables], rhs, modules='numpy')


times = np.linspace(start_time, final_time, 300)

start = numeric.forward({t: start_time, z1: initial_state[0], z2: initial_state[1]})
solution = solve_ivp(as_function(reduced), (start_time, final_time),
                     [float(start[y0]), float(start[y1])],
                     t_eval=times, rtol=1e-11, atol=1e-13, dense_output=True)

# This system has no constants, so the only thing we can know is a point value: z1 at the
# start.  That determines the auxiliary at that one point ...
at_start = numeric.auxiliaries_from(start, known_values={z1: initial_state[0]})

# ... and from there the series is rebuilt, integrating the auxiliary along the way.
reduced_values = {reduced.indep_var: solution.t, y0: solution.y[0], y1: solution.y[1]}
recovered = numeric.reverse(reduced_values, auxiliaries=at_start, invariants_at=solution.sol)

# The auxiliary was never solved for.  This is the quadrature that reverse() runs
# internally, asked for on its own so it can be drawn, against the paper's exact x.
quadrature_auxiliary = numeric.recover_auxiliaries(
    reduced_values, auxiliaries=at_start, invariants_at=solution.sol)[0]
true_auxiliary = exact_x0(times)


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

# What is actually solved: two equations instead of three.
left.plot(solution.t, solution.y[0], color=FIRST, linewidth=2)
left.plot(solution.t, solution.y[1], color=SECOND, linewidth=2)
label_end(left, solution.t, solution.y[0], r'$y_0$', FIRST)
label_end(left, solution.t, solution.y[1], r'$y_1$', SECOND)
style(left, 'What is solved\ntwo invariants, no auxiliary', r'$t$', 'invariant')

# What was dropped, put back.
marks = slice(None, None, 20)
middle.plot(times, quadrature_auxiliary, color=AUXILIARY, linewidth=2)
middle.plot(times[marks], true_auxiliary[marks], linestyle='none', marker='o',
            markersize=8, markerfacecolor='none', markeredgewidth=1.6, color=AUXILIARY)
label_end(middle, times, quadrature_auxiliary, r'$x_0$', AUXILIARY)
style(middle, 'What was dropped, put back\n$x_0$ by quadrature; circles are its exact value',
      r'$t$', 'auxiliary')

# The original system, reconstructed from both, against the paper's exact solution.
right.plot(times, exact_z1(times), color=FIRST, linewidth=2, label='exact (Hubert & Labahn)')
right.plot(times, exact_z2(times), color=SECOND, linewidth=2)
right.plot(times[marks], recovered[z1][marks], linestyle='none', marker='o', markersize=8,
           markerfacecolor='none', markeredgewidth=1.6, color=FIRST,
           label='recovered from reduced')
right.plot(times[marks], recovered[z2][marks], linestyle='none', marker='o', markersize=8,
           markerfacecolor='none', markeredgewidth=1.6, color=SECOND)
label_end(right, times, exact_z1(times), r'$z_1$', FIRST)
label_end(right, times, exact_z2(times), r'$z_2$', SECOND)
style(right, 'The original system\nrebuilt from the invariants and $x_0$', r'$t$', 'value')
legend = right.legend(frameon=False, fontsize=9, loc='upper left', labelcolor=INK_SOFT)
for handle in legend.legend_handles:
    handle.set_color(INK_SOFT)
    handle.set_markeredgecolor(INK_SOFT)

figure.tight_layout()


if __name__ == '__main__':
    print('reduced system:')
    print(reduced)
    print('growth rate of the auxiliary:', numeric.auxiliary_growth_rates())
    print('largest difference from the exact solution:')
    print('  z1: {:.2e}'.format(np.abs(recovered[z1] - exact_z1(times)).max()))
    print('  z2: {:.2e}'.format(np.abs(recovered[z2] - exact_z2(times)).max()))
    print('  x0: {:.2e}'.format(np.abs(quadrature_auxiliary - true_auxiliary).max()))
    plt.show()
