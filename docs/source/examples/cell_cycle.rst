.. _cell_cycle:

A cell cycle control network
==============================

The model of the cell cycle machinery regulating mitosis in frog egg extracts that
Sible and Tyson use as their example in :cite:`sible-tyson`.  The equations are their
Figure 2.  With seven species and thirty parameters it is the largest model among the
examples; it is read directly from its LaTeX, and its reduction removes seven
parameters.

The model
-----------

Cyclin is made at the constant rate :math:`k_1`, is degraded at the rate
:math:`k_2`, and binds Cdk at the rate :math:`k_3` to form MPF.  MPF and preMPF convert
into one another at the rates :math:`k_{\mathrm{wee}}` and :math:`k_{25}`.  Four
further species, Cdc25P, Wee1P, IEP and APC, each follow an equation of one form: a
Michaelis–Menten rate forward and one backward, the backward one driven by a
phosphatase.  MPF drives the forward rate of the first three.  IEP, the phosphorylated
form of an intermediate enzyme, drives APC's.  Three of the four feed back on the
first three equations: :math:`k_{25}` depends on Cdc25P, :math:`k_{\mathrm{wee}}` on
Wee1P, and :math:`k_2`, the rate of cyclin degradation, on APC.

The system is seven differential equations followed by four definitions.  The
definitions give :math:`[\mathrm{Cdk}]`, :math:`k_{25}`, :math:`k_{\mathrm{wee}}` and
:math:`k_2` in terms of species:

    >>> tex = r'''\frac{d}{dt}\textrm{[Cyclin]} &= k_1 - k_2 \textrm{[Cyclin]} - k_3 \textrm{[Cyclin]} \textrm{[Cdk]} \\
    ... \frac{d}{dt}\textrm{[MPF]} &= k_3 \textrm{[Cyclin]} \textrm{[Cdk]} - k_2 \textrm{[MPF]} - k_{\textrm{wee}} \textrm{[MPF]} + k_{25} \textrm{[preMPF]} \\
    ... \frac{d}{dt}\textrm{[preMPF]} &= -k_2 \textrm{[preMPF]} + k_{\textrm{wee}} \textrm{[MPF]} - k_{25} \textrm{[preMPF]} \\
    ... \frac{d}{dt}\textrm{[Cdc25P]} &= \frac{k_a \textrm{[MPF]}\left( \textrm{[total Cdc25]} - \textrm{[Cdc25P]} \right)}{K_a + \textrm{[total Cdc25]} - \textrm{[Cdc25P]}} - \frac{k_b \textrm{[PPase]} \textrm{[Cdc25P]}}{K_b + \textrm{[Cdc25P]}} \\
    ... \frac{d}{dt}\textrm{[Wee1P]} &= \frac{k_e \textrm{[MPF]}\left( \textrm{[total Wee1]} - \textrm{[Wee1P]} \right)}{K_e + \textrm{[total Wee1]} - \textrm{[Wee1P]}} - \frac{k_f \textrm{[PPase]} \textrm{[Wee1P]}}{K_f + \textrm{[Wee1P]}} \\
    ... \frac{d}{dt}\textrm{[IEP]} &= \frac{k_g \textrm{[MPF]}\left( \textrm{[total IE]} - \textrm{[IEP]} \right)}{K_g + \textrm{[total IE]} - \textrm{[IEP]}} - \frac{k_h \textrm{[PPase]} \textrm{[IEP]}}{K_h + \textrm{[IEP]}} \\
    ... \frac{d}{dt}\textrm{[APC]} &= \frac{k_c \textrm{[IEP]}\left( \textrm{[total APC]} - \textrm{[APC]} \right)}{K_c + \textrm{[total APC]} - \textrm{[APC]}} - \frac{k_d \textrm{[PPase]} \textrm{[APC]}}{K_d + \textrm{[APC]}} \\
    ... \textrm{[Cdk]} &= \textrm{[total Cdk]} - \textrm{[MPF]} - \textrm{[preMPF]} \\
    ... k_{25} &= V_{25}' \left( \textrm{[total Cdc25]} - \textrm{[Cdc25P]} \right) + V_{25}'' \textrm{[Cdc25P]} \\
    ... k_{\textrm{wee}} &= V_{\textrm{wee}}' \textrm{[Wee1P]} + V_{\textrm{wee}}'' \left( \textrm{[total Wee1]} - \textrm{[Wee1P]} \right) \\
    ... k_2 &= V_2' \left( \textrm{[total APC]} - \textrm{[APC]} \right) + V_2'' \textrm{[APC]}'''

:meth:`~desr.ode_system.ODESystem.from_tex` reads it.  The four definitions involve
species that change in time, so they are not parameters and not constraints;
``sub_defns=True`` substitutes them into the equations.  A prime is read as part of a
name, so :math:`V_{25}'` and :math:`V_{25}''` are two parameters.

    >>> system = ODESystem.from_tex(tex, sub_defns=True)
    >>> len(system.dependents), len(system.parameters)
    (7, 30)

The order of the variables
----------------------------

The reduction uses as units the parameters at the front of the variable order; see
:meth:`~desr.ode_system.ODESystem.reorder_variables`.  The order is parameters, then
the independent variable, then the dependent variables.

Here the five total amounts come first, so that each species is measured in its total,
followed by :math:`V_2''` and :math:`k_h`.  The other twenty-three parameters follow in
a fixed order, which is the order their reduced counterparts are numbered in.

    >>> used_as_units = ['totalWee1', 'totalCdk', 'totalAPC', 'totalIE', 'totalCdc25',
    ...                  'V_2primeprime', 'k_h']
    >>> remaining_parameters = ['PPase',
    ...                         'K_a', 'K_b', 'K_c', 'K_d', 'K_e', 'K_f', 'K_g', 'K_h',
    ...                         'k_1', 'k_3',
    ...                         'k_a', 'k_b', 'k_c', 'k_d', 'k_e', 'k_f', 'k_g',
    ...                         'V_25prime', 'V_25primeprime',
    ...                         'V_weeprime', 'V_weeprimeprime', 'V_2prime']
    >>> species = ['Cyclin', 'MPF', 'preMPF', 'Cdc25P', 'Wee1P', 'IEP', 'APC']
    >>> system.reorder_variables(used_as_units + remaining_parameters + ['t'] + species)

Seven scaling symmetries
--------------------------

    >>> translator = ODETranslator.from_ode_system(system)
    >>> translator.num_scaling_symmetries
    7

The units the reduction infers show what they are.  For the species:

    >>> dimensions = system.dimensions_of()
    >>> for symbol in system.dependents:
    ...     print(symbol, dimensions[symbol])
    Cyclin [totalCdk]
    MPF [totalCdk]
    preMPF [totalCdk]
    Cdc25P [totalCdc25]
    Wee1P [totalWee1]
    IEP [totalIE]
    APC [totalAPC]

Five units of concentration.  Cyclin, MPF and preMPF share one, because they turn into
one another.  Each of the other four species has its own.  Outside its own
Michaelis–Menten terms it appears in one place only: Cdc25P, Wee1P and APC in one of
the rates :math:`k_{25}`, :math:`k_{\mathrm{wee}}`, :math:`k_2`, and IEP in APC's
equation.  In each of those places it is multiplied by a parameter, so its amount can
be measured in any unit if that parameter is adjusted to match.

The sixth unit is time.  The seventh is not a unit of anything measured.  It belongs to
the phosphatase:

    >>> PPase = sympy.Symbol('PPase')
    >>> dimensions[PPase]
    [totalWee1]/([k_f]*[t])
    >>> for name in ('k_b', 'k_d', 'k_f', 'k_h'):
    ...     print(name, dimensions[PPase] * dimensions[sympy.Symbol(name)])
    k_b [totalCdc25]/[t]
    k_d [totalAPC]/[t]
    k_f [totalWee1]/[t]
    k_h [totalIE]/[t]

The phosphatase concentration appears in the model only multiplied by one of the four
rate constants :math:`k_b`, :math:`k_d`, :math:`k_f`, :math:`k_h`.  Each product is a
concentration per unit time; the phosphatase alone carries a unit, written
:math:`[k_f]` above, that cancels in every product.  Scaling the phosphatase up and
those four down by the same factor changes nothing, so the model cannot tell the
phosphatase's amount from its activity.  The reduction finds that as a symmetry.

Sible and Tyson's own parameter file, their Table 1, sets five parameters to 1 under
the comment "These parameters are arbitrarily set to 1": the four totals
:math:`[\mathrm{total\ Cdc25}]`, :math:`[\mathrm{total\ Wee1}]`,
:math:`[\mathrm{total\ IE}]`, :math:`[\mathrm{total\ APC}]`, and the phosphatase.
Those are five of the seven symmetries used by hand: a unit for each of the four
species, and the phosphatase's amount.  The order chosen above makes the same choice
for the four species.  For the phosphatase's symmetry it sets :math:`k_h` to 1 where
the table sets the phosphatase.

The reduction
---------------

    >>> reduced = translator.translate(system)
    >>> translator.scheme()
    ReductionScheme.PARAMETER
    >>> len(reduced.parameters), len(reduced.dependents)
    (23, 7)

Seven symmetries, and seven parameters fewer: twenty-three in place of thirty.  The
clock is rescaled time, and each species is measured in its total:

    >>> names = reduced.named_symbols
    >>> names[reduced.independent_variable]
    V_2primeprime*t*totalAPC
    >>> for symbol in reduced.dependents:
    ...     print(symbol, '=', names[symbol])
    nu0 = Cyclin/totalCdk
    nu1 = MPF/totalCdk
    nu2 = preMPF/totalCdk
    nu3 = Cdc25P/totalCdc25
    nu4 = Wee1P/totalWee1
    nu5 = IEP/totalIE
    nu6 = APC/totalAPC

Among the twenty-three parameters, the first is what is left of the phosphatase, and
the next eight are the Michaelis constants, each as a fraction of the total of its
species:

    >>> for index in range(9):
    ...     symbol = sympy.Symbol('kappa{}'.format(index))
    ...     print(symbol, '=', names[symbol])
    kappa0 = PPase*k_h/(V_2primeprime*totalAPC*totalIE)
    kappa1 = K_a/totalCdc25
    kappa2 = K_b/totalCdc25
    kappa3 = K_c/totalAPC
    kappa4 = K_d/totalAPC
    kappa5 = K_e/totalWee1
    kappa6 = K_f/totalWee1
    kappa7 = K_g/totalIE
    kappa8 = K_h/totalIE

The equation of the reduced system for cyclin:

    >>> reduced.derivative_dict[sympy.Symbol('nu0')]
    kappa10*nu0*nu1 + kappa10*nu0*nu2 - kappa10*nu0 + kappa22*nu0*nu6 - kappa22*nu0 + kappa9 - nu0*nu6
    >>> for name in ('kappa9', 'kappa10', 'kappa22'):
    ...     print(name, '=', names[sympy.Symbol(name)])
    kappa9 = k_1/(V_2primeprime*totalAPC*totalCdk)
    kappa10 = k_3*totalCdk/(V_2primeprime*totalAPC)
    kappa22 = V_2prime/V_2primeprime

It is the original equation with :math:`1 - \nu_1 - \nu_2` for the free Cdk and
:math:`\kappa_{22} (1 - \nu_6) + \nu_6` for the rate of degradation.  The whole
reduced system can be printed with ``reduced`` or typeset with
:meth:`~desr.ode_system.ODESystem.to_tex`.

Another order
---------------

The number of symmetries, and so the number of parameters removed, does not depend on
the order.  The choice of units does.  Left in its default order, which is
alphabetical, the system reduces to twenty-three parameters as well:

    >>> default = ODESystem.from_tex(tex, sub_defns=True)
    >>> other = ODETranslator.from_ode_system(default)
    >>> other.num_scaling_symmetries
    7
    >>> otherwise = other.translate(default)
    >>> len(otherwise.parameters)
    23

The four species that have a unit of their own are then measured in Michaelis
constants, which come before the totals in the alphabet:

    >>> images = otherwise.named_symbols
    >>> for symbol in otherwise.dependents:
    ...     if len(images[symbol].free_symbols) == 2:
    ...         print(symbol, '=', images[symbol])
    nu0 = APC/K_c
    nu1 = Cdc25P/K_a
    nu3 = IEP/K_g
    nu5 = Wee1P/K_e
