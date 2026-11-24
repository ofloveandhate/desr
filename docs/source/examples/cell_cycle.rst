.. _cell_cycle:

A cell cycle control network
==============================

A model of the control of entry into mitosis, with seven species and thirty
parameters, from *Mathematical modeling as a tool for investigating cell cycle control
networks* (Tyson and co-author, 2007).  It is the largest model among the examples, it
is read directly from its LaTeX, and its reduction removes seven parameters.

The model
-----------

The equations, as transcribed, say this.  Cyclin is made at the constant rate
:math:`k_1`, is degraded at the rate :math:`k_2`, and binds Cdk at the rate
:math:`k_3` to form MPF.  MPF and preMPF convert into one another at the rates
:math:`k_{\mathrm{wee}}` and :math:`k_{25}`.  Four further species, Cdc25P, Wee1P, IEP
and APC, each follow an equation of one form: a Michaelis–Menten rate forward, driven
by MPF, and one backward, driven by a phosphatase.  Three of them feed back:
:math:`k_{25}` depends on Cdc25P, :math:`k_{\mathrm{wee}}` on Wee1P, and :math:`k_2`
on APC.

The system is seven differential equations followed by four definitions.  The
definitions give :math:`[\mathrm{Cdk}]`, :math:`k_{25}`, :math:`k_{\mathrm{wee}}` and
:math:`k_2` in terms of species:

    >>> tex = r'''\frac{d}{dt}\textrm{[Cyclin]} &= k_1 - k_2 \textrm{[Cyclin]} - k_3 \textrm{[Cyclin]} \textrm{[Cdk]} \\
    ... \frac{d}{dt}\textrm{[MPF]} &= k_3 \textrm{[Cyclin]} \textrm{[Cdk]} - k_2 \textrm{[MPF]} - k_{\textrm{wee}} \textrm{[MPF]} + k_{25} \textrm{[preMPF]} \\
    ... \frac{d}{dt}\textrm{[preMPF]} &= -k_2 \textrm{[preMPF]} + k_{\textrm{wee}} \textrm{[MPF]} - k_{25} \textrm{[preMPF]} \\
    ... \frac{d}{dt}\textrm{[Cdc25P]} &= \frac{k_a \textrm{[MPF]}\left( \textrm{[total Cdc25]} - \textrm{[Cdc25P]} \right)}{K_a + \textrm{[total Cdc25]} - \textrm{[Cdc25P]}} - \frac{k_b \textrm{[PPase]} \textrm{[Cdc25P]}}{K_b + \textrm{[Cdc25P]}} \\
    ... \frac{d}{dt}\textrm{[Wee1P]} &= \frac{k_e \textrm{[MPF]}\left( \textrm{[total Wee1]} - \textrm{[Wee1P]} \right)}{K_e + \textrm{[total Wee1]} - \textrm{[Wee1P]}} - \frac{k_f \textrm{[PPase]} \textrm{[Wee1P]}}{K_f + \textrm{[Wee1P]}} \\
    ... \frac{d}{dt}\textrm{[IEP]} &= \frac{k_g \textrm{[MPF]}\left( \textrm{[total IE]} - \textrm{[IEP]} \right)}{K_g + \textrm{[total IE]} - \textrm{[IEP]}} - \frac{k_h \textrm{[PPase]} \textrm{[IEP]}}{K_h + \textrm{[IEP]}} \\
    ... \frac{d}{dt}\textrm{[APC]} &= \frac{k_c \textrm{[MPF]}\left( \textrm{[total APC]} - \textrm{[APC]} \right)}{K_c + \textrm{[total APC]} - \textrm{[APC]}} - \frac{k_d \textrm{[PPase]} \textrm{[APC]}}{K_d + \textrm{[APC]}} \\
    ... \textrm{[Cdk]} &= \textrm{[total Cdk]} - \textrm{[MPF]} - \textrm{[preMPF]} \\
    ... k_{25} &= V_{25}' \left( \textrm{[total Cdc25]} - \textrm{[Cdc25P]} \right) + V_{25}'' \textrm{[Cdc25P]} \\
    ... k_{\textrm{wee}} &= V_{\textrm{wee}}' \textrm{[Wee1P]} + V_{\textrm{wee}}'' \left( \textrm{[total Wee1]} - \textrm{[Wee1P]} \right) \\
    ... k_2 &= V_2' \left( \textrm{[total APC]} - \textrm{[APC]} \right) + V_2'' \textrm{[APC]}'''

:meth:`~desr.ode_system.ODESystem.from_tex` reads it.  The four definitions involve
species that change in time, so they are not parameters and not constraints;
``sub_defns=True`` substitutes them into the equations.  A prime is read as part of a
name, so :math:`V_{25}'` and :math:`V_{25}''` are two parameters.

    >>> system = ODESystem.from_tex(tex, sub_defns=True)
    >>> system.dependents
    (APC, Cdc25P, Cyclin, IEP, MPF, Wee1P, preMPF)
    >>> len(system.parameters)
    30

Seven scaling symmetries
--------------------------

    >>> translator = ODETranslator.from_ode_system(system)
    >>> translator.num_scaling_symmetries
    7

The units the reduction infers show what they are.  For the species:

    >>> dimensions = system.dimensions_of()
    >>> for species in system.dependents:
    ...     print(species, dimensions[species])
    APC [K_c]
    Cdc25P [K_a]
    Cyclin [totalCdk]
    IEP [K_g]
    MPF [totalCdk]
    Wee1P [K_e]
    preMPF [totalCdk]

Five units of concentration.  Cyclin, MPF and preMPF share one, because they turn into
one another.  Each of the other four species has its own: it appears only inside its
own Michaelis–Menten terms and, for three of them, in one of the rates
:math:`k_{25}`, :math:`k_{\mathrm{wee}}`, :math:`k_2`, where it is multiplied by a
parameter.  So its amount can be measured in any unit if that parameter is adjusted to
match.

The sixth unit is time.  The seventh is not a unit of anything measured:

    >>> PPase, k_b, k_h = sympy.symbols('PPase k_b k_h')
    >>> dimensions[k_h], dimensions[PPase] * dimensions[k_h]
    ([k_h], [K_g]/[t])

The phosphatase concentration appears in the model only multiplied by one of the four
rate constants :math:`k_b`, :math:`k_d`, :math:`k_f`, :math:`k_h`.  Scaling it up and
those four down by the same factor changes nothing, so the model cannot tell the
phosphatase's amount from its activity.  The reduction finds that as a symmetry.

The reduction
---------------

    >>> reduced = translator.translate(system)
    >>> translator.scheme()
    ReductionScheme.PARAMETER
    >>> len(reduced.parameters), len(reduced.dependents)
    (23, 7)

Seven symmetries, and seven parameters fewer: twenty-three in place of thirty.  The
clock is rescaled time, and each species is measured in a unit built from the
parameters:

    >>> names = reduced.named_symbols
    >>> names[reduced.independent_variable]
    PPase*k_b*t/K_a
    >>> for symbol in reduced.dependents:
    ...     print(symbol, '=', names[symbol])
    nu0 = APC/K_c
    nu1 = Cdc25P/K_a
    nu2 = Cyclin*PPase*k_b/(K_a*k_1)
    nu3 = IEP/K_g
    nu4 = MPF*PPase*k_b/(K_a*k_1)
    nu5 = Wee1P/K_e
    nu6 = PPase*k_b*preMPF/(K_a*k_1)

Each of the four is measured in its own Michaelis constant.  Four of the remaining
parameters are then ratios of Michaelis constants, and five are the total amounts,
each in the unit of its species:

    >>> for name in ('kappa0', 'kappa1', 'kappa2', 'kappa3'):
    ...     print(name, '=', names[sympy.Symbol(name)])
    kappa0 = K_b/K_a
    kappa1 = K_d/K_c
    kappa2 = K_f/K_e
    kappa3 = K_h/K_g
    >>> totals = [symbol for symbol in reduced.parameters
    ...           if any(str(s).startswith('total') for s in names[symbol].free_symbols)]
    >>> for symbol in totals:
    ...     print(symbol, '=', names[symbol])
    kappa12 = totalAPC/K_c
    kappa13 = totalCdc25/K_a
    kappa14 = totalIE/K_g
    kappa15 = totalWee1/K_e
    kappa22 = PPase*k_b*totalCdk/(K_a*k_1)

One equation of the reduced system, for the species that was cyclin:

    >>> Cyclin = sympy.Symbol('Cyclin')
    >>> cyclin = next(s for s in reduced.dependents if names[s].has(Cyclin))
    >>> reduced.derivative_dict[cyclin]
    -kappa12*kappa18*nu2/kappa5 + kappa18*nu0*nu2/kappa5 - kappa19*nu0*nu2/kappa5 - kappa22*nu2/kappa5 + 1 + nu2*nu4/kappa5 + nu2*nu6/kappa5

The whole reduced system can be printed with ``reduced`` or typeset with
:meth:`~desr.ode_system.ODESystem.to_tex`.

The order of the variables
----------------------------

Which parameters are used as units depends on the order of the variables; see
:meth:`~desr.ode_system.ODESystem.reorder_variables`.  The number of symmetries, and so
the number of parameters removed, does not.  Here the parameters are put in the
reverse of their default order:

    >>> reordered = ODESystem.from_tex(tex, sub_defns=True)
    >>> order = sorted(reordered.parameters, key=str, reverse=True)
    >>> order += [reordered.independent_variable] + list(reordered.dependents)
    >>> reordered.reorder_variables(order)
    >>> other = ODETranslator.from_ode_system(reordered)
    >>> other.num_scaling_symmetries
    7
    >>> len(other.translate(reordered).parameters)
    23
    >>> set(other.auxiliary_monomials()) == set(translator.auxiliary_monomials())
    False

Seven symmetries and twenty-three parameters again, with a different choice of units.
