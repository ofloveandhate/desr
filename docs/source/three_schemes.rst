.. _three-schemes:

The three reduction schemes
=============================

`desr` implements three reductions from Hubert & Labahn, and
:meth:`~desr.ode_translation.ODETranslator.translate` picks between them.  They are not
three algorithms for the same job: each applies to a different shape of system, and
they disagree about what the reduced system's clock *is*.  That last point is the one
worth reading carefully, because the same letter can mean two different quantities.

===================== ================ ============================================
Scheme                Hubert & Labahn  Applies when
===================== ================ ============================================
``PARAMETER``         Section 7        every auxiliary is a monomial in the
                                       parameters alone
``DEP_VAR``           Section 6.1      the scaling leaves the independent variable
                                       alone
``GENERAL``           Section 6.2      always
===================== ================ ============================================

:meth:`~desr.ode_translation.ODETranslator.translate` tries them in that order and uses
the first that fits.  :meth:`~desr.ode_translation.ODETranslator.scheme` says which one
that is, without running it.


The parameter scheme
----------------------

The best case, and the one that deserves the word *nondimensionalisation*.  Every
auxiliary is constant, so every one can be normalised, and what is left has fewer
parameters and no auxiliaries at all:

    >>> mm = ODESystem.from_equations(
    ...     ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
    ...      'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
    >>> ODETranslator.from_ode_system(mm).scheme()
    ReductionScheme.PARAMETER
    >>> ODETranslator.from_ode_system(mm).translate(mm)
    dtau/dtau = 1
    dnu0/dtau = -kappa0*nu0 - kappa1*nu0 - nu0*nu1 + nu1
    dnu1/dtau = kappa1*nu0 + nu0*nu1 - nu1
    parameters: kappa0, kappa1

Four parameters became two.  The clock is called :math:`\tau` because it is not
:math:`t`:

    >>> reduced = ODETranslator.from_ode_system(mm).translate(mm)
    >>> reduced.named_symbols[reduced.independent_variable]
    e_0*k_1*t

Time is now counted in units of :math:`1/(e_0 k_1)`.  Calling that :math:`t` would make
the clock the one symbol in the output whose name lied, so this scheme has always
renamed it.


The dependent-variable scheme
-------------------------------

This one applies when the scaling does not move the independent variable — when
:math:`A`'s column for it is zero.  The clock is then invariant exactly as it stands,
so the reduction keeps it rather than coining a symbol for it.  Hubert & Labahn 6.4:

    >>> hl64 = ODESystem.from_equations(['dz1/dt = z1*(1+z1*z2)',
    ...                                  'dz2/dt = z2*(1/t - z1*z2)'])
    >>> translator = ODETranslator.from_ode_system(hl64)
    >>> translator.scheme()
    ReductionScheme.DEP_VAR
    >>> translator.scaling_matrix[:, hl64.independent_variable_index]
    Matrix([[0]])
    >>> translator.translate(hl64)
    dt/dt = 1
    dx0/dt = x0*(-y0 - 1)
    dy0/dt = y0*(1 + 1/t)

The zero column is the whole criterion, and the consequence is worth stating plainly:
**this reduced system's** :math:`t` **is the original** :math:`t`.  Not a rescaled time
that happens to be spelled the same — the same quantity, whatever else the reduction
normalises:

    >>> translator.substitutions()[hl64.independent_variable]
    t

So a solution of the reduced system is already on the original clock, and a time series
needs no conversion.  That is the thing the dependent-variable scheme is *for*, and it
is easy to miss, because Hubert & Labahn do not set it beside the general scheme and
say so.


The general scheme
--------------------

Always applicable, and the fallback.  It carries the auxiliaries as variables of their
own — they are the quadratures you integrate afterwards to get back — and keeps the
system's own independent variable to integrate against.

Hubert & Labahn 6.6 is a general reduction that keeps :math:`t`:

    >>> hl66 = ODESystem.from_equations(
    ...     ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
    ...      'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])
    >>> ODETranslator.from_ode_system(hl66).scheme()
    ReductionScheme.GENERAL
    >>> ODETranslator.from_ode_system(hl66).translate(hl66)
    dt/dt = 1
    dx0/dt = x0*(2*y0*y1/3 - 1/3)/t
    dy0/dt = y0*(y0*y1 - 1)/t
    dy1/dt = y1*(y1 + 1)/t

Note the scheme is ``GENERAL``, not ``DEP_VAR``: the scaling *does* move :math:`t`
here, since :math:`A`'s column for it is :math:`3`, not zero.  The clock survives as
itself for a different reason, given below.


The clock, and when it is still ``t``
---------------------------------------

The general scheme runs two further steps by default, and both can be turned off.
Turning both off reproduces Hubert & Labahn exactly.

**Normalising** sets every constant auxiliary to 1, which is the choice of units
described in :ref:`dimensions`.

**Identifying** is the one that touches the clock.  Once the units are chosen, an
invariant such as :math:`k_{23} t` may *be* :math:`t` — and then it is a second name
for a variable the reduced system already carries, whose equation reads
:math:`dt/dt = 1`, which the system already has.  The two become one quantity.

From that moment the reduced system is integrated against :math:`k_{23} t` rather than
:math:`t`, so it is a rescaled time and takes the name a rescaled time gets:

    >>> enzyme = ChemicalReactionNetwork.from_diagram(
    ...     'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
    >>> translator = ODETranslator.from_ode_system(enzyme)
    >>> reduced = translator.translate(enzyme)
    >>> reduced.independent_variable
    tau
    >>> reduced.named_symbols[reduced.independent_variable]
    k_2_3*t

Ask for no normalisation and nothing becomes :math:`t`, so nothing is identified, and
the clock is the system's own — with :math:`k_{23} t` standing beside it as an
invariant with a name of its own:

    >>> untouched = translator.translate_general(enzyme, normalise_parameters=False)
    >>> untouched.independent_variable
    t
    >>> untouched.named_symbols[untouched.independent_variable]
    t

This is why the rule is about denotation and not about which scheme ran:

    **The reduced clock is called** :math:`t` **only when it denotes** :math:`t`.

Under the dependent-variable scheme that is always, because its multiplier is written
in the symbols other than the clock, so none of its invariants can turn out to be one.
Under the general scheme it holds whenever nothing was identified — which includes
every reduction that normalises nothing, every one run with
``identify_indep_var=False``, and systems like 6.6 whose clock invariants all involve a
dependent variable.  Under the parameter scheme it never holds, and that scheme has
always said so.

So if you are moving numbers between a reduced system and the original, the question to
ask is not which scheme ran but what the clock denotes.
:attr:`~desr.ode_system.ODESystem.named_symbols` answers it for every symbol including
that one, and :mod:`desr.numerics` asks it rather than trusting the spelling.
