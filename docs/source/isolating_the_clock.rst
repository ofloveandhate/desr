.. _isolating_the_clock:

Isolating the clock
=====================

A reduction writes a system in its scaling invariants.  The invariants a translator
lists are one basis among many (infinitely many),
and the choice of basis affects how the reduced system reads.
This page is about a choice the user has in desr:
whether the independent variable appears in one
invariant or in several.

Example system
----------------

Consider the autonomous ordinary differential equation

.. math::

    \frac{du}{dt} = 2 b^2 u + \frac{1}{a^2}.


In desr, the way the user orders their symbols determines the form of the reduce
system.  With symbols in the order
:math:`(t, u, a, b)` the scaling has two symmetries, and the normal Hermite multiplier
gives these invariants:

    >>> system = ODESystem.from_equations(['du/dt = 2*b**2*u + 1/a**2'])
    >>> system.reorder_variables(['t', 'u', 'a', 'b'])
    >>> translator = ODETranslator(system)
    >>> translator.r
    2
    >>> translator.invariant_monomials()
    Matrix([[a**2*u/t, b**2*t]])

Both invariants contain :math:`t`.  The second, :math:`b^2 t`, is the independent variable times a
parameter, so the reduction takes it as the clock :math:`\tau`.  The first is
:math:`y_0 = a^2 u / t`, and its equation has :math:`\tau` in it:

    >>> translator.translate()
    dtau/dtau = 1
    dy0/dtau = 2*y0 - y0/tau + 1/tau

This is a correct reduction.  Sadly the resulting system is not autonomous, although the original is.
The :math:`1/\tau` terms come from the :math:`1/t` inside :math:`y_0`: differentiating
:math:`a^2 u / t` differentiates the :math:`1/t` as well.

Another basis of the same invariants
---------------------------------------

By choosing a different set of invariants, one might obtain an autonomous system.
Observe first, the product of two invariants is an invariant.  For example

.. math::

    y_0 \cdot \tau = \frac{a^2 u}{t} \cdot b^2 t = a^2 b^2 u ,

which has no :math:`t`.  The pair :math:`(a^2 b^2 u,\ b^2 t)` generates the same
invariants as the pair :math:`(a^2 u / t,\ b^2 t)`: each can be written as a product of
integer powers of the other.  In this basis only one invariant contains :math:`t`.

:meth:`~desr.ode_translation.ODETranslator.isolate_clock` puts a translator in that
basis.  It changes the translator it is called on, as the ``multiplier_*`` operations
do, so a second translator is built here to keep the first for comparison:

    >>> translator_with_isolate_enabled = ODETranslator(system)
    >>> translator_with_isolate_enabled.isolate_clock() # internally manipulates the Hermite multiplier
    >>> translator_with_isolate_enabled.invariant_monomials()
    Matrix([[a**2*b**2*u, b**2*t]])
    >>> translator_with_isolate_enabled.translate()
    dtau/dtau = 1
    dy0/dtau = 2*y0 + 1

Writing :math:`\nu = a^2 b^2 u` and :math:`\tau = b^2 t`, this is
:math:`d\nu/d\tau = 2\nu + 1`, which is autonomous.

Everything the translator reports afterwards is in the new basis: its invariants, the
names it gives them, a reduction, a figure.  There is one multiplier, and it is the one
the reduction uses.

Isolating the clock is something you have to ask for.
A translator on which `isolate_clock` has not been called simply uses the normal
Hermite multiplier:

    >>> translator.invariant_monomials()
    Matrix([[a**2*u/t, b**2*t]])

How isolating the clock works
------------------------------

This part is pretty mathy, and gets into notation details and names of things.  I do
apologize for the technical details here.  Just if you're curious, read on.
Each column of the multiplier :math:`V` is the exponent vector of one monomial, and the
row belonging to :math:`t` holds the power of :math:`t` in each.  For the example:

    >>> translator.herm_mult_n
    Matrix([
    [-1, 1],
    [ 1, 0],
    [ 2, 0],
    [ 0, 2]])

The first row :math:`(-1, 1)` corresponds to :math:`t`.

Adding an integer multiple of one invariant column to another multiplies one invariant
by a power of the other.  The result is still an invariant, and the multiplier is still
unimodular.  The clock-isolating operation uses such steps to reduce :math:`t`'s row to a single
nonzero entry:

1. Among the invariants that contain :math:`t`, take the one with the smallest power
   of :math:`t`.  Where two have the same power, take one with no dependent variable in
   it.
2. Subtract multiples of its column from the others, to reduce their powers of
   :math:`t`.
3. Repeat until one invariant contains :math:`t`.

This is Euclid's algorithm applied to the row.  The invariant left holding :math:`t`
is negated if its power is negative, and is then used in the same way to remove
:math:`t` from the auxiliaries.  Here one step is enough: the first column plus the
second.

    >>> translator_with_isolate_enabled.herm_mult_n
    Matrix([
    [0, 1],
    [1, 0],
    [2, 0],
    [2, 2]])

The scaling matrix is the same, and so is its Hermite form:

    >>> translator_with_isolate_enabled.scaling_matrix == translator.scaling_matrix
    True
    >>> translator_with_isolate_enabled.herm_form == translator.herm_form
    True

Will the result be autonomous?
-------------------------------

Suppose the original system is autonomous, and after the operation the invariant
holding :math:`t` is :math:`\tau = k t` with :math:`k` a monomial in the parameters.
The reduction takes :math:`\tau` as its clock.  Every other invariant, and every
auxiliary, is then free of :math:`t`, so nothing in the reduced system depends on the
clock, and the reduced system is autonomous.

The invariant holding :math:`t` may instead involve a dependent variable.  Hubert &
Labahn's example 6.6 is such a system.  Its invariants are :math:`t z_1^3` and
:math:`z_1^2 z_2 / t`; the operation replaces the second by the product of the two:

    >>> hl66 = ODESystem.from_equations(
    ...     ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
    ...      'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])
    >>> translator = ODETranslator(hl66)
    >>> translator.invariant_monomials()
    Matrix([[t*z1**3, z1**2*z2/t]])
    >>> translator.isolate_clock()
    >>> translator.invariant_monomials()
    Matrix([[t*z1**3, z1**5*z2]])

No invariant here is :math:`t` multiplied by parameters, so the reduction keeps :math:`t` as
its clock in either basis, and the system is not autonomous to begin with.

Notes
---------------------

Clock isolation works with the invariants the multiplier lists.  It does not look for an invariant
:math:`k t` that is a product of the listed ones without being one of them.

It does not change which scalings the system has, the auxiliaries' scaling, or the
number of symbols in the reduced system.

The invariants it produces are not the same as the normal Hermite multiplier gives, and so
not the same Hubert & Labahn print for their examples.  A translator as-built uses the
normal multiplier, and what
:meth:`~desr.ode_translation.ODETranslator.invariant_monomials` reports is always what
a reduction is written in.
