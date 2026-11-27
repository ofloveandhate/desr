.. _no-symmetry:

Systems with no scaling symmetry
==================================

Not every system has a scaling symmetry.  This page describes what `desr` reports when
your system has none.


When a system has none
------------------------

A scaling is a symmetry of a system when it leaves the map :math:`F = t\,z^{-1} \star G`
invariant.  Hubert & Labahn :cite:`Hubert2013c` collect the exponents of :math:`F` into a
matrix :math:`K`, with one row per symbol, and the condition becomes :math:`A K = 0`: the
rows of the scaling matrix :math:`A` lie in the left kernel of :math:`K`.

Their Proposition 5.1 computes that kernel.  With :math:`U K` in row Hermite normal
form, having exactly :math:`r` zero rows, :math:`A` is the last :math:`r` rows of
:math:`U`.  If :math:`K` has full row rank then :math:`U K` has no zero rows, so
:math:`r = 0` and :math:`A` has no rows.

Here is such a system:

    >>> system = ODESystem.from_equations(['dz/dt = z + 1 + t'])
    >>> system.exponent_matrix()
    Matrix([
    [ 1,  2, 1],
    [-1, -1, 0]])
    >>> system.exponent_matrix().rank()
    2

Two symbols and rank two, so the left kernel is trivial.  The three terms :math:`z`,
:math:`1` and :math:`t` cannot all be made to scale alike.


The scaling matrix has no rows.  It still has one column per symbol:

    >>> system.maximal_scaling_matrix().shape
    (0, 2)
    >>> system.scaling_symmetries()
    []

A structural dimension is a symbol's weight under the scaling, so with no scaling every
symbol has dimension 1:

    >>> system.dimensions_of()
    {t: 1, z: 1}

With nothing scaled, every symbol is invariant and there are no auxiliaries:

    >>> translator = ODETranslator(system)
    >>> translator.num_scaling_symmetries
    0
    >>> translator.invariant_monomials()
    Matrix([[t, z]])
    >>> translator.auxiliary_monomials()
    Matrix(0, 0, [])


Reducing by no symmetries returns the system:

    >>> reduced = translator.translate()
    >>> reduced
    dt/dt = 1
    dz/dt = t + z + 1
    >>> reduced.named_symbols
    {t: t, z: z}

Each symbol denotes itself.  An invariant that is one of the original symbols keeps that
symbol's name, since the reduction did nothing to it, and here that is every symbol.

The scheme that applies is the parameter scheme.  Its condition is that every auxiliary
is a monomial in the parameters, and with no auxiliaries that holds.

To have generated names throughout, set ``rename_unchanged_symbols`` in the naming
scheme.  The equations are the same ones in different letters:

    >>> renaming = ODETranslator(system, naming_scheme={'rename_unchanged_symbols': True})
    >>> renaming.translate()
    dtau/dtau = 1
    dnu0/dtau = nu0 + tau + 1
    >>> renaming.translate().named_symbols
    {tau: t, nu0: z}

:attr:`~desr.ode_translation.ODETranslator.num_scaling_symmetries` is the thing to check
before a reduction if it matters whether anything will be removed.


A scaling matrix of your own must have full row rank
------------------------------------------------------

This subsection is a very particular note, the casual user is not likely to
provide their own scaling matrix...

A matrix with no rows is not the same as a row of zeros.  A row of the scaling matrix is
one scaling parameter :math:`\lambda_i`, and its entries are the powers of
:math:`\lambda_i` that multiply each symbol.  A zero row is a parameter that multiplies
nothing.  It adds no scaling, but it would be counted as a symmetry, and a reduction
that quotients by a symmetry which is not there returns equations that are not true.

The same holds for any row that is a combination of the others.  Hubert & Labahn assume
a scaling matrix of full row rank throughout, from their section 2.1, so that the number
of rows is the number of independent scalings.  A translator built from a matrix that is
not of full row rank is refused:

    >>> ODETranslator(system, scaling_matrix=sympy.Matrix([[0, 0]]))
    Traceback (most recent call last):
        ...
    ValueError: The scaling matrix has 1 row but rank 0, and it must have full row rank: each row is taken as one scaling symmetry, and 1 of these adds no scaling the others do not already give.  It defines no scaling at all.  The scaling matrix of a system with no scaling symmetry has no rows: sympy.zeros(0, 2).  A translator built without a scaling matrix computes one of full row rank.

A matrix computed by :meth:`~desr.ode_system.ODESystem.maximal_scaling_matrix` always
has full row rank, by part (i) of the proposition, so this only concerns a matrix
written by hand.

See :ref:`idempotent` for the most common way to meet a system with no symmetry: by
reducing one.
