Extending invariants in four-variable Michaelis-Menten
=======================================================

Raw Michaelis-Menten Equation Analysis
---------------------------------------

We perform an example analysis of the Michael-Mentis equations :cite:`Segel1989`:

.. math::
    :nowrap:

    \begin{align}
    \frac{dE}{dt} &= - k_1 E S + k_{-1} C + k_2 C \\
    \frac{dS}{dt} &= - k_1 E S + k_{-1} C \\
    \frac{dC}{dt} &= k_1 E S - k_{-1} C - k_2 C \\
    \frac{dP}{dt} &= k_2 C
    \end{align}

First we initiate the system from LaTeX and find the maximal scaling matrix such that the system is invariant.
Note that negative subscripts are turned into 'm' so that they comply with :py:mod:`sympy`. The 'm's are turned back into
negatives when printint to LaTeX using :func:`desr.tex_tools`.

    >>> import sympy
    >>> from desr import ODESystem, ODETranslator, scale_action
    >>> system_tex = '''\frac{dE}{dt} &= - k_1 E S + k_{-1} C + k_2 C \\\\
    ...                 \frac{dS}{dt} &= - k_1 E S + k_{-1} C \\\\
    ...                 \frac{dC}{dt} &= k_1 E S - k_{-1} C - k_2 C \\\\
    ...                 \frac{dP}{dt} &= k_2 C'''
    >>> original_system = ODESystem.from_tex(system_tex)
    >>> original_system.symbols
    (k_1, k_2, k_m1, t, C, E, P, S)
    >>> max_scal1 = ODETranslator.from_ode_system(original_system, naming_scheme={'parameters': 'c', 'independent': 't', 'dependents': ['C','E','P','S']})
    >>> print('Variable order: ', max_scal1.symbols_domain)
    Variable order:  (k_1, k_2, k_m1, t, C, E, P, S)
    >>> print(f'Scaling Matrix:\n{max_scal1.scaling_matrix.__repr__()}')
    Scaling Matrix:
    Matrix([
    [-1, -1, -1, 1, 0, 0, 0, 0],
    [-1,  0,  0, 0, 1, 1, 1, 1]])

Inspect the scaling invariants:

    >>> print('Invariants: ', max_scal1.invariant_monomials())
    Invariants:  Matrix([[k_m1/k_2, k_2*t, C*k_1/k_2, E*k_1/k_2, P*k_1/k_2, S*k_1/k_2]])

Finding the reduced system is also easy.
Since the Hermite multiplier and inverse are compatible with the simplest parameter reduction scheme,
:meth:`~desr.ode_translation.ODETranslator.translate` will automatically perform this reduction.

    >>> print(f'Reduced system:\n{max_scal1.translate(original_system)}')
    Reduced system:
    dt/dt = 1
    dC/dt = -C*c0 - C + E*S
    dE/dt = C*c0 + C - E*S
    dP/dt = C
    dS/dt = C*c0 - E*S
    parameters: c0



Changing the variable order
---------------------------

Parameters come first in the variable order, and the leading ones are the ones the algorithm normalises.
In our previous example the parameter order was :math:`k_1, k_2, k_{-1}`, so :math:`k_1` and :math:`k_2` were
normalised and the one surviving parameter was the ratio :math:`k_{-1}/k_2`.
Instead, we can normalise by :math:`k_{-1}` and leave :math:`k_2` in the ratio, by swapping the last two parameters.
Note that we need to recalculate the :class:`~desr.ode_translation.ODETranslator` instance.

    >>> original_system_reorder = original_system.copy()
    >>> original_system_reorder.reorder_variables(['k_1', 'k_m1', 'k_2', 't', 'C', 'E', 'P', 'S'])
    >>> original_system_reorder.symbols
    (k_1, k_m1, k_2, t, C, E, P, S)
    >>> max_scal1_reorder = ODETranslator.from_ode_system(original_system_reorder, naming_scheme={'parameters': 'c', 'independent': 't', 'dependents': ['C','E','P','S']})
    >>> print('Invariants:', ', '.join(map(str, max_scal1_reorder.invariant_monomials())))
    Invariants: k_2/k_m1, k_m1*t, C*k_1/k_m1, E*k_1/k_m1, P*k_1/k_m1, S*k_1/k_m1

Now we can reduce to find another, equivalent system.

    >>> reduced_system = max_scal1_reorder.translate(original_system_reorder)
    >>> reduced_system
    dt/dt = 1
    dC/dt = -C*c0 - C + E*S
    dE/dt = C*c0 + C - E*S
    dP/dt = C*c0
    dS/dt = C - E*S
    parameters: c0

Extending a choice of invariants
--------------------------------

We return to our original variable order: :math:`k_1, k_2, k_{-1}, t, C, E, P, S`.

Suppose we wish to study the invariants :math:`\frac{k_1}{k_{2}}C` and :math:`\frac{k_1}{k_{-1}}P`.
We say so by writing them down:

    >>> invariant_choice = ['C*k_1/k_2', 'P*k_1/k_m1']

Internally these become a matrix of exponents, one invariant per column, read against the
variable order:

.. math::

    P = \left[\begin{matrix}1 & 1\\-1 & 0\\0 & -1\\0 & 0\\1 & 0\\0 & 0\\0 & 1\\0 & 0\end{matrix}\right].

That matrix may be passed instead, but it is not the thing to type: its *rows* are
positions in the variable order, so reordering the variables leaves the same matrix meaning
different monomials, and nothing in it records which invariants were intended.  Writing the
monomials says what is meant and survives a reordering.

    >>> scale_action(max_scal1.symbols_domain, sympy.Matrix([[1, -1, 0, 0, 1, 0, 0, 0],
    ...                                                      [1, 0, -1, 0, 0, 0, 1, 0]]).T)
    Matrix([[C*k_1/k_2, P*k_1/k_m1]])

Finding a maximal scaling matrix that can be used to rewrite the system in terms of these invariants is also simple.

    >>> max_scal2 = max_scal1.extend_from_invariants(invariant_choice=invariant_choice)
    >>> max_scal2
    A=
    Matrix([
    [-1, -1, -1, 1, 0, 0, 0, 0],
    [-1,  0,  0, 0, 1, 1, 1, 1]])
    V=
    Matrix([
    [ 0, -1,  1,  1, 0,  1,  1,  1],
    [-1,  1, -1,  0, 1, -1, -1, -1],
    [ 0,  0,  0, -1, 0,  0,  0,  0],
    [ 0,  0,  0,  0, 1,  0,  0,  0],
    [ 0,  0,  1,  0, 0,  0,  0,  0],
    [ 0,  0,  0,  0, 0,  1,  0,  0],
    [ 0,  0,  0,  1, 0,  0,  1,  0],
    [ 0,  0,  0,  0, 0,  0,  0,  1]])
    W=
    Matrix([
    [-1, -1, -1, 1, 0, 0, 0, 0],
    [-1,  0,  0, 0, 1, 1, 1, 1],
    [ 0,  0,  0, 0, 1, 0, 0, 0],
    [ 0,  0, -1, 0, 0, 0, 0, 0],
    [ 0,  0,  0, 1, 0, 0, 0, 0],
    [ 0,  0,  0, 0, 0, 1, 0, 0],
    [ 0,  0,  1, 0, 0, 0, 1, 0],
    [ 0,  0,  0, 0, 0, 0, 0, 1]])

For Python code that steps through this procedure, see :py:mod:`desr.examples.michaelis_menten`.

Now, this transformation doesn't satisfy the conditions of the parameter reduction scheme, so if we try to reduce it
:meth:`~desr.ode_translation.ODETranslator.translate` will use the dependent reduction scheme implemented in
:meth:`~desr.ode_translation.ODETranslator.translate_dep_var`.

    >>> max_scal2.invariant_monomials()
    Matrix([[C*k_1/k_2, P*k_1/k_m1, k_2*t, E*k_1/k_2, P*k_1/k_2, S*k_1/k_2]])
    >>> max_scal2.translate(original_system)
    dt/dt = 1
    dy0/dt = y0*(-y2 - y2*y4/y1 + y2*y3*y5/y0)/t
    dy1/dt = y0*y1*y2/(t*y4)
    dy2/dt = y2/t
    dy3/dt = y3*(y0*y2/y3 + y0*y2*y4/(y1*y3) - y2*y5)/t
    dy4/dt = y0*y2/t
    dy5/dt = y5*(y0*y2*y4/(y1*y5) - y2*y3)/t
    parameters: x0, x1

Here, :code:`x0` and :code:`x1` are auxiliary variables, which can be fixed at any value at all.
They are reported as *parameters* because each is a monomial in the original parameters alone
and so cannot vary; being auxiliaries is a separate fact, recorded in
:attr:`~desr.ode_system.ODESystem.auxiliary_names`.
The remaining variables :code:`(y0, y1, y2, y3, y4, y5) = (C*k_1/k_2, P*k_1/k_m1, k_2*t, E*k_1/k_2, P*k_1/k_2, S*k_1/k_2)`
are our new dependent invariants; the first two are exactly the invariants we chose to study.

However, we can rearrange these invariants into an order that *does* satisfy the parameter reduction scheme.
While this isn't implemented yet, we can do it by hand for the moment.

None of the six invariants above is a monomial in the parameters alone, so the parameter
scheme has nothing to reduce yet.  One column operation makes one: dividing
:math:`\frac{k_1}{k_{-1}} P` by :math:`\frac{k_1}{k_2} P` cancels the :math:`P` and leaves
:math:`k_2 / k_{-1}`.  Both column numbers come from
:meth:`~desr.ode_translation.ODETranslator.column_of`, which takes the invariant a column
carries and gives back its index.

    >>> ratio = max_scal2.column_of('P*k_1/k_m1')
    >>> product = max_scal2.column_of('P*k_1/k_2')
    >>> ratio, product
    (3, 6)
    >>> max_scal2.multiplier_add_columns(ratio, product, power=-1)
    >>> max_scal2.invariant_monomials()
    Matrix([[C*k_1/k_2, k_2/k_m1, k_2*t, E*k_1/k_2, P*k_1/k_2, S*k_1/k_2]])

We then reorder the last :math:`n-r` columns into the layout the scheme wants: the
parameter invariant first, then the independent variable, then the dependent invariants.
Which column is which is decided from what each one is.  A permutation written as a
literal would be read against the variable order, and would mean something different if
that order changed.

    >>> variables, indep = set(original_system.variables), original_system.independent_variable
    >>> def column_kind(monomial):
    ...     if indep in monomial.free_symbols:
    ...         return 1                                   # the independent variable
    ...     return 2 if monomial.free_symbols & variables else 0   # dependent, else parameter
    >>> monomials = list(max_scal2.invariant_monomials())
    >>> permutation = sorted(range(len(monomials)), key=lambda i: (column_kind(monomials[i]), i))
    >>> permutation
    [1, 2, 0, 3, 4, 5]

    >>> max_scal3 = max_scal2.herm_mult_n[:, permutation]
    >>> print(f'Permuted Vn:\n{repr(max_scal3)}')
    Permuted Vn:
    Matrix([
    [ 0, 0,  1,  1,  1,  1],
    [ 1, 1, -1, -1, -1, -1],
    [-1, 0,  0,  0,  0,  0],
    [ 0, 1,  0,  0,  0,  0],
    [ 0, 0,  1,  0,  0,  0],
    [ 0, 0,  0,  1,  0,  0],
    [ 0, 0,  0,  0,  1,  0],
    [ 0, 0,  0,  0,  0,  1]])
    >>> max_scal3 = sympy.Matrix.hstack(max_scal2.herm_mult_i, max_scal3)
    >>> max_scal3 = ODETranslator(max_scal2.scaling_matrix, symbols_domain=original_system.symbols, hermite_multiplier=max_scal3, naming_scheme={'parameters': 'c', 'independent': 't', 'dependents': ['C','E','P','S']})

:math:`V_i` is unchanged, since only the invariant block was touched.  ``symbols_domain``
tells the new translator which symbols the scaling acts on, so the invariants read in them:

    >>> max_scal3.invariant_monomials()
    Matrix([[k_2/k_m1, k_2*t, C*k_1/k_2, E*k_1/k_2, P*k_1/k_2, S*k_1/k_2]])

    >>> print(max_scal3.translate(original_system).subs({'c0':'alpha'}))
    dt/dt = 1
    dC/dt = -C - C/alpha + E*S
    dE/dt = C + C/alpha - E*S
    dP/dt = C
    dS/dt = C/alpha - E*S
    parameters: alpha

So we have found a third different reparametrization of the Michaelis-Menten equations.


Initial conditions, two of them zero
------------------------------------

The four-species form takes the initial conditions of :cite:`murray` (6.4) directly: all
the substrate and all the enzyme are present at the start, and none of the complex or the
product.

.. math::

    s(0) = s_0, \quad e(0) = e_0, \quad c(0) = 0, \quad p(0) = 0.

In `desr`:

    >>> system = ODESystem.from_tex(system_tex)
    >>> system.update_initial_conditions(['S(0) = S_0', 'E(0) = E_0', 'C(0) = 0', 'P(0) = 0'])
    >>> system.initial_conditions
    {S: S_0, E: E_0, C: C_0, P: P_0}

Each zero is recorded as the value of the symbol that carries the units, and the names are
inferred from the variables':

    >>> system.initial_condition_values
    {C: 0, P: 0}

:meth:`~desr.ode_system.ODESystem.reorder_parameters` moves only the parameter block, so
the layout the parameter scheme needs is preserved.  :math:`S_0` and :math:`k_1` go at the
front, which is the end that gets scaled by, and the two zero-valued symbols go at the
back, which is the end that is kept.

    >>> system.reorder_parameters(['S_0', 'k_1', 'k_2', 'k_m1', 'E_0', 'C_0', 'P_0'])
    >>> system.symbols
    (S_0, k_1, k_2, k_m1, E_0, C_0, P_0, t, C, E, P, S)

    >>> translator = ODETranslator.from_ode_system(system)
    >>> translator.scheme(system)
    ReductionScheme.PARAMETER

The reduction divides by these symbols, so none of them may be zero.  Neither :math:`C_0`
nor :math:`P_0` is among them: each appears only in a numerator, so the invariant holding
it is zero.

    >>> translator.symbols_that_must_not_be_zero()
    (S_0, k_1)

    >>> translator.invariant_monomials()
    Matrix([[k_2/(S_0*k_1), k_m1/(S_0*k_1), E_0/S_0, C_0/S_0, P_0/S_0, S_0*k_1*t, C/S_0, E/S_0, P/S_0, S/S_0]])

    >>> translator.translate(system)
    dtau/dtau = 1
    dnu0/dtau = -kappa0*nu0 - kappa1*nu0 + nu1*nu3
    dnu1/dtau = kappa0*nu0 + kappa1*nu0 - nu1*nu3
    dnu2/dtau = kappa0*nu0
    dnu3/dtau = kappa1*nu0 - nu1*nu3
    nu0(0) = kappa3
    nu1(0) = kappa2
    nu2(0) = kappa4
    nu3(0) = 1
    parameters: kappa0, kappa1, kappa2, kappa3, kappa4
    parameter values: kappa3 = 0, kappa4 = 0

:math:`\kappa_3` is :math:`C_0/S_0` and :math:`\kappa_4` is :math:`P_0/S_0`, so both are
zero.  Read as values, the substrate starts at one and the complex and the product at
nothing:

    >>> translator.translate(system).initial_condition_values
    {nu3: 1, nu0: 0, nu2: 0}

.. todo::

    Add a method to :class:`~desr.ode_translation.ODETranslator` that will try and re-order the last :math:`n-r` columns so
    that the parameter reduction scheme can be applied.