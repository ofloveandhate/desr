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
Then we must create a matrix representing these invariants:

.. math::

    P = \left[\begin{matrix}1 & 1\\-1 & 0\\0 & -1\\0 & 0\\1 & 0\\0 & 0\\0 & 1\\0 & 0\end{matrix}\right].

We can easily check we have correct matrix:

    >>> invariant_choice = sympy.Matrix([[1, -1, 0, 0, 1, 0, 0, 0],
    ...                                  [1, 0, -1, 0, 0, 0, 1, 0]]).T
    >>> scale_action(max_scal1.symbols_domain, invariant_choice)
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
While this isn't implemented yet, we can do it by hand for the moment. First we expose the parameter
:math:`k_2 / k_{-1}` as a single column, by subtracting the :math:`\frac{k_1}{k_2} P` column from the
:math:`\frac{k_1}{k_{-1}} P` column. We then reorder the last :math:`n-r` columns into the layout the
scheme wants: the parameter invariant first, then the independent variable, then the dependent invariants.

    >>> max_scal3 = max_scal2.herm_mult_n
    >>> max_scal3[:, 1] = max_scal3[:, 1] - max_scal3[:, 4]  # P*k_1/k_m1 - P*k_1/k_2 = k_2/k_m1
    >>> max_scal3 = max_scal3[:, [1, 2, 0, 3, 4, 5]]  # parameter first, then time and the dependents
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
    >>> max_scal3 = sympy.Matrix.hstack(max_scal1.herm_mult_i, max_scal3)
    >>> max_scal3 = ODETranslator(max_scal1.scaling_matrix, hermite_multiplier=max_scal3, naming_scheme={'parameters': 'c', 'independent': 't', 'dependents': ['C','E','P','S']})
    >>> print(max_scal3.translate(original_system).subs({'c0':'alpha'}))
    dt/dt = 1
    dC/dt = -C - C/alpha + E*S
    dE/dt = C + C/alpha - E*S
    dP/dt = C
    dS/dt = C/alpha - E*S
    parameters: alpha

So we have found a third different reparametrization of the Michaelis-Menten equations.

.. todo::

    Add a method to :class:`~desr.ode_translation.ODETranslator` that will try and re-order the last :math:`n-r` columns so
    that the parameter reduction scheme can be applied.