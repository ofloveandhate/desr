    

Michaelis-Menten basic
==========================




Michaelis-Menten model, two variables
-----------------------------------------

There are two chemical species :math:`E` and :math:`S`, which reversibly combine to form :math:`ES`, which decomposes into :math:`E` and :math:`P`.

.. math::
    :nowrap:

    \[
    E+S\leftrightarrow ES \rightarrow E+P
    \]

After using the law of mass action, we are able to translate to the ODE system

.. math::
    :nowrap:

    \begin{align}
    \frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\
    \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c
    \end{align}


Model reduction in `desr`
----------------------------


First we make the system in Python from LaTeX code.  

    >>> system_tex = '''\frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\\\
    ...          \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c'''
    >>> system = ODESystem.from_tex(system_tex)
    >>> system.reorder_variables(['e_0', 'k_1', 'k_2', 'k_m1', 't', 's', 'c'])
    >>> system
    dt/dt = 1
    ds/dt = c*k_1*s + c*k_m1 - e_0*k_1*s
    dc/dt = -c*k_1*s - c*k_2 - c*k_m1 + e_0*k_1*s
    parameters: e_0, k_1, k_2, k_m1

We can see the system's properties

    >>> system.exponent_matrix()
    Matrix([
    [1, 0,  0, 0, 0, 0,  1],
    [1, 1,  0, 1, 0, 0,  1],
    [0, 0,  0, 0, 1, 0,  0],
    [0, 0,  1, 0, 0, 1,  0],
    [1, 1,  1, 1, 1, 1,  1],
    [0, 0, -1, 1, 0, 0,  1],
    [0, 1,  1, 0, 0, 0, -1]])

One reduces by first making a `ODETranslator`.  One can query the properties/results of the translator, like `scaling_matrix` and `invariants`.

    >>> translator = ODETranslator.from_ode_system(system, naming_scheme={'parameters': 'c', 'independent': 'tau', 'dependents': ['u','v']})
    >>> translator.scaling_matrix
    Matrix([
    [0, -1, -1, -1, 1, 0, 0],
    [1, -1,  0,  0, 0, 1, 1]])
    >>> translator.invariant_monomials()
    Matrix([[k_2/(e_0*k_1), k_m1/(e_0*k_1), e_0*k_1*t, s/e_0, c/e_0]])

Actually do the reduction by calling the `translate` function on the `ODETranslator` instance.

    >>> translator.translate(system)
    dtau/dtau = 1
    du/dtau = c1*v + u*v - u
    dv/dtau = -c0*v - c1*v - u*v + u
    parameters: c0, c1










