    

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
    >>> system.reorder_variables(['t', 's', 'c', 'k_m1', 'k_2', 'k_1', 'e_0'])
    >>> system
    dt/dt = 1
    ds/dt = c*k_1*s + c*k_m1 - e_0*k_1*s
    dc/dt = -c*k_1*s - c*k_2 - c*k_m1 + e_0*k_1*s
    dk_m1/dt = 0
    dk_2/dt = 0
    dk_1/dt = 0
    de_0/dt = 0

We can see the system's properties

    >>> system.exponent_matrix()
    Matrix([
    [1, 1,  1, 1, 1, 1,  1],
    [0, 0, -1, 1, 0, 0,  1],
    [0, 1,  1, 0, 0, 0, -1],
    [0, 0,  1, 0, 0, 1,  0],
    [0, 0,  0, 0, 1, 0,  0],
    [1, 1,  0, 1, 0, 0,  1],
    [1, 0,  0, 0, 0, 0,  1]])

One reduces by first making a `ODETranslation`.  One can query the properties/results of the translation, like `scaling_matrix` and `invariants`.

    >>> translation = ODETranslation.from_ode_system(system, naming_scheme=('tau',['u','v'], 'c'))
    >>> translation.scaling_matrix
    Matrix([
    [1, 0, 0, -1, -1, -1, 0],
    [0, 1, 1,  0,  0, -1, 1]])
    >>> translation.invariants()
    Matrix([[e_0*k_1*t, s/e_0, c/e_0, k_m1/(e_0*k_1), k_2/(e_0*k_1)]])

Actually do the reduction by calling the `translate` function on the `ODETranslation` instance.

    >>> translation.translate(system)
    dtau/dtau = 1
    du/dtau = c0*v + u*v - u
    dv/dtau = -c0*v - c1*v - u*v + u
    dc0/dtau = 0
    dc1/dtau = 0










