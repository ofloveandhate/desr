import itertools
import re

import sympy

from desr.matrix_normal_forms import hnf_col, hnf_row, normal_hnf_col
from desr.sympy_helper import (expressions_to_variables, unique_array_stable,
                               monomial_to_powers, sympify_expression, to_symbol)
from desr.tex_tools import TexParseError, expr_to_tex, var_to_tex, tex_to_sympy


def _as_substitution_pairs(to_sub, extra_args):
    '''
    Normalise the argument shapes :meth:`sympy.core.basic.Basic.subs` accepts into a list of
    `(old, new)` pairs, so that :meth:`ODESystem.subs` takes the same shapes sympy does.

    Args:
        to_sub: A dict, an iterable of `(old, new)` pairs, or the `old` of an `old, new` pair.
        extra_args (tuple): The positional arguments following `to_sub`; either empty, or the
            single `new` of an `old, new` pair.

    Returns:
        list: The `(old, new)` pairs, in the order given.  A dict keeps its insertion order.
    '''
    if extra_args:
        if len(extra_args) > 1:
            raise ValueError('subs takes either one argument -- a dict or an iterable of '
                             '(old, new) pairs -- or two, an old and a new.  Got {} '
                             'arguments.'.format(1 + len(extra_args)))
        return [(to_sub, extra_args[0])]

    if isinstance(to_sub, dict):
        return list(to_sub.items())

    if isinstance(to_sub, (str, sympy.Basic)):
        raise ValueError('subs was given a single {} with nothing to replace it by.  Pass a '
                         'dict, an iterable of (old, new) pairs, or both an old and a '
                         'new.'.format(type(to_sub).__name__))

    pairs = [tuple(pair) for pair in to_sub]
    for pair in pairs:
        if len(pair) != 2:
            raise ValueError('subs expects (old, new) pairs; got a group of {}: '
                             '{}.'.format(len(pair), pair))
    return pairs


class ODESystem(object):
    '''
    A system of differential equations.

    **The vocabulary**, since it is the thing that most often trips people up -- and did trip
    up this library, which inherited the muddle from Hubert & Labahn:

    - :attr:`~desr.ode_system.ODESystem.symbols` is the collective: every name the system
      declares.  It is what :math:`n` counts.
    - :attr:`~desr.ode_system.ODESystem.parameters` are the :math:`p` symbols with no
      equation of their own.
    - the *variables* are the rest -- the ones that actually vary -- namely
      :attr:`~desr.ode_system.ODESystem.independent_variable` and the :math:`q`
      :attr:`~desr.ode_system.ODESystem.dependents`.

    So :math:`\\textrm{symbols} = \\textrm{parameters} + \\textrm{variables}` and
    :math:`n = p + 1 + q`.  Hubert & Labahn have no word for the collective, which is how
    "variable" came to do that job in their paper and lose its own meaning along the way:
    at the top of their section 7 they separate "state variables" from "constant
    parameters", and a page later call :math:`\\alpha = k_3/(k_1 a)` a variable.

    :attr:`~desr.ode_system.ODESystem.derivatives` is an ordered tuple the same length as
    `symbols`, holding each derivative with respect to the independent variable.  The
    classification is read straight off it, which is why it is exhaustive: the independent
    variable, a :const:`None` for a parameter, or an expression for a dependent.

    An equation that says zero is still an equation.  :math:`dw/dt = 0` states that
    :math:`w` is *conserved*, which is a fact about a dependent; to say instead that it does
    not vary at all, name it in `parameters` or call
    :meth:`~desr.ode_system.ODESystem.make_parameter`.  Both print as `0`, so only the
    classification tells them apart -- see
    :meth:`~desr.ode_translation.ODETranslator.symbol_kinds`.

    Args:
        variables (tuple of sympy.Symbol): Ordered tuple of symbols.  The argument keeps its
            old name for now; it is the `symbols` of the constructed system.
        derivatives (tuple of sympy.Expression): Ordered tuple of derivatives, :const:`None` for a parameter.
        indep_var (sympy.Symbol, optional): Independent variable we are differentiating with respect to.
        initial_conditions (tuple of sympy.Symbol): The initial values of dependent variables
        same_units (dict of sympy.Symbol : sympy.Expression): lists of things that must have the same units.
        parameters (iter of str or sympy.Symbol, optional): Variables to treat as parameters,
            discarding the zero equations they were given.
        scheme (desr.ode_translation.ReductionScheme, optional): For a reduced system, which reduction produced
            it.  **Stored, not inferred**: a reduced system handed on by itself cannot work
            out where it came from, and inferring is only right when
            :meth:`~desr.ode_translation.ODETranslator.translate` did the choosing.
        num_scaling_symmetries (int, optional): For a reduced system, the :math:`r` scaling
            symmetries that were reduced away.
        auxiliary_names (iter of sympy.Symbol, optional): For a reduced system, which of
            its symbols are auxiliaries -- normalised directions, free to be set to
            anything.  Also stored rather than inferred: auxiliary-ness is a fact about
            which block of the Hermite multiplier a symbol came from, and nothing in the
            reduced system's own equations records it.
    '''

    def __init__(self, variables, derivatives, indep_var=None, initial_conditions=None, same_units=None, constraints=None, is_reduced=False, parameters=None, scheme=None, num_scaling_symmetries=None, auxiliary_names=None):
        self._variables = tuple(variables)
        self._derivatives = tuple(derivatives)

        self._indep_var = sympy.var('t') if indep_var is None else indep_var

        # initialize some class member storage
        self._initial_conditions = {}
        self._constraints = []
        self._requires_same_unit = {} # a dict of {Symbol : Expression}
        self._is_reduced = is_reduced

        # What a reduced system cannot re-derive about itself, and so must carry.
        self._scheme = scheme
        self._num_scaling_symmetries = num_scaling_symmetries
        self._auxiliary_names = tuple(auxiliary_names or ())

        # sanity checks
        assert len(self._variables) == len(self._derivatives)
        assert self._derivatives[self.independent_variable_index] == sympy.sympify(1)

        # Reclassify before anything else reads the derivatives, so that the rest of
        # construction sees the system the caller asked for.
        if parameters is not None:
            for parameter in parameters:
                self.make_parameter(parameter)

        # carry over additional arguments
        if initial_conditions is not None:
            self.update_initial_conditions(initial_conditions=initial_conditions)

        if same_units is not None:
            self.require_same_unit(requirements=same_units)

        if constraints is not None:
            if type(constraints) is not list:
                raise TypeError(f'constraints must be a list of sympy equalities')

            for c in constraints: 
                if not isinstance(c,sympy.core.relational.Equality):
                    raise TypeError(f'{c} is not a sympy equality.  Constraints must be sympy equalities.')
            self._constraints = constraints

        self._variable_sanity_check()


    def __eq__(self, other):
        if type(self) is not type(other):
            return False

        # Compare variables
        self_var = sorted(self.symbols, key=str)
        other_var = sorted(other.symbols, key=str)
        if self_var != other_var:
            return False

        # Compare derivatives
        self_der, other_der = self.derivative_dict, other.derivative_dict
        for var1, var2 in zip(self_var, other_var):
            der1 = self_der.get(var1)
            der2 = other_der.get(var2)
            if der1 is None:
                if der2 is not None:
                    return False
            else:
                if der2 is None:
                    return False
                if der1.expand() != der2.expand():
                    return False
        # Compare independent variables
        if self._indep_var != other._indep_var:
            return False

        return True

    def copy(self):
        '''
        Returns:
            ODESystem: A copy of the system.
        '''
        system = ODESystem(self._variables, self._derivatives, indep_var=self._indep_var,
                           is_reduced=self.is_reduced, same_units=self._requires_same_unit,
                           scheme=self._scheme,
                           num_scaling_symmetries=self._num_scaling_symmetries,
                           auxiliary_names=self._auxiliary_names)
        system.update_initial_conditions(self.initial_conditions)
        for eqn in self.constraints:
            system.add_constraint(eqn.lhs, eqn.rhs)

        
        return system

    @property
    def is_reduced(self):
        return self._is_reduced


    def rename_indep_var(self, new_indep):
        """
        >>> _input = {'x': 'c_0*x*y', 'y': 'c_1*(1-x)*(1-y)*t'}
        >>> _input = {sympy.Symbol(k): sympy.sympify(v) for k, v in _input.items()}
        >>> system = ODESystem.from_dict(_input)
        >>> system.update_initial_conditions({'x': 'x_0'})
        >>> system.initial_conditions
        {x: x_0}

        >>> system.rename_indep_var('tau')
        >>> system
        dtau/dtau = 1
        dx/dtau = c_0*x*y
        dy/dtau = c_1*tau*(1 - x)*(1 - y)
        x(0) = x_0
        parameters: c_0, c_1, x_0

        >>> system.symbols
        (c_0, c_1, x_0, tau, x, y)
        """
        if not isinstance(new_indep, sympy.Symbol):
            new_indep = sympify_expression(new_indep)

            expressions_to_variables

        to_sub = {self._indep_var:new_indep}
        
        self._derivatives = [d.subs(to_sub) if d else d for d in self._derivatives]
        self._constraints = [c.subs(to_sub) for c in self.constraints]
        self._initial_conditions = {k:v.subs(to_sub) for k,v in self._initial_conditions.items()}

        ind = self.independent_variable_index
        self._variables = self._variables[:ind] + (new_indep,) + self._variables[ind+1:]
        self._indep_var = new_indep

    @property
    def independent_variable(self):
        """
        Return the independent variable: the one everything is differentiated with respect to.

        Hubert & Labahn call it *time*, which it usually is.  The library does not insist on
        that -- it could be arclength or a spatial coordinate -- so the general name is used
        here, with :py:attr:`~indep_var` as a short alias for typing at a prompt.

        Returns:
            sympy.Symbol
        """
        return self._indep_var

    #: Short alias for :py:attr:`~independent_variable`.  The library itself spells the name
    #: out, so that a search for the long form finds every use within desr.
    indep_var = independent_variable

    @property
    def independent_variable_index(self):
        """
        Return the index of the independent variable in :py:attr:`~symbols`.

        Under the canonical order this equals :py:attr:`~num_parameters`, since the
        parameters come first.

        Return:
             int
        """
        return self.symbols.index(self.independent_variable)

    #: Short alias for :py:attr:`~independent_variable_index`.
    indep_var_index = independent_variable_index

    @property
    def symbols(self):
        '''
        Every symbol the system declares: the parameters, the independent variable and the
        dependent variables, in order.

        **The vocabulary.**  `symbols` is the collective; it is what :math:`n` counts.  The
        *variables* are the ones that actually vary -- the independent variable and the
        dependents -- and the *parameters* are the rest.  So

        .. math::

            \\textrm{symbols} = \\textrm{parameters} + \\textrm{variables},
            \\qquad n = p + 1 + q.

        Hubert & Labahn have no term for the collective, which is why "variable" drifted
        into the job in their paper and in earlier versions of this library: at the top of
        their section 7 they distinguish "state variables" from "constant parameters", and a
        page later call :math:`\\alpha = k_3/(k_1 a)` -- an expression in the parameters
        alone -- a variable.  Naming the collective is what keeps "variable" meaning
        something.

        Returns:
            tuple: Ordered tuple of every symbol in the system.
        '''
        return self._variables

    @property
    def parameters(self):
        '''
        The parameters of the system: the symbols with no equation of their own.

        A symbol whose equation says zero is *conserved*, which is a different thing, and so
        is a dependent rather than a parameter.  See :attr:`~desr.ode_system.ODESystem.dependents`.

        Returns:
            tuple: The parameters, in symbol order.

        >>> system = ODESystem.from_equations(['dz/dt = a*z', 'dw/dt = 0'])
        >>> system.parameters
        (a,)
        >>> system.dependents
        (w, z)
        '''
        return tuple(var for var, deriv in zip(self.symbols, self._derivatives) if deriv is None)

    @property
    def num_parameters(self):
        '''
        Returns:
            int: The number of parameters.
        '''
        return len(self.parameters)

    @property
    def dependents(self):
        '''
        The dependents: the symbols with an equation of their own, excluding the independent
        variable.

        Hubert & Labahn's "state variables", and the right name for them in prose -- the
        state is what a solver integrates.  The accessor is the short form so that it pairs
        with :attr:`~parameters` without either needing a noun bolted on.

        A conserved symbol, whose equation says zero, is a dependent.

        Returns:
            tuple: The dependents, in symbol order.

        >>> _input = {'x': 'c_0*x*y', 'y': 'c_1*(1-x)*(1-y)*t'}
        >>> _input = {sympy.Symbol(k): sympy.sympify(v) for k, v in _input.items()}
        >>> system = ODESystem.from_dict(_input)
        >>> system.dependents
        (x, y)

        A dependent may move at unit rate; being told apart from the independent variable is
        a matter of which symbol it is, not of what its equation happens to say.

        >>> ODESystem.from_equations(['dx/dt = 1']).dependents
        (x,)
        '''
        return tuple(var for var, deriv in zip(self.symbols, self._derivatives)
                     if deriv is not None and var != self.independent_variable)

    @property
    def num_dependents(self):
        '''
        Returns:
            int: The number of dependents, excluding the independent variable.
        '''
        return len(self.dependents)

    @property
    def variables(self):
        '''
        The symbols that actually vary: the independent variable, then the dependents.

        The complement of :attr:`~parameters` within :attr:`~symbols`, and the set the
        classification rule asks about -- an invariant is a parameter exactly when its
        monomial involves none of these.

        Note that this is **not** what `variables` meant before 2.0, when it was the
        collective and is now :attr:`~symbols`.  The rename went in first and this came
        back afterwards, so that nothing could silently receive the shorter tuple.

        Returns:
            tuple: The independent variable followed by the dependents.

        >>> system = ODESystem.from_equations(['dz/dt = a*z', 'dw/dt = 0'])
        >>> system.symbols
        (a, t, w, z)
        >>> system.parameters, system.variables
        ((a,), (t, w, z))
        >>> system.num_parameters + system.num_variables == len(system.symbols)
        True
        '''
        return (self.independent_variable,) + self.dependents

    @property
    def num_variables(self):
        '''
        Returns:
            int: The number of symbols that vary, which is ``1 + num_dependents``.
        '''
        return 1 + self.num_dependents

    @property
    def conserved_variables(self):
        '''
        The dependents whose equation says zero.

        *Conserved* is the right word only for a dependent: a quantity that does not change
        along the trajectory.  A parameter does not change either, but it is not conserved
        *along* anything -- it has no equation at all.  The umbrella word for both is
        "constant", which is useful in prose and deliberately has no accessor: it would be
        implied for every parameter, impossible for the independent variable, and
        informative only here, where it is this.

        A conserved dependent and a parameter both print as `0` or are absent from the
        equations, so nothing in the printed system tells them apart.  What does is the
        classification: see :attr:`~parameters` and
        :meth:`~desr.ode_translation.ODETranslator.symbol_kinds`.

        Returns:
            tuple: The conserved dependents, in symbol order.

        >>> system = ODESystem.from_equations(['dx/dt = a*x - b*x*y',
        ...                                    'dy/dt = b*x*y - c*y',
        ...                                    'dw/dt = 0'])
        >>> system.conserved_variables
        (w,)
        >>> system.parameters
        (a, b, c)

        `w` is not a parameter, and saying `dw/dt = 0` is how you state that:

        >>> system.derivative_dict[sympy.Symbol('w')]
        0
        '''
        return tuple(var for var in self.dependents
                     if self.derivative_dict[var] == 0)

    @property
    def scheme(self):
        '''
        Which reduction produced this system, or :const:`None` if it is not a reduction.

        :rtype: desr.ode_translation.ReductionScheme or None
        '''
        return self._scheme

    @property
    def num_scaling_symmetries(self):
        '''
        The :math:`r` scaling symmetries reduced away, or :const:`None` for an original system.

        :rtype: int or None
        '''
        return self._num_scaling_symmetries

    @property
    def auxiliary_names(self):
        '''
        The symbols that are auxiliaries: normalised directions, free to take any value.

        Empty for an original system, and for one reduced by the parameter scheme, which
        normalises its auxiliaries to 1 and discards them.

        **Orthogonal to the classification, not a fourth kind.**  An auxiliary is a
        parameter or a dependent like any other symbol -- in prey-predator's general
        reduction the auxiliaries ``h/k``, ``K`` and ``K/h`` are monomials in the original
        parameters and so are parameters, while Hubert & Labahn 6.6's auxiliary
        :math:`t z_1^2` involves a dependent and evolves.  Auxiliary-ness says you may
        choose the value; the kind says whether it varies.

        Returns:
            tuple: The auxiliaries, in symbol order.
        '''
        return self._auxiliary_names

    def is_auxiliary(self, variable):
        '''
        Whether `variable` is an auxiliary of this system.

        Args:
            variable (str or sympy.Symbol): The symbol to ask about.

        :rtype: bool
        '''
        if isinstance(variable, str):
            variable = sympy.Symbol(variable)
        return variable in self._auxiliary_names

    @property
    def derivatives(self):
        ''' Getter for an ordered tuple of expressions representing the derivatives of self.symbols.

        A parameter has no equation, and so is reported as :const:`None`.  That is not the
        same answer as ``0``, which is the equation of a conserved variable.

        Returns:
            tuple: Ordered tuple of sympy.Expression, or :const:`None` for a parameter.
        '''
        return self._derivatives

    @property
    def derivative_dict(self):
        '''
        Return a variable: expr mapping, over the variables that have an equation.

        Returns:
            dict: Keys are the independent and dependent variables, value is the derivative with respect to the independent variable.
        '''
        return dict(filter(lambda x: x[1] is not None, zip(self.symbols, self._derivatives)))

    def make_parameter(self, variable):
        '''
        Declare a variable to be a parameter, discarding the equation it was given.

        Only a variable that does not move can become a parameter: one that is already a
        parameter, in which case this does nothing, or one whose equation says zero.  An
        explicit :math:`dw/dt = 0` says :math:`w` is conserved, which is a statement about
        a dependent variable; this is how to say instead that it is a constant of the model.

        The variable keeps its place in the ordering.  Only
        :meth:`~desr.ode_system.ODESystem.default_order_variables` moves it to sit with the
        other parameters, because the ordering decides which parameters are normalised away.

        Args:
            variable (str, sympy.Symbol): The variable to reclassify.

        Raises:
            ValueError: If the variable is unknown to the system, is the independent
                variable, has an equation saying it moves, or has an initial condition.

        >>> system = ODESystem.from_equations(['dz/dt = a*z', 'dw/dt = 0'])
        >>> system.make_parameter('w')
        >>> system.parameters
        (a, w)
        >>> system.dependents
        (z,)

        >>> system.make_parameter('z')
        Traceback (most recent call last):
            ...
        ValueError: Cannot make a parameter of z, which moves: dz/dt = a*z.
        '''
        if isinstance(variable, str):
            variable = to_symbol(variable)

        if variable not in self.symbols:
            raise ValueError('Cannot make a parameter of {}, which is not a variable of the '
                             'system.'.format(variable))

        if variable == self.independent_variable:
            raise ValueError('Cannot make a parameter of {}, the independent '
                             'variable.'.format(variable))

        index = self.symbols.index(variable)
        derivative = self._derivatives[index]

        if derivative is None:
            return

        if derivative != sympy.sympify(0):
            raise ValueError('Cannot make a parameter of {}, which moves: d{}/d{} = {}.'.format(
                variable, variable, self.independent_variable, derivative))

        initial_condition = self._initial_conditions.get(variable)
        if initial_condition is not None:
            raise ValueError('Cannot make a parameter of {}, which has an initial condition '
                             '{}(0) = {}.  Drop the initial condition first.'.format(
                                 variable, variable, initial_condition))

        derivatives = list(self._derivatives)
        derivatives[index] = None
        self._derivatives = tuple(derivatives)

    @property
    def initial_conditions(self):
        '''
        Return a variable: initial-value mapping.

        Returns:
            dict: Keys are non-constant variables, value is the constant representing their initial condition.
        '''
        return self._initial_conditions.copy()

    @property
    def constraints(self):
        '''
        Returns the constraints imposed on the system.

        Returns:

        '''
        return self._constraints[:]
    
    @property
    def dimension_requirements(self):
        '''
        Returns: the list of lists of sympy.Expression that must have the same units.  Only the user knows these, these are NOT computed
        '''
        return self._requires_same_unit


    def simplify_derivatives(self):
        '''
        use sympy to collect the derivatives, in an effort to make them 
        simpler to print to screen / into tex for a paper, etc

        Returns:
            None
        '''

        new_derivs = []

        for d in self._derivatives:
            if d is not None:
                for v in self.symbols:
                    d = sympy.collect(d,v)
            new_derivs.append(d)

        self._derivatives = tuple(new_derivs)


    def update_initial_conditions(self, initial_conditions):
        '''
        Update the internal record of initial conditions.

        Args:
            initial_conditions (dict): non-constant variable: initial value constant.

        >>> _input = {'x': 'c_0*x*y', 'y': 'c_1*(1-x)*(1-y)*t'}
        >>> _input = {sympy.Symbol(k): sympy.sympify(v) for k, v in _input.items()}
        >>> system = ODESystem.from_dict(_input)
        >>> system.update_initial_conditions({'x': 'x_0'})
        >>> system.initial_conditions
        {x: x_0}
    
        >>> system
        dt/dt = 1
        dx/dt = c_0*x*y
        dy/dt = c_1*t*(1 - x)*(1 - y)
        x(0) = x_0
        parameters: c_0, c_1, x_0

        >>> system.update_initial_conditions({'c_0': 'k'})
        Traceback (most recent call last):
            ...
        ValueError: Cannot set initial condition k for c_0, which is not a dependent variable of the system.

        >>> system
        dt/dt = 1
        dx/dt = c_0*x*y
        dy/dt = c_1*t*(1 - x)*(1 - y)
        x(0) = x_0
        parameters: c_0, c_1, x_0

        >>> system.is_reduced
        False

        >>> system.update_initial_conditions({'x': '1'})
        Traceback (most recent call last):
            ...
        ValueError: an initial condition for a non-reduced system should just be a single variable.  you have set()
        
        >>> system.update_initial_conditions({'x': '0'})
        Traceback (most recent call last):
            ...
        ValueError: an initial condition for a non-reduced system should just be a single variable.  you have set()

        >>> system
        dt/dt = 1
        dx/dt = c_0*x*y
        dy/dt = c_1*t*(1 - x)*(1 - y)
        x(0) = x_0
        parameters: c_0, c_1, x_0
        '''
        for variable, init_cond in initial_conditions.items():

            # convert through sympy
            if not isinstance(variable, sympy.Symbol):
                variable = sympy.Symbol(variable)

            if isinstance(init_cond, str):
                init_cond = sympify_expression(init_cond)

            # We can only set initial conditions for dependent variables we already know about.
            if variable not in self.dependents:
                raise ValueError('Cannot set initial condition {} for {}, which is not a dependent variable of the system.'.format(init_cond,
                                                                                                              variable))



            ics_vars = expressions_to_variables([init_cond]) # convert. 



            # reduced systems are allowed to have 1's (and potentially 0's) as initial conditions
            # but non-reduced systems are NOT.  Otherwise, there might be 1's with different units, and this is nonsense.
            if not self.is_reduced:

                if len(ics_vars) != 1:
                    raise ValueError(f'an initial condition for a non-reduced system should just be a single variable.  you have {ics_vars}')

                if init_cond.is_constant():
                    raise ValueError(f"initial condition `{init_cond}` doesn't appear to be a variable expression.  Non-reduced systems can only have variables as initial conditions.  If you want {variable}(0) to be a numeric constant like 1 or 0, accomplish this via substitution (after reduction)")

            # silviana: i think we should guard against having 1 be a variable...
            if (init_cond not in self.symbols) and not init_cond.is_constant(): # checking emptiness on ics_vars
                self._add_parameters([init_cond])

            self._initial_conditions[variable] = init_cond
    

    def require_same_unit(self,requirements):
        """
        Add requirements to the system that some variables have the same units as some corresponding expressions

        Args:
            requirements (dict of sympy.Symbol:sympy.Expression): 

        Keys: var (sympy.Symbol): a variable to require to have the same unit as the expression expr
        Values: expr (sympy.Expression): an expression that must have the same net unit as the variable var 

        This has the effect of adding columns to the Power Matrix / Exponent Matrix.
        """


        result = {}

        for var, expr in requirements.items():

            # convert to symbols and expressions
            if not isinstance(var, sympy.Symbol):
                var = sympy.Symbol(var)

            if isinstance(expr, str):
                expr = sympify_expression(expr)

            expr_vars = expressions_to_variables([expr])

            if var not in self._variables:
                raise ValueError(f'The variable {var} for a same-dimension requirement is not already a variable in this system.')

            for e in expr_vars:
                if e not in self._variables:
                    raise ValueError(f'The symbol {e} used in expression {expr} for a same-dimension requirement is not already a variable in this system.')

            result[var] = expr

        self._requires_same_unit = result


    def add_constraint(self, lhs, rhs):
        '''
        Add constraint that must be obeyed by the system.

        Args:
            lhs (sympy.Expr): The left hand side of the constraint.
            rhs (sympy.Expr): The right hand side of the constraint.

        Todo:
            * Finish docstring and tests, here and for: finding scaling symmetries and also translation
            * Check for 0 case

        >>> eqns = ['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)']

        >>> system = ODESystem.from_equations(eqns)
        >>> system
        dt/dt = 1
        dx/dt = c_0*x*y
        dy/dt = c_1*(1 - x)*(1 - y)
        parameters: c_0, c_1

        >>> system.add_constraint('c_2', 'c_0 + c_1')
        >>> system
        dt/dt = 1
        dx/dt = c_0*x*y
        dy/dt = c_1*(1 - x)*(1 - y)
        c_2 == c_0 + c_1
        parameters: c_0, c_1, c_2

        >>> system.add_constraint('c_2', 'c_0 + x')
        Traceback (most recent call last):
            ...
        ValueError: Cannot add constraints involving the dependent variables {x}.  Try to use substitution to incorporate the constraints directly into the system before construction.  (To do this automatically would be small project...

        >>> system.add_constraint('c_0', 0)
        Traceback (most recent call last):
            ...
        ValueError: Cannot express equality constraint with 0. (c_0 == 0)
        '''
        if isinstance(lhs, str):
            lhs = sympify_expression(lhs)
        if isinstance(rhs, str):
            rhs = sympify_expression(rhs)
        if (lhs == 0) or (rhs == 0):
            raise ValueError(f'Cannot express equality constraint with 0. ({lhs} == {rhs})')

        # get the variables from the constraint
        all_constraint_vars, dependent_vars, new_vars = self.equation_variables(lhs, rhs)

        # A constraint may only relate parameters.  One involving a dependent variable is a
        # statement about the trajectory, which this library cannot act on.
        if dependent_vars:
            raise ValueError('Cannot add constraints involving the dependent variables {}. '.format(dependent_vars) +
                             ' Try to use substitution to incorporate the constraints directly into the system before construction.  (To do this automatically would be small project...')

        # Only parameters involved, no problem.  update the recorded symbols
        self._add_parameters(new_vars)
        self._constraints.append(sympy.Eq(lhs, rhs))

    def _add_parameters(self, new_vars):
        '''
        Record `new_vars` as parameters of the system, at the end of the parameter block.

        New parameters must land *before* the independent variable, not at the end of the
        variable list.  The canonical order is parameters, independent variable, dependents,
        and a parameter sitting after the dependents breaks it: the parameter reduction
        scheme needs the independent variable at index `num_parameters`, so appending to the
        end silently costs the scheme.  It is easy to miss, because :meth:`__repr__` groups
        the variables by classification and so prints a sane-looking `parameters:` line
        either way.

        The insertion point is the independent variable's index rather than
        `num_parameters`, so that a system the caller has deliberately put in a
        non-canonical order -- :meth:`reorder_variables` allows any permutation -- is not
        silently rearranged around the new parameter.  For a canonical system the two
        agree.

        Args:
            new_vars (list of sympy.Symbol): The parameters to add.  May be empty.
        '''
        if not new_vars:
            return
        at = self.independent_variable_index
        variables = list(self._variables)
        derivatives = list(self._derivatives)
        self._variables = tuple(variables[:at] + list(new_vars) + variables[at:])
        self._derivatives = tuple(derivatives[:at] + [None] * len(new_vars) + derivatives[at:])

    def equation_variables(self, lhs, rhs):
        """
        compute the variables in the left and right hand side of a constraint:
        all, the non-constant ones that exist in the system, and the new ones.

        Args:
            lhs (sympy.Expr): The left hand side of the constraint.
            rhs (sympy.Expr): The right hand side of the constraint.

        Returns: length-three tuple of containers of variables: all, existing dependent, new
        """
        all_vars = expressions_to_variables([lhs, rhs])
        dependent_vars = all_vars.intersection(self.dependents)
        new_vars = sorted(all_vars.difference(set(self.symbols)), key=str)

        return all_vars, dependent_vars, new_vars

    def subs(self, to_sub, *args, expand_before=False, expand_after=True, factor_after=False, subs_constraints=True, new_symbols_are_parameters=True):
        '''
        Make substitutions throughout the system, returning a new system.

        Named and shaped after :meth:`sympy.core.basic.Basic.subs`, which it extends from one
        expression to a whole system: the substitutions are made in the derivatives, and --
        unless `subs_constraints` says otherwise -- in the constraints and the initial
        conditions too.  Like sympy's, it does not mutate; it returns a new system.

        Args:
            to_sub: The substitutions, in any of the forms :meth:`sympy.core.basic.Basic.subs`
                takes -- a dict `{old: new}`, an iterable of `(old, new)` pairs, or the two
                arguments `old, new`.  A dict shows the old-new pairs in one object and is the
                form the examples use.
            expand_before (bool): Expand the sympy expression for each derivative before substitution.
            expand_after (bool): Expand the sympy expression for each derivative after substitution.
            factor_after (bool): Factorise the sympy expression for each derivative after substitution.
            subs_constraints (bool): Perform the substitutions into the initial constraints.
            new_symbols_are_parameters (bool): Treat newly encountered symbols as parameters. If a substitution is a variable renaming, this does not apply.  I have no idea what to do if this is false, so it will generate exceptions..

        Returns:
            ODESystem: System with substitutions carried out.

        >>> eqns = ['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)']
        >>> system = ODESystem.from_equations(eqns)
        >>> system.subs({'1-x': 'z'}, expand_before=False, expand_after=False, factor_after=False)
        Traceback (most recent call last):
        ...
        ValueError: unable to make substitution 1 - x -> z.
        >>> system.subs({'1-x': 'z'}, expand_before=True, expand_after=False, factor_after=False)
        Traceback (most recent call last):
        ...
        ValueError: unable to make substitution 1 - x -> z.

        >>> system.subs({'x': '1-z'}, expand_before=True, expand_after=True, factor_after=False)
        Traceback (most recent call last):
        ...
        ValueError: unable to make substitution x -> 1 - z.

        >>> system.add_constraint('c_0', 'c_1**2')
        >>> system.subs({'c_0': '1'}, subs_constraints=False)
        dt/dt = 1
        dx/dt = x*y
        dy/dt = c_1*x*y - c_1*x - c_1*y + c_1
        c_0 == c_1**2
        parameters: c_0, c_1
        >>> system.subs({'c_0': '1'}, subs_constraints=True)
        dt/dt = 1
        dx/dt = x*y
        dy/dt = c_1*x*y - c_1*x - c_1*y + c_1
        1 == c_1**2
        parameters: c_1

        The sympy call shapes work too, and mean the same thing:

        >>> system.subs('c_0', '1') == system.subs({'c_0': '1'})
        True
        >>> system.subs([('c_0', '1')]) == system.subs({'c_0': '1'})
        True
        '''

        # sympify the dict of things to substitute
        to_sub = {sympify_expression(lhs): sympify_expression(rhs)
                  for lhs, rhs in _as_substitution_pairs(to_sub, args)}

        variable_renamings = []
        
        for lhs,rhs in to_sub.items():
            lhs_vars = set(expressions_to_variables([lhs]))
            rhs_vars = set(expressions_to_variables([rhs]))

            
            lhs_dependent_vars = lhs_vars.intersection(set(self.dependents))

            # make sure we can actually do the substitution
            if lhs_dependent_vars:
                # in this conditional, we're trying to substitute away a dependent variable.  can only do, if the substitution is just a renaming.  anything else must be done by the user

                # if the right hand side is just the only variable on the right hand side, then we CAN do it.  check by getting the first variable and see if equal.
                if (rhs!=list(rhs_vars)[0] or lhs != list(lhs_vars)[0]) and not self.is_reduced:
                    # this is not a parameter substitution, nor is it a simple variable renaming.  
                    raise ValueError(f'unable to make substitution {lhs} -> {rhs}.')

                # ok, we made it here so this substitution is a simple renaming.  
                # we will need to not only do the substitution in the rhs of the derivative, but also in the lhs.
                if (rhs==list(rhs_vars)[0] and lhs == list(lhs_vars)[0]):
                    variable_renamings.append(  (lhs, rhs) ) 


        # copy so can do work on the derivatives, variables, initial conditions
        new_derivs = self._derivatives
        variables = list(self.symbols)
        initial_conditions = self.initial_conditions
        constraints = self.constraints

        # expand, if desired by caller
        if expand_before:
            new_derivs = [d.expand() if d is not None else None for d in new_derivs]



        # this is the main action of this function, really.
        # do the substitutions!
        new_derivs = [d.subs(to_sub) if d is not None else None for d in new_derivs]

        # Now need to substitute in the variable order and initial conditions

        
        for renaming in variable_renamings: # i am sure there is a more elegant way of doing this
            lhs, rhs = renaming # unpack

            # replace in variables.  the derivatives are implicitly stored in this order, so this also renames the derivative
            variables = [v if v!=lhs else rhs for v in variables]
            
            # replace in initial conditions.  the key in the dict is the name of the variable that has an ICS
            initial_conditions = { (v if v!=lhs else rhs) : ics.subs(to_sub) for v,ics in initial_conditions.items()}


        # silviana says: why would we ever not want to do this?  if you make a sub, shouldn't it affect the entire system?!
        if subs_constraints: 
            constraints = [eqn.subs(to_sub) for eqn in constraints]

            # filter out the trivial constraints after doing the substitution
            constraints = [c for c in constraints if not c==True]


        # post-substitution actions

        # expand, if desired by caller
        if expand_after:
            new_derivs = [d.expand() if d is not None else None for d in new_derivs]

        # factor, if desired by caller
        if factor_after:
            new_derivs = [sympy.factor(d) if d is not None else None for d in new_derivs]



        # next, we add the new variables from the substitution
        known_vars = set(expressions_to_variables(variables + new_derivs + constraints + list(initial_conditions.values())))
        missing_vars = sorted(list(set(known_vars) - set(variables)),key=str)
        
        if not new_symbols_are_parameters:
            raise NotImplementedError(f"don't want new symbols from substitution to be parameters, but don't know what else to do.  the new symbols: {missing_vars}")

        # A symbol arriving from a substitution has no equation, so it is a parameter, and it
        # belongs at the end of the parameter block -- before the independent variable, not at
        # the end of the list, which would put it after the dependents and cost the parameter
        # reduction scheme.  None rather than 0: giving it `0` would declare it a conserved
        # dependent.  The derivatives are stored in the variable order, so both lists are cut
        # at the same point.
        at = variables.index(self.independent_variable)
        new_vars = variables[:at] + missing_vars + variables[at:]
        new_derivs = new_derivs[:at] + [None] * len(missing_vars) + new_derivs[at:]


        # remove variables that don't appear in any expressions

        variables_to_inspect = new_vars # because we're about to modify `new_vars` in the loop below
        extant_vars = set(expressions_to_variables( [d for d in new_derivs if d is not None] + constraints + list(initial_conditions.values()) + [self.independent_variable] + list(self.dependents)))
        # An auxiliary is kept even when it appears in no equation.  Under the general
        # scheme an auxiliary is often a monomial in the original parameters, so it has no
        # equation of its own, and pruning it here would throw away part of the record of
        # how this system was reduced -- which the system cannot reconstruct.
        #
        # Both names, because a substitution naming such a parameter does not rename it:
        # `variable_renamings` only picks up a symbol that occurs in an expression, so a
        # parameter with no equation keeps its old name.  Keeping both is what stops it
        # being dropped instead; see ODESystem.subs' note on renaming an equation-less
        # parameter.
        extant_vars |= set(self._auxiliary_names)
        extant_vars |= {to_sub.get(aux, aux) for aux in self._auxiliary_names}

        for v in variables_to_inspect:

            if v not in extant_vars:
                var_loc = new_vars.index(v)
                
                # omit using list slicing.  must omit from both the variable ordering and from the derivatives (which are stored impilictly in the order of the variables)
                new_vars = new_vars[:var_loc] + new_vars[var_loc+1:]
                new_derivs = new_derivs[:var_loc] + new_derivs[var_loc+1:]

        # we're finally ready to construct the new post-substitution System. Huzzah.
        # A substitution renames and rewrites; it does not change where the system came from,
        # so the stored provenance travels with it.  Auxiliaries are filtered to those that
        # survived: a substitution may rename one, or remove it by making it vanish from
        # every expression.
        subs_system = ODESystem(new_vars, new_derivs,
                                initial_conditions=initial_conditions,
                                indep_var=self.independent_variable,
                                constraints=constraints,
                                is_reduced = self.is_reduced,
                                scheme=self._scheme,
                                num_scaling_symmetries=self._num_scaling_symmetries,
                                auxiliary_names=[
                                    aux if aux in new_vars else to_sub.get(aux)
                                    for aux in self._auxiliary_names
                                    if aux in new_vars or to_sub.get(aux) in new_vars])

        subs_system._variable_sanity_check()
        return subs_system

    @classmethod
    def from_equations(cls, equations, indep_var=sympy.var('t'), initial_conditions=None, is_reduced=False, parameters=None):
        '''
        Instantiate from multiple equations.

        Args:
            equations (str, iter of str): Equations of the form "dx/dt = expr", optionally seperated by :code:`\\n`.
            indep_var (sympy.Symbol): The independent variable, usually :code:`t`.
            initial_conditions (tuple of sympy.Symbol): The initial values of dependent variables
            parameters (iter of str or sympy.Symbol, optional): Variables to treat as parameters,
                discarding the zero equations they were given.

        Returns:
            ODESystem: System of equations.

        >>> eqns = ['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)']
        >>> ODESystem.from_equations(eqns)
        dt/dt = 1
        dx/dt = c_0*x*y
        dy/dt = c_1*(1 - x)*(1 - y)
        parameters: c_0, c_1
        >>> eqns = '\\n'.join(['dy/dx = c_0*x*y', 'dz/dx = c_1*(1-y)*z**2'])
        >>> ODESystem.from_equations(eqns, indep_var=sympy.Symbol('x'))
        dx/dx = 1
        dy/dx = c_0*x*y
        dz/dx = c_1*z**2*(1 - y)
        parameters: c_0, c_1

        A variable written down with a zero equation is conserved, and stays a dependent
        variable; `parameters` is how to ask for it to be a constant of the model instead.

        >>> ODESystem.from_equations(['dz/dt = a*z', 'dw/dt = 0'])
        dt/dt = 1
        dw/dt = 0
        dz/dt = a*z
        parameters: a
        >>> ODESystem.from_equations(['dz/dt = a*z', 'dw/dt = 0'], parameters=['w'])
        dt/dt = 1
        dz/dt = a*z
        parameters: a, w
        '''
        if isinstance(equations, str):
            equations = equations.strip().split('\n')

        deriv_dict = dict(map(lambda x: parse_de(x, indep_var=str(indep_var)), equations))
        system = cls.from_dict(deriv_dict=deriv_dict, indep_var=indep_var, initial_conditions=initial_conditions, is_reduced=is_reduced, parameters=parameters)
        system.default_order_variables()

        system._variable_sanity_check()

        return system

    @classmethod
    def from_dict(cls, deriv_dict, indep_var=sympy.var('t'), initial_conditions=None, is_reduced=False, parameters=None):
        '''
        Instantiate from a text of equations.

        Args:
            deriv_dict (dict): {variable: derivative} mapping.  A variable absent from the
                mapping, but appearing in the expressions, is a parameter.
            indep_var (sympy.Symbol): Independent variable, that the derivatives are with respect to.
            initial_conditions (tuple of sympy.Symbol): The initial values of dependent variables
            parameters (iter of str or sympy.Symbol, optional): Variables to treat as parameters,
                discarding the zero equations they were given.


        Returns:
            ODESystem: System of ODEs.

        >>> _input = {'x': 'c_0*x*y', 'y': 'c_1*(1-x)*(1-y)'}
        >>> _input = {sympy.Symbol(k): sympy.sympify(v) for k, v in _input.items()}
        >>> ODESystem.from_dict(_input)
        dt/dt = 1
        dx/dt = c_0*x*y
        dy/dt = c_1*(1 - x)*(1 - y)
        parameters: c_0, c_1

        >>> _input = {'y':  'c_0*x*y', 'z': 'c_1*(1-y)*z**2'}
        >>> _input = {sympy.Symbol(k): sympy.sympify(v) for k, v in _input.items()}
        >>> ODESystem.from_dict(_input, indep_var=sympy.Symbol('x'))
        dx/dx = 1
        dy/dx = c_0*x*y
        dz/dx = c_1*z**2*(1 - y)
        parameters: c_0, c_1
        '''
        # Make a tuple of all variables.
        variables = set(expressions_to_variables(deriv_dict.values())).union(set(deriv_dict.keys()))
        if initial_conditions is not None:
            variables.update(map(expressions_to_variables, initial_conditions.values()))
        variables = tuple(variables.union(set([indep_var])))

        assert ((deriv_dict.get(indep_var) is None) or (deriv_dict.get(indep_var) == 1))
        deriv_dict[indep_var] = sympy.sympify(1)

        system = cls(variables,
                     tuple([deriv_dict.get(var) for var in variables]),
                     indep_var=indep_var,
                     initial_conditions=initial_conditions,
                     is_reduced=is_reduced,
                     parameters=parameters)

        system.default_order_variables()

        system._variable_sanity_check()

        return system

    def __repr__(self):
        # None means "no equation", so the None's are exactly the parameters; they are
        # named below instead of being given a line claiming they do not move.  That leaves
        # `d_/dt = 0` free to mean what it says -- a conserved dependent variable.
        lines = ['d{}/d{} = {}'.format(var, self.independent_variable, expr)
                 for var, expr in zip(self.symbols, self._derivatives) if expr is not None]
        for v in self.dependents:
            init_cond = self.initial_conditions.get(v)
            if init_cond is not None:
                lines.append('{}(0) = {}'.format(v, init_cond))
        for eqn in self.constraints:
            lines.append('{} == {}'.format(eqn.lhs, eqn.rhs))
        if self.parameters:
            lines.append('parameters: {}'.format(', '.join(map(str, self.parameters))))
        return '\n'.join(lines)

    def to_tex(self):
        '''
        Returns:
            str: TeX representation.

        Parameters are left out.  :math:`dk_1/dt = 0` is not part of the mathematics of the
        system; it is how this library records that :math:`k_1` has no equation, and it was
        being deleted by hand from the generated tex anyway.  A conserved variable, which
        does have an equation, keeps its line.

        >>> eqns = ['dC/dt = -C*k_2 - C*k_m1 + E*S*k_1',
        ... 'dE/dt = C*k_2 + C*k_m1 - E*S*k_1',
        ... 'dP/dt = C*k_2',
        ... 'dS/dt = C*k_m1 - E*S*k_1']
        >>> system = ODESystem.from_equations('\\n'.join(eqns))
        >>> print(system.to_tex())
        \\frac{dt}{dt} &= 1 \\\\
        \\frac{dC}{dt} &= - C k_{2} - C k_{-1} + E S k_{1} \\\\
        \\frac{dE}{dt} &= C k_{2} + C k_{-1} - E S k_{1} \\\\
        \\frac{dP}{dt} &= C k_{2} \\\\
        \\frac{dS}{dt} &= C k_{-1} - E S k_{1}

        >>> system.update_initial_conditions({'C': 'C_0'})
        >>> print(system.to_tex())
        \\frac{dt}{dt} &= 1 \\\\
        \\frac{dC}{dt} &= - C k_{2} - C k_{-1} + E S k_{1} \\\\
        \\frac{dE}{dt} &= C k_{2} + C k_{-1} - E S k_{1} \\\\
        \\frac{dP}{dt} &= C k_{2} \\\\
        \\frac{dS}{dt} &= C k_{-1} - E S k_{1} \\\\
        C\\left(0\\right) &= C_{0}

        >>> system.add_constraint('K_m', '(k_m1 + k_2) / k_1')
        >>> print(system.to_tex())
        \\frac{dt}{dt} &= 1 \\\\
        \\frac{dC}{dt} &= - C k_{2} - C k_{-1} + E S k_{1} \\\\
        \\frac{dE}{dt} &= C k_{2} + C k_{-1} - E S k_{1} \\\\
        \\frac{dP}{dt} &= C k_{2} \\\\
        \\frac{dS}{dt} &= C k_{-1} - E S k_{1} \\\\
        C\\left(0\\right) &= C_{0} \\\\
        K_{m} &= \\frac{k_{2} + k_{-1}}{k_{1}}

        '''
        line_template = '\\frac{{d{}}}{{d{}}} &= {}'
        # `expr is not None` is the test for "this variable has an equation".  A parameter
        # is stored as None precisely because it has none, so skipping the None's drops the
        # parameters and keeps every real equation -- including a conserved variable's
        # `= 0`, which is an equation and must appear.  Testing `expr != 0` instead would
        # delete the conservation law along with the parameters.
        lines = [line_template.format(var_to_tex(var), var_to_tex(self.independent_variable), expr_to_tex(expr))
                 for var, expr in zip(self.symbols, self._derivatives) if expr is not None]
        # Here, by contrast, None only means that `self.initial_conditions` has no entry
        # for this variable.
        for v in self.dependents:
            init_cond = self.initial_conditions.get(v)
            if init_cond is not None:
                lines.append('{}\\left(0\\right) &= {}'.format(var_to_tex(v), expr_to_tex(init_cond)))
        for eqn in self.constraints:
            lines.append('{} &= {}'.format(expr_to_tex(eqn.lhs), expr_to_tex(eqn.rhs)))
        return ' \\\\\n'.join(lines)

    @classmethod
    def from_tex(cls, tex, is_reduced=False, parameters=None, primes='name',
                 sub_defns=False):
        """
        Given the LaTeX of a system of differential equations, return a ODESystem of it.

        Every line of the tex is one of three things, and which one is decided by what the
        line mentions rather than by where it sits:

        - a **derivative**, which is an equation of the system;
        - a relation among **parameters** only, which becomes a constraint;
        - a **definition** of something in terms of a dependent variable, which cannot be
          a constraint -- a constraint is a statement about parameters -- and so is either
          substituted away or refused.

        Args:
            tex (str): LaTeX
            parameters (iter of str or sympy.Symbol, optional): Variables to treat as
                parameters, discarding the zero equations they were given.
            primes (str): ``'name'`` (the default) or ``'derivative'``; see
                :func:`desr.tex_tools.tex_to_sympy_one_line`.
            sub_defns (bool): Substitute definitions away rather than refusing them.
                :const:`False` by default, so that a system quietly different from the one
                written down is not built by accident.

        Returns:
            ODESystem: System of ODEs.

        Raises:
            TexParseError: If a line is a definition and `sub_defns` was not asked for.
            NotImplementedError: If a definition refers to another definition.

        >>> eqns = ['\\frac{dE}{dt} &= - k_1 E S + k_{-1} C + k_2 C \\\\',
        ... '\\frac{dS}{dt} &= - k_1 E S + k_{-1} C \\\\',
        ... '\\frac{dC}{dt} &= k_1 E S - k_{-1} C - k_2 C \\\\',
        ... '\\frac{dP}{dt} &= k_2 C']
        >>> ODESystem.from_tex('\\n'.join(eqns))
        dt/dt = 1
        dC/dt = -C*k_2 - C*k_m1 + E*S*k_1
        dE/dt = C*k_2 + C*k_m1 - E*S*k_1
        dP/dt = C*k_2
        dS/dt = C*k_m1 - E*S*k_1
        parameters: k_1, k_2, k_m1

        Todo:
            * Allow initial conditions to be set from tex.
        """
        sympification = tex_to_sympy(tex, primes=primes)
        derivative_dict = {}
        indep_var = None

        relations = {}

        for eqn in sympification:
            if not isinstance(eqn.lhs, sympy.Derivative):
                # Not a derivative, so it relates existing quantities.  Which kind of
                # relation it is cannot be known until every derivative line has been
                # read, because that is what says which variables are dependent.
                relations[eqn.lhs] = eqn.rhs

            else:
                # add to the dict from which we will construct the system
                derivative_dict[eqn.lhs.args[0]] = eqn.rhs

                # Check we always have the same independent variable.
                if indep_var is None:
                    indep_var = eqn.lhs.args[1]
                else:
                    if indep_var != eqn.lhs.args[1]:
                        raise ValueError('Must be ordinary differential equation, but found two independent variables {} and {}.'.format(indep_var, eqn.lhs.args[1]))

        # make the system which we will eventually return.
        system = cls.from_dict(deriv_dict=derivative_dict, is_reduced=is_reduced,
                               parameters=parameters)

        constraints, definitions = cls._classify_relations(system, relations, sub_defns)

        if definitions:
            system = system.subs(definitions)

        # Add the constraints.  These mention parameters only, so they cannot be refused.
        for lhs, rhs in constraints.items():
            system.add_constraint(lhs, rhs)

        # Add the initial conditions
        # todo: code here


        # cuz i'd prefer to NOT be insane.
        system._variable_sanity_check()

        return system

    @staticmethod
    def _classify_relations(system, relations, sub_defns):
        '''
        Split the non-derivative lines of some tex into constraints and definitions.

        A relation mentioning only parameters is a **constraint**: it says something true
        of the constants of the model, and the system can carry it as one.  A relation
        mentioning a dependent variable is a **definition** -- it gives a quantity in terms
        of things that move -- and :meth:`add_constraint` rightly refuses those, so it is
        either substituted away or reported.

        Args:
            system (ODESystem): The system built from the derivative lines, which is what
                knows which variables are dependent.
            relations (dict): {lhs: rhs} for each non-derivative line.
            sub_defns (bool): Whether a definition may be substituted away.

        Returns:
            tuple: (constraints, definitions), both {lhs: rhs} dicts.

        Raises:
            TexParseError: On a definition when `sub_defns` is False.
            NotImplementedError: On a definition that refers to another definition.
        '''
        dependent = set(system.dependents)

        constraints, definitions = {}, {}
        for lhs, rhs in relations.items():
            if expressions_to_variables([lhs, rhs]) & dependent:
                definitions[lhs] = rhs
            else:
                constraints[lhs] = rhs

        if definitions and not sub_defns:
            raise TexParseError(
                'The relation{} {} mention{} a dependent variable, so {} a definition '
                'rather than a constraint -- a constraint is a statement about '
                'parameters.  Pass sub_defns=True to substitute {} away, which is what a '
                'definition is for.'.format(
                    '' if len(definitions) == 1 else 's',
                    ', '.join('{} = {}'.format(lhs, rhs)
                              for lhs, rhs in sorted(definitions.items(), key=str)),
                    's' if len(definitions) == 1 else '',
                    'it is' if len(definitions) == 1 else 'they are',
                    'it' if len(definitions) == 1 else 'them'))

        # A definition whose right-hand side names another definition would have an answer
        # that depended on the order the two were applied, and a pair that named each other
        # would not have one at all.  Neither is guessed at.
        defined = set(definitions)
        for lhs, rhs in sorted(definitions.items(), key=str):
            chained = expressions_to_variables([rhs]) & defined
            if chained:
                raise NotImplementedError(
                    'The definition of {} refers to {}, which {} also defined.  '
                    'Substituting definitions into one another is not implemented: the '
                    'result would depend on the order they were applied in, and mutually '
                    'referring definitions have no result at all.  Substitute by hand '
                    'with subs, in the order you intend.'.format(
                        lhs, ', '.join(sorted(map(str, chained))),
                        'is' if len(chained) == 1 else 'are'))

        return constraints, definitions

    def _expressions(self):
        """
        A helper function that produces a list of expressions combining all of the expression types

        Args:

        Returns:
            List of sympy expressions
        """

        # first, the derivatives.  The independent variable is excluded by being the
        # independent variable, not by its derivative happening to equal 1: a dependent
        # variable moving at unit rate has a real equation, and `t * 1 / x` is a real
        # constraint on the scaling.
        exprs = [self._indep_var * expr / var for var, expr in self.derivative_dict.items()
                 if var != self.independent_variable]

        # next, initial conditions
        exprs.extend([var / init_cond for var, init_cond in self.initial_conditions.items()])

        # then, the constraints on constant variables
        exprs.extend([eq.lhs / eq.rhs for eq in self.constraints])

        # finally, the things that are required to have the same unit
        exprs.extend([expr/var for var,expr in self._requires_same_unit.items()])
        
        return exprs

    def _variable_sanity_check(self):
        """
        Makes sure that all the variables in the system are accounted for
        """

        exprs = self._expressions()

        expr_vars = expressions_to_variables(exprs)
        
        unlisted_vars = expr_vars - set(self.symbols)
        

        if unlisted_vars:
            raise RuntimeError(f"there are variables in the expressions that are not in the list of variables: {unlisted_vars}")

        for v in self.symbols:
            found_vars = expressions_to_variables([v]) # convert.  constants don't make it through.

            if not found_vars:
                raise RuntimeError(f"variable `{v}` doesn't appear to be a variable")

            if list(found_vars)[0] != sympy.sympify(v):
                raise RuntimeError(f'variable `{v}` appears to not actually be a single variable.  have {found_vars}')


    def exponent_matrix(self):
        '''
        Determine the 'exponent' or 'power' matrix of the system, denoted by :math:`K` in the literature,
        by gluing together the power matrices of each derivative, initial condition, constraint, and require-same-variable.

        In particular, it concatenates :math:`K_{\\left(\\frac{t}{x} \\cdot \\frac{dx}{dt}\\right)}` for :math:`x` in :attr:`~variables`,
        where :math:`t` is the independent variable.

        >>> eqns = '\\n'.join(['ds/dt = -k_1*e_0*s + (k_1*s + k_m1)*c',
        ... 'dc/dt = k_1*e_0*s - (k_1*s + k_m1 + k_2)*c'])
        >>> system = ODESystem.from_equations(eqns)
        >>> system.symbols
        (e_0, k_1, k_2, k_m1, t, c, s)
        >>> system.exponent_matrix()
        Matrix([
        [0, 0, 0,  1, 1, 0,  0],
        [1, 0, 0,  1, 1, 1,  0],
        [0, 1, 0,  0, 0, 0,  0],
        [0, 0, 1,  0, 0, 0,  1],
        [1, 1, 1,  1, 1, 1,  1],
        [0, 0, 0, -1, 0, 1,  1],
        [1, 0, 0,  1, 0, 0, -1]])

        While we get a different answer to the example in the paper, this is just due to choosing our reference exponent in a different way.

        Todo:
            * Change the code to agree with the paper.

        The initial-condition constant is appended, so it lands after the dependent
        variables rather than among the parameters where it belongs; see
        :meth:`~desr.ode_system.ODESystem.default_order_variables`.

        >>> system.update_initial_conditions({'s': 's_0'})
        >>> system.exponent_matrix()
        Matrix([
        [0, 0, 0,  1, 1, 0,  0,  0],
        [1, 0, 0,  1, 1, 1,  0,  0],
        [0, 1, 0,  0, 0, 0,  0,  0],
        [0, 0, 1,  0, 0, 0,  1,  0],
        [0, 0, 0,  0, 0, 0,  0, -1],
        [1, 1, 1,  1, 1, 1,  1,  0],
        [0, 0, 0, -1, 0, 1,  1,  0],
        [1, 0, 0,  1, 0, 0, -1,  1]])
        '''

        exprs = self._expressions()

        matrices = [rational_expr_to_exponent_matrix(expr, self.symbols) for expr in exprs]
        out = sympy.Matrix.hstack(*matrices)
        assert out.shape[0] == len(self.symbols)
        return out

    def maximal_scaling_matrix(self):
        '''
        Determine the maximal scaling matrix leaving this system invariant.

        Returns:
            sympy.Matrix: Maximal scaling matrix.


        >>> eqns = '\\n'.join(['ds/dt = -k_1*e_0*s + (k_1*s + k_m1)*c',
        ... 'dc/dt = k_1*e_0*s - (k_1*s + k_m1 + k_2)*c'])
        >>> system = ODESystem.from_equations(eqns)
        >>> system.maximal_scaling_matrix()
        Matrix([
        [0, -1, -1, -1, 1, 0, 0],
        [1, -1,  0,  0, 0, 1, 1]])
        '''
        exprs = self._expressions()

        return maximal_scaling_matrix(exprs, variables=self.symbols)

    def reorder_variables(self, variables):
        '''
        Store the new variable order, and reorder the equations according to the new order of variables.

        The order within the parameters is the one that changes the answer: the parameters
        placed **first** tend to be the ones normalised to 1, and the rest become the
        reduced constants.  "Tend to" rather than "are": :math:`r` of them are normalised,
        the leading ones are the natural candidates, but the Hermite structure decides which
        set actually can be, so a parameter that cannot participate is skipped over.  Put
        the ones you would rather see disappear at the front, then check
        :meth:`~desr.ode_translation.ODETranslator.translate_parameter_substitutions` to see
        which went.

        Args:
            variables (str, iter):
                Another ordering of the variables.

        >>> eqns = ['dz_1/dt = z_1*z_3', 'dz_2/dt = z_1*z_2 / (z_3 ** 2)']
        >>> system = ODESystem.from_equations('\\n'.join(eqns))
        >>> system.symbols
        (z_3, t, z_1, z_2)

        `z_3` is a parameter, so it comes first, has no equation, and its entry is
        :const:`None`.  The classification travels with the variable through a reordering.

        >>> system.derivatives
        (None, 1, z_1*z_3, z_1*z_2/z_3**2)

        >>> system.reorder_variables(['z_2', 'z_3', 't', 'z_1'])
        >>> system.symbols
        (z_2, z_3, t, z_1)
        >>> system.derivatives
        (z_1*z_2/z_3**2, None, 1, z_1*z_3)
        >>> system.parameters
        (z_3,)
        '''
        if isinstance(variables, str):
            if ' ' in variables:
                variables = variables.split(' ')
            else:
                variables = tuple(variables)

        variables = [sympy.Symbol(v) if isinstance(v, str) else v for v in variables]

        if set(variables) != set(self.symbols):
            raise ValueError('Mismatching variables:\n{} vs\n{}'.format(sorted(self.symbols, key=str), sorted(variables, key=str)))

        var_to_index = {var: i for i, var in enumerate(self.symbols)}
        column_shuffle = [var_to_index[new_var] for new_var in variables]

        self._variables = tuple( [self._variables[i] for i in column_shuffle])
        self._derivatives = tuple( [self._derivatives[i] for i in column_shuffle])

    def _describe_variable(self, variable):
        '''
        What kind of variable this is, as a phrase for an error message.
        '''
        if variable == self.independent_variable:
            return 'the independent variable'
        if variable in self.parameters:
            return 'a parameter'
        if variable in self.dependents:
            return 'a dependent variable'
        return 'not a variable of this system'

    def _block_reorder(self, order, current, what):
        '''
        Work out the new order of one block of variables from a partial list.

        The named ones lead, in the order given; the rest follow in the order they are in
        now.  That is the useful default, because the reason to name a variable at all is to
        move it towards the front.

        Args:
            order (str, iter): The variables to bring to the front of the block.  May be
                partial, and may be empty.
            current (tuple): The block as it stands.
            what (str): What this block is, for error messages.

        Returns:
            list: The whole block, reordered.
        '''
        if isinstance(order, str):
            order = order.split(' ') if ' ' in order else (order,)

        named = [sympy.Symbol(v) if isinstance(v, str) else v for v in order]

        seen = set()
        for var in named:
            if var in seen:
                raise ValueError('{} is named twice in the {} order.'.format(var, what))
            seen.add(var)
            if var not in current:
                if var in self.symbols:
                    raise ValueError(
                        '{} is not one of the {}; it is {}.  Use reorder_variables to move a '
                        'variable between blocks, or make_parameter to change what it '
                        'is.'.format(var, what, self._describe_variable(var)))
                raise ValueError('{} is not a variable of this system.  The {} are {}.'.format(
                    var, what, ', '.join(map(str, current))))

        return named + [var for var in current if var not in seen]

    def reorder_parameters(self, order):
        '''
        Reorder the parameters, leaving everything else where it is.

        The order within the parameters is the one that changes the answer -- the leading
        ones tend to be normalised to 1 -- and this is the way to set it without having to
        write out the independent variable and the dependents as well.  See
        :meth:`reorder_variables` for the "tend to" caveat.

        The list may be partial.  The parameters named come first, in the order given, and
        the rest follow in the order they are already in, so naming one parameter is enough
        to say "normalise this one if you can".

        Args:
            order (str, iter): The parameters to bring to the front.

        The layout is left canonical: parameters, independent variable, dependents.  A
        system deliberately put in some other order by :meth:`reorder_variables` is
        canonicalised by this call, which is the point of it -- say what is meaningful, and
        let the library place the blocks.

        >>> eqns = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
        ...         'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
        >>> system = ODESystem.from_equations(eqns)
        >>> system.parameters
        (e_0, k_1, k_2, k_m1)

        Name one, and it moves to the front; the others keep their order.

        >>> system.reorder_parameters('k_2')
        >>> system.symbols
        (k_2, e_0, k_1, k_m1, t, c, s)

        Name several, and they lead in the order given.

        >>> system.reorder_parameters(['k_m1', 'e_0'])
        >>> system.parameters
        (k_m1, e_0, k_2, k_1)

        A dependent variable is refused, and the message says what to use instead:

        >>> system.reorder_parameters('s')
        Traceback (most recent call last):
            ...
        ValueError: s is not one of the parameters; it is a dependent variable.  Use reorder_variables to move a variable between blocks, or make_parameter to change what it is.
        '''
        parameters = self._block_reorder(order, self.parameters, 'parameters')
        self.reorder_variables(parameters + [self.independent_variable]
                               + list(self.dependents))

    def reorder_dependents(self, order):
        '''
        Reorder the dependent variables, leaving everything else where it is.

        This sets which dependent gets which name in the reduced system -- the `nu`
        numbering, or whatever :data:`~desr.ode_translation.NAMING_SCHEME_KEYS` says -- and
        the order the reduced equations are printed in.  It does not change which parameters
        are normalised; that is :meth:`reorder_parameters`.

        The list may be partial, with the named ones leading.

        Args:
            order (str, iter): The dependent variables to bring to the front.

        >>> eqns = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
        ...         'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']
        >>> system = ODESystem.from_equations(eqns)
        >>> system.dependents
        (c, s)
        >>> system.reorder_dependents('s')
        >>> system.symbols
        (e_0, k_1, k_2, k_m1, t, s, c)

        The independent variable is not a dependent variable, and saying so is more use
        than a mismatched-variables error:

        >>> system.reorder_dependents('t')
        Traceback (most recent call last):
            ...
        ValueError: t is not one of the dependent variables; it is the independent variable.  Use reorder_variables to move a variable between blocks, or make_parameter to change what it is.
        '''
        dependents = self._block_reorder(order, self.dependents,
                                         'dependent variables')
        self.reorder_variables(list(self.parameters) + [self.independent_variable] + dependents)

    def default_order_variables(self):
        '''
        Reorder the variables into (parameters, independent variable, dependent variables),
        which is the order the parameter scheme of translation requires.
        Variables of the same type are sorted by their string representations.

        A conserved variable is a dependent variable, so it sorts with the dependents rather
        than with the parameters.

        This is Hubert & Labahn's ordering: they write the system as :math:`F(c, t, z)` with
        :math:`p` parameters, the one independent variable, and :math:`q` dependent
        variables, so :math:`n = p + 1 + q`.  It goes with the Hermite normal form
        convention -- see :data:`desr.matrix_normal_forms._CONVENTION` -- rather than being
        a free choice beside it: the identity block of :math:`V_n` sits on the **last**
        :math:`1 + q` rows, and only this order puts :math:`t` and the dependent variables
        there.  Order the parameters last instead and the identity lands on *parameters*,
        which satisfies the block conditions while giving each parameter its own invariant
        and smearing the independent variable across several -- not a parameter reduction
        at all.

        desr used to order these as (independent, dependent, parameters), alongside the
        mirrored normal form.


        >>> eqns = ['dz_1/dt = z_1*z_3', 'dz_2/dt = z_1*z_2 / (z_3 ** 2)']
        >>> system = ODESystem.from_equations('\\n'.join(eqns))
        >>> system.symbols
        (z_3, t, z_1, z_2)

        >>> system.reorder_variables(['z_2', 'z_3', 't', 'z_1'])
        >>> system.symbols
        (z_2, z_3, t, z_1)

        >>> system.default_order_variables()
        >>> system.symbols
        (z_3, t, z_1, z_2)
        '''
        dep_var = sorted(self.dependents, key=str)
        parameters = sorted(self.parameters, key=str)

        # Order variables as parameters, independent, dependent
        variables = parameters + [self.independent_variable] + dep_var
        assert len(variables) == len(set(variables)) # assures no duplicates
        self.reorder_variables(variables=variables)

def parse_de(diff_eq, indep_var='t'):
    ''' Parse a first order ordinary differential equation and return (variable of derivative, rational function

        >>> parse_de('dn/dt = n( r(1 - n/K) - kp/(n+d) )')
        (n, n(-kp/(d + n) + r(1 - n/K)))

        >>> parse_de('dp/dt==sp(1 - hp / n)')
        (p, sp(-hp/n + 1))
    '''
    diff_eq = diff_eq.strip()
    match = re.match(r'd([a-zA-Z0-9_]*)/d([a-zA-Z0-9_]*)\s*=*\s*(.*)', diff_eq)
    if match is None:
        raise ValueError("Invalid differential equation: {}".format(diff_eq))
    if match.group(2) != indep_var:
        raise ValueError('We only work in ordinary DEs in {}'.format(indep_var))
    # Read every name as a symbol, so that a species called S, E, I, beta or pi is an
    # ordinary variable rather than one of sympy's own objects.
    return to_symbol(match.group(1)), sympify_expression(match.group(3))

def rational_expr_to_exponent_matrix(expr, variables):
    '''
    Take a rational expression and determine the power matrix wrt an ordering on the variables, as on page 497 of
    Hubert-Labahn.

    >>> exprs = list(map(sympy.sympify, "n*( r*(1 - n/K) - k*p/(n+d) );s*p*(1 - h*p / n)".split(';')))
    >>> variables = sorted(expressions_to_variables(exprs), key=str)
    >>> variables
    [K, d, h, k, n, p, r, s]
    >>> rational_expr_to_exponent_matrix(exprs[0], variables)
    Matrix([
    [0, -1, -1, 0, 0,  0],
    [0,  1,  0, 1, 0,  1],
    [0,  0,  0, 0, 0,  0],
    [1,  0,  0, 0, 0,  0],
    [0,  1,  2, 0, 1, -1],
    [1,  0,  0, 0, 0,  0],
    [0,  1,  1, 1, 1,  0],
    [0,  0,  0, 0, 0,  0]])

    >>> rational_expr_to_exponent_matrix(exprs[1], variables)
    Matrix([
    [ 0, 0],
    [ 0, 0],
    [ 1, 0],
    [ 0, 0],
    [-1, 0],
    [ 2, 1],
    [ 0, 0],
    [ 1, 1]])

    A zero expression -- the equation of a conserved variable, such as a catalyst in a
    reaction network -- has no monomials, and so constrains no scaling.  It contributes a
    column of zeros, one entry per variable, rather than nothing at all:

    >>> rational_expr_to_exponent_matrix(sympy.sympify(0), variables).T
    Matrix([[0, 0, 0, 0, 0, 0, 0, 0]])
    '''
    expr = expr.cancel()
    num, denom = expr.as_numer_denom()
    num_const, num_terms = num.as_coeff_add()
    denom_const, denom_terms = denom.as_coeff_add()
    num_terms = sorted(num_terms, key=str)
    denom_terms = sorted(denom_terms, key=str)

    if denom_const != 0:
        ref_power = 1
        # If we have another constant in the numerator, add it onto the terms for processing.
        if num_const != 0:
            num_terms = list(num_terms)
            num_terms.append(num_const)
    else:
        if num_const != 0:
            ref_power = 1
        else:
            denom_terms = list(denom_terms)

            # Find the lowest power
            ref_power = min(denom_terms, key=lambda x: list(map(abs, monomial_to_powers(x, variables))))

            denom_terms.remove(ref_power)  # Use the last term of the denominator as our reference power

    powers = []
    for mon in itertools.chain(num_terms, denom_terms):
        powers.append(monomial_to_powers(mon / ref_power, variables))

    if not powers:
        # A zero expression has no monomials.  It constrains no scaling, so the honest
        # answer is a column of zeros.  Building a matrix from the empty list instead would
        # give a 0x0 one, which ignores `variables` and cannot be glued to the others:
        # sympy's hstack absorbs a zero-row block on the left but raises on the right, so
        # whether a system reduced at all would depend on where its conserved variables
        # fell in the ordering -- that is, on what they were named.
        return sympy.zeros(len(variables), 1)

    powers = sympy.Matrix(powers).T
    return powers

def maximal_scaling_matrix(exprs, variables=None):
    ''' Determine the maximal scaling matrix leaving this system invariant, in row Hermite normal form.

    Args:
        exprs (iter): Iterable of sympy.Expressions.
        variables: An ordering on the variables. If None, sort according to the string representation.
    Returns:
        sympy.Matrix

    >>> exprs = ['z_1*z_3', 'z_1*z_2 / (z_3 ** 2)']
    >>> exprs = list(map(sympy.sympify, exprs))
    >>> maximal_scaling_matrix(exprs)
    Matrix([[-1, 3, 1]])

    >>> exprs = ['(z_1 + z_2**2) / z_3']
    >>> exprs = list(map(sympy.sympify, exprs))
    >>> maximal_scaling_matrix(exprs)
    Matrix([[2, 1, 2]])
    '''
    if variables is None:
        variables = sorted(expressions_to_variables(exprs), key=str)
    matrices = [rational_expr_to_exponent_matrix(expr, variables) for expr in exprs]
    exponent_matrix = sympy.Matrix.hstack(*matrices)
    assert exponent_matrix.shape[0] == len(variables)

    hermite_rform, multiplier_rform = hnf_row(exponent_matrix)

    # Find the non-zero rows at the bottom
    row_is_zero = [all([i == 0 for i in row]) for row in hermite_rform.tolist()]
    # Make sure they all come at the end
    num_nonzero = sum(map(int, row_is_zero))
    if num_nonzero == 0:
        return sympy.zeros(1, len(variables))
    assert hermite_rform[-num_nonzero:, :].is_zero_matrix

    # Make sure we have the right number of columns
    assert multiplier_rform.shape[1] == len(variables)
    # Return the last num_nonzero rows of the Hermite multiplier
    return hnf_row(multiplier_rform[-num_nonzero:, :])[0]

__all__ = ['ODESystem','parse_de','rational_expr_to_exponent_matrix','maximal_scaling_matrix']

if __name__ == '__main__':
    import doctest
    doctest.testmod()