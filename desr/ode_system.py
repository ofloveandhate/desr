import itertools
import re

import sympy
from sympy.abc import _clash1

from desr.matrix_normal_forms import hnf_col, hnf_row, normal_hnf_col
from desr.sympy_helper import expressions_to_variables, unique_array_stable, monomial_to_powers
from desr.tex_tools import expr_to_tex, var_to_tex, tex_to_sympy

class ODESystem(object):
    '''
    A system of differential equations.

    The main attributes are :attr:`~desr.ode_system.ODESystem.variables` and :attr:`~desr.ode_system.ODESystem.derivatives`.
    :attr:`~desr.ode_system.ODESystem.variables` is an ordered tuple of variables, which includes the independent variable.
    :attr:`~desr.ode_system.ODESystem.derivatives` is an ordered tuple of the same length that contains the derivatives with respect to :attr:`~desr.ode_system.ODESystem.indep_var`.

    Args:
        variables (tuple of sympy.Symbol): Ordered tuple of variables.
        derivatives (tuple of sympy.Expression): Ordered tuple of derivatives.
        indep_var (sympy.Symbol, optional): Independent variable we are differentiating with respect to.
        initial_conditions (tuple of sympy.Symbol): The initial values of non-constant variables
        same_units (dict of sympy.Symbol : sympy.Expression): lists of things that must have the same units.
    '''

    def __init__(self, variables, derivatives, indep_var=None, initial_conditions=None, same_units=None, constraints=None, is_reduced=False):
        self._variables = tuple(variables)
        self._derivatives = tuple(derivatives)

        self._indep_var = sympy.var('t') if indep_var is None else indep_var

        # initialize some class member storage
        self._initial_conditions = {}
        self._constraints = []
        self._requires_same_unit = {} # a dict of {Symbol : Expression}
        self._is_reduced = is_reduced

        # sanity checks
        assert len(self._variables) == len(self._derivatives)
        assert self.derivatives[self.indep_var_index] == sympy.sympify(1)

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
        self_var = sorted(self.variables, key=str)
        other_var = sorted(other.variables, key=str)
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
        system = ODESystem(self._variables, self._derivatives, indep_var=self._indep_var, is_reduced = self.is_reduced, same_units=self._requires_same_unit)
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
        dc_0/dtau = 0
        dc_1/dtau = 0
        dx_0/dtau = 0
        x(0) = x_0

        >>> system.variables
        (tau, x, y, c_0, c_1, x_0)
        """
        if not isinstance(new_indep, sympy.Symbol):
            new_indep = sympy.sympify(new_indep, locals=_clash1)

            expressions_to_variables

        to_sub = {self._indep_var:new_indep}
        
        self._derivatives = [d.subs(to_sub) if d else d for d in self._derivatives]
        self._constraints = [c.subs(to_sub) for c in self.constraints]
        self._initial_conditions = {k:v.subs(to_sub) for k,v in self._initial_conditions.items()}

        ind = self.indep_var_index
        self._variables = self._variables[:ind] + (new_indep,) + self._variables[ind+1:]
        self._indep_var = new_indep

    @property
    def indep_var(self):
        """
        Return the independent variable.

        Returns:
            sympy.Symbol: The independent variable, which we are differentiating with respect to.
        """
        return self._indep_var


    @property
    def indep_var_index(self):
        """
        Return the independent variable index.

        Return:
             int: The index of :py:attr:`~indep_var` in :py:attr:`~self.variables`.
        """
        return self.variables.index(self.indep_var)

    @property
    def variables(self):
        '''
        Return ALL the variables appearing in the system.

        Returns:
            tuple: Ordered tuple of variables appearing in the system.
        '''
        return self._variables

    @property
    def constant_variables(self):
        '''
        Return the constant variables - specifically those which have a None derivative.

        Returns:
            tuple: The constant variables.
        '''
        return tuple(var for var, deriv in zip(self.variables, self._derivatives) if deriv is None)

    @property
    def non_constant_variables(self):
        '''
        Return the non-constant non-independent variables - specifically those which have a derivative that isn't None or 1.

        Returns:
            tuple: The non-constant non-independent variables.

        >>> _input = {'x': 'c_0*x*y', 'y': 'c_1*(1-x)*(1-y)*t'}
        >>> _input = {sympy.Symbol(k): sympy.sympify(v) for k, v in _input.items()}
        >>> system = ODESystem.from_dict(_input)
        >>> system.non_constant_variables
        (x, y)
        '''
        return tuple(var for var, deriv in zip(self.variables, self._derivatives) if
                     ((deriv is not None) and (deriv != 1)))

    @property
    def num_nonconstants(self):
        return len(self.non_constant_variables)

    @property
    def num_constants(self):
        '''
        Return the number of constant variables - specifically those which have a :const:`None` derivative

        Returns:
            int: Number of non-constant variables.
        '''
        return len(self.constant_variables)
        


    @property
    def derivatives(self):
        ''' Getter for an ordered tuple of expressions representing the derivatives of self.variables.

        Returns:
            tuple: Ordered tuple of sympy.Expressions.
        '''
        return [expr if expr is not None else sympy.sympify(0) for expr in self._derivatives]

    @property
    def derivative_dict(self):
        '''
        Return a variable: expr mapping, filtering out the :const:`None`'s in expr.

        Returns:
            dict: Keys are non-constant variables, value is the derivative with respect to the independent variable.
        '''
        return dict(filter(lambda x: x[1] is not None, zip(self.variables, self._derivatives)))

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
                for v in self.variables:
                    d = sympy.collect(d,v)
            new_derivs.append(d)

        self._derivatives = new_derivs


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
        dc_0/dt = 0
        dc_1/dt = 0
        dx_0/dt = 0
        x(0) = x_0

        >>> system.update_initial_conditions({'c_0': 'k'})
        Traceback (most recent call last):
            ...
        ValueError: Cannot set initial condition k for variable c_0 with derivative None.
        
        >>> system
        dt/dt = 1
        dx/dt = c_0*x*y
        dy/dt = c_1*t*(1 - x)*(1 - y)
        dc_0/dt = 0
        dc_1/dt = 0
        dx_0/dt = 0
        x(0) = x_0
        
        >>> system.is_reduced
        False

        >>> system.update_initial_conditions({'x': '1'})
        Traceback (most recent call last):
            ...
        ValueError: initial condition `1` doesn't appear to be a variable expression.  if you want it to be a numeric constant like 1 or 0, accomplish this via substitution (after reduction)
        
        >>> system.update_initial_conditions({'x': '0'})
        Traceback (most recent call last):
            ...
        ValueError: initial condition `0` doesn't appear to be a variable expression.  if you want it to be a numeric constant like 1 or 0, accomplish this via substitution (after reduction)

        >>> system
        dt/dt = 1
        dx/dt = c_0*x*y
        dy/dt = c_1*t*(1 - x)*(1 - y)
        dc_0/dt = 0
        dc_1/dt = 0
        dx_0/dt = 0
        x(0) = x_0
        '''
        for variable, init_cond in initial_conditions.items():

            # convert through sympy
            if not isinstance(variable, sympy.Symbol):
                variable = sympy.Symbol(variable)

            if isinstance(init_cond, str):
                init_cond = sympy.sympify(init_cond, locals=_clash1)

            # We can only set initial conditions for non-constant variables we already know about.
            if variable not in self.non_constant_variables:
                raise ValueError('Cannot set initial condition {} for constant variable {} with derivative {}.'.format(init_cond,
                                                                                                              variable,
                                                                                                              self.derivative_dict.get(variable)))



            ics_vars = expressions_to_variables([init_cond]) # convert. 



            # reduced systems are allowed to have 1's (and potentially 0's) as initial conditions
            # but non-reduced systems are NOT.  Otherwise, there might be 1's with different units, and this is nonsense.
            if not self.is_reduced:

                if len(ics_vars) != 1:
                    raise ValueError(f'an initial condition for a non-reduced system should just be a single variable.  you have {ics_vars}')

                if init_cond.is_constant():
                    raise ValueError(f"initial condition `{init_cond}` doesn't appear to be a variable expression.  Non-reduced systems can only have variables as initial conditions.  If you want {variable}(0) to be a numeric constant like 1 or 0, accomplish this via substitution (after reduction)")

            # silviana: i think we should guard against having 1 be a variable...
            if (init_cond not in self.variables) and not init_cond.is_constant(): # checking emptiness on ics_vars
                self._variables = tuple(list(self._variables) + [init_cond])
                self._derivatives = tuple(list(self._derivatives) + [None])

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
                expr = sympy.sympify(expr, locals=_clash1)

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
        dc_0/dt = 0
        dc_1/dt = 0

        >>> system.add_constraint('c_2', 'c_0 + c_1')
        >>> system
        dt/dt = 1
        dx/dt = c_0*x*y
        dy/dt = c_1*(1 - x)*(1 - y)
        dc_0/dt = 0
        dc_1/dt = 0
        dc_2/dt = 0
        c_2 == c_0 + c_1

        >>> system.add_constraint('c_2', 'c_0 + x')
        Traceback (most recent call last):
            ...
        ValueError: Cannot add constraints on non-constant parameters set([x]). This would make an interesting project though...

        >>> system.add_constraint('c_0', 0)
        Traceback (most recent call last):
            ...
        ValueError: Cannot express equality with 0.
        '''
        if isinstance(lhs, str):
            lhs = sympy.sympify(lhs, locals=_clash1)
        if isinstance(rhs, str):
            rhs = sympy.sympify(rhs, locals=_clash1)
        if (lhs == 0) or (rhs == 0):
            raise ValueError(f'Cannot express equality constraint with 0. ({lhs} == {rhs})')

        # get the variables from the constraint
        all_constraint_vars, nonconst_vars, new_vars = self.equation_variables(lhs, rhs)

        # check if have non-constant variables present.  if so, raise with a helpful error message
        if nonconst_vars:
            raise ValueError('Cannot add constraints on non-constant parameters {}. '.format(nonconst_vars) +
                             '.  Try to use substitution to incorporate the constraints directly into the system before construction.  (To do this automatically would be small project...')

        # no non-constants, no problem.  update the recorded symbols
        self._variables = tuple(list(self._variables) + new_vars)
        self._derivatives = tuple(list(self._derivatives) + [None for _ in new_vars])
        self._constraints.append(sympy.Eq(lhs, rhs))

    def equation_variables(self, lhs, rhs):
        """
        compute the variables in the left and right hand side of a constraint:
        all, the non-constant ones that exist in the system, and the new ones.

        Args:
            lhs (sympy.Expr): The left hand side of the constraint.
            rhs (sympy.Expr): The right hand side of the constraint.

        Returns: length-three tuple of containers of variables: all, non-constant existing, new
        """
        all_vars = expressions_to_variables([lhs, rhs])
        nonconst_vars = all_vars.intersection(self.non_constant_variables)
        new_vars = sorted(all_vars.difference(set(self.variables)), key=str)

        return all_vars, nonconst_vars, new_vars

    def diff_subs(self, to_sub, expand_before=False, expand_after=True, factor_after=False, subs_constraints=False, new_symbols_are_constants=True):
        '''
        Make substitutions into the derivatives, returning a new system.

        Args:
            to_sub (dict): Dictionary of substitutions to make.
            expand_before (bool): Expand the sympy expression for each derivative before substitution.
            expand_after (bool): Expand the sympy expression for each derivative after substitution.
            factor_after (bool): Factorise the sympy expression for each derivative after substitution.
            subs_constraints (bool): Perform the substitutions into the initial constraints.
            new_symbols_are_constants (bool): Treat newly encountered symbols as constants. If a substitution is a variable renaming, this does not apply.  I have no idea what to do if this is false, so it will generate exceptions..
        
        Returns:
            ODESystem: System with substitutions carried out.

        >>> eqns = ['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)']
        >>> system = ODESystem.from_equations(eqns)
        >>> system.diff_subs({'1-x': 'z'}, expand_before=False, expand_after=False, factor_after=False)
        Traceback (most recent call last):
        ...
        ValueError: unable to make substitution 1 - x -> z.
        >>> system.diff_subs({'1-x': 'z'}, expand_before=True, expand_after=False, factor_after=False)
        Traceback (most recent call last):
        ...
        ValueError: unable to make substitution 1 - x -> z.

        >>> system.diff_subs({'x': '1-z'}, expand_before=True, expand_after=True, factor_after=False)
        Traceback (most recent call last):
        ...
        ValueError: unable to make substitution x -> 1 - z.

        >>> system.add_constraint('c_0', 'c_1**2')
        >>> system.diff_subs({'c_0': '1'}, subs_constraints=False)
        dt/dt = 1
        dx/dt = x*y
        dy/dt = c_1*x*y - c_1*x - c_1*y + c_1
        dc_0/dt = 0
        dc_1/dt = 0
        c_0 == c_1**2
        >>> system.diff_subs({'c_0': '1'}, subs_constraints=True)
        dt/dt = 1
        dx/dt = x*y
        dy/dt = c_1*x*y - c_1*x - c_1*y + c_1
        dc_1/dt = 0
        1 == c_1**2
        '''

        # sympify the dict of things to substitute
        to_sub = {sympy.sympify(lhs, locals=_clash1): sympy.sympify(rhs, locals=_clash1) for lhs, rhs in to_sub.items()}

        variable_renamings = []
        
        for lhs,rhs in to_sub.items():
            lhs_vars = set(expressions_to_variables([lhs]))
            rhs_vars = set(expressions_to_variables([rhs]))

            
            lhs_nonconst_vars = lhs_vars.intersection(set(self.non_constant_variables))

            # make sure we can actually do the substitution
            if lhs_nonconst_vars:
                # in this conditional, we're trying to substitute away a non-constant variable.  can only do, if the substitution is just a renaming.  anything else must be done by the user

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
        variables = list(self.variables)
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
        
        # record that we have them, for the new post-substitution system
        new_vars = variables + missing_vars # this is a list extension operation

        if new_symbols_are_constants:
            for x in missing_vars:
                new_derivs.append(None)
        else:
            raise NotImplementedError(f"don't want new symbols from substitution to be constant, but don't know what else to do.  the new symbols: {missing_vars}")


        # remove variables that don't appear in any expressions

        variables_to_inspect = new_vars # because we're about to modify `new_vars` in the loop below
        extant_vars = set(expressions_to_variables( new_derivs + constraints + list(initial_conditions.values()) + [self.indep_var] + list(self.non_constant_variables)))

        for v in variables_to_inspect:

            if v not in extant_vars:
                var_loc = new_vars.index(v)
                
                # omit using list slicing.  must omit from both the variable ordering and from the derivatives (which are stored impilictly in the order of the variables)
                new_vars = new_vars[:var_loc] + new_vars[var_loc+1:]
                new_derivs = new_derivs[:var_loc] + new_derivs[var_loc+1:]

        # we're finally ready to construct the new post-substitution System. Huzzah.
        subs_system = ODESystem(new_vars, new_derivs,
                                initial_conditions=initial_conditions,
                                indep_var=self.indep_var,
                                constraints=constraints,
                                is_reduced = self.is_reduced)

        subs_system._variable_sanity_check()
        return subs_system

    @classmethod
    def from_equations(cls, equations, indep_var=sympy.var('t'), initial_conditions=None, is_reduced=False):
        '''
        Instantiate from multiple equations.

        Args:
            equations (str, iter of str): Equations of the form "dx/dt = expr", optionally seperated by :code:`\\n`.
            indep_var (sympy.Symbol): The independent variable, usually :code:`t`.
            initial_conditions (tuple of sympy.Symbol): The initial values of non-constant variables

        Returns:
            ODESystem: System of equations.

        >>> eqns = ['dx/dt = c_0*x*y', 'dy/dt = c_1*(1-x)*(1-y)']
        >>> ODESystem.from_equations(eqns)
        dt/dt = 1
        dx/dt = c_0*x*y
        dy/dt = c_1*(1 - x)*(1 - y)
        dc_0/dt = 0
        dc_1/dt = 0
        >>> eqns = '\\n'.join(['dy/dx = c_0*x*y', 'dz/dx = c_1*(1-y)*z**2'])
        >>> ODESystem.from_equations(eqns, indep_var=sympy.Symbol('x'))
        dx/dx = 1
        dy/dx = c_0*x*y
        dz/dx = c_1*z**2*(1 - y)
        dc_0/dx = 0
        dc_1/dx = 0
        '''
        if isinstance(equations, str):
            equations = equations.strip().split('\n')

        deriv_dict = dict(map(lambda x: parse_de(x, indep_var=str(indep_var)), equations))
        system = cls.from_dict(deriv_dict=deriv_dict, indep_var=indep_var, initial_conditions=initial_conditions, is_reduced=is_reduced)
        system.default_order_variables()

        system._variable_sanity_check()

        return system

    @classmethod
    def from_dict(cls, deriv_dict, indep_var=sympy.var('t'), initial_conditions=None, is_reduced=False):
        '''
        Instantiate from a text of equations.

        Args:
            deriv_dict (dict): {variable: derivative} mapping.
            indep_var (sympy.Symbol): Independent variable, that the derivatives are with respect to.
            initial_conditions (tuple of sympy.Symbol): The initial values of non-constant variables


        Returns:
            ODESystem: System of ODEs.

        >>> _input = {'x': 'c_0*x*y', 'y': 'c_1*(1-x)*(1-y)'}
        >>> _input = {sympy.Symbol(k): sympy.sympify(v) for k, v in _input.items()}
        >>> ODESystem.from_dict(_input)
        dt/dt = 1
        dx/dt = c_0*x*y
        dy/dt = c_1*(1 - x)*(1 - y)
        dc_0/dt = 0
        dc_1/dt = 0

        >>> _input = {'y':  'c_0*x*y', 'z': 'c_1*(1-y)*z**2'}
        >>> _input = {sympy.Symbol(k): sympy.sympify(v) for k, v in _input.items()}
        >>> ODESystem.from_dict(_input, indep_var=sympy.Symbol('x'))
        dx/dx = 1
        dy/dx = c_0*x*y
        dz/dx = c_1*z**2*(1 - y)
        dc_0/dx = 0
        dc_1/dx = 0
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
                     is_reduced=is_reduced)

        system.default_order_variables()

        system._variable_sanity_check()

        return system

    def __repr__(self):
        lines = ['d{}/d{} = {}'.format(var, self.indep_var, expr) for var, expr in zip(self.variables, self.derivatives)]
        for v in self.non_constant_variables:
            init_cond = self.initial_conditions.get(v)
            if init_cond is not None:
                lines.append('{}(0) = {}'.format(v, init_cond))
        for eqn in self.constraints:
            lines.append('{} == {}'.format(eqn.lhs, eqn.rhs))
        return '\n'.join(lines)

    def to_tex(self):
        '''
        Returns:
            str: TeX representation.


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
        \\frac{dS}{dt} &= C k_{-1} - E S k_{1} \\\\
        \\frac{dk_{1}}{dt} &= 0 \\\\
        \\frac{dk_{2}}{dt} &= 0 \\\\
        \\frac{dk_{-1}}{dt} &= 0

        >>> system.update_initial_conditions({'C': 'C_0'})
        >>> print(system.to_tex())
        \\frac{dt}{dt} &= 1 \\\\
        \\frac{dC}{dt} &= - C k_{2} - C k_{-1} + E S k_{1} \\\\
        \\frac{dE}{dt} &= C k_{2} + C k_{-1} - E S k_{1} \\\\
        \\frac{dP}{dt} &= C k_{2} \\\\
        \\frac{dS}{dt} &= C k_{-1} - E S k_{1} \\\\
        \\frac{dk_{1}}{dt} &= 0 \\\\
        \\frac{dk_{2}}{dt} &= 0 \\\\
        \\frac{dk_{-1}}{dt} &= 0 \\\\
        \\frac{dC_{0}}{dt} &= 0 \\\\
        C\\left(0\\right) &= C_{0}

        >>> system.add_constraint('K_m', '(k_m1 + k_2) / k_1')
        >>> print(system.to_tex())
        \\frac{dt}{dt} &= 1 \\\\
        \\frac{dC}{dt} &= - C k_{2} - C k_{-1} + E S k_{1} \\\\
        \\frac{dE}{dt} &= C k_{2} + C k_{-1} - E S k_{1} \\\\
        \\frac{dP}{dt} &= C k_{2} \\\\
        \\frac{dS}{dt} &= C k_{-1} - E S k_{1} \\\\
        \\frac{dk_{1}}{dt} &= 0 \\\\
        \\frac{dk_{2}}{dt} &= 0 \\\\
        \\frac{dk_{-1}}{dt} &= 0 \\\\
        \\frac{dC_{0}}{dt} &= 0 \\\\
        \\frac{dK_{m}}{dt} &= 0 \\\\
        C\\left(0\\right) &= C_{0} \\\\
        K_{m} &= \\frac{k_{2} + k_{-1}}{k_{1}}

        '''
        line_template = '\\frac{{d{}}}{{d{}}} &= {}'
        lines = [line_template.format(var_to_tex(var), var_to_tex(self.indep_var), expr_to_tex(expr))
                 for var, expr in zip(self.variables, self.derivatives)]
        for v in self.non_constant_variables:
            init_cond = self.initial_conditions.get(v)
            if init_cond is not None:
                lines.append('{}\\left(0\\right) &= {}'.format(var_to_tex(v), expr_to_tex(init_cond)))
        for eqn in self.constraints:
            lines.append('{} &= {}'.format(expr_to_tex(eqn.lhs), expr_to_tex(eqn.rhs)))
        return ' \\\\\n'.join(lines)

    @classmethod
    def from_tex(cls, tex, is_reduced=False):
        """
        Given the LaTeX of a system of differential equations, return a ODESystem of it.

        Args:
            tex (str): LaTeX

        Returns:
            ODESystem: System of ODEs.

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
        dk_1/dt = 0
        dk_2/dt = 0
        dk_m1/dt = 0

        Todo:
            * Allow initial conditions to be set from tex.
        """
        sympification = tex_to_sympy(tex)
        derivative_dict = {}
        indep_var = None

        constraints = {}

        for eqn in sympification:
            if not isinstance(eqn.lhs, sympy.Derivative):
                # since it's not a derivative, it must be a constraint.

                # add to a local container, and then we'll process later
                constraints[eqn.lhs] = eqn.rhs

            else:
                # add to the dict from which we will construct the system
                derivative_dict[eqn.lhs.args[0]] = eqn.rhs

                # Check we always have the same independent variable.
                if indep_var is None:
                    indep_var = eqn.lhs.args[1]
                else:
                    if indep_var != eqn.lhs.args[1]:
                        raise ValueError('Must be ordinary differential equation, but found two independent variables {} and {}.'.format(indep_var, eq.lhs.args[1]))

        # make the system which we will eventually return.
        system = cls.from_dict(deriv_dict=derivative_dict, is_reduced=is_reduced)

        # Add the constraints.  This might `raise` if the constraints involve non-constant variables
        for lhs, rhs in constraints.items():
            system.add_constraint(lhs,rhs)

        # Add the initial conditions
        # todo: code here


        # cuz i'd prefer to NOT be insane.
        system._variable_sanity_check()

        return system

    def _expressions(self):
        """
        A helper function that produces a list of expressions combining all of the expression types

        Args:

        Returns:
            List of sympy expressions
        """

        # first, the derivatives
        exprs = [self._indep_var * expr / var for var, expr in self.derivative_dict.items() if expr != 1]

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
        
        unlisted_vars = expr_vars - set(self.variables)
        

        if unlisted_vars:
            raise RuntimeError(f"there are variables in the expressions that are not in the list of variables: {unlisted_vars}")

        for v in self.variables:
            found_vars = expressions_to_variables([v]) # convert.  constants don't make it through.

            if not found_vars:
                raise RuntimeError(f"variable `{v}` doesn't appear to be a variable")

            if list(found_vars)[0] != sympy.sympify(v):
                raise RuntimeError(f'variable `{v}` appears to not actually be a single variable.  have {found_vars}')


    def power_matrix(self):
        '''
        Determine the 'exponent' or 'power' matrix of the system, denoted by :math:`K` in the literature,
        by gluing together the power matrices of each derivative, initial condition, constraint, and require-same-variable.

        In particular, it concatenates :math:`K_{\\left(\\frac{t}{x} \\cdot \\frac{dx}{dt}\\right)}` for :math:`x` in :attr:`~variables`,
        where :math:`t` is the independent variable.

        >>> eqns = '\\n'.join(['ds/dt = -k_1*e_0*s + (k_1*s + k_m1)*c',
        ... 'dc/dt = k_1*e_0*s - (k_1*s + k_m1 + k_2)*c'])
        >>> system = ODESystem.from_equations(eqns)
        >>> system.variables
        (t, c, s, e_0, k_1, k_2, k_m1)
        >>> system.power_matrix()
        Matrix([
        [1, 1, 1,  1, 1, 1,  1],
        [0, 0, 0, -1, 0, 1,  1],
        [1, 0, 0,  1, 0, 0, -1],
        [0, 0, 0,  1, 1, 0,  0],
        [1, 0, 0,  1, 1, 1,  0],
        [0, 1, 0,  0, 0, 0,  0],
        [0, 0, 1,  0, 0, 0,  1]])

        While we get a different answer to the example in the paper, this is just due to choosing our reference exponent in a different way.

        Todo:
            * Change the code to agree with the paper.

        >>> system.update_initial_conditions({'s': 's_0'})
        >>> system.power_matrix()
        Matrix([
        [1, 1, 1,  1, 1, 1,  1,  0],
        [0, 0, 0, -1, 0, 1,  1,  0],
        [1, 0, 0,  1, 0, 0, -1,  1],
        [0, 0, 0,  1, 1, 0,  0,  0],
        [1, 0, 0,  1, 1, 1,  0,  0],
        [0, 1, 0,  0, 0, 0,  0,  0],
        [0, 0, 1,  0, 0, 0,  1,  0],
        [0, 0, 0,  0, 0, 0,  0, -1]])
        '''

        exprs = self._expressions()

        matrices = [rational_expr_to_power_matrix(expr, self.variables) for expr in exprs]
        out = sympy.Matrix.hstack(*matrices)
        assert out.shape[0] == len(self.variables)
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
        [1, 0, 0, 0, -1, -1, -1],
        [0, 1, 1, 1, -1,  0,  0]])
        '''
        exprs = self._expressions()

        return maximal_scaling_matrix(exprs, variables=self.variables)

    def reorder_variables(self, variables):
        '''
        Store the new variable order, and reorder the equations according to the new order of variables.

        Args:
            variables (str, iter):
                Another ordering of the variables.

        >>> eqns = ['dz_1/dt = z_1*z_3', 'dz_2/dt = z_1*z_2 / (z_3 ** 2)']
        >>> system = ODESystem.from_equations('\\n'.join(eqns))
        >>> system.variables
        (t, z_1, z_2, z_3)
        >>> system.derivatives
        [1, z_1*z_3, z_1*z_2/z_3**2, 0]

        >>> system.reorder_variables(['z_2', 'z_3', 't', 'z_1'])
        >>> system.variables
        (z_2, z_3, t, z_1)
        >>> system.derivatives
        [z_1*z_2/z_3**2, 0, 1, z_1*z_3]
        '''
        if isinstance(variables, str):
            if ' ' in variables:
                variables = variables.split(' ')
            else:
                variables = tuple(variables)

        # test that the set of variables is the same in the new order.
        # silvaina suspects the code would be clearer if we just set-ified them.
        if not sorted(list(map(str, variables))) == sorted(list(map(str, self.variables))):
            raise ValueError('Mismatching variables:\n{} vs\n{}'.format(sorted(list(map(str, self.variables))), sorted(list(map(str, variables)))))
        
        column_shuffle = []
        for new_var in variables:
            for i, var in enumerate(self.variables):
                if str(var) == str(new_var):
                    column_shuffle.append(i)

        self._variables = tuple( [self._variables[i] for i in column_shuffle])
        self._derivatives = tuple( [self._derivatives[i] for i in column_shuffle])

    def default_order_variables(self):
        '''
        Reorder the variables into (independent variable, dependent variables, constant variables),
        which generally gives the simplest reductions.
        Variables of the same type are sorted by their string representations.


        >>> eqns = ['dz_1/dt = z_1*z_3', 'dz_2/dt = z_1*z_2 / (z_3 ** 2)']
        >>> system = ODESystem.from_equations('\\n'.join(eqns))
        >>> system.variables
        (t, z_1, z_2, z_3)

        >>> system.reorder_variables(['z_2', 'z_3', 't', 'z_1'])
        >>> system.variables
        (z_2, z_3, t, z_1)

        >>> system.default_order_variables()
        >>> system.variables
        (t, z_1, z_2, z_3)
        '''
        all_var = self.variables
        dep_var = sorted(self.derivative_dict.keys(), key=str)
        dep_var.remove(self.indep_var)
        const_var = sorted(set(all_var).difference(dep_var).difference(set([self.indep_var])), key=str)

        # Order variables as independent, dependent, parameters
        variables = [self.indep_var] + dep_var + const_var
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
    # Feed in _clash1 so that we can use variables S, C, etc., which are special characters in sympy.
    return sympy.var(match.group(1)), sympy.sympify(match.group(3), _clash1)

def rational_expr_to_power_matrix(expr, variables):
    '''
    Take a rational expression and determine the power matrix wrt an ordering on the variables, as on page 497 of
    Hubert-Labahn.

    >>> exprs = list(map(sympy.sympify, "n*( r*(1 - n/K) - k*p/(n+d) );s*p*(1 - h*p / n)".split(';')))
    >>> variables = sorted(expressions_to_variables(exprs), key=str)
    >>> variables
    [K, d, h, k, n, p, r, s]
    >>> rational_expr_to_power_matrix(exprs[0], variables)
    Matrix([
    [0, -1, -1, 0, 0,  0],
    [0,  1,  0, 1, 0,  1],
    [0,  0,  0, 0, 0,  0],
    [1,  0,  0, 0, 0,  0],
    [0,  1,  2, 0, 1, -1],
    [1,  0,  0, 0, 0,  0],
    [0,  1,  1, 1, 1,  0],
    [0,  0,  0, 0, 0,  0]])

    >>> rational_expr_to_power_matrix(exprs[1], variables)
    Matrix([
    [ 0, 0],
    [ 0, 0],
    [ 1, 0],
    [ 0, 0],
    [-1, 0],
    [ 2, 1],
    [ 0, 0],
    [ 1, 1]])
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
    Matrix([[1, -3, -1]])

    >>> exprs = ['(z_1 + z_2**2) / z_3']
    >>> exprs = list(map(sympy.sympify, exprs))
    >>> maximal_scaling_matrix(exprs)
    Matrix([[2, 1, 2]])
    '''
    if variables is None:
        variables = sorted(expressions_to_variables(exprs), key=str)
    matrices = [rational_expr_to_power_matrix(expr, variables) for expr in exprs]
    power_matrix = sympy.Matrix.hstack(*matrices)
    assert power_matrix.shape[0] == len(variables)

    hermite_rform, multiplier_rform = hnf_row(power_matrix)

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


if __name__ == '__main__':
    import doctest
    doctest.testmod()