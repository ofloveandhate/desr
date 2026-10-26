"""
Created on Fri Dec 26 12:35:16 2014

Helper functions to deal with sympy expressions and equations

Author: Richard Tanburn (richard.tanburn@gmail.com)
"""

import fractions
import re
import sympy
from sympy.abc import _clash1, _clash2

#: Names sympy would otherwise read as something other than a symbol.  `_clash1` holds the
#: single letters -- `E I N O Q S` -- and `_clash2` the multi-letter ones -- `beta gamma
#: pi zeta`.  A modeller writing `dS/dt = -beta*S*I` means five symbols, not Euler's
#: number, the imaginary unit and two sympy functions, so desr reads every name as a symbol
#: unless told otherwise.
SYMBOL_CLASHES = dict(_clash1, **_clash2)


def sympify_expression(expr, numeric=()):
    '''
    Sympify a user-supplied expression, reading every name as a symbol.

    This is the one place desr turns text into sympy, so that a name means the same thing
    whichever way it entered: through :meth:`~desr.ode_system.ODESystem.from_tex`, through
    an initial condition, through a constraint, or through
    :mod:`desr.numerics`.

    Args:
        expr: A string, or anything else sympify accepts.  Non-strings are passed straight
            through, so an expression that is already sympy is untouched.
        numeric (iter of str): Names to leave with their usual sympy meaning instead of
            reading them as symbols.  `numeric=['pi']` is the way to write 3.14159 as
            `pi`; without it, `pi` is a parameter like any other.

    Returns:
        sympy.Expr

    >>> sympify_expression('beta*S*I')
    I*S*beta
    >>> sorted(sympify_expression('beta*S*I').free_symbols, key=str)
    [I, S, beta]

    Left to itself sympy reads `beta` as a function and the expression does not parse at
    all.  `pi` is worse, because it parses and quietly loses the symbol:

    >>> sympy.sympify('pi*x').free_symbols
    {x}

    So `pi` is a symbol here too, unless it is asked for as a number:

    >>> sympify_expression('pi*x')
    pi*x
    >>> sympify_expression('pi*x').free_symbols == {sympy.Symbol('pi'), sympy.Symbol('x')}
    True
    >>> sympify_expression('pi*x', numeric=['pi'])
    pi*x
    >>> sympify_expression('pi*x', numeric=['pi']).free_symbols
    {x}
    '''
    if not isinstance(expr, str):
        return sympy.sympify(expr)
    locals_ = {name: symbol for name, symbol in SYMBOL_CLASHES.items()
               if name not in set(numeric)}
    return sympy.sympify(expr, locals=locals_)


def to_symbol(name, numeric=()):
    '''
    Turn a single user-supplied name into a symbol, by the same convention as
    :func:`sympify_expression`.

    Args:
        name: A string, or a sympy object, which is returned unchanged.
        numeric (iter of str): As :func:`sympify_expression`.

    Returns:
        sympy.Symbol

    >>> to_symbol('S')
    S
    >>> to_symbol('S') == sympy.Symbol('S')
    True
    >>> to_symbol('beta') == sympy.Symbol('beta')
    True
    '''
    return sympify_expression(name, numeric=numeric)


def is_monomial(expr):
   ''' Determine whether expr is a monomial

       >>> is_monomial(sympy.sympify('a*b**2/c'))
       True
       >>> is_monomial(sympy.sympify('a*b**2/c + d/e'))
       False
       >>> is_monomial(sympy.sympify('a*b**2/c + 1'))
       False
       >>> is_monomial(sympy.sympify('a*(b**2/c + 1)'))
       False

   '''
   _const, _expr = expr.expand().as_coeff_add()
   if (_const != 0 and len(_expr)) or (len(_expr) > 1):
       return False
   return True

def monomial_to_powers(monomial, variables):
   ''' Given a monomial, return powers wrt some variables

       >>> variables = sympy.var('a b c d e')
       >>> monomial_to_powers(sympy.sympify('a*b'), variables)
       [1, 1, 0, 0, 0]

       >>> monomial_to_powers(sympy.sympify('a*b**2/c'), variables)
       [1, 2, -1, 0, 0]

       >>> monomial_to_powers(sympy.sympify('a*b**2/c + d/e'), variables)
       Traceback (most recent call last):
           ...
       ValueError: a*b**2/c + d/e is not a monomial

       >>> monomial_to_powers(sympy.sympify('a*b**2/c + 1'), variables)
       Traceback (most recent call last):
           ...
       ValueError: a*b**2/c + 1 is not a monomial
   '''
   # Check we have a monomial
   if not is_monomial(monomial):
       raise ValueError('{} is not a monomial'.format(monomial))

   powers = []
   power_dict = monomial.as_powers_dict()
   for var in variables:
       powers.append(power_dict.get(var, 0))
   return powers


def unique_array_stable(array):
    ''' Given a list of things, return a new list with unique elements with
        original order preserved (by first occurence)

        >>> print(unique_array_stable([1, 3, 5, 4, 7, 4, 2, 1, 9]))
        [1, 3, 5, 4, 7, 2, 9]
    '''
    seen = set()
    seen_add = seen.add
    return [x for x in array if not (x in seen or seen_add(x))]

## Helper functions

def degree(expr):
    ''' Return the degree of a sympy expression. I.e. the largest number of
        variables multiplied together.
        NOTE DOES take into account idempotency of binary variables

        >>> str_eqns = ['x + y',
        ...             'x*y*z - 1',
        ...             'x ** 2 + a*b*c',
        ...             'x**2 + y',
        ...             'x',
        ...             'x*y',]
        >>> eqns = str_exprs_to_sympy_eqns(str_eqns)
        >>> for e in eqns: print(degree(e.lhs - e.rhs))
        1
        3
        3
        1
        1
        2

        Check we deal with constants correctly
        >>> (degree(0), degree(1), degree(4),
        ... degree(sympy.S.Zero), degree(sympy.S.One), degree(sympy.sympify(4)))
        (0, 0, 0, 0, 0, 0)
    '''
    if is_constant(expr):
        return 0
    degree = 0
    for term in expr.as_coefficients_dict().keys():
        degree = max(degree, len(term.free_symbols))
    return degree

def is_constant(expr):
    ''' Determine whether an expression is constant
        >>> expr = 'x + 2*y'
        >>> is_constant(sympy.sympify(expr))
        False
        >>> expr = 'x + 5'
        >>> is_constant(sympy.sympify(expr))
        False
        >>> expr = '3'
        >>> is_constant(sympy.sympify(expr))
        True
        >>> expr = '2*x - 4'
        >>> is_constant(sympy.sympify(expr))
        False
    '''
    if isinstance(expr, (int, float)):
        return True
    return not expr.free_symbols

def is_equation(eqn, check_true=True):
    ''' Return True if it is an equation rather than a boolean value.
        If it is False, raise a ContradictionException. We never want anything
        that might be False.

        Optionally, we can turn the check off, but THE DEFAULT VALUE SHOULD
        ALWAYS BE TRUE. Otherwise bad things will happen.

        >>> x, y = sympy.symbols('x y')
        >>> eq1 = sympy.Eq(x, y)
        >>> eq2 = sympy.Eq(x, x)
        >>> eq3 = sympy.Eq(x, y).subs(y, x)
        >>> eq4 = sympy.Eq(2*x*y, 2)

        >>> is_equation(eq1)
        True
        >>> is_equation(eq2)
        False
        >>> is_equation(eq3)
        False
        >>> is_equation(eq4)
        True

        Now check that it raises exceptions for the right things

        >>> is_equation(0)
        False
    '''
    if sympy.__version__ == '0.7.5':
        return isinstance(eqn, sympy.Equality)
    elif re.match(r'1\..*', sympy.__version__):
        return isinstance(eqn, sympy.Equality)
    else:
        return eqn is True

def standardise_equation(eqn):
    ''' Remove binary squares etc '''
    if not is_equation(eqn):
        return eqn
    eqn = remove_binary_squares_eqn(eqn.expand())
    eqn = balance_terms(eqn)
    eqn = cancel_constant_factor(eqn)
    return eqn

def expressions_to_variables(exprs):
    ''' Take a list of equations or expressions and return a set of variables

        >>> eqn = sympy.Eq(sympy.sympify('x*a + 1'),0)
        >>> expr = sympy.sympify('x + y*z + 2*a^b')
        >>> to_test = [expr, eqn]
        >>> vars = expressions_to_variables(to_test)
        >>> print(sorted(list(map(str,vars))))
        ['a', 'b', 'x', 'y', 'z']
    '''
    if len(exprs) == 0:
        return set()
    if sympy.__version__ == '0.7.5':
        assert all(map(lambda x: isinstance(x, sympy.Basic), exprs))
    return set.union(*[expr.atoms(sympy.Symbol) for expr in exprs if expr is not None])

def eqns_with_variables(eqns, variables, strict=False):
    ''' Given a set of atoms, return only equations that have something in
        common

        >>> x, y, z1, z2 = sympy.symbols('x y z1 z2')
        >>> eqns = ['x + y == 1', '2*z1 + 1 == z2', 'x*z1 == 0']
        >>> eqns = str_eqns_to_sympy_eqns(eqns)
        >>> eqns_with_variables(eqns, [x])
        [Eq(x + y - 1, 0), Eq(x*z1, 0)]
        >>> eqns_with_variables(eqns, [z1])
        [Eq(2*z1 - z2 + 1, 0), Eq(x*z1, 0)]
        >>> eqns_with_variables(eqns, [y])
        [Eq(x + y - 1, 0)]

        >>> eqns_with_variables(eqns, [x], strict=True)
        []
        >>> eqns_with_variables(eqns, [x, z1], strict=True)
        [Eq(x*z1, 0)]
        >>> eqns_with_variables(eqns, [x, y, z1], strict=True)
        [Eq(x + y - 1, 0), Eq(x*z1, 0)]
    '''
    if strict:
        return [eqn for eqn in eqns if eqn.atoms(sympy.Symbol).issubset(variables)]
    else:
        return [eqn for eqn in eqns if len(eqn.atoms(sympy.Symbol).intersection(variables))]

def dict_as_eqns(dict_):
    ''' Given a dictionary of lhs: rhs, return the sympy Equations in a list

        >>> x, y, z = sympy.symbols('x y z')
        >>> dict_as_eqns({x: 1, y: z, x*y: 1 - z})
        [Eq(x, 1), Eq(y, z), Eq(x*y, 1 - z)]
    '''
    return [sympy.Eq(lhs, rhs) for lhs, rhs in dict_.items()]

def str_eqns_to_sympy_eqns(str_eqns):
    ''' Take string equations and sympify

        >>> str_eqns = ['x + y == 1', 'x*y*z - 3*a == -3']
        >>> eqns = str_eqns_to_sympy_eqns(str_eqns)
        >>> for e in eqns: print(e)
        Eq(x + y - 1, 0)
        Eq(-3*a + x*y*z + 3, 0)
    '''
    str_exprs = []
    for str_eqn in str_eqns:
        str_exprs.append('{} - ({})'.format(*str_eqn.split('==')))
    return str_exprs_to_sympy_eqns(str_exprs)

def str_exprs_to_sympy_eqns(str_exprs):
    ''' Take some strings and return the sympy expressions

        >>> str_eqns = ['x + y - 1', 'x*y*z - 3*a + 3', '2*a - 4*b']
        >>> eqns = str_exprs_to_sympy_eqns(str_eqns)
        >>> for e in eqns: print(e)
        Eq(x + y - 1, 0)
        Eq(-3*a + x*y*z + 3, 0)
        Eq(2*a - 4*b, 0)
    '''
    exprs = map(sympy.sympify, str_exprs)
    exprs = list(map(lambda e: sympy.Eq(e,0), exprs))
    return exprs


if __name__ == "__main__":
    import doctest
    doctest.testmod()
