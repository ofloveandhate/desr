'''
Statements a reduction must satisfy, checked without reference to how it was computed.

Shared by the test modules that need them.  Not a test module itself.
'''
import sympy

from desr.ode_system import ODESystem


def false_equations(system, reduced):
    '''
    The reduced equations that are not true, by the chain rule.

    Each reduced symbol denotes a monomial in the original symbols, and the reduced
    system gives its rate with respect to the reduced clock.  The original system gives
    the rate of the same monomial directly.  The two must agree.

    Returns:
        list of sympy.Symbol: the reduced symbols whose equations fail; empty when the
        reduction is true.
    '''
    denotes = reduced.named_symbols
    rates = {symbol: (sympy.Integer(0) if rate is None else rate)
             for symbol, rate in zip(system.symbols, system.derivatives)}
    rates[system.independent_variable] = sympy.Integer(1)

    def rate_of(monomial):
        return sum(sympy.diff(monomial, symbol) * rate for symbol, rate in rates.items())

    clock = reduced.independent_variable
    clock_rate = rate_of(denotes[clock])
    wrong = []
    for symbol, rate in zip(reduced.symbols, reduced.derivatives):
        if symbol == clock:
            continue
        claimed = (sympy.Integer(0) if rate is None
                   else rate.subs(denotes, simultaneous=True) * clock_rate)
        if sympy.simplify(rate_of(denotes[symbol]) - claimed) != 0:
            wrong.append(symbol)
    return wrong


def forget(reduced):
    ''' The same equations as a plain system, with no record of being a reduction. '''
    return ODESystem._from_variables(reduced.symbols, reduced.derivatives,
                     indep_var=reduced.independent_variable)
