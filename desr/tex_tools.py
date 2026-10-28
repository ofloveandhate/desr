r"""
Reading and writing the LaTeX that papers are written in.

Created on Wed Aug 12 01:37:16 2015

Author: Richard Tanburn (richard.tanburn@gmail.com)

The reading half is a **convenience**, not a LaTeX implementation: it exists because
pasting a system out of a paper's source is easier than typing it out in sympy by hand.
So it covers the notation that papers actually use -- text-mode species names,
concentration brackets, ``\left(``, the Leibniz ``\frac{d}{dt}x``, primes, multi-character
subscripts -- and raises :class:`TexParseError` on anything it cannot read or could read
two ways.  It will not cover every paper, and a raise is the right answer when it does
not: the alternative is a system that is quietly not the one that was written down.

General LaTeX math parsers were considered and do not fit.  :func:`sympy.parsing.latex.parse_latex`
fails outright on ``\textrm``, square brackets, ``&`` alignment and single primes; keeps
the braces inside symbol names, so ``k_a`` arrives as ``k_{a}`` and collides with the
``k_{-1}`` / ``k_m1`` convention here; returns an ambiguity tree rather than an error for
``k_3 (x - y)``; and needs antlr4 or lark as a hard dependency.  None of them classify the
lines of a system, which is the part that matters most.
"""

import re
import sympy

from desr.sympy_helper import sympify_expression, to_symbol

VAR_RE = r'[A-Za-z*][\d_]*'

def matrix_to_tex(matrix_):
    '''
    Given a matrix, write out the TeX.

    Args:
        matrix_ (sympy.Matrix): Matrix to turn into TeX
    Returns:
        str

    Printing is the correct way to use this function, but the docstring looks a bit odd.
    >>> print(matrix_to_tex(sympy.eye(2)))
    1 & 0 \\\\
    0 & 1 \\\\
    '''
    lines = []
    for line in matrix_.tolist():
        lines.append(' & '.join(map(str, line)) + ' \\\\')
    return '\n'.join(lines)

def _var_repler(var):
    var = var.group()
    if len(var) == 1:
        return var[0]
    var_letter, subscript = var[0], var[1:]
    if subscript[0] == '_':
        subscript = subscript[1:]
    subscript = subscript.replace('_', '')
    return '{}_{{{}}}'.format(var_letter, subscript)


def var_to_tex(var):
    """
    Given a sympy variable, write out the TeX.

    Args:
        var (sympy.Symbol): Variable to turn into TeX
    Returns:
        str

    >>> print(list(map(var_to_tex, sympy.symbols('x y_1 Kw_3 z_{3} k_m1'))))
    ['x', 'y_{1}', 'Kw_{3}', 'z_{3}', 'k_{-1}']
    """
    return expr_to_tex(var)

def expr_to_tex(expr):
    """
    Given a sympy expression, write out the TeX.

    Args:
        expr (sympy.Expression): Expression to turn into TeX
    Returns:
        str

    >>> print(list(map(expr_to_tex, map(lambda x: sympy.sympify(x,rational=True), ['(x + y - 1.5)**2', '(x + y_m1)**1', 'k_m1*t']))))
    ['\\\\left(x + y - \\\\frac{3}{2}\\\\right)^{2}', 'x + y_{-1}', 'k_{-1} t']
    """
    tex = sympy.latex(expr)
    # Substitute _{m...} for _{-...}
    tex = re.sub(r'\_\{?m([^\}]+)\}?', r'_{-\1}', tex)
    return tex

def eqn_to_tex(eqn):
    eqn = str(eqn).replace(' ', '')

    expr1, expr2 = eqn.split('==')

    tex = '{} &= {}'.format(expr_to_tex(expr1), expr_to_tex(expr2))
    return tex


def eqns_to_tex(eqns):
    ''' To convert to array environment, copy the output into a lyx LaTeX cell,
        then copy this entire cell into an eqnarray of sufficient size
    '''
    return '\\\\'.join(map(eqn_to_tex, eqns))

class TexParseError(ValueError):
    '''
    The tex could not be read, or could be read in more than one way.

    A subclass of :class:`ValueError`, so ``except ValueError`` still catches it.
    '''


#: The commands whose contents are a literal name rather than mathematics.  Inside one of
#: these a space is part of the name, which is the only reason ``[total Cdc25]`` can be
#: read at all -- everywhere else in the tex, a space between names means multiplication.
TEXT_MODE_COMMANDS = ('textrm', 'text', 'mathrm')

_NAME = re.compile(r'[A-Za-z_][A-Za-z0-9_]*\Z')

_TEXT_MODE = re.compile(r'\\(?:' + '|'.join(TEXT_MODE_COMMANDS) + r')\s*\{([^{}]*)\}')

_LEIBNIZ = re.compile(r'\\frac\s*\{\s*d\s*\}\s*\{\s*d\s*([A-Za-z][A-Za-z0-9_]*)\s*\}'
                      r'\s*(?:\{\s*([A-Za-z][A-Za-z0-9_]*)\s*\}|([A-Za-z][A-Za-z0-9_]*))')

_BRACKETED = re.compile(r'\[([^\[\]]*)\]')

_SUBSCRIPT = re.compile(r'_\{([^{}]*)\}')

#: Anything that makes bracket contents mathematics rather than a name.
_NOT_A_NAME = re.compile(r'[-+*/^()\\]')

#: The independent variable assumed when a prime is read as a derivative.  ``x'`` does not
#: say what it is differentiated with respect to, so there is nothing to infer.
PRIME_INDEP_VAR = 't'


def _as_name(text, source):
    '''Close up the spaces in a text-mode or bracketed group and insist it is a name.'''
    name = re.sub(r'\s+', '', text.strip())
    if not name or not _NAME.match(name):
        raise TexParseError(
            '{!r} does not name a variable, so {!r} cannot be read.  A name starts with a '
            'letter or underscore and continues with letters, digits or '
            'underscores.'.format(name, source.strip()))
    return name


def _resolve_text_mode(tex):
    r'''
    Replace every text-mode group by the name it holds.

    Concentration brackets inside the group are part of the notation, not the name, so
    ``\textrm{[total Cdc25]}`` is the variable ``totalCdc25``.  Runs to a fixed point, so a
    text-mode group inside a subscript -- ``k_{\textrm{wee}}`` -- is resolved too.
    '''
    def replace(match):
        content = match.group(1).strip()
        if content.startswith('[') and content.endswith(']'):
            content = content[1:-1]
        return _as_name(content, match.group(0))

    previous = None
    while previous != tex:
        previous = tex
        tex = _TEXT_MODE.sub(replace, tex)
    return tex


def _resolve_brackets(tex):
    '''
    Decide what each remaining ``[...]`` means.

    Square brackets are usually a grouping paren and occasionally a concentration.  The
    contents say which: anything with an operator in it is mathematics, a single token is a
    name (and is the same expression read either way).  Several bare names separated by
    spaces could be one species or a product, and that is refused rather than guessed --
    :data:`TEXT_MODE_COMMANDS` is how to say which was meant.
    '''
    def replace(match):
        content = match.group(1).strip()
        if _NOT_A_NAME.search(content):
            return '({})'.format(content)
        if len(content.split()) > 1:
            raise TexParseError(
                '{!r} is ambiguous: it could be the single name {!r} or the product of '
                '{}.  Write it as \\textrm{{{}}} if a name was meant.'.format(
                    match.group(0), re.sub(r'\s+', '', content),
                    ' and '.join(repr(part) for part in content.split()), content))
        return _as_name(content, match.group(0))

    previous = None
    while previous != tex:
        previous = tex
        tex = _BRACKETED.sub(replace, tex)
    return tex


def _resolve_subscripts(tex):
    '''
    Flatten ``_{...}`` to ``_...``, mapping a minus to an ``m`` on the way.

    The minus has to be dealt with here rather than afterwards: ``k_{-1}`` is the symbol
    ``k_m1``, and flattening first would give ``k_-1``, which is not a name at all.
    '''
    def replace(match):
        return '_' + re.sub(r'\s+', '', match.group(1)).replace('-', 'm')

    previous = None
    while previous != tex:
        previous = tex
        tex = _SUBSCRIPT.sub(replace, tex)
    return tex


def _normalise_tex(tex):
    r'''
    Turn the tex a paper is written in into the tex the parser can read.

    The order matters, and each step depends on the ones before it:

    1. text-mode groups become names, so that later steps see ``Cyclin`` and not
       ``\textrm{[Cyclin]}``;
    2. ``\left`` and ``\right`` go, being pure noise -- and they must go before the
       brackets are judged, so that ``\left[ x + y \right]`` is seen as ``[x + y]``;
    3. ``\frac{d}{dt}X`` becomes ``\frac{dX}{dt}``, which needs step 1 to have turned a
       text-mode variable into a plain name;
    4. square brackets are resolved;
    5. subscripts are flattened.

    Everything after this point may safely discard backslashes, which is what the rest of
    the pipeline does.
    '''
    tex = _resolve_text_mode(tex)
    tex = re.sub(r'\\(?:left|right)\s*', '', tex)
    tex = _LEIBNIZ.sub(lambda m: r'\frac{{d{}}}{{d{}}}'.format(m.group(2) or m.group(3),
                                                               m.group(1)), tex)
    tex = _resolve_brackets(tex)
    tex = _resolve_subscripts(tex)
    return tex


def tex_to_sympy(tex, primes='name'):
    """
    Given some possibly multi-line TeX, turn it into sympy expressions and equations.
    Each line is parsed seperately.

    Lines may be separated by newlines, by the TeX line break ``\\\\``, or both, since a
    paper's tex is one long string broken by ``\\\\``.

    Args:
        tex (str): LaTeX
        primes (str): ``'name'`` (the default) reads a prime as part of the variable's
            name; ``'derivative'`` reads one on the left of an equation as a derivative.
            See :func:`tex_to_sympy_one_line`.

    Returns:
        list

    >>> lines = [r'\\frac{dE}{dt} &= - k_1 E S + k_{-1} C + k_2 C \\\\',
    ... r'\\frac{dS}{dt} &= - k_1 E S + k_{-1} C \\\\',
    ... r'\\frac{dC}{dt} &= k_1 E S - k_{-1} C - k_2 C \\\\',
    ... r'\\frac{dP}{dt} &= k_2 C']
    >>> sym = tex_to_sympy('\\n'.join(lines))
    >>> for s in sym: print(s)
    Eq(Derivative(E, t), C*k_2 + C*k_m1 - E*S*k_1)
    Eq(Derivative(S, t), C*k_m1 - E*S*k_1)
    Eq(Derivative(C, t), -C*k_2 - C*k_m1 + E*S*k_1)
    Eq(Derivative(P, t), C*k_2)


    >>> print(tex_to_sympy('k_2 &= V_2d ( APCT - APCs ) + V_2dd APCs'))
    [Eq(k_2, APCs*V_2dd + V_2d*(APCT - APCs))]
    """
    # Parse each line individually.  A paper's tex is one long string broken by `\\`, which
    # may or may not coincide with newlines.
    lines = re.split(r'\\\\|\n', tex)
    return [tex_to_sympy_one_line(line, primes=primes)
            for line in lines if line.strip()]


def tex_to_sympy_one_line(tex, primes='name'):
    r"""
    Process just one line of tex.

    Args:
        tex (str): LaTeX, just ONE line of it
        primes (str): How to read a prime.

            ``'name'``, the default, makes it part of the variable's name, so ``V_{25}'``
            is the symbol ``V_25prime``.  Papers use primes this way for two forms of the
            same quantity, and it is the reading that cannot be wrong by accident.

            ``'derivative'`` reads a prime on the **left** of an equation as a derivative
            with respect to :data:`PRIME_INDEP_VAR`.  A prime anywhere else raises, because
            there is no way to represent it; so does a double prime, which is a second
            derivative and outside what this library describes.

    Returns:
        sympyfication of that line of code

    Raises:
        TexParseError: If the line cannot be read, or could be read two ways.

    >>> tex_to_sympy_one_line(r'\frac{d}{dt}\textrm{[total Cdc25]}')
    Derivative(totalCdc25, t)
    >>> tex_to_sympy_one_line(r"V_{25}'")
    V_25prime
    """
    if primes not in ('name', 'derivative'):
        raise ValueError("primes must be 'name' or 'derivative', not {!r}".format(primes))

    return _parse_normalised(_normalise_tex(tex), primes)


def _parse_normalised(tex, primes):
    '''Parse tex that :func:`_normalise_tex` has already been over.'''
    # Remove alignment characters
    tex = tex.strip().replace('&', '')

    sides = tex.split('=')
    if len(sides) == 2:
        lhs, rhs = sides
        if primes == 'derivative':
            return _parse_prime_equation(lhs, rhs)
        return sympy.Eq(_parse_side(lhs, primes), _parse_side(rhs, primes))
    if len(sides) != 1:
        raise TexParseError('Too many = in {}.'.format(tex))

    return _parse_side(tex, primes)


def _parse_prime_equation(lhs, rhs):
    '''An equation under ``primes='derivative'``, where a prime may only be on the left.'''
    if "'" in rhs:
        raise TexParseError(
            "primes='derivative' found a prime on the right of {!r}, where there is no "
            "way to represent it.  Use primes='name' if the prime is part of a variable's "
            "name.".format('='.join((lhs, rhs)).strip()))

    match = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*('+)\s*\Z", lhs)
    if match is None:
        if "'" in lhs:
            raise TexParseError(
                "primes='derivative' can only read a prime on a bare variable, not in "
                "{!r}.".format(lhs.strip()))
        return sympy.Eq(_parse_side(lhs, 'derivative'), _parse_side(rhs, 'derivative'))

    name, marks = match.groups()
    if len(marks) > 1:
        raise TexParseError(
            "{}{} is a second derivative, and desr describes first-order systems only.  "
            "If the primes are part of the variable's name, parse with "
            "primes='name'.".format(name, marks))

    return sympy.Eq(sympy.Derivative(to_symbol(name), to_symbol(PRIME_INDEP_VAR)),
                    _parse_side(rhs, 'derivative'))


def _parse_side(tex, primes):
    '''One side of an equation, or a whole line that is not one.'''
    if primes == 'name':
        tex = tex.replace("'", 'prime')
    elif "'" in tex:
        raise TexParseError(
            "primes='derivative' found a prime in {!r}, where it can only mean a "
            "derivative on the left of an equation.".format(tex.strip()))

    # Past normalisation, the only backslashes left belong to \frac, whose name the
    # patterns below match with or without it.
    tex = tex.replace('\\', '')

    # Turn \frac{d }{d } into sympy.Derivatives
    diff_match = re.match(r'\s*[\\f\ff]?rac\{d(.+)\}\{d(.+)\}', tex)
    if diff_match:
        return sympy.Derivative(*[to_symbol(group.strip()) for group in diff_match.groups()])

    # Turn \frac into ratios. Consume the shortest amount possible
    tex = re.sub(r'[\\f\ff]?rac{(.*?)}{(.*?)}', '((\\1) / (\\2))', tex)

    # Turn spaces between variables into *. Do this by matching anything that isn't an operation
    # Use a lookahead assertion to get overlapping instances.
    tex = re.sub(r'([^+\-*\s/(]+)\s+(?=[^+\-*\s/)]+)', '\\1 * ', tex)

    # A name butted against a parenthesis is multiplication: desr has no function-valued
    # variables, so `k_3(x - y)` can only be a product.
    tex = re.sub(r'([A-Za-z0-9_]+)\s*\(', '\\1*(', tex)

    # Read every name as a symbol: a species called S, E, I, beta or pi is an ordinary
    # variable, not one of sympy's own objects.
    return sympify_expression(tex)

if __name__ == '__main__':
    import doctest
    doctest.testmod()