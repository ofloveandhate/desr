'''
Reading the LaTeX that papers are actually written in.

The constructs here were all taken from one place: the Tyson cell-cycle system as
transcribed in `examples/cell_cycle_control_network.py`, which needed sixty lines of
hand-rolled regex preprocessing before `from_tex` would look at it.  That preprocessing is
the specification, and this is it in test form.

Two rules decided in conversation, both chosen so that a wrong reading is loud rather than
silent:

  - A prime is a **name** by default (`V_{25}'` is one of two forms of an enzyme, not a
    derivative).  `primes='derivative'` opts into the other reading, where a prime is
    meaningful only on the left-hand side, and a double prime -- a second derivative --
    raises, because this library is first-order only.
  - A non-derivative line relating only parameters is a constraint.  One involving a
    dependent variable is a definition, and raises unless `sub_defns=True` asks for it to
    be substituted away.
'''
from unittest import TestCase

import sympy

from desr.ode_system import ODESystem
from desr.tex_tools import TexParseError, tex_to_sympy, tex_to_sympy_one_line


class TestTextMode(TestCase):
    '''
    `\\textrm{...}`, `\\text{...}` and `\\mathrm{...}` hold a literal name.  Inside them a
    space is part of the name, which is what makes `[total Cdc25]` readable at all: a space
    means multiplication everywhere else in the tex.
    '''

    def test_textrm_is_stripped(self):
        self.assertEqual(tex_to_sympy_one_line(r'\textrm{Cyclin}'), sympy.Symbol('Cyclin'))

    def test_text_and_mathrm_too(self):
        for command in ('text', 'mathrm'):
            with self.subTest(command=command):
                self.assertEqual(tex_to_sympy_one_line(r'\%s{Cyclin}' % command),
                                 sympy.Symbol('Cyclin'))

    def test_concentration_brackets_inside_text_mode(self):
        self.assertEqual(tex_to_sympy_one_line(r'\textrm{[Cyclin]}'),
                         sympy.Symbol('Cyclin'))

    def test_a_space_inside_text_mode_is_part_of_the_name(self):
        self.assertEqual(tex_to_sympy_one_line(r'\textrm{[total Cdc25]}'),
                         sympy.Symbol('totalCdc25'))

    def test_text_mode_inside_a_subscript(self):
        self.assertEqual(tex_to_sympy_one_line(r'k_{\textrm{wee}}'),
                         sympy.Symbol('k_wee'))

    def test_adjacent_species_still_multiply(self):
        self.assertEqual(tex_to_sympy_one_line(r'\textrm{[Cyclin]} \textrm{[Cdk]}'),
                         sympy.Symbol('Cyclin') * sympy.Symbol('Cdk'))


class TestBrackets(TestCase):
    '''
    A bare `[...]` is usually a grouping paren and occasionally a concentration.  The
    contents say which, and when they genuinely do not, that is an error rather than a
    guess.
    '''

    def test_brackets_around_an_expression_are_a_group(self):
        x, y = sympy.symbols('x y')
        self.assertEqual(tex_to_sympy_one_line(r'k [x + y]'),
                         sympy.Symbol('k') * (x + y))

    def test_brackets_around_one_name_are_harmless(self):
        # Read either way, `[x]` is the same expression.
        self.assertEqual(tex_to_sympy_one_line(r'[x]'), sympy.Symbol('x'))

    def test_several_bare_names_in_brackets_is_ambiguous(self):
        # `AB` the species, or `A*B` the product?  Only `\textrm` can say.
        with self.assertRaises(TexParseError) as caught:
            tex_to_sympy_one_line(r'[A B]')
        self.assertIn('textrm', str(caught.exception))

    def test_bracket_contents_that_cannot_be_a_name(self):
        with self.assertRaises(TexParseError):
            tex_to_sympy_one_line(r'\textrm{[2x]}')


class TestDelimitersAndFractions(TestCase):
    ''' The noise that papers are full of. '''

    def test_left_and_right_are_stripped(self):
        x, y = sympy.symbols('x y')
        self.assertEqual(tex_to_sympy_one_line(r'k_a \left( x - y \right)'),
                         sympy.Symbol('k_a') * (x - y))

    def test_implicit_multiplication_before_a_parenthesis(self):
        x, y = sympy.symbols('x y')
        self.assertEqual(tex_to_sympy_one_line(r'k_3 (x - y)'),
                         sympy.Symbol('k_3') * (x - y))

    def test_a_fraction_still_works(self):
        x = sympy.Symbol('x')
        self.assertEqual(tex_to_sympy_one_line(r'\frac{k_a x}{K_a + x}'),
                         sympy.Symbol('k_a') * x / (sympy.Symbol('K_a') + x))

    def test_a_minus_in_a_subscript_still_becomes_m(self):
        # Long-standing convention: k_{-1} is the symbol k_m1.  Flattening `_{...}` naively
        # would give `k_-1`, which is not a name at all.
        self.assertEqual(tex_to_sympy_one_line(r'k_{-1} C'),
                         sympy.Symbol('k_m1') * sympy.Symbol('C'))

    def test_a_multi_character_subscript_is_flattened(self):
        self.assertEqual(tex_to_sympy_one_line(r'V_{25}'), sympy.Symbol('V_25'))


class TestLinesThatCollapse(TestCase):
    '''
    An equation whose two sides are already comparable is not an equation: sympy.Eq
    evaluates it to a boolean.  Such a line states nothing about any system, and returning
    the boolean loses it -- `from_tex` then asked the boolean for its `.lhs` and died with
    `'BooleanTrue' object has no attribute 'lhs'`.

    This is what `x(0) = 0` used to do, because the implicit-multiplication rule read it
    as `x*0 = 0`.  That line is now recognised as an initial condition before that rule
    runs, so it no longer reaches here; see `TestInitialConditionLines`.  The guard stays
    for the lines that really do collapse.
    '''

    DYNAMIC = '\n'.join([r'\frac{d}{dt}x &= k_1 - k_2 x \\',
                         r'\frac{d}{dt}y &= k_2 x - k_3 y'])

    def test_a_tautology_is_refused(self):
        with self.assertRaises(TexParseError) as caught:
            tex_to_sympy_one_line(r'k_1 &= k_1')
        self.assertIn('nothing about the system', str(caught.exception))

    def test_a_falsehood_says_so(self):
        with self.assertRaises(TexParseError) as caught:
            tex_to_sympy_one_line(r'0 &= 1')
        self.assertIn('no system satisfies it', str(caught.exception))

    def test_from_tex_reports_it_instead_of_an_attribute_error(self):
        with self.assertRaises(TexParseError):
            ODESystem.from_tex(self.DYNAMIC + '\\\\\n' + r'k_1 &= k_1')

    def test_a_variable_applied_to_something_other_than_zero_still_multiplies(self):
        # The rule that a name against a parenthesis is multiplication is unchanged for
        # everything but the literal zero, so this collapses as it always did.
        with self.assertRaises(TexParseError):
            tex_to_sympy_one_line(r'x\left(0 + 0\right) &= 0')

    def test_a_real_equation_is_untouched(self):
        x = sympy.Symbol('x')
        self.assertEqual(tex_to_sympy_one_line(r'k_4 &= k_1 + x'),
                         sympy.Eq(sympy.Symbol('k_4'), sympy.Symbol('k_1') + x))

    def test_a_line_that_is_not_an_equation_is_untouched(self):
        self.assertEqual(tex_to_sympy_one_line(r'k_{-1} C'),
                         sympy.Symbol('k_m1') * sympy.Symbol('C'))


class TestDerivativeNotation(TestCase):
    ''' Both ways a paper writes a derivative. '''

    def test_the_named_form(self):
        t, x = sympy.symbols('t x')
        self.assertEqual(tex_to_sympy_one_line(r'\frac{dx}{dt}'), sympy.Derivative(x, t))

    def test_the_leibniz_form(self):
        # `\frac{d}{dt}x` is how most papers write it, and from_tex used to refuse it.
        t, x = sympy.symbols('t x')
        self.assertEqual(tex_to_sympy_one_line(r'\frac{d}{dt}x'), sympy.Derivative(x, t))

    def test_the_leibniz_form_with_a_text_mode_variable(self):
        t = sympy.Symbol('t')
        self.assertEqual(tex_to_sympy_one_line(r'\frac{d}{dt}\textrm{[Cyclin]}'),
                         sympy.Derivative(sympy.Symbol('Cyclin'), t))


class TestPrimes(TestCase):
    ''' A prime is a name unless asked otherwise, and never a second derivative. '''

    def test_a_prime_is_part_of_the_name_by_default(self):
        self.assertEqual(tex_to_sympy_one_line(r"V_{25}'"), sympy.Symbol('V_25prime'))

    def test_a_double_prime_is_part_of_the_name_by_default(self):
        self.assertEqual(tex_to_sympy_one_line(r"V_{2}''"),
                         sympy.Symbol('V_2primeprime'))

    def test_a_prime_on_the_left_can_mean_a_derivative(self):
        t, x = sympy.symbols('t x')
        parsed = tex_to_sympy_one_line(r"x' = -k x", primes='derivative')
        self.assertEqual(parsed.lhs, sympy.Derivative(x, t))

    def test_a_second_derivative_is_refused(self):
        # This library is first-order only, so a double prime read as a derivative is a
        # violation of its assumptions and is caught at parse time.
        with self.assertRaises(TexParseError) as caught:
            tex_to_sympy_one_line(r"x'' = -k x", primes='derivative')
        self.assertIn('second', str(caught.exception).lower())

    def test_a_prime_on_the_right_is_refused(self):
        # Under primes='derivative' there is no way to represent it, so say so eagerly
        # rather than read it as something it is not.
        with self.assertRaises(TexParseError):
            tex_to_sympy_one_line(r"x = y'", primes='derivative')

    def test_an_unknown_primes_option_is_refused(self):
        with self.assertRaises(ValueError):
            tex_to_sympy_one_line(r'x', primes='whatever')


class TestLineClassification(TestCase):
    '''
    A line is one of three things: an equation of the system, a constraint among
    parameters, or a definition of something in terms of dependent variables.
    '''

    DYNAMIC = '\n'.join([r'\frac{d}{dt}x &= k_1 - k_2 x \\',
                         r'\frac{d}{dt}y &= k_2 x - k_3 y'])

    def test_derivative_lines_become_the_system(self):
        system = ODESystem.from_tex(self.DYNAMIC)
        self.assertEqual(set(map(str, system.dependents)), {'x', 'y'})
        self.assertEqual(set(map(str, system.parameters)), {'k_1', 'k_2', 'k_3'})

    def test_a_relation_among_parameters_only_is_a_constraint(self):
        tex = self.DYNAMIC + '\\\\\n' + r'k_4 &= k_1 + k_2'
        system = ODESystem.from_tex(tex)
        self.assertEqual(len(system.constraints), 1)
        self.assertIn(sympy.Symbol('k_4'), system.parameters)

    def test_a_definition_involving_a_dependent_raises_by_default(self):
        tex = self.DYNAMIC + '\\\\\n' + r'z &= k_1 - x'
        with self.assertRaises(TexParseError) as caught:
            ODESystem.from_tex(tex)
        message = str(caught.exception)
        self.assertIn('z', message)
        self.assertIn('sub_defns', message)

    def test_sub_defns_substitutes_the_definition_away(self):
        tex = self.DYNAMIC + '\\\\\n' + r'z &= k_1 - x'
        system = ODESystem.from_tex(tex, sub_defns=True)
        self.assertNotIn(sympy.Symbol('z'), system.symbols)

    def test_a_definition_referring_to_another_definition_is_not_implemented(self):
        # Substitution order would decide the answer, and cycles would have to be caught.
        # Refused for now rather than answered wrongly.
        tex = '\n'.join([self.DYNAMIC + r' \\',
                         r'z &= k_1 - x \\',
                         r'w &= z + y'])
        with self.assertRaises(NotImplementedError):
            ODESystem.from_tex(tex, sub_defns=True)


class TestTheTysonCellCycleSystem(TestCase):
    '''
    The system that motivated all of this: Tyson's cell-cycle control network, transcribed
    from the 2007 paper.  Seven derivative lines and four definitions, every one of which
    involves a dependent variable.
    '''

    DERIVATIVE_LINES = r'''\frac{d}{dt}\textrm{[Cyclin]} &= k_1 - k_2 \textrm{[Cyclin]} - k_3 \textrm{[Cyclin]} \textrm{[Cdk]} \\
\frac{d}{dt}\textrm{[MPF]} &= k_3 \textrm{[Cyclin]} \textrm{[Cdk]} - k_2 \textrm{[MPF]} - k_{\textrm{wee}} \textrm{[MPF]} + k_{25} \textrm{[preMPF]} \\
\frac{d}{dt}\textrm{[preMPF]} &= -k_2 \textrm{[preMPF]} + k_{\textrm{wee}} \textrm{[MPF]} - k_{25} \textrm{[preMPF]} \\
\frac{d}{dt}\textrm{[Cdc25P]} &= \frac{k_a \textrm{[MPF]}\left( \textrm{[total Cdc25]} - \textrm{[Cdc25P]} \right)}{K_a + \textrm{[total Cdc25]} - \textrm{[Cdc25P]}} - \frac{k_b \textrm{[PPase]} \textrm{[Cdc25P]}}{K_b + \textrm{[Cdc25P]}} \\
\frac{d}{dt}\textrm{[Wee1P]} &= \frac{k_e \textrm{[MPF]}\left( \textrm{[total Wee1]} - \textrm{[Wee1P]} \right)}{K_e + \textrm{[total Wee1]} - \textrm{[Wee1P]}} - \frac{k_f \textrm{[PPase]} \textrm{[Wee1P]}}{K_f + \textrm{[Wee1P]}} \\
\frac{d}{dt}\textrm{[IEP]} &= \frac{k_g \textrm{[MPF]}\left( \textrm{[total IE]} - \textrm{[IEP]} \right)}{K_g + \textrm{[total IE]} - \textrm{[IEP]}} - \frac{k_h \textrm{[PPase]} \textrm{[IEP]}}{K_h + \textrm{[IEP]}} \\
\frac{d}{dt}\textrm{[APC]} &= \frac{k_c \textrm{[MPF]}\left( \textrm{[total APC]} - \textrm{[APC]} \right)}{K_c + \textrm{[total APC]} - \textrm{[APC]}} - \frac{k_d \textrm{[PPase]} \textrm{[APC]}}{K_d + \textrm{[APC]}}'''

    DEFINITION_LINES = r'''\textrm{[Cdk]} &= \textrm{[total Cdk]} - \textrm{[MPF]} - \textrm{[preMPF]} \\
k_{25} &= V_{25}' \left( \textrm{[total Cdc25]} - \textrm{[Cdc25P]} \right) + V_{25}'' \textrm{[Cdc25P]} \\
k_{\textrm{wee}} &= V_{\textrm{wee}}' \textrm{[Wee1P]} + V_{\textrm{wee}}'' \left( \textrm{[total Wee1]} - \textrm{[Wee1P]} \right) \\
k_2 &= V_2' \left( \textrm{[total APC]} - \textrm{[APC]} \right) + V_2'' \textrm{[APC]}'''

    TEX = DERIVATIVE_LINES + ' \\\\\n' + DEFINITION_LINES

    SPECIES = ('Cyclin', 'MPF', 'preMPF', 'Cdc25P', 'Wee1P', 'IEP', 'APC')

    def test_the_dynamic_lines_parse(self):
        system = ODESystem.from_tex(self.DERIVATIVE_LINES)
        self.assertEqual(set(map(str, system.dependents)), set(self.SPECIES))

    def test_a_fraction_inside_a_fraction_survives(self):
        # The Cdc25P line is a difference of two quotients, each with a `\left(...\right)`
        # and text-mode names inside.  Normalising the commands away first is what leaves
        # the \frac with no nested braces for the fraction regex to trip over.
        system = ODESystem.from_tex(self.DERIVATIVE_LINES)
        rate = system.derivative_dict[sympy.Symbol('Cdc25P')]
        for expected in ('k_a', 'K_a', 'k_b', 'K_b', 'PPase', 'totalCdc25'):
            with self.subTest(name=expected):
                self.assertIn(sympy.Symbol(expected), rate.free_symbols)

    def test_the_species_names_come_through_intact(self):
        system = ODESystem.from_tex(self.DERIVATIVE_LINES)
        names = set(map(str, system.symbols))
        for expected in self.SPECIES + ('Cdk', 'k_wee', 'k_25', 'totalCdc25'):
            with self.subTest(name=expected):
                self.assertIn(expected, names)

    def test_every_definition_is_recognised_as_one(self):
        # All four mention a dependent variable, so none of them can be a constraint.
        with self.assertRaises(TexParseError) as caught:
            ODESystem.from_tex(self.TEX)
        message = str(caught.exception)
        for defined in ('Cdk', 'k_25', 'k_wee', 'k_2'):
            with self.subTest(name=defined):
                self.assertIn(defined, message)

    def test_the_whole_thing_with_definitions_substituted(self):
        system = ODESystem.from_tex(self.TEX, sub_defns=True)
        gone = {sympy.Symbol(n) for n in ('Cdk', 'k_25', 'k_wee', 'k_2')}
        self.assertFalse(gone & set(system.symbols))
        self.assertEqual(set(map(str, system.dependents)), set(self.SPECIES))
        # The primes are names here, so both forms of each enzyme survive as parameters.
        for expected in ('V_25prime', 'V_25primeprime', 'V_weeprime', 'V_2primeprime'):
            with self.subTest(name=expected):
                self.assertIn(sympy.Symbol(expected), system.parameters)

    def test_it_reduces(self):
        from desr import ODETranslator
        system = ODESystem.from_tex(self.TEX, sub_defns=True)
        translator = ODETranslator.from_ode_system(system)
        self.assertTrue(translator.translate(system).symbols)


class TestMultiLine(TestCase):
    ''' tex_to_sympy over a whole block. '''

    def test_lines_separated_by_double_backslash_only(self):
        # A paper's tex is one long string broken by `\\`, not necessarily by newlines.
        parsed = tex_to_sympy(r'\frac{d}{dt}x &= -k x \\ \frac{d}{dt}y &= k x')
        self.assertEqual(len(parsed), 2)
        for equation in parsed:
            self.assertIsInstance(equation.lhs, sympy.Derivative)


class TestInitialConditionLines(TestCase):
    '''
    `C(0) = C_0` is how a paper writes an initial condition, and how `to_tex` has always
    written one.  It has to be recognised before implicit multiplication runs, which would
    read `C(0)` as `C*0`.  It comes back as an applied undefined function, the way a
    derivative line comes back as a sympy.Derivative, so the caller classifies the line by
    what it is.
    '''

    DYNAMIC = '\n'.join([r'\frac{d}{dt}x &= k_1 - k_2 x \\',
                         r'\frac{d}{dt}y &= k_2 x - k_3 y'])

    def test_a_variable_at_zero_parses_as_a_function_at_zero(self):
        parsed = tex_to_sympy_one_line(r'x\left(0\right) &= x_0')
        self.assertEqual(parsed, sympy.Eq(sympy.Function('x')(0), sympy.Symbol('x_0')))

    def test_the_bare_parenthesis_form_works_too(self):
        self.assertEqual(tex_to_sympy_one_line(r'x(0) = x_0'),
                         sympy.Eq(sympy.Function('x')(0), sympy.Symbol('x_0')))

    def test_a_zero_initial_condition_parses(self):
        self.assertEqual(tex_to_sympy_one_line(r'x\left(0\right) &= 0'),
                         sympy.Eq(sympy.Function('x')(0), sympy.Integer(0)))

    def test_the_number_may_come_first(self):
        # A paper may write `x_0 = x(0)`.
        self.assertEqual(tex_to_sympy_one_line(r'x_0 = x\left(0\right)'),
                         sympy.Eq(sympy.Function('x')(0), sympy.Symbol('x_0')))

    def test_both_sides_at_zero_is_refused(self):
        with self.assertRaises(TexParseError) as caught:
            tex_to_sympy_one_line(r'x\left(0\right) = y\left(0\right)')
        self.assertIn('both sides', str(caught.exception))

    def test_a_product_with_a_parenthesis_is_still_a_product(self):
        x, y = sympy.symbols('x y')
        self.assertEqual(tex_to_sympy_one_line(r'k_3 (x - y)'),
                         sympy.Symbol('k_3') * (x - y))

    def test_a_variable_at_t_is_not_an_initial_condition(self):
        # Only the literal zero, so nothing else changes meaning.
        parsed = tex_to_sympy_one_line(r'k_4 = x(t)')
        self.assertEqual(parsed, sympy.Eq(sympy.Symbol('k_4'),
                                          sympy.Symbol('x') * sympy.Symbol('t')))

    def test_a_text_mode_species_name_works(self):
        parsed = tex_to_sympy_one_line(r'\textrm{[ES]}\left(0\right) &= 0')
        self.assertEqual(parsed, sympy.Eq(sympy.Function('ES')(0), sympy.Integer(0)))

    def test_from_tex_records_it_as_an_initial_condition(self):
        system = ODESystem.from_tex(
            self.DYNAMIC + '\\\\\n' + r'x\left(0\right) &= x_0')
        self.assertEqual(system.initial_conditions,
                         {sympy.Symbol('x'): sympy.Symbol('x_0')})

    def test_from_tex_takes_a_zero(self):
        system = ODESystem.from_tex(
            self.DYNAMIC + '\\\\\n' + r'x\left(0\right) &= 0')
        self.assertEqual(system.initial_conditions,
                         {sympy.Symbol('x'): sympy.Symbol('x_0')})
        self.assertEqual(system.parameter_values[sympy.Symbol('x_0')], sympy.Integer(0))

    def test_an_initial_condition_for_something_that_is_not_dependent_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            ODESystem.from_tex(self.DYNAMIC + '\\\\\n' + r'k_1\left(0\right) &= 1')
        self.assertIn('not a dependent variable', str(caught.exception))


class TestParameterValueLines(TestCase):
    '''
    A symbol given a number is a value, which is what a paper means by `k_1 = 2`.  A
    constraint cannot hold one anyway: `add_constraint` refuses a side that is zero, so
    `C_0 = 0` had nowhere to go at all.
    '''

    DYNAMIC = '\n'.join([r'\frac{d}{dt}x &= k_1 - k_2 x \\',
                         r'\frac{d}{dt}y &= k_2 x - k_3 y'])

    def test_a_number_is_a_value_not_a_constraint(self):
        system = ODESystem.from_tex(self.DYNAMIC + '\\\\\n' + r'k_1 &= 2')
        self.assertEqual(system.parameter_values[sympy.Symbol('k_1')], sympy.Integer(2))
        self.assertEqual(system.constraints, [])

    def test_a_zero_value_is_taken(self):
        system = ODESystem.from_tex(self.DYNAMIC + '\\\\\n' + r'k_1 &= 0')
        self.assertEqual(system.parameter_values[sympy.Symbol('k_1')], sympy.Integer(0))

    def test_the_number_may_come_first(self):
        system = ODESystem.from_tex(self.DYNAMIC + '\\\\\n' + r'2 &= k_1')
        self.assertEqual(system.parameter_values[sympy.Symbol('k_1')], sympy.Integer(2))

    def test_a_relation_among_parameters_is_still_a_constraint(self):
        system = ODESystem.from_tex(self.DYNAMIC + '\\\\\n' + r'k_4 &= k_1 + k_2')
        self.assertEqual(len(system.constraints), 1)
        self.assertEqual(system.parameter_values, {})

    def test_a_dependent_given_a_number_is_still_a_definition(self):
        # It says something about the trajectory, not about a constant, so it keeps going
        # down the definition path and its message.
        with self.assertRaises(TexParseError) as caught:
            ODESystem.from_tex(self.DYNAMIC + '\\\\\n' + r'x &= 2')
        self.assertIn('sub_defns', str(caught.exception))

    def test_values_are_checked_against_the_constraints_in_the_same_tex(self):
        tex = self.DYNAMIC + '\\\\\n' + '\\\\\n'.join(
            [r'k_4 &= k_1 + k_2', r'k_1 &= 1', r'k_2 &= 1', r'k_4 &= 3'])
        with self.assertRaises(ValueError) as caught:
            ODESystem.from_tex(tex)
        self.assertIn('do not satisfy it', str(caught.exception))


class TestTexRoundTrip(TestCase):
    '''
    `to_tex` and `from_tex` should be inverse.  Initial conditions were being written and
    not read, and values were not being written at all.
    '''

    TEX = '\n'.join(['\\frac{dS}{dt} &= - k_1 E S + k_{-1} C \\\\',
                     '\\frac{dE}{dt} &= - k_1 E S + k_{-1} C + k_2 C \\\\',
                     '\\frac{dC}{dt} &= k_1 E S - k_{-1} C - k_2 C \\\\',
                     '\\frac{dP}{dt} &= k_2 C'])

    def system(self):
        system = ODESystem.from_tex(self.TEX)
        system.update_initial_conditions({'S': 'S_0', 'E': 'E_0', 'C': 0, 'P': 0})
        system.add_constraint('K_m = (k_m1 + k_2) / k_1')
        return system

    def test_the_equations_survive(self):
        system = self.system()
        self.assertEqual(ODESystem.from_tex(system.to_tex()).derivative_dict,
                         system.derivative_dict)

    def test_the_initial_conditions_survive(self):
        system = self.system()
        self.assertEqual(ODESystem.from_tex(system.to_tex()).initial_conditions,
                         system.initial_conditions)

    def test_the_constraints_survive(self):
        system = self.system()
        self.assertEqual(ODESystem.from_tex(system.to_tex()).constraints,
                         system.constraints)

    def test_the_values_survive(self):
        system = self.system()
        self.assertEqual(ODESystem.from_tex(system.to_tex()).parameter_values,
                         system.parameter_values)

    def test_to_tex_writes_the_values(self):
        self.assertIn('C_{0} &= 0', self.system().to_tex())

    def test_a_second_pass_changes_nothing(self):
        once = self.system().to_tex()
        self.assertEqual(ODESystem.from_tex(once).to_tex(), once)
