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
        self.assertEqual(set(map(str, system.dependent_variables)), {'x', 'y'})
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
        self.assertNotIn(sympy.Symbol('z'), system.variables)

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
        self.assertEqual(set(map(str, system.dependent_variables)), set(self.SPECIES))

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
        names = set(map(str, system.variables))
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
        self.assertFalse(gone & set(system.variables))
        self.assertEqual(set(map(str, system.dependent_variables)), set(self.SPECIES))
        # The primes are names here, so both forms of each enzyme survive as parameters.
        for expected in ('V_25prime', 'V_25primeprime', 'V_weeprime', 'V_2primeprime'):
            with self.subTest(name=expected):
                self.assertIn(sympy.Symbol(expected), system.parameters)

    def test_it_reduces(self):
        from desr import ODETranslator
        system = ODESystem.from_tex(self.TEX, sub_defns=True)
        translator = ODETranslator.from_ode_system(system)
        self.assertTrue(translator.translate(system).variables)


class TestMultiLine(TestCase):
    ''' tex_to_sympy over a whole block. '''

    def test_lines_separated_by_double_backslash_only(self):
        # A paper's tex is one long string broken by `\\`, not necessarily by newlines.
        parsed = tex_to_sympy(r'\frac{d}{dt}x &= -k x \\ \frac{d}{dt}y &= k x')
        self.assertEqual(len(parsed), 2)
        for equation in parsed:
            self.assertIsInstance(equation.lhs, sympy.Derivative)
