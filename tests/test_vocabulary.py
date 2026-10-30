'''
The 2.0 vocabulary, pinned exhaustively.

This is the interface another project is to be built on, so every name, every count
relation and every classification boundary is asserted here rather than left to be
inferred from the examples.

The vocabulary, and why it is this:

    symbols              n = p + 1 + q   the collective; what n counts
    parameters           p               no equation
    variables            1 + q           independent_variable + dependents
    independent_variable 1               `indep_var` is a short alias
    dependents           q               Hubert & Labahn's "state variables"

Hubert & Labahn have no word for the collective, so "variable" was borrowed for the job
and lost its own meaning -- at the top of their section 7 they separate "state variables"
from "constant parameters", and a page later call alpha = k_3/(k_1 a) a variable.  Naming
the collective is what lets "variable" mean something again.

Orthogonal to the partition, and stored rather than inferred because a reduced system
cannot re-derive it: `auxiliary_names`.  Derived, and deliberately without an
accessor: "constant", which is implied for every parameter, impossible for the
independent variable, and informative only for a dependent -- where it is
`conserved_variables`.
'''
from unittest import TestCase

import sympy

from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator, ReductionScheme
from desr.sympy_helper import sympify_expression


class TestThePartition(TestCase):
    ''' Every symbol is exactly one of independent variable, parameter, dependent. '''

    def setUp(self):
        self.system = ODESystem.from_equations(
            ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)'])

    def test_the_three_groups_partition_the_symbols(self):
        parts = ((self.system.independent_variable,) + self.system.parameters
                 + self.system.dependents)
        self.assertEqual(len(self.system.symbols), len(parts))
        self.assertEqual(set(self.system.symbols), set(parts))

    def test_the_groups_do_not_overlap(self):
        self.assertEqual(set(), set(self.system.parameters) & set(self.system.dependents))
        self.assertNotIn(self.system.independent_variable, self.system.parameters)
        self.assertNotIn(self.system.independent_variable, self.system.dependents)

    def test_the_partition_is_the_trichotomy_on_the_stored_derivative(self):
        '''
        Completeness is not a judgement call: a stored derivative is the independent
        variable's, or None, or an expression.  Nothing else is representable.
        '''
        for symbol, derivative in zip(self.system.symbols, self.system.derivatives):
            with self.subTest(symbol=str(symbol)):
                if symbol == self.system.independent_variable:
                    self.assertEqual(1, derivative)
                elif derivative is None:
                    self.assertIn(symbol, self.system.parameters)
                else:
                    self.assertIn(symbol, self.system.dependents)

    def test_n_equals_p_plus_one_plus_q(self):
        self.assertEqual(len(self.system.symbols),
                         self.system.num_parameters + 1 + self.system.num_dependents)

    def test_variables_are_the_symbols_that_vary(self):
        self.assertEqual((self.system.independent_variable,) + self.system.dependents,
                         self.system.variables)

    def test_symbols_are_parameters_plus_variables(self):
        ''' The sentence that goes in the paper, as an assertion. '''
        self.assertEqual(len(self.system.symbols),
                         self.system.num_parameters + self.system.num_variables)
        self.assertEqual(set(self.system.symbols),
                         set(self.system.parameters) | set(self.system.variables))

    def test_num_variables_is_one_plus_num_dependents(self):
        self.assertEqual(1 + self.system.num_dependents, self.system.num_variables)

    def test_every_count_matches_its_collection(self):
        for count, collection in (('num_parameters', 'parameters'),
                                  ('num_dependents', 'dependents'),
                                  ('num_variables', 'variables')):
            with self.subTest(accessor=count):
                self.assertEqual(len(getattr(self.system, collection)),
                                 getattr(self.system, count))

    def test_the_symbols_are_in_canonical_order(self):
        ''' Parameters, then the independent variable, then the dependents. '''
        self.assertEqual(tuple(self.system.parameters)
                         + (self.system.independent_variable,)
                         + tuple(self.system.dependents),
                         self.system.symbols)
        self.assertEqual(self.system.num_parameters,
                         self.system.independent_variable_index)


class TestTheIndependentVariableAlias(TestCase):
    ''' `indep_var` is a short alias for `independent_variable`, not a second thing. '''

    def setUp(self):
        self.system = ODESystem.from_equations(['dx/dt = a*x'])

    def test_the_alias_is_the_same_value(self):
        self.assertEqual(self.system.independent_variable, self.system.indep_var)

    def test_the_index_alias_is_the_same_value(self):
        self.assertEqual(self.system.independent_variable_index,
                         self.system.indep_var_index)

    def test_the_index_is_the_position_in_symbols(self):
        self.assertEqual(self.system.symbols[self.system.independent_variable_index],
                         self.system.independent_variable)

    def test_the_alias_follows_a_rename(self):
        self.system.rename_indep_var('tau')
        self.assertEqual(sympy.Symbol('tau'), self.system.independent_variable)
        self.assertEqual(sympy.Symbol('tau'), self.system.indep_var)

    def test_the_library_reads_the_attribute_by_its_full_name(self):
        '''
        The alias exists for typing at a prompt.  The library itself reads the attribute by
        its full name, so that a search for `independent_variable` finds every place the
        library asks a system for it -- checked here rather than trusted.

        Only *attribute access* is covered.  The private `_indep_var` field, the
        `indep_var=` constructor keyword, parameters of functions that take one, and method
        names like `rename_indep_var` and `indep_var_invariant` keep the short spelling;
        they are not reads of this attribute.
        '''
        import pathlib
        import re
        source = pathlib.Path(__file__).parent.parent / 'desr'
        attribute_read = re.compile(r'\.indep_var(_index)?(?!\w)')
        offenders = []
        for path in sorted(source.glob('*.py')):
            for number, line in enumerate(path.read_text().splitlines(), start=1):
                if attribute_read.search(line):
                    offenders.append('{}:{}: {}'.format(path.name, number, line.strip()))
        self.assertEqual([], offenders)

    def test_both_aliases_are_the_same_descriptor(self):
        ''' Not two implementations that could drift, one property under two names. '''
        self.assertIs(type(self.system).independent_variable,
                      type(self.system).indep_var)
        self.assertIs(type(self.system).independent_variable_index,
                      type(self.system).indep_var_index)


class TestConservedAndConstant(TestCase):
    '''
    `conserved_variables` is the dependents whose equation says zero.  "Constant" is the
    umbrella for those *and* the parameters, and has no accessor on purpose.
    '''

    def setUp(self):
        self.system = ODESystem.from_equations(['dx/dt = a*x - b*x*y',
                                                'dy/dt = b*x*y - c*y',
                                                'dw/dt = 0'])

    def test_a_zero_equation_is_a_conserved_dependent(self):
        self.assertEqual((sympy.Symbol('w'),), self.system.conserved_variables)
        self.assertIn(sympy.Symbol('w'), self.system.dependents)
        self.assertNotIn(sympy.Symbol('w'), self.system.parameters)

    def test_conserved_is_a_subset_of_the_dependents(self):
        self.assertLessEqual(set(self.system.conserved_variables),
                             set(self.system.dependents))

    def test_a_parameter_is_not_conserved(self):
        ''' A parameter does not vary, but it is not conserved *along* anything. '''
        for parameter in self.system.parameters:
            with self.subTest(parameter=str(parameter)):
                self.assertNotIn(parameter, self.system.conserved_variables)

    def test_a_system_with_nothing_conserved_reports_nothing(self):
        system = ODESystem.from_equations(['dx/dt = a*x'])
        self.assertEqual((), system.conserved_variables)

    def test_make_parameter_moves_a_conserved_dependent_across(self):
        ''' Saying `dw/dt = 0` is a statement; this is how you say the other thing. '''
        self.system.make_parameter('w')
        self.assertIn(sympy.Symbol('w'), self.system.parameters)
        self.assertNotIn(sympy.Symbol('w'), self.system.dependents)
        self.assertEqual((), self.system.conserved_variables)

    def test_constant_has_no_accessor(self):
        '''
        Deliberate: it would be implied for every parameter, impossible for the independent
        variable, and informative only where it coincides with `conserved_variables`.  The
        1.x name meant "parameters", so re-using it would silently answer a wider question.
        '''
        for name in ('constant_variables', 'num_constants', 'non_constant_variables',
                     'constants'):
            with self.subTest(name=name):
                self.assertFalse(hasattr(self.system, name))


class TestMonomialsVersusNames(TestCase):
    '''
    An invariant has two faces, and the accessors say which one they give.

    `invariant_monomials()` gives the expressions -- `d/K`, `k*t/h` -- and
    `invariant_names()` gives the names those get -- `kappa0`, `tau`.  Neither is spelled
    with the bare word, so there is no accessor whose aspect you have to guess.
    '''

    def setUp(self):
        self.system = ODESystem.from_equations(
            ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)'])
        self.translator = ODETranslator.from_ode_system(self.system)

    def test_monomials_are_expressions_and_names_are_symbols(self):
        monomials = list(self.translator.invariant_monomials())
        names = self.translator.invariant_names(self.system)
        self.assertTrue(all(not m.is_Symbol for m in monomials[:3]))
        self.assertTrue(all(n.is_Symbol for n in names))

    def test_they_are_in_the_same_order_and_pair_up(self):
        monomials = list(self.translator.invariant_monomials())
        names = self.translator.invariant_names(self.system)
        self.assertEqual(len(monomials), len(names))
        self.assertEqual(dict(zip(names, monomials)),
                         self.translator.named_symbols(self.system))

    def test_every_invariant_is_a_monomial(self):
        '''
        Each column of V_n is an exponent vector, so no invariant is a sum.  That is what
        makes the reduction work -- these generate, and any rational invariant is a
        rational function of them -- and it is why the accessor says `monomials`.
        '''
        for monomial in self.translator.invariant_monomials():
            with self.subTest(monomial=str(monomial)):
                self.assertFalse(monomial.is_Add,
                                 '{} is a sum, so it is not a monomial'.format(monomial))
                # A Laurent monomial: a product of powers, nothing else.
                self.assertTrue(monomial.is_Mul or monomial.is_Pow or monomial.is_Symbol)

    def test_the_auxiliary_pair_matches_too(self):
        monomials = list(self.translator.auxiliary_monomials())
        names = self.translator.auxiliary_names(self.system, scheme='general')
        self.assertEqual(len(monomials), len(names))
        self.assertEqual(('x0', 'x1', 'x2'), tuple(map(str, names)))

    def test_the_parameter_scheme_names_no_auxiliaries(self):
        ''' It normalises them to 1 and keeps none, so there is nothing to name. '''
        self.assertEqual((), self.translator.auxiliary_names(self.system,
                                                             scheme='parameter'))
        # The monomials still exist -- they are a property of the multiplier, not of a
        # reduction that was run.
        self.assertEqual(3, len(self.translator.auxiliary_monomials()))

    def test_the_bare_words_are_gone(self):
        '''
        `invariants` and `auxiliaries` were ambiguous between the two faces, and the
        ambiguity was across a workflow rather than within one class.  Neither spelling
        survives.
        '''
        for name in ('invariants', 'auxiliaries'):
            with self.subTest(name=name):
                self.assertFalse(hasattr(self.translator, name))

    def test_names_match_what_translate_actually_produces(self):
        for scheme in ('parameter', 'general'):
            with self.subTest(scheme=scheme):
                if scheme == 'parameter':
                    reduced = self.translator.translate_parameter(self.system)
                    expected = self.translator.invariant_names(self.system, scheme)
                else:
                    reduced = self.translator.translate_general(self.system)
                    expected = ((reduced.independent_variable,)
                                + self.translator.auxiliary_names(self.system, scheme)
                                + self.translator.invariant_names(self.system, scheme))
                self.assertEqual(tuple(reduced.symbols), tuple(expected))


class TestStoredProvenance(TestCase):
    '''
    What a reduced system cannot re-derive about itself, and so carries: which scheme made
    it, how many symmetries went, and which of its symbols are auxiliaries.
    '''

    def setUp(self):
        self.system = ODESystem.from_equations(
            ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)'])
        self.translator = ODETranslator.from_ode_system(self.system)

    def test_an_original_system_has_no_provenance(self):
        self.assertIsNone(self.system.scheme)
        self.assertIsNone(self.system.num_scaling_symmetries)
        self.assertEqual((), self.system.auxiliary_names)

    def test_the_parameter_scheme_records_itself(self):
        reduced = self.translator.translate_parameter(self.system)
        self.assertEqual(ReductionScheme.PARAMETER, reduced.scheme)
        self.assertEqual(self.translator.r, reduced.num_scaling_symmetries)

    def test_the_parameter_scheme_has_no_auxiliaries(self):
        ''' It normalises them to 1 and discards them, so there are none to record. '''
        reduced = self.translator.translate_parameter(self.system)
        self.assertEqual((), reduced.auxiliary_names)

    def test_the_general_scheme_records_itself_and_its_auxiliaries(self):
        reduced = self.translator.translate_general(self.system)
        self.assertEqual(ReductionScheme.GENERAL, reduced.scheme)
        self.assertEqual(self.translator.r, reduced.num_scaling_symmetries)
        self.assertEqual(('x0', 'x1', 'x2'),
                         tuple(map(str, reduced.auxiliary_names)))

    def test_dropping_the_auxiliaries_records_none(self):
        reduced = self.translator.translate_general(self.system, include_aux_vars=False)
        self.assertEqual((), reduced.auxiliary_names)
        self.assertEqual(ReductionScheme.GENERAL, reduced.scheme)

    def test_the_dep_var_scheme_records_itself(self):
        system = ODESystem.from_equations(['dz1/dt = z1*(1+z1*z2)',
                                           'dz2/dt = z2*(1/t - z1*z2)'])
        translator = ODETranslator.from_ode_system(system)
        reduced = translator.translate_dep_var(system)
        self.assertEqual(ReductionScheme.DEP_VAR, reduced.scheme)
        self.assertEqual(translator.r, reduced.num_scaling_symmetries)

    def test_is_auxiliary_takes_a_string_or_a_symbol(self):
        reduced = self.translator.translate_general(self.system)
        self.assertTrue(reduced.is_auxiliary('x0'))
        self.assertTrue(reduced.is_auxiliary(sympy.Symbol('x0')))
        self.assertFalse(reduced.is_auxiliary('y3'))

    def test_the_provenance_survives_a_copy(self):
        reduced = self.translator.translate_general(self.system)
        duplicate = reduced.copy()
        self.assertEqual(reduced.scheme, duplicate.scheme)
        self.assertEqual(reduced.num_scaling_symmetries,
                         duplicate.num_scaling_symmetries)
        self.assertEqual(reduced.auxiliary_names, duplicate.auxiliary_names)

    def test_the_provenance_survives_a_substitution(self):
        reduced = self.translator.translate_general(self.system)
        substituted = reduced.subs({'y3': 'Y'})
        self.assertEqual(reduced.scheme, substituted.scheme)
        self.assertEqual(reduced.num_scaling_symmetries,
                         substituted.num_scaling_symmetries)
        self.assertEqual(reduced.auxiliary_names, substituted.auxiliary_names)

    def test_an_auxiliary_is_not_pruned_for_having_no_equation(self):
        '''
        `subs` drops symbols that appear in no expression.  An auxiliary that is a monomial
        in the original parameters has no equation, so it would be pruned -- losing part of
        the record of how the system was reduced, which it cannot reconstruct.
        '''
        reduced = self.translator.translate_general(self.system)
        substituted = reduced.subs({'y3': 'Y'})
        self.assertEqual(('x0', 'x1', 'x2'),
                         tuple(map(str, substituted.auxiliary_names)))
        for auxiliary in reduced.auxiliary_names:
            self.assertIn(auxiliary, substituted.symbols)

    def test_renaming_an_equation_less_parameter_leaves_the_name_alone(self):
        '''
        A known limitation, pinned so it is not mistaken for this work: `subs` detects a
        renaming only for a symbol that occurs in an expression, so naming a parameter with
        no equation does nothing.  It used to *delete* it; now it is left as it was, which
        is at least not lossy.
        '''
        reduced = self.translator.translate_general(self.system)
        substituted = reduced.subs({'x0': 'w0'})
        self.assertIn(sympy.Symbol('x0'), substituted.symbols)
        self.assertNotIn(sympy.Symbol('w0'), substituted.symbols)


class TestAuxiliaryNessIsOrthogonal(TestCase):
    '''
    An auxiliary is a parameter or a dependent like any other symbol.  Auxiliary-ness says
    you may choose the value; the kind says whether it varies.
    '''

    def test_an_auxiliary_can_be_a_parameter(self):
        system = ODESystem.from_equations(
            ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)'])
        reduced = ODETranslator.from_ode_system(system).translate_general(system)
        for auxiliary in reduced.auxiliary_names:
            with self.subTest(auxiliary=str(auxiliary)):
                self.assertIn(auxiliary, reduced.parameters)

    def test_an_auxiliary_can_be_a_dependent(self):
        ''' Hubert & Labahn 6.6: the auxiliary is t*z1**2, which evolves. '''
        system = ODESystem.from_equations(
            ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
             'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])
        reduced = ODETranslator.from_ode_system(system).translate_general(system)
        auxiliary, = reduced.auxiliary_names
        self.assertIn(auxiliary, reduced.dependents)
        self.assertNotEqual(0, reduced.derivative_dict[auxiliary])

    def test_both_kinds_of_auxiliary_are_still_auxiliaries(self):
        for equations in (['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )',
                           'dp/dt = s*p*(1 - h*p / n)'],
                          ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
                           'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)']):
            with self.subTest(equations=equations[0]):
                system = ODESystem.from_equations(equations)
                translator = ODETranslator.from_ode_system(system)
                reduced = translator.translate_general(system)
                self.assertEqual(translator.r, len(reduced.auxiliary_names))


class TestClassificationSurvivesTheReduction(TestCase):
    '''
    The point of the whole exercise: a monomial in the original parameters alone is a
    parameter of the reduced system, not a dependent whose equation says zero.
    '''

    def setUp(self):
        self.system = ODESystem.from_equations(
            ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)'])
        self.translator = ODETranslator.from_ode_system(self.system)

    def test_a_parameter_monomial_becomes_a_parameter(self):
        reduced = self.translator.translate_general(self.system)
        named = self.translator.named_symbols(self.system, scheme='general')
        dynamic = set(self.system.variables)
        for symbol, monomial in named.items():
            if symbol == reduced.independent_variable:
                continue
            with self.subTest(symbol=str(symbol), monomial=str(monomial)):
                if monomial.free_symbols & dynamic:
                    self.assertIn(symbol, reduced.dependents)
                else:
                    self.assertIn(symbol, reduced.parameters)

    def test_no_reduced_parameter_has_an_equation(self):
        reduced = self.translator.translate_general(self.system)
        for parameter in reduced.parameters:
            with self.subTest(parameter=str(parameter)):
                self.assertNotIn(parameter, reduced.derivative_dict)

    def test_a_reduction_no_longer_destroys_parameter_ness(self):
        '''
        Before this, the general reduction of a six-parameter system reported zero
        parameters and nine dependents, six of them with rate zero.
        '''
        reduced = self.translator.translate_general(self.system)
        self.assertEqual(6, self.system.num_parameters)
        self.assertEqual(6, reduced.num_parameters)
        self.assertEqual(3, reduced.num_dependents)

    def test_only_the_genuinely_moving_symbols_need_integrating(self):
        reduced = self.translator.translate_general(self.system)
        self.assertEqual((), reduced.conserved_variables)
        self.assertEqual(3, len(reduced.dependents))

    def test_a_conserved_dependent_stays_a_dependent(self):
        '''
        The catalyst case, which is what the structural rule exists for.  E*k_0_1/k_2_3 is
        conserved because the catalyst is, and it must not become a parameter.
        '''
        system = ChemicalReactionNetwork.from_diagram(
            'E + S -> E + P \n P -> ').to_ode_system()
        translator = ODETranslator.from_ode_system(system)
        reduced = translator.translate_general(system)

        conserved = sympify_expression('E*k_0_1/k_2_3')
        named = translator.named_symbols(system, scheme='general')
        symbol = next(s for s, monomial in named.items() if monomial == conserved)

        self.assertIn(symbol, reduced.dependents)
        self.assertNotIn(symbol, reduced.parameters)
        self.assertIn(symbol, reduced.conserved_variables)
        self.assertEqual(0, reduced.derivative_dict[symbol])

    def test_the_two_zero_rate_cases_are_told_apart(self):
        '''
        Both print as nothing-that-moves, and only the classification separates them: a
        parameter has no equation, a conserved dependent has one that says zero.
        '''
        catalyst = ChemicalReactionNetwork.from_diagram(
            'E + S -> E + P \n P -> ').to_ode_system()
        reduced = ODETranslator.from_ode_system(catalyst).translate_general(catalyst)
        self.assertTrue(reduced.parameters)
        self.assertTrue(reduced.conserved_variables)
        self.assertEqual(set(), set(reduced.parameters) & set(reduced.conserved_variables))


class TestTheTranslatorAgreesWithTheReducedSystem(TestCase):
    ''' Asking the translator and asking the reduced system must give the same answer. '''

    def setUp(self):
        self.system = ODESystem.from_equations(
            ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)'])
        self.translator = ODETranslator.from_ode_system(self.system)

    def test_the_kinds_agree_under_the_parameter_scheme(self):
        reduced = self.translator.translate_parameter(self.system)
        kinds = self.translator.symbol_kinds(self.system, scheme='parameter')
        self.assertEqual(set(reduced.parameters),
                         {s for s, kind in kinds.items() if kind == 'parameter'})
        self.assertEqual(set(reduced.dependents),
                         {s for s, kind in kinds.items() if kind == 'dependent'})
        self.assertEqual({reduced.independent_variable},
                         {s for s, kind in kinds.items() if kind == 'independent'})

    def test_the_auxiliaries_agree_under_the_general_scheme(self):
        reduced = self.translator.translate_general(self.system)
        kinds = self.translator.symbol_kinds(self.system, scheme='general')
        self.assertEqual(set(reduced.auxiliary_names),
                         {s for s, kind in kinds.items() if kind == 'auxiliary'})

    def test_the_expected_count_matches_what_comes_out(self):
        for scheme, include in (('parameter', True), ('general', True), ('general', False)):
            with self.subTest(scheme=scheme, include_aux_vars=include):
                if scheme == 'parameter':
                    reduced = self.translator.translate_parameter(self.system)
                else:
                    reduced = self.translator.translate_general(
                        self.system, include_aux_vars=include)
                self.assertEqual(
                    self.translator.expected_variable_count(
                        self.system, scheme=scheme, include_aux_vars=include),
                    len(reduced.symbols))

    def test_the_names_match_the_reduced_systems_symbols(self):
        for scheme in ('parameter', 'general'):
            with self.subTest(scheme=scheme):
                if scheme == 'parameter':
                    reduced = self.translator.translate_parameter(self.system)
                else:
                    reduced = self.translator.translate_general(self.system)
                named = self.translator.named_symbols(self.system, scheme=scheme)
                self.assertEqual(set(reduced.symbols), set(named))
