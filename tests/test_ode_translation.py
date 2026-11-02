from unittest import TestCase

import sympy

from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator, ReductionScheme, scale_action
from desr.sympy_helper import monomial_to_powers, sympify_expression


class TestODETranslatorTraps(TestCase):
    ''' Things that used to fail silently, or somewhere other than where the mistake was. '''

    def test_non_unimodular_multiplier_is_refused_at_construction(self):
        scaling = sympy.Matrix([[1, 0, 1, 1, -1], [0, 1, 0, -1, 1]])
        doubled = sympy.Matrix([[0, 0, 2, 0, 0], [0, 0, 0, 1, 0], [1, 1, -2, -1, 0],
                                [0, 0, 0, 0, 1], [0, 1, 0, -1, 1]])
        with self.assertRaises(ValueError) as caught:
            ODETranslator(scaling, hermite_multiplier=doubled)
        self.assertIn('not unimodular', str(caught.exception))
        self.assertIn('roots of degree 2', str(caught.exception))

    def test_a_conserved_variable_does_not_depend_on_its_name(self):
        # A catalyst gets a derivative of identically zero straight from the stoichiometry.
        # That expression has no monomials, so it used to contribute a 0x0 block to the
        # exponent matrix; sympy's hstack absorbs such a block on the left but raises on
        # the right, so the same network reduced or crashed according to where the catalyst
        # fell in the variable ordering -- that is, according to what it was called.
        scalings = {}
        for catalyst in ('e', 'z'):
            network = ChemicalReactionNetwork.from_diagram(
                'x + {c} -> y + {c}'.format(c=catalyst))
            system = network.to_ode_system()
            self.assertEqual(system.derivative_dict[sympy.Symbol(catalyst)],
                             sympy.sympify(0))
            scalings[catalyst] = (list(system.symbols), system.maximal_scaling_matrix())

        # The two differ only by the name of the catalyst, so permuting the columns of one
        # into the other's variable order must give the same scaling action.
        first_variables, first = scalings['e']
        other_variables, other = scalings['z']
        renamed = [sympy.Symbol('z') if v == sympy.Symbol('e') else v
                   for v in first_variables]
        order = [other_variables.index(v) for v in renamed]
        permuted = other[:, order]
        self.assertEqual(first.rref()[0], permuted.rref()[0])

    def test_invariant_constant_that_cancels_is_still_a_variable(self):
        # Hubert & Labahn 6.6, with initial-condition constants adjoined so that each z_i0
        # scales like its z_i.  kappa0 = t0*z10**3 is the reduced initial time: it appears
        # in no right-hand side, but it is one of the n - r coordinates of the reduced
        # system, and translate_parameter_substitutions names it.  Building the reduced
        # system from the symbols surviving into the equations used to lose it.
        eq1 = 'dz1/dt = z1/t * (-2/3 + 1/3 *z1**5 *z2 )'
        eq2 = 'dz2/dt = z2/t * ( 10/3 - 2/3 *z1**5 *z2 + z1**2*z2/t )'
        system = ODESystem.from_equations([eq1, eq2])
        t0, z10, z20 = sympy.symbols('t0 z10 z20')
        extended = ODESystem(list(system.symbols) + [t0, z10, z20],
                             list(system.derivatives) + [None, None, None],
                             indep_var=system.indep_var)
        # Parameters, then the independent variable, then the dependents.  The scaling
        # matrix below is unchanged by the move because its pattern repeats with period 3.
        extended.reorder_variables(['t0', 'z10', 'z20', 't', 'z1', 'z2'])
        translator = ODETranslator(sympy.Matrix([[3, -1, 5, 3, -1, 5]]),
                                     symbols_domain=extended.symbols)
        reduced = translator.translate(extended)

        self.assertEqual(len(reduced.symbols),
                         len(extended.symbols) - translator.r)
        # Every symbol the substitutions mention must exist in the reduced system.
        mentioned = set()
        for value in translator.translate_parameter_substitutions(extended).values():
            mentioned |= value.free_symbols
        self.assertFalse(mentioned - set(reduced.symbols))

    def test_rational_section_returns_its_result(self):
        equations = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)']
        system = ODESystem.from_equations(equations)
        translator = ODETranslator.from_ode_system(system)
        section = translator.rational_section()
        self.assertIsNotNone(section)
        self.assertEqual(section.shape, (1, translator.r))

    def test_reverse_translate_dispatches_by_count(self):
        translator = ODETranslator(scaling_matrix=sympy.Matrix([[0, 1, -1]]),
                                     hermite_multiplier=sympy.Matrix([[0, 1, 0],
                                                                      [1, 0, 1],
                                                                      [0, 0, 1]]))
        t, x0, y0 = sympy.symbols('t x0 y0')
        self.assertEqual(translator.reverse_translate((x0, y0), 0), (x0, y0 / x0))
        self.assertEqual(translator.reverse_translate((x0, t, y0), 0), (t, x0, y0 / x0))
        # The general branch checks its inputs are genuine solutions, so use the exact one.
        general = ODETranslator(scaling_matrix=sympy.Matrix([[3, -1, 5]]),
                                 hermite_multiplier=sympy.Matrix([[1, 1, -1],
                                                                  [2, 3, 2],
                                                                  [0, 0, 1]]))
        exact = (t,
                 sympy.sympify('c3/(t**(1/3)*(ln(t-c1)-ln(t)+c2)**(2/3))'),
                 sympy.sympify('c1/(t*(ln(t-c1)-ln(t)+c2))'),
                 sympy.sympify('t/(c1 - t)'))
        self.assertEqual(general.reverse_translate(exact),
                         general.reverse_translate_general(exact))


class TestColumnLookup(TestCase):
    '''
    The column numbers a reduction depends on, derived rather than hardcoded.

    They stay numbers -- the author reports them in the paper -- but where each one came
    from is now in the code that uses it, instead of being a magic constant.
    '''

    def setUp(self):
        self.system = ODESystem.from_equations(
            ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
             'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
        self.translator = ODETranslator.from_ode_system(self.system)

    def test_the_index_is_an_int(self):
        ''' Reportable, and what the multiplier operations take. '''
        column = self.translator.column_of('e_0*k_1*t')
        self.assertIsInstance(column, int)

    def test_column_monomials_is_indexed_by_column(self):
        monomials = self.translator.column_monomials()
        self.assertEqual(self.translator.herm_mult.cols, len(monomials))
        for column, monomial in enumerate(monomials):
            with self.subTest(column=column):
                self.assertEqual(column, self.translator.column_of(monomial))

    def test_column_names_is_indexed_the_same_way(self):
        names = self.translator.column_names(self.system)
        self.assertEqual(self.translator.herm_mult.cols, len(names))
        for column, name in enumerate(names):
            with self.subTest(column=column):
                self.assertEqual(column, self.translator.column_of(name, self.system))

    def test_the_auxiliary_columns_come_first(self):
        ''' V_i then V_n, which is the multiplier's own column order. '''
        names = self.translator.column_names(self.system)
        self.assertEqual(('x0', 'x1'), tuple(map(str, names[:self.translator.r])))

    def test_auxiliary_columns_are_named_even_under_the_parameter_scheme(self):
        '''
        The reduced system keeps none of them, so `auxiliary_names` is empty -- but the
        columns exist in the multiplier, and naming them is what lets an operation refer
        to one.
        '''
        self.assertEqual((), self.translator.auxiliary_names(self.system))
        names = self.translator.column_names(self.system)
        self.assertEqual(self.translator.r,
                         len([n for n in names if str(n).startswith('x')]))

    def test_a_monomial_lookup_needs_no_system(self):
        ''' The monomial is what the mathematics points at; it resolves on its own. '''
        self.assertEqual(4, self.translator.column_of('e_0*k_1*t'))

    def test_a_name_lookup_without_a_system_is_refused_rather_than_guessed(self):
        with self.assertRaises(ValueError):
            self.translator.column_of('tau')
        self.assertEqual(4, self.translator.column_of('tau', self.system))

    def test_an_unknown_column_says_what_the_columns_are(self):
        with self.assertRaises(ValueError) as caught:
            self.translator.column_of('nonsense')
        message = str(caught.exception)
        self.assertIn('nonsense', message)
        self.assertIn('e_0*k_1*t', message)

    def test_a_sympy_expression_works_as_well_as_a_string(self):
        self.assertEqual(self.translator.column_of('c/e_0'),
                         self.translator.column_of(sympify_expression('c/e_0')))

    def test_derived_indices_reproduce_a_hand_indexed_reduction(self):
        '''
        The Segel (21a-d) reduction from `docs/source/examples/segel_slemrod.rst`, whose
        five column operations carry the indices 7, 2, 6, 9, 3, 5 and 4.  Each is recovered
        by naming the monomial it refers to, and the resulting invariants match the page.
        '''
        system = ODESystem.from_tex(
            r'''\frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\
                \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c''')
        system.update_initial_conditions({'s': 's_0'})
        system.add_constraint('K_m', '(k_2 + k_m1) / k_1')
        system.add_constraint('epsilon', 'e_0 / (s_0 + K_m)')
        system.reorder_variables(['s_0', 'e_0', 'K_m', 'k_1', 'k_2', 'k_m1', 'epsilon',
                                  't', 's', 'c'])
        translator = ODETranslator.from_ode_system(system, naming_scheme={
            'parameters': 'c', 'independent': 'tau', 'dependents': ['s', 'c']})

        columns = {text: translator.column_of(text) for text in
                   ('k_1*s_0*t', 'e_0/s_0', 'epsilon', 'c/s_0', 'K_m/s_0',
                    'k_m1/(k_1*s_0)', 'k_2/(k_1*s_0)')}
        self.assertEqual({'k_1*s_0*t': 7, 'e_0/s_0': 2, 'epsilon': 6, 'c/s_0': 9,
                          'K_m/s_0': 3, 'k_m1/(k_1*s_0)': 5, 'k_2/(k_1*s_0)': 4},
                         columns)

        translator.multiplier_add_columns(columns['k_1*s_0*t'], columns['e_0/s_0'])
        translator.multiplier_add_columns(columns['k_1*s_0*t'], columns['epsilon'],
                                          power=-1)
        translator.multiplier_add_columns(columns['c/s_0'], columns['epsilon'], power=-1)
        translator.multiplier_negate_column(columns['K_m/s_0'])
        translator.multiplier_add_columns(columns['k_m1/(k_1*s_0)'],
                                          columns['k_2/(k_1*s_0)'], power=-1)

        self.assertEqual(
            ['e_0/s_0', 's_0/K_m', 'k_2/(k_1*s_0)', 'k_m1/k_2', 'epsilon',
             'e_0*k_1*t/epsilon', 's/s_0', 'c/(epsilon*s_0)'],
            [str(m) for m in translator.invariant_monomials()])


class TestInvariantClassification(TestCase):
    '''
    Each invariant is a parameter or a dependent variable, decided **structurally** from
    the exponents in the Hermite multiplier and never from its rate of change.

    The rate cannot decide it.  A reduced parameter always has rate zero, and a reduced
    dependent may have rate zero because it is conserved, so both print as `dy/dt = 0`.
    '''

    def setUp(self):
        self.prey_predator = ODESystem.from_equations(
            'dn/dt = n*( r*(1 - n/K) - k*p/(n+d) );dp/dt = s*p*(1 - h*p / n)'.split(';'))
        self.hl66 = ODESystem.from_equations(
            ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
             'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])

    def test_a_parameter_monomial_is_a_parameter(self):
        translator = ODETranslator.from_ode_system(self.prey_predator)
        self.assertEqual((sympy.sympify('d/K'), sympy.sympify('h*r/k'),
                          sympy.sympify('h*s/k')),
                         translator.parameter_invariants(self.prey_predator))

    def test_an_invariant_touching_a_dependent_row_is_a_dependent(self):
        translator = ODETranslator.from_ode_system(self.prey_predator)
        self.assertEqual((sympy.sympify('n/K'), sympy.sympify('h*p/K')),
                         translator.dependent_invariants(self.prey_predator))

    def test_there_are_as_many_dependent_invariants_as_dependents(self):
        '''
        The clock is reported separately, so this count is exactly the user's remark that
        the loss of symbols is among the parameters only.
        '''
        translator = ODETranslator.from_ode_system(self.prey_predator)
        self.assertEqual(self.prey_predator.num_dependents,
                         len(translator.dependent_invariants(self.prey_predator)))

    def test_the_kinds_partition_the_invariants(self):
        for system in (self.prey_predator, self.hl66):
            with self.subTest(system=repr(system).splitlines()[1]):
                translator = ODETranslator.from_ode_system(system)
                groups = [translator.parameter_invariants(system),
                          translator.dependent_invariants(system)]
                clock = translator.indep_var_invariant(system)
                if clock is not None:
                    groups.append((clock,))
                found = [invariant for group in groups for invariant in group]
                self.assertEqual(len(translator.invariant_monomials()), len(found))
                self.assertEqual(set(translator.invariant_monomials()), set(found))

    def test_a_system_with_no_parameters_has_no_parameter_invariants(self):
        translator = ODETranslator.from_ode_system(self.hl66)
        self.assertEqual((), translator.parameter_invariants(self.hl66))

    def test_an_entangled_clock_is_reported_as_none(self):
        '''
        HL 6.6's invariants are t*z1**3 and z1**2*z2/t.  Both move, so there is no
        invariant that is a function of t and the parameters alone, and saying None is
        better than handing back one that also moves.
        '''
        translator = ODETranslator.from_ode_system(self.hl66)
        self.assertIsNone(translator.indep_var_invariant(self.hl66))

    def test_the_clock_is_found_when_there_is_one(self):
        translator = ODETranslator.from_ode_system(self.prey_predator)
        self.assertEqual(sympy.sympify('k*t/h'),
                         translator.indep_var_invariant(self.prey_predator))

    def test_the_counting_law_holds_exactly_for_the_parameter_scheme(self):
        '''
        The user's remark, made testable: "the loss of symbols is in the parameters only,
        and the invariant dependents involve the new invariant parameters but there's the
        same number of them."
        '''
        catalyst = ChemicalReactionNetwork.from_diagram(
            'E + S -> E + P \n P -> ').to_ode_system()
        for system in (self.prey_predator, self.hl66, catalyst):
            with self.subTest(system=tuple(map(str, system.symbols))):
                translator = ODETranslator.from_ode_system(system)
                law = (len(translator.parameter_invariants(system))
                       == system.num_parameters - translator.r
                       and translator.indep_var_invariant(system) is not None
                       and len(translator.dependent_invariants(system))
                       == system.num_dependents)
                self.assertEqual(translator.scheme(system) == 'parameter', law)

    def test_a_catalyst_invariant_is_conserved_and_still_a_dependent(self):
        '''
        The case this rule exists for, and the one hand-written ODEs do not produce.  A
        catalyst in a reaction network gives an invariant E*k_0_1/k_2_3 whose rate is zero
        -- because the catalyst is conserved -- and which must stay a dependent variable.
        Calling it a parameter would say it is a constant of the model rather than a
        quantity conserved along the trajectory.
        '''
        system = ChemicalReactionNetwork.from_diagram(
            'E + S -> E + P \n P -> ').to_ode_system()
        translator = ODETranslator.from_ode_system(system)
        reduced = translator.translate_general(system)

        # sympify_expression, not sympy.sympify: E is Exp1 to the latter.
        conserved = sympify_expression('E*k_0_1/k_2_3')
        self.assertIn(conserved, list(translator.invariant_monomials()))
        self.assertIn(conserved, translator.dependent_invariants(system))
        self.assertNotIn(conserved, translator.parameter_invariants(system))

        # And its rate really is zero, so the rate could not have told us.
        names = reduced.symbols[1 + translator.r:]
        rates = {invariant: reduced.derivative_dict[name]
                 for name, invariant in zip(names, translator.invariant_monomials())}
        self.assertEqual(0, rates[conserved])

    def test_the_two_faces_of_an_invariant_share_their_keys(self):
        '''
        An invariant is both a new symbol and the monomial that symbol stands for, so the
        map saying what it means and the map saying what kind it is are keyed alike and
        join without reconstruction.
        '''
        translator = ODETranslator.from_ode_system(self.prey_predator)
        meanings = translator.named_symbols(self.prey_predator)
        kinds = translator.symbol_kinds(self.prey_predator)
        self.assertEqual(list(meanings), list(kinds))
        self.assertEqual({'kappa0': ('d/K', 'parameter'),
                          'kappa1': ('h*r/k', 'parameter'),
                          'kappa2': ('h*s/k', 'parameter'),
                          'tau': ('k*t/h', 'independent'),
                          'nu0': ('n/K', 'dependent'),
                          'nu1': ('h*p/K', 'dependent')},
                         {str(name): (str(meanings[name]), kinds[name])
                          for name in meanings})

    def test_the_names_are_the_reduced_systems_own_variables(self):
        ''' The keys must be lookup-friendly: the symbols the reduced system actually has. '''
        for scheme in ('parameter', 'general'):
            with self.subTest(scheme=scheme):
                translator = ODETranslator.from_ode_system(self.prey_predator)
                if scheme == 'parameter':
                    reduced = translator.translate_parameter(self.prey_predator)
                else:
                    reduced = translator.translate_general(self.prey_predator)
                meanings = translator.named_symbols(self.prey_predator, scheme=scheme)
                self.assertEqual(tuple(reduced.symbols), tuple(meanings))

    def test_exactly_one_variable_is_the_one_to_integrate_against(self):
        '''
        The reason the 'independent' kind exists: it says which variable to hand a solver.
        There is one under every scheme, even where no invariant is a rescaled clock.
        '''
        for system in (self.prey_predator, self.hl66):
            for scheme in ('parameter', 'general'):
                if scheme == 'parameter' and system is self.hl66:
                    continue  # the parameter scheme does not apply to HL 6.6
                with self.subTest(scheme=scheme, variables=tuple(map(str, system.symbols))):
                    translator = ODETranslator.from_ode_system(system)
                    kinds = translator.symbol_kinds(system, scheme=scheme)
                    clocks = [name for name, kind in kinds.items() if kind == 'independent']
                    self.assertEqual(1, len(clocks))
                    if scheme == 'general':
                        self.assertEqual(system.indep_var, clocks[0])

    def test_a_rescaled_clock_is_a_dependent_under_the_general_scheme(self):
        '''
        A kind is the variable's role in the reduced system, not the shape of its monomial.
        Under the general scheme k*t/h has an equation and evolves, while t is what the
        system is integrated against -- so k*t/h is a dependent there and the clock under
        the parameter scheme.  The monomial is the same either way.
        '''
        translator = ODETranslator.from_ode_system(self.prey_predator)
        clock = translator.indep_var_invariant(self.prey_predator)
        self.assertEqual(sympify_expression('k*t/h'), clock)

        for scheme, expected in (('parameter', 'independent'), ('general', 'dependent')):
            with self.subTest(scheme=scheme):
                meanings = translator.named_symbols(self.prey_predator, scheme=scheme)
                kinds = translator.symbol_kinds(self.prey_predator, scheme=scheme)
                name = next(n for n, monomial in meanings.items() if monomial == clock)
                self.assertEqual(expected, kinds[name])

    def test_the_general_schemes_extra_variable_is_the_retained_clock(self):
        '''
        n + 1 = 1 (the original independent variable, kept) + r (auxiliaries)
                  + (n - r) (invariants).

        The parameter scheme has no such extra, because it keeps nothing: the auxiliaries
        are normalised away and the clock invariant becomes the new independent variable
        instead of joining one.
        '''
        catalyst = ChemicalReactionNetwork.from_diagram(
            'E + S -> E + P \n P -> ').to_ode_system()
        for system in (self.prey_predator, self.hl66, catalyst):
            with self.subTest(variables=tuple(map(str, system.symbols))):
                translator = ODETranslator.from_ode_system(system)
                n, r = len(system.symbols), translator.r
                kinds = translator.symbol_kinds(system, scheme='general')
                tally = {kind: list(kinds.values()).count(kind)
                         for kind in set(kinds.values())}
                self.assertEqual(1, tally['independent'])
                self.assertEqual(r, tally.get('auxiliary', 0))
                self.assertEqual(n - r, tally.get('parameter', 0) + tally.get('dependent', 0))
                self.assertEqual(n + 1, len(kinds))
                self.assertEqual(n + 1, translator.expected_variable_count(
                    system, scheme='general'))

    def test_auxiliaries_are_their_own_kind(self):
        translator = ODETranslator.from_ode_system(self.prey_predator)
        kinds = translator.symbol_kinds(self.prey_predator, scheme='general')
        auxiliaries = [name for name, kind in kinds.items() if kind == 'auxiliary']
        self.assertEqual(translator.r, len(auxiliaries))
        self.assertEqual(('x0', 'x1', 'x2'), tuple(map(str, auxiliaries)))
        # And the parameter scheme makes none.
        kinds = translator.symbol_kinds(self.prey_predator, scheme='parameter')
        self.assertNotIn('auxiliary', kinds.values())

    def test_dropping_the_auxiliaries_drops_them_from_both_maps(self):
        translator = ODETranslator.from_ode_system(self.prey_predator)
        reduced = translator.translate_general(self.prey_predator, include_aux_vars=False)
        meanings = translator.named_symbols(self.prey_predator, scheme='general',
                                              include_aux_vars=False)
        kinds = translator.symbol_kinds(self.prey_predator, scheme='general',
                                          include_aux_vars=False)
        self.assertEqual(tuple(reduced.symbols), tuple(meanings))
        self.assertNotIn('auxiliary', kinds.values())

    def test_the_general_scheme_names_follow_hubert_labahn(self):
        translator = ODETranslator.from_ode_system(self.hl66)
        self.assertEqual({'t': 't', 'x0': 't*z1**2', 'y0': 't*z1**3', 'y1': 'z1**2*z2/t'},
                         {str(name): str(monomial) for name, monomial
                          in translator.named_symbols(self.hl66).items()})

    def test_the_scheme_may_be_said_rather_than_inferred(self):
        '''
        The names depend on the scheme, so a caller who will run `translate_general` on a
        system `translate` would have given to the parameter scheme must be able to say so.
        '''
        translator = ODETranslator.from_ode_system(self.prey_predator)
        inferred = translator.named_symbols(self.prey_predator)
        general = translator.named_symbols(self.prey_predator, scheme='general')
        self.assertEqual(('kappa0', 'kappa1', 'kappa2', 'tau', 'nu0', 'nu1'),
                         tuple(str(n) for n in inferred))
        self.assertEqual(('t', 'x0', 'x1', 'x2', 'y0', 'y1', 'y2', 'y3', 'y4', 'y5'),
                         tuple(str(n) for n in general))


class TestNamingScheme(TestCase):
    '''
    The naming scheme is a dict, and every reduction scheme honours it.

    Until 2.0 it was a positional triple, and the general scheme accepted it and threw it
    away -- a scheme passed to `translate_general` changed nothing at all, not even the
    independent variable's name.  Both halves are pinned here.
    '''

    def setUp(self):
        # Hubert & Labahn 6.6: the general scheme, with one auxiliary and two invariants.
        self.general_system = ODESystem.from_equations(
            ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
             'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])
        self.general_system.reorder_variables(['t', 'z1', 'z2'])
        # Michaelis-Menten: the parameter scheme.
        self.parameter_system = ODESystem.from_equations(
            ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
             'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
        self.parameter_system.reorder_variables(['e_0', 'k_1', 'k_2', 'k_m1', 't', 's', 'c'])

    def names(self, system):
        return tuple(map(str, system.symbols))

    def test_a_tuple_is_refused_and_says_what_to_write(self):
        with self.assertRaises(ValueError) as caught:
            ODETranslator.from_ode_system(self.parameter_system,
                                          naming_scheme=('tau', 'nu', 'kappa'))
        message = str(caught.exception)
        self.assertIn('dict', message)
        for key in ('parameters', 'independent', 'dependents'):
            self.assertIn(key, message)

    def test_an_unknown_key_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            ODETranslator.from_ode_system(self.parameter_system,
                                          naming_scheme={'indep': 'tau'})
        self.assertIn('indep', str(caught.exception))

    def test_a_list_where_a_pattern_is_needed_is_refused(self):
        for key in ('parameters', 'auxiliaries'):
            with self.assertRaises(ValueError):
                ODETranslator.from_ode_system(self.parameter_system,
                                              naming_scheme={key: ['a', 'b']})

    def test_parameter_scheme_defaults_are_unchanged(self):
        translator = ODETranslator.from_ode_system(self.parameter_system)
        reduced = translator.translate(self.parameter_system)
        self.assertEqual(('kappa0', 'kappa1', 'tau', 'nu0', 'nu1'), self.names(reduced))

    def test_general_scheme_defaults_follow_hubert_labahn(self):
        translator = ODETranslator.from_ode_system(self.general_system)
        reduced = translator.translate_general(self.general_system)
        self.assertEqual(('t', 'x0', 'y0', 'y1'), self.names(reduced))

    def test_general_scheme_honours_the_scheme(self):
        ''' The regression: this used to be discarded in full. '''
        translator = ODETranslator.from_ode_system(
            self.general_system,
            naming_scheme={'independent': 'T', 'dependents': ['alpha', 'beta'],
                           'auxiliaries': 'aux'})
        reduced = translator.translate_general(self.general_system)
        self.assertEqual(('T', 'aux0', 'alpha', 'beta'), self.names(reduced))

    def test_general_scheme_renames_only_the_auxiliaries_if_that_is_all_that_is_set(self):
        translator = ODETranslator.from_ode_system(self.general_system,
                                                   naming_scheme={'auxiliaries': 'w'})
        reduced = translator.translate_general(self.general_system)
        self.assertEqual(('t', 'w0', 'y0', 'y1'), self.names(reduced))

    def test_parameter_scheme_takes_a_partial_scheme(self):
        translator = ODETranslator.from_ode_system(self.parameter_system,
                                                   naming_scheme={'parameters': 'c'})
        reduced = translator.translate(self.parameter_system)
        self.assertEqual(('c0', 'c1', 'tau', 'nu0', 'nu1'), self.names(reduced))

    def test_one_scheme_serves_both_reductions(self):
        '''
        Each reduction falls back to its own defaults for what the scheme leaves unset, so
        the general one keeps y for its invariants while both take T and aux.
        '''
        shared = {'parameters': 'p', 'independent': 'T', 'auxiliaries': 'aux'}
        reduced = ODETranslator.from_ode_system(
            self.parameter_system, naming_scheme=shared).translate(self.parameter_system)
        self.assertEqual(('p0', 'p1', 'T', 'nu0', 'nu1'), self.names(reduced))
        reduced = ODETranslator.from_ode_system(
            self.general_system, naming_scheme=shared).translate_general(self.general_system)
        self.assertEqual(('T', 'aux0', 'y0', 'y1'), self.names(reduced))

    def test_a_wrong_length_dependents_list_is_refused_with_both_counts(self):
        translator = ODETranslator.from_ode_system(
            self.general_system, naming_scheme={'dependents': ['only_one']})
        with self.assertRaises(ValueError) as caught:
            translator.translate_general(self.general_system)
        message = str(caught.exception)
        self.assertIn('1 names', message)
        self.assertIn('2 of', message)

    def test_the_scheme_as_set_is_readable(self):
        translator = ODETranslator.from_ode_system(self.general_system,
                                                   naming_scheme={'auxiliaries': 'w'})
        self.assertEqual({'auxiliaries': 'w'}, translator.naming_scheme)

    def test_setting_no_scheme_leaves_it_empty(self):
        translator = ODETranslator.from_ode_system(self.general_system)
        self.assertEqual({}, translator.naming_scheme)


class TestSymbolsThatMustNotBeZero(TestCase):
    '''
    The scaling group acts where every coordinate is nonzero, so a reduction is only
    defined where the symbols it divides by are.  Which ones those are is a sign condition
    on one row of the multiplier, and it is a *report*: it follows from the columns that
    were chosen, and a column operation changes it.
    '''

    EQUATIONS = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                 'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']

    def system(self, parameter_order=None):
        system = ODESystem.from_equations(self.EQUATIONS)
        if parameter_order is not None:
            system.reorder_parameters(parameter_order)
        return system

    def test_the_symbols_scaled_by_are_the_ones_listed(self):
        system = self.system()
        translator = ODETranslator.from_ode_system(system)
        # The auxiliaries are 1/(e_0*k_1) and e_0, so those two are what it divides by.
        self.assertEqual(translator.symbols_that_must_not_be_zero(),
                         tuple(sympy.symbols('e_0 k_1')))

    def test_a_symbol_only_ever_in_a_numerator_is_not_listed(self):
        # The good case: such an invariant is simply 0, which is a reduced parameter whose
        # value is zero.  That is how a zero initial condition survives a reduction.
        system = self.system(['k_2', 'k_m1', 'k_1', 'e_0'])
        translator = ODETranslator.from_ode_system(system)
        self.assertNotIn(sympy.Symbol('e_0'), translator.symbols_that_must_not_be_zero())

    def test_every_listed_symbol_is_in_an_auxiliary_or_a_denominator(self):
        system = self.system()
        translator = ODETranslator.from_ode_system(system)
        listed = translator.symbols_that_must_not_be_zero()
        for row, symbol in enumerate(system.symbols):
            in_auxiliary = any(e != 0 for e in translator.herm_mult_i[row, :])
            in_denominator = any(e < 0 for e in translator.herm_mult_n[row, :])
            with self.subTest(symbol=symbol):
                self.assertEqual(symbol in listed, in_auxiliary or in_denominator)

    def test_a_column_operation_changes_the_answer(self):
        # Which is why it is a method and not a cached property, and why the guard cannot
        # sit at construction time: the choice is still being made.
        system = self.system()
        translator = ODETranslator.from_ode_system(system)
        before = translator.symbols_that_must_not_be_zero()
        # Column 2 is the first invariant column; an auxiliary column cannot be the source
        # of the operation, since adding it would change the Hermite form.
        translator.multiplier_add_columns(0, 2)
        self.assertNotEqual(before, translator.symbols_that_must_not_be_zero())

    def test_the_system_may_be_passed_instead_of_relying_on_the_domain(self):
        system = self.system()
        translator = ODETranslator.from_ode_system(system)
        self.assertEqual(translator.symbols_that_must_not_be_zero(system),
                         translator.symbols_that_must_not_be_zero())

    def test_the_order_is_the_systems_order(self):
        system = self.system()
        translator = ODETranslator.from_ode_system(system)
        listed = translator.symbols_that_must_not_be_zero()
        positions = [list(system.symbols).index(symbol) for symbol in listed]
        self.assertEqual(positions, sorted(positions))


class TestTranslatingWithAZeroValue(TestCase):
    '''
    A value of zero for a symbol the reduction divides by makes the reduced system
    undefined, and the reduction is where that can first be known: the translator does not
    hold the system, so it has no values, and column operations keep changing which
    symbols are involved.
    '''

    EQUATIONS = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                 'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']

    def system(self, parameter_order=None):
        system = ODESystem.from_equations(self.EQUATIONS)
        if parameter_order is not None:
            system.reorder_parameters(parameter_order)
        return system

    def test_a_zero_scale_is_refused_at_translate_time(self):
        system = self.system()
        system.update_parameter_values({'e_0': 0})
        translator = ODETranslator.from_ode_system(system)
        with self.assertRaises(ValueError) as caught:
            translator.translate(system)
        message = str(caught.exception)
        self.assertIn('e_0', message)
        self.assertIn('divides by', message)
        self.assertIn('reorder_parameters', message)

    def test_building_the_translator_is_not_refused(self):
        # Exploring orderings and column operations must stay free of the values; the
        # complaint belongs at the point a *system* is asked for.
        system = self.system()
        system.update_parameter_values({'e_0': 0})
        translator = ODETranslator.from_ode_system(system)
        self.assertIn(sympy.Symbol('e_0'), translator.symbols_that_must_not_be_zero())

    def test_a_zero_value_on_a_symbol_it_does_not_divide_by_is_fine(self):
        system = self.system(['k_2', 'k_m1', 'k_1', 'e_0'])
        system.update_parameter_values({'e_0': 0})
        translator = ODETranslator.from_ode_system(system)
        reduced = translator.translate(system)
        self.assertTrue(reduced.is_reduced)

    def test_the_zero_reaches_the_reduced_system_as_a_parameter_value(self):
        # e_0 appears only in the numerator of one invariant, so that reduced parameter is
        # zero.  The value is derived, not restated.
        system = self.system(['k_2', 'k_m1', 'k_1', 'e_0'])
        system.update_parameter_values({'e_0': 0})
        translator = ODETranslator.from_ode_system(system)
        reduced = translator.translate(system)
        self.assertIn(sympy.Integer(0), set(reduced.parameter_values.values()))

    def test_a_nonzero_value_is_never_refused(self):
        system = self.system()
        system.update_parameter_values({'e_0': 2})
        translator = ODETranslator.from_ode_system(system)
        self.assertTrue(translator.translate(system).is_reduced)

    def test_a_system_with_no_values_behaves_exactly_as_before(self):
        system = self.system()
        translator = ODETranslator.from_ode_system(system)
        self.assertTrue(translator.translate(system).is_reduced)


class TestBothDirections(TestCase):
    '''
    The translator is where both sets of symbols are visible, so it is where both
    directions can be asked about.  Forward, a reduced symbol is a monomial in the
    originals, so only an original can be in a denominator; in reverse an original is a
    monomial in the reduced ones, so only a reduced symbol can.  One symbol set per
    direction, which is why neither method needs a `direction` argument.
    '''

    EQUATIONS = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                 'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']

    def setUp(self):
        self.system = ODESystem.from_equations(self.EQUATIONS)
        self.translator = ODETranslator.from_ode_system(self.system)

    def test_the_two_directions_give_different_answers(self):
        self.assertNotEqual(
            set(self.translator.symbols_that_must_not_be_zero()),
            set(self.translator.reduced_symbols_that_must_not_be_zero(self.system)))

    def test_the_reverse_answer_is_over_the_reduced_symbols(self):
        reverse = self.translator.reduced_symbols_that_must_not_be_zero(self.system)
        self.assertEqual(reverse, tuple(sympy.symbols('x1 kappa0 kappa1')))

    def test_the_forward_answer_is_over_the_original_symbols(self):
        forward = self.translator.symbols_that_must_not_be_zero()
        self.assertTrue(set(forward) <= set(self.system.symbols))

    def test_no_original_symbol_appears_in_the_reverse_answer(self):
        reverse = self.translator.reduced_symbols_that_must_not_be_zero(self.system)
        self.assertEqual(set(reverse) & set(self.system.symbols), set())

    def test_the_reverse_answer_is_exactly_the_denominators_of_the_recovery(self):
        # Each original symbol is recovered as a monomial in the reduced ones; collect the
        # reduced symbols that show up in a denominator and the two must agree.
        names = self.translator.column_names(self.system)
        recovered = scale_action(names, self.translator.inv_herm_mult.T)
        in_a_denominator = set()
        for expression in recovered:
            for symbol, power in expression.as_powers_dict().items():
                if power < 0:
                    in_a_denominator.add(symbol)
        self.assertEqual(
            set(self.translator.reduced_symbols_that_must_not_be_zero(self.system)),
            in_a_denominator)

    def test_the_reverse_answer_reads_against_the_column_names(self):
        names = list(self.translator.column_names(self.system))
        reverse = self.translator.reduced_symbols_that_must_not_be_zero(self.system)
        positions = [names.index(name) for name in reverse]
        self.assertEqual(positions, sorted(positions))

    def test_a_reduction_can_be_sound_forward_and_unrecoverable_backward(self):
        # k_2 = 0 makes kappa0 zero.  The forward reduction never divides by k_2, so it
        # goes through; but e_0 = tau/(kappa0*kappa1*x1) cannot then be recovered.
        self.system.update_parameter_values({'k_2': 0})
        reduced = self.translator.translate(self.system)
        self.assertNotIn(sympy.Symbol('k_2'),
                         self.translator.symbols_that_must_not_be_zero())
        self.assertIn(sympy.Symbol('kappa0'),
                      self.translator.reduced_symbols_that_must_not_be_zero(self.system))
        self.assertEqual(reduced.parameter_values[sympy.Symbol('kappa0')],
                         sympy.Integer(0))


class TestParameterValuesThroughTranslation(TestCase):
    ''' A reduced parameter is an invariant, so its value is derived, not restated. '''

    EQUATIONS = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                 'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']

    def system(self, parameter_order=None):
        system = ODESystem.from_equations(self.EQUATIONS)
        if parameter_order is not None:
            system.reorder_parameters(parameter_order)
        return system

    def test_a_fully_known_invariant_gets_its_value(self):
        system = self.system()
        # kappa0 = k_2/(e_0*k_1), so all three are needed before it is known.
        system.update_parameter_values({'k_2': 6, 'e_0': 2, 'k_1': 3})
        reduced = ODETranslator.from_ode_system(system).translate(system)
        self.assertEqual(reduced.parameter_values[sympy.Symbol('kappa0')],
                         sympy.Integer(1))

    def test_a_partly_known_invariant_gets_none(self):
        system = self.system()
        system.update_parameter_values({'k_2': 6})
        reduced = ODETranslator.from_ode_system(system).translate(system)
        self.assertNotIn(sympy.Symbol('kappa0'), reduced.parameter_values)

    def test_no_values_means_no_values(self):
        system = self.system()
        reduced = ODETranslator.from_ode_system(system).translate(system)
        self.assertEqual(reduced.parameter_values, {})

    def test_only_parameters_of_the_reduced_system_are_given_values(self):
        system = self.system()
        system.update_parameter_values({'k_2': 6, 'e_0': 2, 'k_1': 3, 'k_m1': 3})
        reduced = ODETranslator.from_ode_system(system).translate(system)
        self.assertTrue(set(reduced.parameter_values) <= set(reduced.parameters))

    def test_the_general_scheme_propagates_too(self):
        # Called directly, since `translate` would pick the parameter scheme here.  Every
        # symbol of a reduced system stands for a monomial whatever the scheme, so there is
        # nothing scheme-specific about knowing its value.
        system = self.system()
        system.update_parameter_values({'k_1': 2, 'k_2': 3, 'k_m1': 4, 'e_0': 5})
        reduced = ODETranslator.from_ode_system(system).translate_general(system)
        # y0 = k_2/(e_0*k_1) = 3/10
        self.assertEqual(reduced.parameter_values[sympy.Symbol('y0')],
                         sympy.Rational(3, 10))

    def test_the_general_scheme_values_its_auxiliaries(self):
        # An auxiliary is free as far as the reduced *dynamics* go, but it still denotes a
        # monomial, so these values are the ones this particular original system has.
        system = self.system()
        system.update_parameter_values({'k_1': 2, 'k_2': 3, 'k_m1': 4, 'e_0': 5})
        reduced = ODETranslator.from_ode_system(system).translate_general(system)
        self.assertEqual(reduced.parameter_values[sympy.Symbol('x1')],   # x1 = e_0
                         sympy.Integer(5))

    def test_dropping_the_auxiliaries_drops_their_values(self):
        system = self.system()
        system.update_parameter_values({'k_1': 2, 'k_2': 3, 'k_m1': 4, 'e_0': 5})
        reduced = ODETranslator.from_ode_system(system).translate_general(
            system, include_aux_vars=False)
        self.assertNotIn(sympy.Symbol('x1'), reduced.parameter_values)
        self.assertIn(sympy.Symbol('y0'), reduced.parameter_values)

    def test_the_dependent_variable_scheme_propagates_too(self):
        system = ODESystem.from_equations(['dz1/dt = a*z1*(1+z1*z2)',
                                           'dz2/dt = z2*(1/t - z1*z2)'])
        system.update_parameter_values({'a': 7})
        translator = ODETranslator.from_ode_system(system)
        self.assertEqual(translator.scheme(system), ReductionScheme.DEP_VAR)
        reduced = translator.translate_dep_var(system)
        self.assertEqual(reduced.parameter_values[sympy.Symbol('y0')],
                         sympy.Integer(7))

    def test_every_scheme_values_exactly_its_known_parameters(self):
        # The invariant across all three: a reduced parameter has a value exactly when the
        # monomial it stands for does, and nothing else is given one.
        system = self.system()
        system.update_parameter_values({'k_1': 2, 'k_2': 3, 'k_m1': 4, 'e_0': 5})
        translator = ODETranslator.from_ode_system(system)
        for name, reduce in [(ReductionScheme.PARAMETER,
                              translator.translate_parameter),
                             (ReductionScheme.GENERAL,
                              translator.translate_general)]:
            with self.subTest(scheme=name):
                reduced = reduce(system)
                named = translator.named_symbols(system, name)
                expected = {symbol for symbol in reduced.parameters
                            if system.value_of(named[symbol]) is not None}
                self.assertEqual(set(reduced.parameter_values), expected)


class TestGeneratedNameClashes(TestCase):
    '''
    The original and the reduced system are different namespaces, and `desr.numerics` maps
    between them, so one spelling meaning two quantities is not cosmetic: a `known_values`
    entry would not say which was meant.  A generated name landing on a symbol the original
    already has is refused.
    '''

    def test_a_parameter_already_called_kappa0_is_refused(self):
        # It used to reduce quietly, and the reduced `kappa0` was an invariant monomial in
        # the original `c_1` and `kappa0` -- same spelling, different quantity.
        system = ODESystem.from_equations(['dx/dt = kappa0*x*y',
                                           'dy/dt = c_1*(1-x)*(1-y)'])
        translator = ODETranslator.from_ode_system(system)
        with self.assertRaises(ValueError) as caught:
            translator.translate(system)
        message = str(caught.exception)
        self.assertIn('kappa0', message)
        self.assertIn('naming_scheme', message)

    def test_a_naming_scheme_colliding_with_a_minted_name_is_refused(self):
        # The route worth worrying about: an initial condition mints `c_0` from the
        # variable `c`, and then a naming scheme of 'c_' generates `c_0` as well.  It
        # needs enough parameters that one survives the reduction to be named.
        system = ODESystem.from_equations(['dc/dt = -k_1*c*s + k_2*p',
                                           'ds/dt = k_1*c*s - k_3*s',
                                           'dp/dt = k_3*s - k_2*p'])
        system.update_initial_conditions('c(0) = c_0')
        system.reorder_parameters(['k_1', 'k_2', 'k_3', 'c_0'])
        translator = ODETranslator.from_ode_system(system)
        with self.assertRaises(ValueError) as caught:
            translator.translate(system, naming_scheme={'parameters': 'c_'})
        self.assertIn('c_0', str(caught.exception))

    def test_a_scheme_whose_names_are_all_normalised_out_cannot_clash(self):
        # Two parameters and two scaling symmetries, so no reduced parameter is generated
        # at all and there is nothing to collide with.
        system = ODESystem.from_equations(['dc/dt = -k_1*c*s', 'ds/dt = k_1*c*s'])
        system.update_initial_conditions('c(0) = c_0')
        translator = ODETranslator.from_ode_system(system)
        reduced = translator.translate(system, naming_scheme={'parameters': 'c_'})
        self.assertEqual(set(system.symbols) & set(reduced.symbols), set())

    def test_a_scheme_that_avoids_the_clash_is_accepted(self):
        system = ODESystem.from_equations(['dx/dt = kappa0*x*y',
                                           'dy/dt = c_1*(1-x)*(1-y)'])
        translator = ODETranslator.from_ode_system(system)
        reduced = translator.translate(system, naming_scheme={'parameters': 'theta'})
        self.assertTrue(reduced.is_reduced)

    def test_keeping_the_systems_own_clock_is_not_a_clash(self):
        # The general and dependent-variable schemes retain the independent variable on
        # purpose, so that one shared name is allowed.
        equations = 'dz1/dt = z1*(1+z1*z2);dz2/dt = z2*(1/t - z1*z2)'.split(';')
        system = ODESystem.from_equations(equations)
        reduced = ODETranslator.from_ode_system(system).translate(system)
        self.assertEqual(reduced.independent_variable, system.independent_variable)

    def test_an_ordinary_reduction_is_unaffected(self):
        equations = 'dn/dt = n*( r*(1 - n/K) - k*p/(n+d) );dp/dt = s*p*(1 - h*p / n)'
        system = ODESystem.from_equations(equations.split(';'))
        reduced = ODETranslator.from_ode_system(system).translate(system)
        self.assertEqual(set(system.symbols) & set(reduced.symbols), set())


class TestExplicitNamesAreTheCallersChoice(TestCase):
    '''
    Only a *generated* name can clash by accident.  Reusing the original letters on purpose
    is a real thing to want: Hubert & Labahn write their reduced predator-prey system in
    the same `n` and `p`.
    '''

    EQUATIONS = ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )', 'dp/dt = s*p*(1 - h*p / n)']

    def system(self):
        system = ODESystem.from_equations(self.EQUATIONS)
        system.reorder_variables(['K', 'h', 'r', 'd', 'k', 's', 't', 'n', 'p'])
        return system

    def test_explicitly_reusing_the_original_dependent_names_is_allowed(self):
        system = self.system()
        translator = ODETranslator.from_ode_system(
            system, naming_scheme={'parameters': 'c', 'independent': 't',
                                   'dependents': ['n', 'p']})
        reduced = translator.translate_parameter(system)
        self.assertEqual(set(map(str, reduced.dependents)), {'n', 'p'})

    def test_a_generated_name_landing_on_an_original_is_still_refused(self):
        # Same system, but `n` arrives from the pattern rather than from a list.
        system = ODESystem.from_equations(['dx/dt = n0*x*y', 'dy/dt = c_1*(1-x)*(1-y)'])
        translator = ODETranslator.from_ode_system(system,
                                                   naming_scheme={'parameters': 'n'})
        with self.assertRaises(ValueError) as caught:
            translator.translate(system)
        self.assertIn('n0', str(caught.exception))

    def test_an_explicit_independent_variable_name_is_allowed(self):
        system = self.system()
        translator = ODETranslator.from_ode_system(
            system, naming_scheme={'parameters': 'c', 'independent': 't',
                                   'dependents': ['n', 'p']})
        reduced = translator.translate_parameter(system)
        self.assertEqual(str(reduced.independent_variable), 't')


class TestExtendFromWrittenInvariants(TestCase):
    '''
    An invariant is said by writing it.  The exponent matrix still works, but its rows are
    positions in the variable order, so the same matrix means different monomials after a
    reorder -- and nothing in it records which invariants were intended.
    '''

    EQUATIONS = 'dn/dt = n*( r*(1 - n/K) - k*p/(n+d) );dp/dt = s*p*(1 - h*p / n)'.split(';')

    def translator(self, parameter_order=None):
        system = ODESystem.from_equations(self.EQUATIONS)
        if parameter_order is not None:
            system.reorder_parameters(parameter_order)
        return system, ODETranslator.from_ode_system(system)

    def test_a_list_of_expressions(self):
        system, translator = self.translator()
        h, p, d, r, t = sympy.symbols('h p d r t')
        extended = translator.extend_from_invariants([r * t, h * p / d])
        self.assertEqual(tuple(extended.invariant_monomials())[:2], (r * t, h * p / d))

    def test_a_list_of_strings(self):
        system, translator = self.translator()
        extended = translator.extend_from_invariants(['r*t', 'h*p/d'])
        self.assertEqual(tuple(extended.invariant_monomials())[:2],
                         (sympy.sympify('r*t'), sympy.sympify('h*p/d')))

    def test_one_expression_need_not_be_wrapped(self):
        system, translator = self.translator()
        extended = translator.extend_from_invariants(sympy.sympify('h*p/d'))
        self.assertEqual(tuple(extended.invariant_monomials())[0],
                         sympy.sympify('h*p/d'))

    def test_one_string_need_not_be_wrapped(self):
        system, translator = self.translator()
        extended = translator.extend_from_invariants('h*p/d')
        self.assertEqual(tuple(extended.invariant_monomials())[0],
                         sympy.sympify('h*p/d'))

    def test_the_exponent_matrix_still_works(self):
        system, translator = self.translator()
        matrix = sympy.Matrix([monomial_to_powers(sympy.sympify('h*p/d'),
                                                  list(system.symbols))]).T
        extended = translator.extend_from_invariants(matrix)
        self.assertEqual(tuple(extended.invariant_monomials())[0],
                         sympy.sympify('h*p/d'))

    def test_the_two_forms_agree(self):
        system, translator = self.translator()
        matrix = sympy.Matrix([monomial_to_powers(sympy.sympify('h*p/d'),
                                                  list(system.symbols))]).T
        self.assertEqual(
            tuple(translator.extend_from_invariants(matrix).invariant_monomials()),
            tuple(translator.extend_from_invariants('h*p/d').invariant_monomials()))

    def test_a_written_invariant_survives_a_reorder(self):
        # The whole point.  The same thing written gives the same thing chosen.
        for order in (['K', 'd', 'h', 'k', 'r', 's'], ['s', 'r', 'k', 'h', 'd', 'K']):
            with self.subTest(order=''.join(order)):
                system, translator = self.translator(order)
                extended = translator.extend_from_invariants('h*p/d')
                self.assertEqual(tuple(extended.invariant_monomials())[0],
                                 sympy.sympify('h*p/d'))

    def test_the_exponent_matrix_does_not_survive_a_reorder(self):
        # Why the monomial form exists: built against one order, the identical matrix
        # denotes `k*p/r` against another.
        system, translator = self.translator(['K', 'd', 'h', 'k', 'r', 's'])
        matrix = sympy.Matrix([monomial_to_powers(sympy.sympify('h*p/d'),
                                                  list(system.symbols))]).T
        reordered_system, reordered = self.translator(['s', 'r', 'k', 'h', 'd', 'K'])
        self.assertEqual(
            scale_action(reordered_system.symbols, matrix)[0],
            sympy.sympify('k*p/r'))

    def test_something_that_is_not_invariant_is_named(self):
        system, translator = self.translator()
        with self.assertRaises(ValueError) as caught:
            translator.extend_from_invariants('h*p')
        message = str(caught.exception)
        self.assertIn('h*p', message)
        self.assertIn('not an invariant', message)

    def test_a_sum_is_refused_as_not_a_monomial(self):
        system, translator = self.translator()
        with self.assertRaises(ValueError) as caught:
            translator.extend_from_invariants('h*p/d + r*t')
        self.assertIn('not a monomial', str(caught.exception))

    def test_an_unknown_symbol_is_refused_rather_than_ignored(self):
        # monomial_to_powers looks up only the variables it is given, so a stranger
        # contributes nothing and would silently denote a different invariant.
        system, translator = self.translator()
        with self.assertRaises(ValueError) as caught:
            translator.extend_from_invariants('h*p/dd')
        message = str(caught.exception)
        self.assertIn('dd', message)
        self.assertIn('does not act on', message)

    def test_a_fractional_power_is_refused_with_a_reason(self):
        system, translator = self.translator()
        with self.assertRaises(ValueError) as caught:
            translator.extend_from_invariants('sqrt(h*p/d)')
        self.assertIn('fractional power', str(caught.exception))

    def test_an_empty_choice_is_refused(self):
        system, translator = self.translator()
        with self.assertRaises(ValueError) as caught:
            translator.extend_from_invariants([])
        self.assertIn('No invariants', str(caught.exception))

    def test_an_explicit_variable_order_may_be_given(self):
        # For a translator built without a system, where symbols_domain is unset.
        translator = ODETranslator(sympy.Matrix([[1, 0, 3, 0, 2, 2],
                                                 [0, 2, 0, 1, 0, 1],
                                                 [2, 0, 0, 3, 0, 0]]))
        variables = sympy.symbols(' '.join('y{}'.format(i) for i in range(6)))
        extended = translator.extend_from_invariants('y4**3/y2**2', variables=variables)
        self.assertEqual(tuple(extended.invariant_monomials(variables=variables))[0],
                         sympy.sympify('y4**3/y2**2'))


class TestPassingAPrebuiltExpression(TestCase):
    '''
    Building the expression first and passing it in works.  The one way it quietly does
    not is `sympy.sympify`, which reads `E`, `I`, `S` and friends as sympy's own objects --
    and the result prints exactly like the right one.
    '''

    EQUATIONS = ['dS/dt = -k_1*E*S + k_m1*C',
                 'dE/dt = -k_1*E*S + k_m1*C + k_2*C',
                 'dC/dt = k_1*E*S - k_m1*C - k_2*C']

    def setUp(self):
        self.system = ODESystem.from_equations(self.EQUATIONS)
        self.translator = ODETranslator.from_ode_system(self.system)

    def test_an_expression_built_with_sympify_expression_works(self):
        chosen = sympify_expression('E*k_1/k_2')
        extended = self.translator.extend_from_invariants(chosen)
        self.assertEqual(tuple(extended.invariant_monomials())[0], chosen)

    def test_an_expression_built_from_symbols_works(self):
        E, k_1, k_2 = sympy.symbols('E k_1 k_2')
        extended = self.translator.extend_from_invariants(E * k_1 / k_2)
        self.assertEqual(tuple(extended.invariant_monomials())[0], E * k_1 / k_2)

    def test_sympys_own_sympify_is_caught_and_the_cause_named(self):
        # `sympy.sympify('E*k_1/k_2')` has Euler's number where the enzyme should be, and
        # prints identically.  The only tell is that `E` is absent from free_symbols.
        wrong = sympy.sympify('E*k_1/k_2')
        self.assertNotIn(sympy.Symbol('E'), wrong.free_symbols)
        with self.assertRaises(ValueError) as caught:
            self.translator.extend_from_invariants(wrong)
        message = str(caught.exception)
        self.assertIn('sympify_expression', message)
        self.assertIn('print the same', message)

    def test_the_two_sympifications_print_the_same(self):
        # Which is why the error has to say what is wrong rather than show the expression.
        self.assertEqual(str(sympy.sympify('E*k_1/k_2')),
                         str(sympify_expression('E*k_1/k_2')))

    def test_a_numeric_factor_is_refused_separately(self):
        with self.assertRaises(ValueError) as caught:
            self.translator.extend_from_invariants(sympify_expression('2*E*k_1/k_2'))
        message = str(caught.exception)
        self.assertIn('without its factor of 2', message)
        self.assertNotIn('sympify_expression', message)

    def test_a_rational_factor_in_a_denominator_is_refused_too(self):
        with self.assertRaises(ValueError) as caught:
            self.translator.extend_from_invariants(sympify_expression('E*k_1/(2*k_2)'))
        self.assertIn('without its factor of', str(caught.exception))
