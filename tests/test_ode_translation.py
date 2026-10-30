from unittest import TestCase

import sympy

from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.ode_system import ODESystem
from desr.ode_translation import ODETranslator
from desr.sympy_helper import sympify_expression


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
