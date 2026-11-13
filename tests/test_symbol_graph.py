'''
The exponent graph.

The graph is a reading of the Hermite multiplier, so most of these compare it against the
multiplier and against the methods that answer the same questions in other words.  The
pair that matter most are the two that would have caught the orientation bug in
`reduced_symbols_that_must_not_be_zero`: the edges must agree with the matrix entry by
entry, and composing the two directions must give the identity.
'''
from unittest import TestCase

import matplotlib
matplotlib.use('Agg')

import networkx
import sympy

from desr import ODESystem, ODETranslator, ReductionScheme
from desr.chemical_reaction_network import ChemicalReactionNetwork
from desr.symbol_graph import (ORIGINAL, REDUCED, _node_label, _positions,
                               draw_symbol_graph, symbol_graph)


MICHAELIS_MENTEN = ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
                    'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c']


def _edges(graph):
    '''
    The edges in a comparable form.  Sorted by string: a node key holds a sympy symbol,
    and sympy answers `<` between two symbols with an unevaluated relational rather than a
    bool, which `sorted` cannot use.
    '''
    return sorted(((str(tail), str(head), exponent)
                   for tail, head, exponent in graph.edges(data='exponent')))


def michaelis_menten_with_initial_conditions():
    '''The two-variable form carrying both initial conditions, one of them zero.'''
    system = ODESystem.from_equations(MICHAELIS_MENTEN)
    system.add_constraint('K_m', '(k_m1 + k_2) / k_1')
    system.update_initial_conditions(['s(0) = s_0', 'c(0) = 0'])
    system.reorder_variables(['s_0', 'e_0', 'k_1', 'k_2', 'k_m1', 'K_m', 'c_0',
                              't', 's', 'c'])
    return system


class TestTheGraphIsTheMultiplier(TestCase):
    ''' Edge for entry, in both directions. '''

    def setUp(self):
        self.system = ODESystem.from_equations(MICHAELIS_MENTEN)
        self.translator = ODETranslator.from_ode_system(self.system)

    def test_both_layers_have_one_vertex_per_symbol(self):
        graph = symbol_graph(self.translator, self.system)
        n = self.translator.n
        self.assertEqual(graph.number_of_nodes(), 2 * n)
        for layer in (ORIGINAL, REDUCED):
            with self.subTest(layer=layer):
                in_layer = [node for node, data in graph.nodes(data=True)
                            if data['layer'] == layer]
                self.assertEqual(len(in_layer), n)

    def test_the_forward_edges_are_the_multiplier_entry_for_entry(self):
        graph = symbol_graph(self.translator, self.system)
        originals = list(self.translator.symbols_domain)
        names = list(self.translator.column_names(self.system))
        multiplier = self.translator.herm_mult

        expected = {}
        for row, original in enumerate(originals):
            for column, name in enumerate(names):
                if multiplier[row, column] != 0:
                    expected[((ORIGINAL, original), (REDUCED, name))] = \
                        int(multiplier[row, column])

        self.assertEqual(
            {(tail, head): exponent
             for tail, head, exponent in graph.edges(data='exponent')},
            expected)

    def test_the_reverse_edges_are_the_inverse_entry_for_entry(self):
        graph = symbol_graph(self.translator, self.system, direction='reverse')
        originals = list(self.translator.symbols_domain)
        names = list(self.translator.column_names(self.system))
        inverse = self.translator.inv_herm_mult

        # The inverse is indexed [reduced, original], so the reduced symbol is the row.
        expected = {}
        for row, name in enumerate(names):
            for column, original in enumerate(originals):
                if inverse[row, column] != 0:
                    expected[((REDUCED, name), (ORIGINAL, original))] = \
                        int(inverse[row, column])

        self.assertEqual(
            {(tail, head): exponent
             for tail, head, exponent in graph.edges(data='exponent')},
            expected)

    def test_the_two_directions_compose_to_the_identity(self):
        # Two hops: the exponents multiply along a path and the parallel paths sum.  This
        # is the assertion that catches a transposed read of either matrix, where a
        # statement about which symbols divide cannot.
        forward = symbol_graph(self.translator, self.system)
        reverse = symbol_graph(self.translator, self.system, direction='reverse')

        for start in self.translator.symbols_domain:
            for finish in self.translator.symbols_domain:
                total = sum(
                    forward[(ORIGINAL, start)][middle]['exponent']
                    * reverse[middle][(ORIGINAL, finish)]['exponent']
                    for middle in forward[(ORIGINAL, start)]
                    if (ORIGINAL, finish) in reverse[middle])
                with self.subTest(start=start, finish=finish):
                    self.assertEqual(total, 1 if start == finish else 0)

    def test_the_graph_is_bipartite_in_the_networkx_sense(self):
        graph = symbol_graph(self.translator, self.system)
        undirected = graph.to_undirected()
        self.assertTrue(networkx.algorithms.bipartite.is_bipartite(undirected))
        left, right = networkx.algorithms.bipartite.sets(undirected)
        self.assertEqual({len(left), len(right)}, {self.translator.n})

    def test_a_column_operation_changes_the_graph(self):
        # Read from the multiplier as it stands, which is what makes the graph useful for
        # choosing a column operation rather than only for reporting one.
        before = symbol_graph(self.translator, self.system)
        self.translator.multiplier_add_columns(
            self.translator.column_of('k_2/(e_0*k_1)'),
            self.translator.column_of('k_m1/(e_0*k_1)'))
        after = symbol_graph(self.translator, self.system)
        self.assertNotEqual(_edges(before), _edges(after))

    def test_an_unknown_direction_is_refused_and_names_the_two(self):
        with self.assertRaises(ValueError) as caught:
            symbol_graph(self.translator, self.system, direction='sideways')
        self.assertIn("'forward'", str(caught.exception))
        self.assertIn("'reverse'", str(caught.exception))


class TestTheTwoLayersDoNotMerge(TestCase):
    '''
    A reduced symbol can *be* an original symbol -- a naming scheme may reuse the original
    letters, and two of the three schemes keep the system's own clock.  Bare symbols as
    nodes would silently collapse those vertices and quietly change the graph.
    '''

    def setUp(self):
        self.system = ODESystem.from_equations(MICHAELIS_MENTEN)
        self.translator = ODETranslator.from_ode_system(
            self.system,
            naming_scheme={'parameters': 'kappa', 'independent': 't',
                           'dependents': ['s', 'c']})

    def test_the_two_layers_really_do_share_symbols(self):
        shared = set(self.translator.column_names(self.system)) & set(self.system.symbols)
        self.assertEqual(shared, {sympy.Symbol('t'), sympy.Symbol('s'),
                                  sympy.Symbol('c')})

    def test_the_shared_symbols_stay_four_distinct_vertices(self):
        graph = symbol_graph(self.translator, self.system)
        self.assertEqual(graph.number_of_nodes(), 2 * self.translator.n)
        for name in ('t', 's', 'c'):
            with self.subTest(symbol=name):
                symbol = sympy.Symbol(name)
                self.assertIn((ORIGINAL, symbol), graph)
                self.assertIn((REDUCED, symbol), graph)
                self.assertEqual(graph.nodes[(ORIGINAL, symbol)]['layer'], ORIGINAL)
                self.assertEqual(graph.nodes[(REDUCED, symbol)]['layer'], REDUCED)


class TestWhatTheNodesCarry(TestCase):

    def setUp(self):
        self.system = michaelis_menten_with_initial_conditions()
        self.translator = ODETranslator.from_ode_system(
            self.system,
            naming_scheme={'parameters': 'kappa', 'independent': 'tau',
                           'dependents': ['u', 'v']})
        self.graph = symbol_graph(self.translator, self.system)

    def reduced(self, name):
        return self.graph.nodes[(REDUCED, sympy.Symbol(name))]

    def original(self, name):
        return self.graph.nodes[(ORIGINAL, sympy.Symbol(name))]

    def test_an_original_symbol_carries_what_it_equals_coming_back(self):
        # Symmetric with the reduced side: every node says what it is in terms of the
        # other layer, so whichever direction is drawn, the arrows point at nodes that
        # name themselves.  Compared against the public method rather than against a
        # second reading of the inverse multiplier.
        names = self.translator.column_names(self.system)
        recoveries = self.translator.reverse_translate_parameter(names)
        for symbol, recovery in zip(self.translator.symbols_domain, recoveries):
            with self.subTest(symbol=str(symbol)):
                self.assertEqual(self.graph.nodes[(ORIGINAL, symbol)]['monomial'],
                                 recovery)

    def test_the_round_trip_is_visible_in_the_two_monomials(self):
        # Substituting the forward monomials into an original symbol's own expression
        # returns that symbol.  The graph carries both halves, so this is checkable on the
        # node data alone.
        substitution = {self.graph.nodes[(REDUCED, name)]['symbol']:
                        self.graph.nodes[(REDUCED, name)]['monomial']
                        for name in self.translator.column_names(self.system)}
        for symbol in self.translator.symbols_domain:
            with self.subTest(symbol=str(symbol)):
                coming_back = self.graph.nodes[(ORIGINAL, symbol)]['monomial']
                self.assertEqual(
                    sympy.simplify(coming_back.subs(substitution) - symbol), 0)

    def test_the_monomials_are_the_columns(self):
        names = list(self.translator.column_names(self.system))
        for name, monomial in zip(names, self.translator.column_monomials()):
            with self.subTest(column=name):
                self.assertEqual(self.graph.nodes[(REDUCED, name)]['monomial'], monomial)

    def test_the_first_r_columns_are_the_auxiliaries(self):
        names = list(self.translator.column_names(self.system))
        flags = [self.graph.nodes[(REDUCED, name)]['auxiliary'] for name in names]
        self.assertEqual(flags, [True] * self.translator.r
                         + [False] * (self.translator.n - self.translator.r))

    def test_being_an_auxiliary_is_not_one_of_the_kinds(self):
        # It is orthogonal to the classification, so it never displaces it.  An auxiliary
        # is a parameter or a dependent like any other symbol.
        for _, data in self.graph.nodes(data=True):
            with self.subTest(symbol=str(data['symbol'])):
                self.assertIn(data['kind'],
                              ('parameter', 'dependent', 'independent'))

    def test_a_rescaled_clock_is_a_variable_not_the_clock(self):
        # The general scheme keeps the system's own `t` and integrates against it, so an
        # invariant that is a rescaled clock -- the enzyme chain's `k_2_3*t` -- is a
        # variable with an equation, not the thing being integrated against.  Calling it
        # independent would say the reduced system has two independent variables.  This is
        # the role, which `symbol_kinds` gives; the *structural* classification says
        # 'independent' here and is a different question.
        system = ChemicalReactionNetwork.from_diagram(
            'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
        translator = ODETranslator.from_ode_system(system)
        graph = symbol_graph(translator, system)

        y_0 = graph.nodes[(REDUCED, sympy.Symbol('y0'))]
        self.assertEqual(y_0['monomial'], sympy.sympify('k_2_3*t'))
        self.assertEqual(y_0['kind'], 'dependent')

        # And the reduced system agrees that there is one clock.  It agrees more strongly
        # than it used to: `y0` is not a second dependent any more, because the reduction
        # identifies it *with* `t` -- so the reduced system's own clock is what `y0`
        # denotes, the rescaled clock, exactly as the parameter scheme's `tau` is.
        reduced = translator.translate(system)
        self.assertEqual(reduced.independent_variable, sympy.Symbol('t'))
        self.assertNotIn(sympy.Symbol('y0'), reduced.symbols)
        self.assertEqual(sympy.sympify('k_2_3*t'),
                         reduced.named_symbols[reduced.independent_variable])

    def test_exactly_one_vertex_per_layer_is_the_independent_variable(self):
        # The invariant the previous test protects: a layer with two independent variables
        # is a figure that cannot be true of any system.
        for label, system in (('enzyme', ChemicalReactionNetwork.from_diagram(
                                   'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()),
                              ('michaelis-menten', ODESystem.from_equations(MICHAELIS_MENTEN))):
            translator = ODETranslator.from_ode_system(system)
            graph = symbol_graph(translator, system)
            for layer in (ORIGINAL, REDUCED):
                with self.subTest(system=label, layer=layer):
                    clocks = [data for (side, _), data in graph.nodes(data=True)
                              if side == layer and data['kind'] == 'independent']
                    self.assertLessEqual(len(clocks), 1)

    def test_the_original_kinds_are_the_systems_own(self):
        self.assertEqual(self.original('t')['kind'], 'independent')
        self.assertEqual(self.original('s')['kind'], 'dependent')
        self.assertEqual(self.original('c')['kind'], 'dependent')
        self.assertEqual(self.original('k_1')['kind'], 'parameter')

    def test_a_known_zero_is_carried_and_propagated(self):
        # The whole point of the annotation: c_0 is zero because the complex starts at
        # nothing, and kappa4 is c_0/s_0, so it is zero for every s_0.  One substitution
        # settles it; nothing had to be translated first.
        self.assertEqual(self.original('c_0')['value'], 0)
        self.assertEqual(self.reduced('kappa4')['value'], 0)

    def test_the_propagated_value_agrees_with_the_reduced_system(self):
        reduced_system = self.translator.translate(self.system)
        for symbol, value in reduced_system.parameter_values.items():
            with self.subTest(symbol=symbol):
                self.assertEqual(self.graph.nodes[(REDUCED, symbol)]['value'], value)

    def test_a_symbol_with_no_known_value_carries_none(self):
        self.assertNotIn('value', self.original('k_1'))
        self.assertNotIn('value', self.reduced('kappa0'))

    def test_must_not_be_zero_is_the_forward_answer_on_the_original_layer(self):
        named = {symbol for (layer, symbol), data in self.graph.nodes(data=True)
                 if layer == ORIGINAL and data['must_not_be_zero']}
        self.assertEqual(named,
                         set(self.translator.symbols_that_must_not_be_zero(self.system)))

    def test_nothing_on_the_reduced_layer_is_flagged_going_forward(self):
        # Going forward only an original symbol can sit in a denominator, so the reduced
        # layer has nothing to say.  The reverse graph is where it does.
        flagged = [data for (layer, _), data in self.graph.nodes(data=True)
                   if layer == REDUCED and data['must_not_be_zero']]
        self.assertEqual(flagged, [])

    def test_must_not_be_zero_is_the_reverse_answer_coming_back(self):
        back = symbol_graph(self.translator, self.system, direction='reverse')
        named = {symbol for (layer, symbol), data in back.nodes(data=True)
                 if layer == REDUCED and data['must_not_be_zero']}
        self.assertEqual(
            named,
            set(self.translator.reduced_symbols_that_must_not_be_zero(self.system)))

    def test_a_zero_initial_condition_is_not_something_the_reduction_divides_by(self):
        # The three-case reading, on the graph: c_0 has one edge, it is positive, and it
        # goes to an invariant rather than an auxiliary.  So the invariant is simply zero.
        c_0 = (ORIGINAL, sympy.Symbol('c_0'))
        edges = list(self.graph.out_edges(c_0, data='exponent'))
        self.assertEqual(len(edges), 1)
        (_, target, exponent), = edges
        self.assertGreater(exponent, 0)
        self.assertFalse(self.graph.nodes[target]['auxiliary'])
        self.assertFalse(self.graph.nodes[c_0]['must_not_be_zero'])


class TestWhatBecomesOfTheAuxiliaries(TestCase):
    '''
    The auxiliary *columns* exist in the multiplier whatever reduction runs.  Whether the
    reduction sets them to 1 and throws them away, or keeps them as variables, is the
    scheme's business, so the graph derives it rather than assuming the commoner case.
    '''

    def test_the_parameter_scheme_normalises_them_going_out(self):
        system = ODESystem.from_equations(MICHAELIS_MENTEN)
        translator = ODETranslator.from_ode_system(system)
        self.assertEqual(translator.scheme(system), ReductionScheme.PARAMETER)
        graph = symbol_graph(translator, system)

        auxiliaries = [data for _, data in graph.nodes(data=True)
                       if data['auxiliary']]
        self.assertEqual(len(auxiliaries), translator.r)
        for data in auxiliaries:
            with self.subTest(symbol=str(data['symbol'])):
                self.assertTrue(data['normalised'])
                self.assertFalse(data['must_be_supplied'])

    def test_the_normalisation_lands_on_the_callers_own_symbols(self):
        # Where the 1s come from.  Setting the auxiliaries to 1 pins exactly r of the
        # original symbols, and which ones is not the graph's opinion: it is what
        # `translate_parameter_substitutions` -- the thing that actually rewrites the
        # system -- replaces with 1.
        system = ODESystem.from_equations(MICHAELIS_MENTEN)
        translator = ODETranslator.from_ode_system(system)
        graph = symbol_graph(translator, system)

        marked = {data['symbol'] for (layer, _), data in graph.nodes(data=True)
                  if layer == ORIGINAL and data['normalised']}
        substitutions = translator.translate_parameter_substitutions(system)
        self.assertEqual(marked,
                         {symbol for symbol, replacement in substitutions.items()
                          if replacement == 1})
        self.assertEqual(len(marked), translator.r)

    def test_a_symbol_the_reduction_spends_says_so_with_an_arrow(self):
        # `s_0 |-> 1` is an action, not an equality: the reduction measures concentration
        # in units of s_0, which is why `u = s/s_0` has no s_0 left in it.
        system = michaelis_menten_with_initial_conditions()
        translator = ODETranslator.from_ode_system(system)
        graph = symbol_graph(translator, system)
        s_0 = graph.nodes[(ORIGINAL, sympy.Symbol('s_0'))]
        self.assertTrue(s_0['normalised'])
        self.assertEqual(_node_label(s_0, latex=False), 's_0 |-> 1')

    def test_the_general_scheme_normalises_the_auxiliaries_it_can(self):
        # Not all or nothing.  An auxiliary is normalised exactly when it is a monomial in
        # the original parameters alone, because only then is it constant -- so the enzyme
        # chain normalises three of its four and keeps `x3 = x`, which moves.
        system = ChemicalReactionNetwork.from_diagram(
            'x + e -> y + e\ny -> z\nz + f -> x + f').to_ode_system()
        translator = ODETranslator.from_ode_system(system)
        graph = symbol_graph(translator, system)

        normalised = sorted(str(data['symbol']) for _, data in graph.nodes(data=True)
                            if data['auxiliary'] and data['normalised'])
        self.assertEqual(['x0', 'x1', 'x2'], normalised)
        kept = sorted(str(data['symbol']) for _, data in graph.nodes(data=True)
                      if data['auxiliary'] and not data['normalised'])
        self.assertEqual(['x3'], kept)

        # and on the original layer, the symbols that pays for
        self.assertEqual(
            ['k_0_1', 'k_2_3', 'k_4_5'],
            sorted(str(data['symbol']) for (layer, _), data in graph.nodes(data=True)
                   if layer == ORIGINAL and data['normalised']))

    def test_a_normalised_auxiliary_is_never_given_the_value_one(self):
        # The reduction reaches that 1 by choosing a group element, which rescales the
        # original symbols.  So `x0 = 1` holds of the rescaled coordinates, not of the
        # caller's, and writing it beside `x0 = 1/(e_0*k_1)` would read as a constraint on
        # `e_0` and `k_1` -- which are the caller's to choose and are under no such
        # constraint.  True in both directions, so neither carries it.
        system = ODESystem.from_equations(MICHAELIS_MENTEN)
        translator = ODETranslator.from_ode_system(system)
        for direction in ('forward', 'reverse'):
            graph = symbol_graph(translator, system, direction=direction)
            for _, data in graph.nodes(data=True):
                if data['normalised']:
                    with self.subTest(direction=direction,
                                      symbol=str(data['symbol'])):
                        self.assertNotIn('value', data)

    def test_coming_back_they_are_unknowns_and_not_one(self):
        # The opposite fact, and the reason this is direction-dependent.  A discarded
        # auxiliary is 1 because the reduction chose it; recovering the original values
        # has to be *told* r of them, which is what InsufficientKnownValues is about.
        # Carrying `value = 1` into the reverse graph would make `s_0 = x1` read as
        # `s_0 = 1` for every system, which is false.
        system = ODESystem.from_equations(MICHAELIS_MENTEN)
        translator = ODETranslator.from_ode_system(system)
        graph = symbol_graph(translator, system, direction='reverse')

        auxiliaries = [data for _, data in graph.nodes(data=True)
                       if data['auxiliary']]
        self.assertEqual(len(auxiliaries), translator.r)
        for data in auxiliaries:
            with self.subTest(symbol=str(data['symbol'])):
                self.assertTrue(data['normalised'])
                self.assertTrue(data['must_be_supplied'])
                self.assertNotIn('value', data)

    def test_how_many_must_be_supplied_is_the_number_of_symmetries(self):
        system = michaelis_menten_with_initial_conditions()
        translator = ODETranslator.from_ode_system(system)
        graph = symbol_graph(translator, system, direction='reverse')
        supplied = [data for _, data in graph.nodes(data=True)
                    if data['must_be_supplied']]
        self.assertEqual(len(supplied), translator.r)

    def test_nothing_must_be_supplied_going_out(self):
        system = michaelis_menten_with_initial_conditions()
        translator = ODETranslator.from_ode_system(system)
        graph = symbol_graph(translator, system)
        self.assertEqual([data for _, data in graph.nodes(data=True)
                          if data['must_be_supplied']], [])

    def test_the_reverse_label_does_not_claim_an_auxiliary_is_one(self):
        system = ODESystem.from_equations(MICHAELIS_MENTEN)
        translator = ODETranslator.from_ode_system(system)
        graph = symbol_graph(translator, system, direction='reverse')
        x1 = graph.nodes[(REDUCED, sympy.Symbol('x1'))]
        self.assertEqual(
            _node_label(x1, with_monomial=False, latex=False,
                        show_normalisation=False),
            'x1')

    def test_the_general_scheme_keeps_them_so_they_are_not_normalised(self):
        # Hubert & Labahn 6.6: no parameters at all and a clock entangled with a dependent
        # variable, so only the general scheme applies and its auxiliaries survive as
        # variables of the reduced system.
        system = ODESystem.from_equations(
            ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
             'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])
        translator = ODETranslator.from_ode_system(system)
        self.assertEqual(translator.scheme(system), ReductionScheme.GENERAL)
        graph = symbol_graph(translator, system)

        auxiliaries = [data for _, data in graph.nodes(data=True)
                       if data['auxiliary']]
        self.assertTrue(auxiliaries)
        for data in auxiliaries:
            with self.subTest(symbol=str(data['symbol'])):
                self.assertFalse(data['normalised'])
        # And they really are variables of the reduction, which is why.
        reduced = translator.translate(system)
        for data in auxiliaries:
            self.assertIn(data['symbol'], reduced.symbols)

    def test_normalisation_marks_both_sides_of_the_same_act(self):
        # Setting the auxiliaries to 1 and spending r of the caller's symbols are not two
        # things: `x1 = s_0`, so normalising `x1` *is* measuring in units of `s_0`.  The
        # figure marks both, and the counts match because they are the same act.
        system = michaelis_menten_with_initial_conditions()
        translator = ODETranslator.from_ode_system(system)
        graph = symbol_graph(translator, system)

        marked = {ORIGINAL: set(), REDUCED: set()}
        for (layer, symbol), data in graph.nodes(data=True):
            if data['normalised']:
                marked[layer].add(symbol)

        self.assertEqual(len(marked[REDUCED]), translator.r)
        self.assertEqual(len(marked[ORIGINAL]), translator.r)
        # On the reduced side it is exactly the auxiliaries, and nothing else.
        for symbol in marked[REDUCED]:
            with self.subTest(symbol=str(symbol)):
                self.assertTrue(graph.nodes[(REDUCED, symbol)]['auxiliary'])
        # On the original side it is the symbols those auxiliaries are built from.
        self.assertEqual(
            marked[ORIGINAL],
            set().union(*(monomial.free_symbols
                          for monomial in translator.auxiliary_monomials())))

    def test_a_normalised_auxiliary_is_labelled_with_what_it_is(self):
        # What the auxiliary *is*, and nothing about what the reduction does to it: that
        # part is the hatch's to say, not an equation's.
        system = ODESystem.from_equations(MICHAELIS_MENTEN)
        translator = ODETranslator.from_ode_system(system)
        graph = symbol_graph(translator, system)
        x0 = graph.nodes[(REDUCED, sympy.Symbol('x0'))]
        self.assertEqual(_node_label(x0, with_monomial=True, latex=False),
                         'x0 = 1/(e_0*k_1) |-> 1')


class TestTheSchemeIsAFeatureOfTheMultiplier(TestCase):
    '''
    Which scheme applies is a condition on the multiplier, so a column operation can change
    it.  That is why `symbol_graph` takes no scheme and remembers none: either would be
    stale the moment a column moved.
    '''

    def setUp(self):
        self.system = ODESystem.from_equations(MICHAELIS_MENTEN)
        self.translator = ODETranslator.from_ode_system(self.system)

    def move_the_clock_invariant(self):
        '''Multiply the clock invariant by a dependent one, which leaves the block form.'''
        self.translator.multiplier_add_columns(
            self.translator.column_of('e_0*k_1*t'),
            self.translator.column_of('c/e_0'))

    def test_a_column_operation_can_change_which_scheme_applies(self):
        self.assertEqual(self.translator.scheme(self.system), ReductionScheme.PARAMETER)
        self.move_the_clock_invariant()
        self.assertEqual(self.translator.scheme(self.system), ReductionScheme.GENERAL)

    def test_the_graph_follows_the_multiplier_rather_than_what_was_translated(self):
        # Translate first, so anything that recorded a scheme would have recorded
        # `parameter`; then move a column and ask again.  The graph must have moved too.
        self.translator.translate(self.system)
        before = {str(symbol) for (layer, symbol) in symbol_graph(
            self.translator, self.system).nodes if layer == REDUCED}
        self.assertIn('kappa0', before)

        self.move_the_clock_invariant()
        after = {str(symbol) for (layer, symbol) in symbol_graph(
            self.translator, self.system).nodes if layer == REDUCED}
        self.assertNotIn('kappa0', after)
        self.assertIn('y0', after)

    def test_normalisation_follows_the_column_and_not_the_scheme(self):
        # It used to be read off the scheme -- all of them under the parameter scheme,
        # none otherwise.  It is per column now: an auxiliary is normalised when its
        # monomial is built from the original parameters alone.  Moving the clock
        # invariant changes the scheme here, and leaves that answer alone, because it does
        # not change what the auxiliary columns are made of.
        before = symbol_graph(self.translator, self.system)
        self.assertTrue(all(data['normalised'] for _, data in before.nodes(data=True)
                            if data['auxiliary']))

        self.move_the_clock_invariant()
        after = symbol_graph(self.translator, self.system)
        auxiliaries = [data for _, data in after.nodes(data=True) if data['auxiliary']]
        self.assertTrue(auxiliaries)
        self.assertTrue(all(data['normalised'] for data in auxiliaries))


class TestAnAuxiliaryIsAParameterOrADependent(TestCase):
    '''
    Being an auxiliary says the scaling moves the symbol.  It says nothing about whether
    the symbol evolves, and the two get confused because `symbol_kinds` reports
    `'auxiliary'` in the slot where the other answer would go.

    The mathematics, which Hubert & Labahn leave implicit: the auxiliary system is a
    quadrature whose integrand is the logarithmic derivatives of the original symbols
    combined by a column of V_i, and a parameter contributes zero to it.  So an auxiliary
    that is a monomial in the parameters alone is constant -- a constant of integration
    rather than an integral.  That is why the parameter scheme can normalise its
    auxiliaries away: its block form puts every column of V_i on parameter rows.
    '''

    CASES = {
        'enzyme chain': ['x + e -> y + e', 'y -> z', 'z + f -> x + f'],
    }

    def enzyme(self):
        return ChemicalReactionNetwork.from_diagram(
            '\n'.join(self.CASES['enzyme chain'])).to_ode_system()

    def test_a_parameter_auxiliary_and_a_dependent_auxiliary_can_coexist(self):
        # The case that motivates all of this: three auxiliaries that sit still and one
        # that moves, in one reduction.
        system = self.enzyme()
        translator = ODETranslator.from_ode_system(system)
        graph = symbol_graph(translator, system)
        kinds = {str(data['symbol']): data['kind']
                 for _, data in graph.nodes(data=True) if data['auxiliary']}
        self.assertEqual(kinds, {'x0': 'parameter', 'x1': 'parameter',
                                 'x2': 'parameter', 'x3': 'dependent'})

    def test_the_kinds_agree_with_the_reduced_system(self):
        # The oracle.  `1/k_2_3` is called a parameter here because its derivative is
        # zero; the reduced system is where that claim can be checked, since it lists its
        # own parameters and dependents.  Derived one way, verified the other.
        for label, system in (('enzyme chain', self.enzyme()),
                              ('HL 6.6', ODESystem.from_equations(
                                  ['dz1/dt = z1*(z1**5*z2 - 2)/(3*t)',
                                   'dz2/dt = z2*(10 - 2*z1**5*z2 + 3*z1**2*z2/t )/(3*t)'])),
                              ('prey-predator', ODESystem.from_equations(
                                  ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )',
                                   'dp/dt = s*p*(1 - h*p / n)']))):
            translator = ODETranslator.from_ode_system(system)
            graph = symbol_graph(translator, system)
            # The figure shows every column, including the ones a reduction normalises
            # away, so the reduction that keeps them all is the oracle that can speak
            # about every node.  A default reduction has nothing to say about `x0`.
            reduced = translator.translate_general(system, normalise_parameters=False)
            for _, data in graph.nodes(data=True):
                if not data['auxiliary']:
                    continue
                symbol = data['symbol']
                with self.subTest(system=label, symbol=str(symbol)):
                    if data['kind'] == 'parameter':
                        self.assertIn(symbol, reduced.parameters)
                    else:
                        self.assertIn(symbol, reduced.dependents)

    def test_a_constant_auxiliary_really_has_no_equation(self):
        # The derivative being zero is the whole argument, so assert it directly rather
        # than only through the classification.
        system = self.enzyme()
        translator = ODETranslator.from_ode_system(system)
        reduced = translator.translate_general(system)
        graph = symbol_graph(translator, system)
        for _, data in graph.nodes(data=True):
            if data['auxiliary'] and data['kind'] == 'parameter':
                with self.subTest(symbol=str(data['symbol'])):
                    self.assertNotIn(data['symbol'], reduced.derivative_dict)
                    self.assertFalse(data['monomial'].free_symbols
                                     & set(system.variables))

    def test_under_the_parameter_scheme_every_auxiliary_is_a_parameter(self):
        # Forced by the block form the scheme requires: V_i is supported on parameter
        # rows, so no auxiliary can contain a variable.
        for label, equations in (('michaelis-menten', MICHAELIS_MENTEN),
                                 ('prey-predator',
                                  ['dn/dt = n*( r*(1 - n/K) - k*p/(n+d) )',
                                   'dp/dt = s*p*(1 - h*p / n)'])):
            system = ODESystem.from_equations(equations)
            translator = ODETranslator.from_ode_system(system)
            self.assertEqual(translator.scheme(system), ReductionScheme.PARAMETER)
            graph = symbol_graph(translator, system)
            for _, data in graph.nodes(data=True):
                if data['auxiliary']:
                    with self.subTest(system=label, symbol=str(data['symbol'])):
                        self.assertEqual(data['kind'], 'parameter')


class TestDrawing(TestCase):

    def setUp(self):
        self.system = michaelis_menten_with_initial_conditions()
        self.translator = ODETranslator.from_ode_system(self.system)

    def test_the_source_layer_is_on_the_left_in_both_directions(self):
        # So the arrows always run left to right.
        for direction, source in (('forward', ORIGINAL), ('reverse', REDUCED)):
            with self.subTest(direction=direction):
                graph = symbol_graph(self.translator, self.system, direction=direction)
                positions = _positions(graph)
                sides = {layer: {positions[node][0]
                                 for node, data in graph.nodes(data=True)
                                 if data['layer'] == layer}
                         for layer in (ORIGINAL, REDUCED)}
                other = REDUCED if source == ORIGINAL else ORIGINAL
                self.assertEqual(sides[source], {0.0})
                self.assertEqual(sides[other], {1.0})

    def test_the_reduced_column_lines_up_across_the_pair(self):
        # Draw forward beside reverse and the reduced column ends the first panel and
        # begins the second.  For a round trip to read as one path across the two, those
        # have to be at the same heights and in the same order.
        forward = symbol_graph(self.translator, self.system)
        reverse = symbol_graph(self.translator, self.system, direction='reverse')
        forward_positions, reverse_positions = _positions(forward), _positions(reverse)
        for node, data in forward.nodes(data=True):
            if data['layer'] == REDUCED:
                with self.subTest(node=str(node[1])):
                    self.assertEqual(forward_positions[node][1],
                                     reverse_positions[node][1])

    def test_a_column_operation_does_not_move_anything(self):
        # The other comparison the fixed layout exists for: before against after.
        before = _positions(symbol_graph(self.translator, self.system))
        self.translator.multiplier_add_columns(
            self.translator.column_of('c_0/s_0'), self.translator.column_of('e_0/s_0'))
        after = _positions(symbol_graph(self.translator, self.system))
        self.assertEqual(before, after)

    def test_the_layout_does_not_depend_on_anything_but_the_order(self):
        twice = [_positions(symbol_graph(self.translator, self.system)) for _ in range(2)]
        self.assertEqual(*twice)

    def test_the_two_layers_are_two_columns(self):
        positions = _positions(symbol_graph(self.translator, self.system))
        columns = {layer: {positions[node][0]
                           for node, data in
                           symbol_graph(self.translator, self.system).nodes(data=True)
                           if data['layer'] == layer}
                   for layer in (ORIGINAL, REDUCED)}
        self.assertEqual(len(columns[ORIGINAL]), 1)
        self.assertEqual(len(columns[REDUCED]), 1)
        self.assertNotEqual(columns[ORIGINAL], columns[REDUCED])

    def test_a_label_chains_the_symbol_its_expression_and_its_value(self):
        # Every part of it is genuinely equal, so the chain is true as written and says
        # the whole zero-initial-condition story in one line.
        graph = symbol_graph(self.translator, self.system)
        kappa4 = graph.nodes[(REDUCED, sympy.Symbol('kappa4'))]
        self.assertEqual(_node_label(kappa4, with_monomial=True, latex=False),
                         'kappa4 = c_0/s_0 = 0')

    def test_a_label_without_a_known_value_stops_at_the_expression(self):
        graph = symbol_graph(self.translator, self.system)
        kappa0 = graph.nodes[(REDUCED, sympy.Symbol('kappa0'))]
        self.assertEqual(_node_label(kappa0, with_monomial=True, latex=False),
                         'kappa0 = e_0/s_0')

    def test_latex_labels_are_mathtext(self):
        graph = symbol_graph(self.translator, self.system)
        label = _node_label(graph.nodes[(REDUCED, sympy.Symbol('kappa4'))],
                            with_monomial=True)
        self.assertTrue(label.startswith('$') and label.endswith('$'))
        self.assertIn(r'\kappa', label)
        # Display size, because a graph is a figure and an inline fraction comes out
        # half height with its subscripts smaller again.
        self.assertIn(r'\dfrac', label)

    def test_a_negative_subscript_is_set_as_a_minus(self):
        # desr writes k_m1 for what a paper writes as k_{-1}, and sympy's own printer
        # sets that m literally.
        from desr.tex_tools import expr_to_tex
        system = ODESystem.from_equations(
            ['ds/dt = -k_1*e_0*s + k_1*c*s + k_m1*c',
             'dc/dt = k_1*e_0*s - k_1*c*s - k_m1*c - k_2*c'])
        translator = ODETranslator.from_ode_system(system)
        graph = symbol_graph(translator, system)
        label = _node_label(graph.nodes[(ORIGINAL, sympy.Symbol('k_m1'))])
        self.assertIn('k_{-1}', label)
        self.assertNotIn('k_{m1}', label)
        self.assertEqual('k_{-1}', expr_to_tex(sympy.Symbol('k_m1')))

    def test_drawing_returns_the_axes_it_drew_on(self):
        import matplotlib.pyplot
        _, ax = matplotlib.pyplot.subplots()
        self.assertIs(draw_symbol_graph(symbol_graph(self.translator, self.system),
                                        ax=ax), ax)
        matplotlib.pyplot.close('all')

    def test_drawing_without_axes_makes_its_own(self):
        import matplotlib.pyplot
        ax = draw_symbol_graph(symbol_graph(self.translator, self.system))
        self.assertIsNotNone(ax.get_figure())
        matplotlib.pyplot.close('all')

    def test_both_directions_draw(self):
        import matplotlib.pyplot
        for direction in ('forward', 'reverse'):
            with self.subTest(direction=direction):
                graph = symbol_graph(self.translator, self.system, direction=direction)
                self.assertIsNotNone(draw_symbol_graph(graph))
        matplotlib.pyplot.close('all')


class TestTheGraphAsTheNarrativePageShowsIt(TestCase):
    '''
    `docs/source/examples/mm_init_cond.rst` reads the zero initial condition off the graph
    after both of its column operations.  Its doctests are a backstop; the assertions are
    here, so a regression surfaces before anyone builds the docs.
    '''

    def setUp(self):
        self.system = michaelis_menten_with_initial_conditions()
        self.translator = ODETranslator.from_ode_system(
            self.system,
            naming_scheme={'parameters': 'kappa', 'independent': 'tau',
                           'dependents': ['u', 'v']})
        self.translator.multiplier_add_columns(
            self.translator.column_of('k_1*s_0*t'),
            self.translator.column_of('e_0/s_0'))
        self.translator.multiplier_add_columns(
            self.translator.column_of('c/s_0'),
            self.translator.column_of('e_0/s_0'), power=-1)
        self.graph = self.translator.symbol_graph(self.system)

    def test_the_size_the_page_quotes(self):
        self.assertEqual(
            (self.graph.number_of_nodes(), self.graph.number_of_edges()), (20, 22))

    def test_the_zero_has_one_positive_edge_and_it_is_not_into_an_auxiliary(self):
        edges = list(self.graph.out_edges((ORIGINAL, sympy.Symbol('c_0')),
                                          data='exponent'))
        self.assertEqual(edges, [((ORIGINAL, sympy.Symbol('c_0')),
                                  (REDUCED, sympy.Symbol('kappa4')), 1)])
        self.assertFalse(
            self.graph.nodes[(REDUCED, sympy.Symbol('kappa4'))]['auxiliary'])

    def test_the_invariant_it_feeds_is_zero(self):
        kappa4 = self.graph.nodes[(REDUCED, sympy.Symbol('kappa4'))]
        self.assertEqual(kappa4['monomial'], sympy.sympify('c_0/s_0'))
        self.assertEqual(kappa4['value'], 0)

    def test_the_two_answers_the_page_quotes(self):
        self.assertEqual(self.translator.symbols_that_must_not_be_zero(self.system),
                         tuple(sympy.symbols('s_0 e_0 k_1')))
        self.assertEqual(
            self.translator.reduced_symbols_that_must_not_be_zero(self.system),
            tuple(sympy.symbols('x0 x1 kappa0')))

    def test_the_zero_returns_the_way_it_came(self):
        back = self.translator.symbol_graph(self.system, direction='reverse')
        self.assertEqual(back.nodes[(ORIGINAL, sympy.Symbol('c_0'))]['monomial'],
                         sympy.sympify('kappa4*x1'))

    def test_two_auxiliaries_must_be_supplied_coming_back(self):
        back = self.translator.symbol_graph(self.system, direction='reverse')
        supplied = [data['symbol'] for _, data in back.nodes(data=True)
                    if data['must_be_supplied']]
        self.assertEqual(sorted(map(str, supplied)), ['x0', 'x1'])


class TestTheTranslatorMethods(TestCase):
    ''' The method is the way in; the free functions are what it forwards to. '''

    def setUp(self):
        self.system = ODESystem.from_equations(MICHAELIS_MENTEN)
        self.translator = ODETranslator.from_ode_system(self.system)

    def test_the_method_gives_the_same_graph_as_the_function(self):
        by_method = self.translator.symbol_graph(self.system)
        by_function = symbol_graph(self.translator, self.system)
        self.assertEqual(_edges(by_method), _edges(by_function))
        self.assertEqual(set(by_method.nodes), set(by_function.nodes))

    def test_the_method_passes_the_direction_through(self):
        back = self.translator.symbol_graph(self.system, direction='reverse')
        self.assertEqual(back.graph['direction'], 'reverse')

    def test_the_drawing_method_forwards_its_keywords(self):
        import matplotlib.pyplot
        _, ax = matplotlib.pyplot.subplots()
        self.assertIs(self.translator.draw_symbol_graph(self.system, ax=ax), ax)
        matplotlib.pyplot.close('all')
