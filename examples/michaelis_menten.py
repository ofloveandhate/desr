"""Michaelis-Menten Equations

This example script walks through an analysis of the Michaelis Menten equations.

The original source for the equations is
L. Michaelis and M. Menten, Die kinetik der Invertinwirkung, Biochem Z., 49 (1913), pp. 333–369.

and for a modern source, see

J. Gunawardena, Time-scale separation–Michaelis and Menten’s old idea, still bearing
fruit., The FEBS journal, 281 (2014), pp. 473–88, https://doi.org/10.1111/febs.12532,
http://www.ncbi.nlm.nih.gov/pubmed/24103070http://www.pubmedcentral.nih.gov/
articlerender.fcgi?artid=PMC3991559.
"""

import sympy
from desr import ODESystem, ODETranslator
from desr.matrix_normal_forms import smith_normal_form
from desr.sympy_helper import monomial_to_powers, sympify_expression




def michaelis_menten_four_variables():
    """
    Example: Michaelis-Menten kinetics in four variables with initial conditions

    We start with ODEs.
    """

    print('\n\n----------------------\nmichaelis_menten_four_variables\n-----------------------\n\n')


    def system_four_variable_michaelis_menten():


        # make string to hold the system in latex form
        # note that, so make valid variable names,
        # minus signs  `-` in subscripts are changed to the character `m`
        system_tex = \
              '''\frac{dS}{dt} &= - k_1 E S + k_{-1} C \\\\
                 \frac{dE}{dt} &= - k_1 E S + k_{-1} C + k_2 C \\\\
                 \frac{dC}{dt} &=   k_1 E S - k_{-1} C - k_2 C \\\\
                 \frac{dP}{dt} &=   k_2 C'''

        system = ODESystem.from_tex(system_tex)

        return system


    # There used to be a `put_at_front` helper here, moving one symbol at a time to the
    # front of the whole symbol list.  `reorder_parameters` with the order written out says
    # the same thing and says it once -- and it cannot put a parameter after the
    # dependents, which is what silently costs the parameter reduction scheme.



    def four_variable_with_ics():
        print('\n\nFour-variable with ICs\n-----------------------------------')

        system = system_four_variable_michaelis_menten()
        # The complex and the product both start at nothing.  Saying so is the point.
        system.update_initial_conditions(['S(0) = S_0', 'E(0) = E_0',
                                          'C(0) = 0', 'P(0) = 0'])
        print('initial conditions:', system.initial_conditions)
        print('as values         :', system.initial_condition_values)

        # Written out, rather than built by moving one symbol at a time: the front of the
        # order is the end that gets scaled by, so S_0 and k_1 are what this normalises --
        # the same clock the two-variable section uses, k_1 S_0 t.  The initial values of
        # the complex and the product go at the back, which is the end that is kept, since
        # a quantity known to be zero must not be scaled by.
        #
        # `reorder_parameters` rather than `reorder_variables`, so the parameter block is
        # the only thing that moves.  Reordering the whole symbol list by hand is how a
        # parameter ends up after the dependents, which silently costs the parameter
        # reduction scheme.
        system.reorder_parameters(['S_0', 'k_1', 'k_2', 'k_m1', 'E_0', 'C_0', 'P_0'])
        print('system\n',system)


        translator = ODETranslator.from_ode_system(system)


        print('variable order:', translator.symbols_domain)
        print('scheme:', translator.scheme(system))

        # Asked, not assumed.  S_0 and k_1 are what the reduction divides by; C_0 and P_0
        # are not, which is what makes a zero harmless for them.
        print('must not be zero:', translator.symbols_that_must_not_be_zero())

        # Print scaling matrices
        print('Scaling matrix:')
        print(translator.scaling_matrix.__repr__())

        # Print invariants
        print('Invariants: ',
              ['{}={}'.format(a, b) for a, b
               in zip(translator.invariant_names(system),
                      translator.invariant_monomials())])

        # Print translated result
        print('Reduced system:')
        reduced_system = translator.translate(system)
        print(reduced_system)

        # The zeros arrive here derived, not restated: each is the invariant holding
        # C_0 or P_0 in its numerator, so the complex and the product start at nothing in
        # the reduced system too.
        print('initial conditions, as values:',
              reduced_system.initial_condition_values)
        assert (set(reduced_system.initial_condition_values.values())
                == {0, 1}), reduced_system.initial_condition_values

        print('Substitutions:')
        print(translator.translate_parameter_substitutions(system=system))

    def four_variable_choosing_invariants():
        """
        Impose a choice of invariants, and rearrange into the parameter scheme by hand.

        The same reduction as `docs/source/examples/mm_four_variable.rst`, with the
        extension step shown rather than taken on trust, because that step is what the
        paper describes.
        """
        print('\n\nFour-variable, imposing a choice of invariants\n-----------------------------------')

        system = system_four_variable_michaelis_menten()
        print('Variable order:', system.symbols)

        translator = ODETranslator.from_ode_system(system)
        print('Invariants:', ', '.join(map(str, translator.invariant_monomials())))

        # The reduction the algorithm picks on its own, for comparison.
        print('Reduced system:', translator.translate(system))


        print('\nImposing a choice of invariants')
        print('--------------------------------')
        # Say which invariants, by writing them.  An exponent matrix would say the same
        # thing, but its rows are *positions* in the variable order, so reordering the
        # variables would leave the same matrix meaning different monomials -- and nothing
        # in it records which invariants were meant.  This is how the section used to be
        # written, and it had rotted into denoting `E*k_2/P, E*t/S`, which are not
        # invariants at all.
        invariant_choice = ['C*k_1/k_2', 'P*k_1/k_m1']
        print('Chosen invariants:', invariant_choice)

        chosen = translator.extend_from_invariants(invariant_choice)
        print('Extended to:', ', '.join(map(str, chosen.invariant_monomials())))


        print('\nThe extension step, for the paper')
        print('--------------------------------')
        # What `extend_from_invariants` does, spelled out.  The exponent matrix is needed
        # here because this is the linear algebra; it is derived from the monomials above,
        # so there is one statement of the choice rather than two that can disagree.
        exponents = sympy.Matrix(
            [monomial_to_powers(sympify_expression(m), list(system.symbols))
             for m in invariant_choice]).T
        print('P =')
        print(exponents.__repr__())

        assert (translator.scaling_matrix * exponents).is_zero_matrix, (
            'the choice is not invariant')

        WdP = translator.inv_herm_mult_d * exponents
        print(f'{WdP=}')
        normal_form, row_ops, snf_col_ops = smith_normal_form(WdP)
        print('Smith normal form of W_d P:')
        print('{}\n*\n{}\n*\n{}\n=\n{}'.format(
            *map(repr, (row_ops, WdP, snf_col_ops, normal_form))))
        print('U^{-1} =', row_ops.inv().__repr__())

        col_ops = sympy.Matrix.hstack(WdP, row_ops.inv()[:, len(invariant_choice):])
        by_hand = ODETranslator(
            translator.scaling_matrix, symbols_domain=system.symbols,
            hermite_multiplier=sympy.Matrix.hstack(
                translator.herm_mult_i, translator.herm_mult_n * col_ops))
        print('invariants this way:', ', '.join(map(str, by_hand.invariant_monomials())))

        # The same invariants, in a different column order.  Worth asserting rather than
        # asserting in prose: the hand-rolled extension and the library's are the same
        # reduction, and the column order is all that separates them.
        assert (sorted(map(str, by_hand.invariant_monomials()))
                == sorted(map(str, chosen.invariant_monomials()))), (
            'the two extensions disagree on the invariants themselves')


        print('\nRearranging into the parameter scheme')
        print('--------------------------------')
        # None of the six invariants is a monomial in the parameters alone, so the
        # parameter scheme has nothing to reduce yet.  One column operation makes one:
        # dividing P k_1/k_{-1} by P k_1/k_2 cancels the P and leaves k_2/k_{-1}.
        ratio = chosen.column_of('P*k_1/k_m1')
        product = chosen.column_of('P*k_1/k_2')
        print(f'{ratio=}, {product=}')
        chosen.multiplier_add_columns(ratio, product, power=-1)
        print('invariants now:', ', '.join(map(str, chosen.invariant_monomials())))

        # The scheme wants the parameters first, then the independent variable, then the
        # dependents.  Derived from what each column *is*, so the permutation cannot go
        # stale when the variable order changes -- which is what a literal
        # `[1, 2, 0, 3, 4, 5]` would do.
        def column_kind(monomial):
            symbols = monomial.free_symbols
            if system.independent_variable in symbols:
                return 1
            return 2 if symbols & set(system.variables) else 0

        monomials = list(chosen.invariant_monomials())
        permutation = sorted(range(len(monomials)),
                             key=lambda i: (column_kind(monomials[i]), i))
        print('permutation derived from the kinds:', permutation)

        rearranged = ODETranslator(
            chosen.scaling_matrix, symbols_domain=system.symbols,
            hermite_multiplier=sympy.Matrix.hstack(
                chosen.herm_mult_i, chosen.herm_mult_n[:, permutation]),
            naming_scheme={'parameters': 'c', 'independent': 't',
                           'dependents': ['C', 'E', 'P', 'S']})
        print('final invariants:',
              ', '.join(map(str, rearranged.invariant_monomials())))

        print('Reduced system:')
        reduced_system = rearranged.translate(system)
        print(reduced_system.subs({'c0': 'alpha'}))
        assert len(reduced_system.parameters) == 1, (
            'expected a one-parameter reduction, got {}'.format(
                reduced_system.parameters))



    four_variable_with_ics()
    four_variable_choosing_invariants()



def michaelis_menten_two_variable(verbose=True):
    """
    Example `sec:michaelis_menten_comparison_with_classical`
    """

    print('\n\n----------------------\nmichaelis_menten_two_variable\n-----------------------\n\n')

    # Enable pretty printing
    sympy.init_printing(pretty_print=True, use_latex=True)

    def system_tex_two_var():
        return '''\frac{ds}{dt} &= - k_1 e_0 s + k_1 c s + k_{-1} c \\\\
             \frac{dc}{dt} &= k_1 e_0 s - k_1 c s - k_{-1} c - k_2 c'''


    def with_initial_conditions_and_a_constant():
        '''
        this code reproduces the michaelis-menten example in the literature.
        '''
        system_tex = system_tex_two_var()


        system = ODESystem.from_tex(system_tex)
        system.add_constraint('Km', '(k_m1 + k_2) / k_1')




        print('\n\nAdding in the initial conditions')
        print('-------------------------------------')
        # The complex starts at nothing.  Saying so used to be impossible, and the
        # literature just glosses over it; `c(0) = 0` records the value against a symbol
        # `c_0`, which is what carries the units -- a zero concentration is still a
        # concentration, and the reduction stays entirely symbolic.
        system.update_initial_conditions(['s(0) = s0', 'c(0) = 0'])
        print('initial conditions:', system.initial_conditions)
        print('as values         :', system.initial_condition_values)

        # `c_0` goes towards the back, which is the end that is *kept*.  A zero must not be
        # scaled by, and the front of the order is the scale-by end.
        system.reorder_variables(['s0','e_0', 'k_1', 'k_2', 'k_m1', 'Km', 'c_0', 't', 'c', 's'])


        print('System:', system)
        print('Variable order: ', system.symbols)

        print('Exponent Matrix:', system.exponent_matrix().__repr__())
        # The dependents are in the order (c, s), so Murray's names go (v, u): his `u` is
        # the substrate and his `v` is the complex.
        #
        # `kappa` rather than `c` for the parameters.  They are substituted to Murray's
        # names below, so the letter is arbitrary -- but `c` is a dependent variable here,
        # whose initial value is the minted `c_0`, and printing `c_0` beside `c0, c1, ...`
        # asks a reader to tell apart two things that look the same.  They are different
        # symbols, so nothing breaks; it just cannot be read.
        translator = ODETranslator.from_ode_system(system, naming_scheme={'parameters': 'kappa', 'independent': 'tau', 'dependents': ['v','u']})
        print('symbols this reduction divides by, so none of them may be zero:',
              translator.symbols_that_must_not_be_zero())


        print('Scaling matrix:',translator.scaling_matrix.__repr__())

        print('Hermite multiplier', translator.herm_mult.__repr__())
        print('Invariants: ', translator.invariant_monomials())


        print('Reduced system:')
        print(translator.translate(system))#.to_tex()


        print('Changing invariants by column operations on the Hermite multiplier')

        # Two groups still differ from Murray's: his clock is k_1 e_0 t rather than
        # k_1 s0 t, and he measures the complex in units of the enzyme, c/e_0, rather than
        # of the substrate.  Both are a factor of epsilon = e_0/s0 away, so both are one
        # column operation.
        #
        # The column numbers come from `column_of`, so they are derived from the invariant
        # each one carries rather than counted off by hand.  They can still be printed and
        # quoted -- and they move when the variable order moves, which is exactly why they
        # should not be written as literals.
        clock = translator.column_of('k_1*s0*t')
        epsilon = translator.column_of('e_0/s0')
        complexed = translator.column_of('c/s0')
        print(f'{clock=}, {epsilon=}, {complexed=}')

        translator.multiplier_add_columns(clock, epsilon)              # tau -> k_1 e_0 t
        translator.multiplier_add_columns(complexed, epsilon, power=-1)  # v -> c/e_0

        print('Invariants: ',
            ['{}={}'.format(a,b) for a,b in zip(translator.invariant_names(system), translator.invariant_monomials())])
        print('Hermite multiplier', translator.herm_mult.__repr__())


        # print('\nMichaelis-Menten Reparametrisation 2')
        # print('Changing invariants by column operations on the Hermite multiplier')
        # # Divide time through by epsilon
        # translator.multiplier_add_columns(i=0, j=8, power=1)

        # print('Hermite multiplier', translator.herm_mult.__repr__())
        # print('Invariants: ', translator.invariant_monomials())

        print('W:', translator.inv_herm_mult.__repr__())
        # Print translated system

        print('substitutions:',translator.translate_parameter_substitutions(system))

        reduced_system = translator.translate(system)

        print('Reduced system, in arbitrary parameter names:')
        print(reduced_system)#.to_tex()

        print('the initial conditions, as values:',
              reduced_system.initial_condition_values)

        # Read off the monomials rather than guessed at, because the order of the kappa_i
        # follows the variable order and so is not something to remember:
        #
        #   kappa0 = e_0/s0        = epsilon
        #   kappa1 = k_2/(k_1 s0)  = lambda, which this script calls L
        #   kappa2 = k_m1/(k_1 s0) = K - L
        #   kappa3 = Km/s0         = K
        #   kappa4 = c_0/s0        = the initial complex, in substrate units
        #
        # Murray has no name for the last one because he sets it to zero, which is exactly
        # what we have recorded.  `gamma` until then.
        #
        # The check is the constraint.  `kappa3 == kappa1 + kappa2` has to become
        # `K == L + (K - L)`, a tautology, and so disappear.  If a constraint survives the
        # substitution, the mapping is wrong -- which is how this was caught: it used to
        # leave the false `epsilon == K` behind.
        to_sub = dict(zip(translator.invariant_names(system),
                          ['epsilon', 'L', 'K - L', 'K', 'gamma']))
        print('substituting Murray\'s names:', to_sub)
        reduced_system = reduced_system.subs(to_sub, subs_constraints=True)
        assert not reduced_system.constraints, (
            'a constraint survived the renaming, so the names do not match the monomials: '
            '{}'.format(reduced_system.constraints))

        print('after subbing into canonical names from Murray')
        print(reduced_system)#.to_tex()

        # `v(0) = gamma/epsilon`, and gamma is zero, so v starts at nothing however small
        # epsilon is.  Both faces of that are available: the expression says *what* the
        # initial value is, and the value says what it *is*.
        print('the initial conditions, as values:',
              reduced_system.initial_condition_values)

        # Murray writes three groups, not four, because he sets the initial complex to
        # zero.  Substituting the value the system already holds does that, and the group
        # goes with it -- derived from the stated initial condition rather than asserted.
        # Substituting anything else here is refused, since it would contradict `c(0) = 0`.
        print('\nSubstituting the initial complex, which is known to be zero')
        print('-----------------------------------------------------------')
        reduced_system = reduced_system.subs({'gamma': 0})
        print(reduced_system)
        assert set(map(str, reduced_system.parameters)) == {'K', 'L', 'epsilon'}, (
            'expected Murray (6.20) three groups, got {}'.format(
                reduced_system.parameters))

        return




    with_initial_conditions_and_a_constant()
    # another_way()

if __name__ == '__main__':
    # Enable pretty printing
    sympy.init_printing(pretty_print=True, use_latex=True)

    michaelis_menten_two_variable()

    michaelis_menten_four_variables()
