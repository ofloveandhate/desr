"""
Michaelis-Menten kinetics in two variables, with initial conditions and the
Michaelis constant, brought to the form in Murray, Mathematical Biology I, (6.20).

The equations are those of Michaelis and Menten, Die Kinetik der Invertinwirkung,
Biochem. Z. 49 (1913) 333-369.  For a modern account see Gunawardena, Time-scale
separation -- Michaelis and Menten's old idea, still bearing fruit, FEBS Journal 281
(2014) 473-488, https://doi.org/10.1111/febs.12532.
"""
import sympy

from desr import ODESystem, ODETranslator


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
        # The complex starts at nothing.  `c(0) = 0` records the value against a symbol
        # `c_0`, and `c_0` carries the units: a zero concentration is a concentration, so
        # the reduction stays symbolic.
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
        # `kappa` for the parameters.  They are substituted to Murray's names below, so the
        # letter is free; `c` is a dependent variable here and its initial value is `c_0`,
        # which a parameter named `c0` would be hard to read beside.
        translator = ODETranslator.from_ode_system(system, naming_scheme={'parameters': 'kappa', 'independent': 'tau', 'dependents': ['v','u']})
        print('symbols this reduction divides by, so none of them may be zero:',
              translator.forward_symbols_that_must_not_be_zero())


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
        # `column_of` takes the invariant a column carries and gives back its index.  The
        # indices are still numbers, printed below, and they follow the variable order.
        clock = translator.column_of('k_1*s0*t')
        epsilon = translator.column_of('e_0/s0')
        complexed = translator.column_of('c/s0')
        print(f'{clock=}, {epsilon=}, {complexed=}')

        translator.multiplier_add_columns(clock, epsilon)              # tau -> k_1 e_0 t
        translator.multiplier_add_columns(complexed, epsilon, power=-1)  # v -> c/e_0

        print('Invariants: ',
            ['{}={}'.format(a,b) for a,b in zip(translator.invariant_names(), translator.invariant_monomials())])
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

        # The names come from the monomials.  The order of the kappa_i follows the variable
        # order:
        #
        #   kappa0 = e_0/s0        = epsilon
        #   kappa1 = k_2/(k_1 s0)  = lambda, which this script calls L
        #   kappa2 = k_m1/(k_1 s0) = K - L
        #   kappa3 = Km/s0         = K
        #   kappa4 = c_0/s0        = the initial complex, in substrate units
        #
        # Murray has no name for the last one, because he sets it to zero.  `gamma` here.
        #
        # The check is the constraint.  `kappa3 == kappa1 + kappa2` has to become
        # `K == L + (K - L)`, a tautology, and so disappear.  If a constraint survives the
        # substitution, the mapping is wrong -- which is how this was caught: it used to
        # leave the false `epsilon == K` behind.
        to_sub = dict(zip(translator.invariant_names(),
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
        # zero.  Substituting the value the system holds for gamma does the same, and the
        # group goes with it.  Any other value is refused, as it would contradict
        # `c(0) = 0`.
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
    michaelis_menten_two_variable()
