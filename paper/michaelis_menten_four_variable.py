"""
Michaelis-Menten kinetics in four variables.

With initial conditions, and then with a choice of invariants imposed and rearranged
into the parameter scheme by hand.

The equations are those of Michaelis and Menten, Die Kinetik der Invertinwirkung,
Biochem. Z. 49 (1913) 333-369.  For a modern account see Gunawardena, Time-scale
separation -- Michaelis and Menten's old idea, still bearing fruit, FEBS Journal 281
(2014) 473-488, https://doi.org/10.1111/febs.12532.
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


    # Variable orders below are written out and passed to `reorder_parameters`, which moves
    # the parameter block and leaves the independent variable and the dependents where they
    # are.  The parameter reduction scheme needs that layout.



    def four_variable_with_ics():
        print('\n\nFour-variable with ICs\n-----------------------------------')

        system = system_four_variable_michaelis_menten()
        # Murray, Mathematical Biology I, (6.4): all the substrate and all the enzyme are
        # there at the start, and none of the complex or the product.
        system.update_initial_conditions(['S(0) = S_0', 'E(0) = E_0',
                                          'C(0) = 0', 'P(0) = 0'])
        print('initial conditions:', system.initial_conditions)
        print('as values         :', system.initial_condition_values)

        # The front of the order is the end that gets scaled by, so S_0 and k_1 are the two
        # this normalises, giving the clock k_1 S_0 t that the two-variable section uses.
        # C_0 and P_0 go at the back, the end that is kept: a symbol known to be zero
        # cannot be scaled by.
        system.reorder_parameters(['S_0', 'k_1', 'k_2', 'k_m1', 'E_0', 'C_0', 'P_0'])
        print('system\n',system)


        translator = ODETranslator.from_ode_system(system)


        print('variable order:', translator.symbols_domain)
        print('scheme:', translator.scheme())

        # S_0 and k_1 are what this reduction divides by.  C_0 and P_0 are not, so their
        # being zero leaves the reduction defined.
        print('must not be zero:', translator.forward_symbols_that_must_not_be_zero())

        # Print scaling matrices
        print('Scaling matrix:')
        print(translator.scaling_matrix.__repr__())

        # Print invariants
        print('Invariants: ',
              ['{}={}'.format(a, b) for a, b
               in zip(translator.invariant_names(),
                      translator.invariant_monomials())])

        # Print translated result
        print('Reduced system:')
        reduced_system = translator.translate(system)
        print(reduced_system)

        # Each zero reaches the reduced system through the invariant holding C_0 or P_0 in
        # its numerator, so the complex and the product start at nothing here too.
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
        extension step written out, since that step is what the paper describes.
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
        # The invariants are written as monomials.  `extend_from_invariants` also takes an
        # exponent matrix, whose rows are positions in the variable order; the monomials
        # say which invariants are meant whatever that order is.
        invariant_choice = ['C*k_1/k_2', 'P*k_1/k_m1']
        print('Chosen invariants:', invariant_choice)

        chosen = translator.extend_from_invariants(invariant_choice)
        print('Extended to:', ', '.join(map(str, chosen.invariant_monomials())))


        print('\nThe extension step, for the paper')
        print('--------------------------------')
        # What `extend_from_invariants` does, written out.  The linear algebra needs the
        # exponent matrix, built here from the monomials above so that the choice is stated
        # in one place.
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

        # This extension and `extend_from_invariants` give the same invariants; only the
        # column order differs.
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
        # dependents.  Which column is which is read off what each one is, so the
        # permutation follows the variable order instead of being fixed to one.
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


if __name__ == '__main__':
    michaelis_menten_four_variables()
