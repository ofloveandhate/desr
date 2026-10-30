'''
What `import desr` gives you.

A whole workflow should be one import.  The constraint pulling the other way is that
numpy and scipy are optional for *using* desr -- only `desr.numerics` needs them -- so
`import desr` must not require them.  Both halves are pinned here, because they are easy
to break by adding one convenient import at the top of a module.
'''
import subprocess
import sys
from unittest import TestCase

import desr


class TestPublicNames(TestCase):
    ''' The names a caller reaches for are on the package. '''

    def test_a_whole_workflow_is_one_import(self):
        from desr import ODESystem, ODETranslator
        system = ODESystem.from_equations(['dn/dt = r*n*(1 - n/k)'])
        translator = ODETranslator.from_ode_system(system)
        reduced = translator.translate(system)
        self.assertTrue(reduced.symbols)

    def test_the_chemical_reaction_network_names_are_exported(self):
        from desr import ChemicalReactionNetwork, ChemicalSpecies, Complex, Reaction
        for name in (ChemicalReactionNetwork, ChemicalSpecies, Complex, Reaction):
            self.assertTrue(callable(name))

    def test_the_parsing_convention_is_exported(self):
        from desr import sympify_expression, to_symbol
        self.assertEqual(to_symbol('S'), sympify_expression('S'))

    def test_the_numeric_names_resolve_through_the_package(self):
        from desr import (ConflictingKnownValues, InsufficientKnownValues,
                          NumericTranslation)
        self.assertTrue(issubclass(InsufficientKnownValues, ValueError))
        self.assertTrue(issubclass(ConflictingKnownValues, ValueError))
        self.assertTrue(callable(NumericTranslation))

    def test_the_reduction_scheme_is_eager_because_it_needs_no_numpy(self):
        '''
        `ReductionScheme` names the three reductions, which is symbolic business, so it
        lives in `ode_translation` and imports eagerly.  It was in `numerics` until
        `ODETranslator` needed to return one, which would have dragged numpy into the core.
        '''
        from desr import ReductionScheme
        self.assertEqual(ReductionScheme.GENERAL, 'general')
        self.assertNotIn('ReductionScheme', desr._LAZY)
        # Still reachable where it used to be, since that is where it was introduced.
        from desr.numerics import ReductionScheme as same
        self.assertIs(ReductionScheme, same)

    def test_submodules_are_reachable_as_attributes(self):
        for name in ('chemical_reaction_network', 'matrix_normal_forms', 'numerics',
                     'ode_system', 'ode_translation', 'sympy_helper', 'tex_tools'):
            with self.subTest(module=name):
                self.assertTrue(hasattr(desr, name))

    def test_everything_promised_in_all_is_actually_there(self):
        for name in desr.__all__:
            with self.subTest(name=name):
                self.assertTrue(hasattr(desr, name))

    def test_an_unknown_name_still_raises_attribute_error(self):
        # The lazy __getattr__ must not turn a typo into an ImportError or a silent None.
        with self.assertRaises(AttributeError):
            desr.no_such_thing


class TestOptionalDependencies(TestCase):
    '''
    numpy and scipy are optional for using desr.  They are installed here -- the test
    suite needs them -- so asking whether `import desr` pulled them in has to happen in a
    fresh interpreter.

    `desr.numerics` needs numpy at import; scipy is needed only for recovering auxiliary
    variables that vary along the solution, and is imported where it is used.
    '''

    @staticmethod
    def modules_after(statement):
        '''Which of numpy and scipy are loaded after running `statement` in a fresh python.'''
        script = (statement + '\n'
                  'import sys\n'
                  'print(",".join(sorted(n for n in ("numpy", "scipy")'
                  ' if n in sys.modules)))\n')
        result = subprocess.run([sys.executable, '-I', '-c', script],
                                capture_output=True, text=True, check=True)
        return result.stdout.strip()

    def test_importing_desr_does_not_require_numpy_or_scipy(self):
        self.assertEqual(self.modules_after('import desr'), '')

    def test_reaching_for_the_numeric_part_imports_numpy(self):
        # The lazy attribute must actually deliver the module, not merely stay quiet.
        self.assertEqual(self.modules_after('from desr import NumericTranslation'),
                         'numpy')

    def test_scipy_stays_unimported_until_a_solver_is_needed(self):
        # desr.numerics needs numpy to exist at all, but scipy only for recovering
        # auxiliary variables that vary along the solution, so it is imported at the point
        # of use with its own ImportError message.
        self.assertNotIn('scipy', self.modules_after('import desr.numerics'))
