'''
desr -- Differential Equation Symmetry Reduction.

Reduce the number of parameters in a dynamical system by finding its scaling symmetries,
after Hubert & Labahn (2013).

The names a caller needs are here, so that a whole workflow is one import::

    from desr import ODESystem, ODETranslator

    system = ODESystem.from_equations(['dn/dt = r*n*(1 - n/k)'])
    translator = ODETranslator.from_ode_system(system)
    reduced = translator.translate(system)

The submodules remain available for the rest -- :mod:`desr.matrix_normal_forms` for the
integer matrix algorithms, :mod:`desr.sympy_helper` for the parsing convention,
:mod:`desr.tex_tools` for LaTeX.
'''
from . import ode_system
from . import matrix_normal_forms
from . import chemical_reaction_network
from . import ode_translation
from . import sympy_helper
from . import tex_tools

from .chemical_reaction_network import (ChemicalReactionNetwork, ChemicalSpecies, Complex,
                                        Reaction)
from .ode_system import ODESystem
from .ode_translation import ODETranslator, scale_action
from .sympy_helper import sympify_expression, to_symbol

# desr.numerics needs numpy and scipy, which are optional for using the rest of the
# library.  Importing it here would make them mandatory for `import desr`, so the numeric
# names are resolved on first access instead: `from desr import NumericTranslation` works,
# and an ImportError for numpy arrives only if you actually reach for the numeric part.
_LAZY = {
    'ConflictingKnownValues': 'numerics',
    'InsufficientKnownValues': 'numerics',
    'NumericTranslation': 'numerics',
    'Scheme': 'numerics',
    'numerics': None,
}


def __getattr__(name):
    module_name = _LAZY.get(name, False)
    if module_name is False:
        raise AttributeError('module {!r} has no attribute {!r}'.format(__name__, name))

    import importlib
    module = importlib.import_module('.{}'.format(module_name or name), __name__)
    return module if module_name is None else getattr(module, name)


def __dir__():
    return sorted(set(__all__) | set(globals()))


__all__ = [
    # Systems of ODEs
    'ODESystem',
    # Reduction
    'ODETranslator',
    'scale_action',
    # Chemical reaction networks
    'ChemicalReactionNetwork',
    'ChemicalSpecies',
    'Complex',
    'Reaction',
    # Numerical translation (needs numpy and scipy; imported on first access)
    'NumericTranslation',
    'InsufficientKnownValues',
    'ConflictingKnownValues',
    'Scheme',
    # The one place text becomes sympy
    'sympify_expression',
    'to_symbol',
    # Submodules
    'chemical_reaction_network',
    'matrix_normal_forms',
    'numerics',
    'ode_system',
    'ode_translation',
    'sympy_helper',
    'tex_tools',
]
