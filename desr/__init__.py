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

#: The one place the version is written.  `setup.py` reads it from here without importing
#: the package, and `docs/source/conf.py` reads it too, so a bump happens once and cannot
#: leave the three of them disagreeing -- which is what had happened: `setup.py` said
#: 1.1.0 while the docs said version 1, release 0.
__version__ = '3.0.0'

from . import ode_system
from . import matrix_normal_forms
from . import chemical_reaction_network
from . import ode_translation
from . import sympy_helper
from . import tex_tools

from .chemical_reaction_network import (ChemicalReactionNetwork, ChemicalSpecies, Complex,
                                        Reaction)
from .ode_system import ConflictingParameterValues, ODESystem
from .ode_translation import ODETranslator, ReductionScheme, scale_action
from .sympy_helper import sympify_expression, to_symbol
from .tex_tools import TexParseError

# desr.numerics loads numpy, and the reduction itself does not need it.  Importing it here
# would make every `import desr` pay for numpy, so the numeric names are resolved on first
# access instead: `from desr import NumericTranslator` works, and numpy is loaded only when
# something reaches for the numeric part.
_LAZY = {
    'ConflictingKnownValues': 'numerics',
    'InsufficientKnownValues': 'numerics',
    'NumericTranslator': 'numerics',
    'numerics': None,
    # desr.symbol_graph loads networkx, and its drawing matplotlib.  Same bargain as
    # numerics: resolved on first access, so `import desr` loads neither.  The
    # submodule only, not its functions: the module and its main function share a name, so
    # a top-level `symbol_graph` would be ambiguous.  `ODETranslator.symbol_graph` is the
    # way in, and `from desr.symbol_graph import draw_symbol_graph` gets the drawing.
    'symbol_graph': None,
    # desr.figures draws with matplotlib, which is not a requirement of desr, and the
    # symbol map reaches networkx through desr.symbol_graph.  The same bargain again:
    # `from desr.figures import flowchart_figure` works where matplotlib is installed,
    # and `import desr` works where it is not.
    'figures': None,
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
    'ReductionScheme',
    'scale_action',
    # Chemical reaction networks
    'ChemicalReactionNetwork',
    'ChemicalSpecies',
    'Complex',
    'Reaction',
    # Numerical translation (needs numpy and scipy; imported on first access)
    'NumericTranslator',
    # Every exception the library defines, so that one can be caught without knowing
    # which module raises it.  The two numeric ones are imported on first access.
    'ConflictingParameterValues',
    'TexParseError',
    'InsufficientKnownValues',
    'ConflictingKnownValues',
    # The one place text becomes sympy
    'sympify_expression',
    'to_symbol',
    # Submodules
    'chemical_reaction_network',
    'matrix_normal_forms',
    'numerics',
    'ode_system',
    'ode_translation',
    'symbol_graph',
    'sympy_helper',
    'tex_tools',
]
