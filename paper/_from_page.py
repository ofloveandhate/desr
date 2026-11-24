'''
Print every computation on one of the documentation's example pages.

Two of the paper's examples have no script of their own, because the documentation page
is where they were worked out.  For those, the page is the source: this runs the
examples on it in order and prints each one with its result, so that what goes into the
paper comes from the same place a reader of the documentation sees.

Not a script itself.
'''
import doctest
import pathlib
import sys

PAGES = pathlib.Path(__file__).resolve().parent.parent / 'docs' / 'source' / 'examples'

#: What the documentation's doctests are given before they start: `docs/source/conf.py`.
SETUP = '''
import sympy
from desr.ode_system import *
from desr.ode_translation import *
from desr.matrix_normal_forms import *
from desr.chemical_reaction_network import *
from desr.tex_tools import *
from desr.sympy_helper import *
'''


def print_page(name):
    '''
    Run the examples of `docs/source/examples/<name>` and print each with its result.

    Args:
        name (str): The page's file name, such as ``'mm_basic.rst'``.
    '''
    text = (PAGES / name).read_text()
    names = {}
    exec(SETUP, names)
    for example in doctest.DocTestParser().get_examples(text, name):
        lines = example.source.rstrip('\n').split('\n')
        print('>>> ' + lines[0])
        for line in lines[1:]:
            print('... ' + line)
        if example.exc_msg is not None:
            # An example that shows a refusal: run it, and print what is raised.
            try:
                exec(compile(example.source, name, 'single'), names)
            except Exception as refused:
                print('{}: {}'.format(type(refused).__name__, refused))
            continue
        exec(compile(example.source, name, 'single'), names)
    sys.stdout.flush()
