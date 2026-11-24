"""
Write the output of every script here into `expected/`, and its figures into `figures/`.

Run this after a change that is meant to alter what a script prints, read the
difference it makes to `expected/` as you would any other change, and commit both.
`tests/test_paper_scripts.py` fails whenever a script's output and its file in
`expected/` differ, so a change to the library that moves a result in the paper is
seen when it is made.

    python paper/regenerate_expected.py
"""
import os
import pathlib
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
EXPECTED = HERE / 'expected'
FIGURES = HERE / 'figures'

#: The figures each script draws, by the script's name.  A script not listed draws none.
DRAWN = {
    'michaelis_menten_initial_conditions': ['michaelis_menten_initial_conditions_flowchart'],
}


def scripts():
    ''' The scripts that reproduce something in the paper, in name order. '''
    return sorted(path for path in HERE.glob('*.py')
                  if not path.name.startswith('_') and path.name != pathlib.Path(__file__).name)


def output_of(script, hash_seed='0', figures=FIGURES):
    '''
    What a script prints.  Any figures it draws are written into `figures`.

    Run in a process of its own, as a reader would run it.  The hash seed is fixed so
    that a script which prints a set in whatever order it happens to iterate in is
    caught: the test runs with another seed and must get the same text.
    '''
    environment = dict(os.environ, PYTHONHASHSEED=hash_seed, MPLBACKEND='Agg',
                       DESR_PAPER_FIGURES=str(figures))
    result = subprocess.run([sys.executable, str(script)], capture_output=True, text=True,
                            cwd=str(HERE.parent), env=environment)
    if result.returncode != 0:
        raise RuntimeError('{} failed:\n{}'.format(script.name, result.stderr))
    return result.stdout


if __name__ == '__main__':
    EXPECTED.mkdir(exist_ok=True)
    for script in scripts():
        text = output_of(script)
        (EXPECTED / (script.stem + '.txt')).write_text(text)
        print('{:45s} {:4d} lines'.format(script.name, len(text.splitlines())))
