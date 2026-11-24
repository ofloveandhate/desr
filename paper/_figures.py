'''
Where the paper's figures are written, and how.

A script that draws a figure calls :func:`save`.  The figures go to `paper/figures/`,
or to the directory named by the environment variable `DESR_PAPER_FIGURES`, which is
how the tests draw them somewhere else.

Not a script itself.
'''
import os
import pathlib

import matplotlib

matplotlib.use('Agg')

#: The default place.
FIGURES = pathlib.Path(__file__).resolve().parent / 'figures'


def directory():
    ''' The directory figures are written to, made if it is not there. '''
    chosen = pathlib.Path(os.environ.get('DESR_PAPER_FIGURES') or FIGURES)
    chosen.mkdir(parents=True, exist_ok=True)
    return chosen


def save(figure, name):
    '''
    Write `figure` as `name.pdf`, for the paper, and `name.png`, for looking at.

    The PDF is written without a creation date, so that drawing the same figure twice
    gives the same file and a difference in the file is a difference in the figure.

    Returns:
        tuple of pathlib.Path: the two files.
    '''
    import matplotlib.pyplot as plt

    target = directory()
    pdf, png = target / (name + '.pdf'), target / (name + '.png')
    figure.savefig(pdf, metadata={'CreationDate': None})
    figure.savefig(png, dpi=200)
    plt.close(figure)
    return pdf, png
