# setuptools, not distutils.core.  distutils was removed from the standard library in
# Python 3.12 (PEP 632) and only resolves today because setuptools shims it, so a bare
# 3.12+ interpreter cannot run a distutils-based setup.py at all.  And distutils.core.setup
# ignores `python_requires` and `extras_require` -- they are setuptools features -- so
# declaring them above a distutils import would have had no effect whatsoever.
from setuptools import setup

import ast
import pathlib


def read_version():
    """The version from desr/__init__.py, which is the one place it is written.

    Parsed rather than imported: importing the package here would require sympy and
    python-flint to be installed before the install that provides them.
    """
    source = pathlib.Path(__file__).parent / 'desr' / '__init__.py'
    tree = ast.parse(source.read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if getattr(target, 'id', None) == '__version__':
                    return ast.literal_eval(node.value)
    raise RuntimeError('no __version__ in {}'.format(source))


setup(
    name='desr',
    version=read_version(),
    packages=['desr', 'desr.figures'],
    # Both authors by name.  One contact address, the maintainer's: an index shows it as
    # where to write, and that is the person who answers.
    author='Richard Tanburn, Silviana Amethyst',
    maintainer='Silviana Amethyst',
    maintainer_email='amethyst@mpi-cbg.de',
    description='Simplify ordinary differential equations by finding scaling symmetries.',
    # An SPDX expression (PEP 639), which is what an index files the licence under.  It
    # takes the place of a `License ::` classifier.
    license_expression='Apache-2.0',
    license_files=['LICENSE'],
    # With the encoding named: left out, it follows the locale of the machine that builds.
    long_description=(pathlib.Path(__file__).parent / 'README.md').read_text(encoding='utf-8'),
    # Without this the Markdown README is served as plain text by any index that renders
    # long_description, since the default is text/plain.
    long_description_content_type='text/markdown',
    # No extras.  Someone who finds a part missing installs it and runs again; nobody goes
    # back to reinstall desr with an extra.  So numpy and scipy (desr.numerics) and networkx
    # (desr.symbol_graph) are required.  matplotlib, which only the drawing needs and which
    # is the heaviest of them, is not: drawing without it stops at "No module named
    # 'matplotlib'", which says what to install.
    install_requires=[
        'sympy',
        'python-flint >= 0.7',
        'numpy',
        'scipy',
        'networkx',
    ],
    # 3.10 because that is python-flint's own floor, not an independent judgement of ours.
    # python-flint is the sole Hermite and Smith normal form backend and is not optional,
    # so desr cannot install anywhere it cannot; following its support window means our
    # claim never drifts from the truth, and when python-flint raises its floor the CI
    # matrix fails and tells us rather than our having to remember.
    #
    # Python 3.9 reached end of life on 2025-10-31 and python-flint requires >= 3.10 in any
    # case, so it is not a candidate.  CI tests 3.10 through 3.14.
    python_requires='>=3.10',
    # The DOI first because it is the citable identifier and does not depend on which
    # forge hosts the code.  Then the repository where the work happens, and last the
    # institutional copy, labelled as a mirror: it is downstream, and an issue or merge
    # request opened there is not read.
    # What an index sorts and filters by.  The Python versions are the ones CI tests.
    classifiers=[
        # Beta while the interface is still being worked on.
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3.14',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
    keywords=['nondimensionalization', 'scaling symmetry', 'ordinary differential equations',
              'symmetry reduction', 'Hermite normal form'],
    project_urls={
        "DOI": "https://doi.org/10.5281/zenodo.17671203",
        "Source": "https://github.com/ofloveandhate/desr",
        "Issues": "https://github.com/ofloveandhate/desr/issues",
        "Documentation": "https://desr.readthedocs.io",
        "Institutional mirror": "https://git.mpi-cbg.de/ASB/desr",
    }
)
