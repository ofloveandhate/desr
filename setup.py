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
    packages=['desr',],
    author='Richard Tanburn, Silviana Amethyst',
    author_email='richard.tanburn@gmail.com, amethyst@mpi-cbg.de',
    description='Simplify ordinary differential equations by finding scaling symmetries.',
    license='Apache License Version 2.0',
    long_description=open('README.md').read(),
    # Without this the Markdown README is served as plain text by any index that renders
    # long_description, since the default is text/plain.
    long_description_content_type='text/markdown',
    install_requires=[
        'sympy',
        'python-flint >= 0.7',
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
    extras_require={
        # desr itself is pure sympy.  Only desr.numerics needs these, and `import desr`
        # must keep working without them -- tests/test_namespace.py pins that, and CI
        # proves it on a bare install.
        'numerics': ['numpy', 'scipy'],
    },
    project_urls={
        "DOI": "https://doi.org/10.5281/zenodo.17671203",
        "repo_MPI": "https://git.mpi-cbg.de/ASB/desr"
    }
)




