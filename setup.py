from distutils.core import setup

setup(
    name='desr',
    version='1.0.2',
    packages=['desr',],
    author='Richard Tanburn, Silviana Amethyst',
    author_email='richard.tanburn@gmail.com, amethyst@mpi-cbg.de',
    description='Simplify ordinary differential equations by finding scaling symmetries.',
    license='Apache License Version 2.0',
    long_description=open('README.md').read(),
    install_requires=[
        'sympy',
    ],
)