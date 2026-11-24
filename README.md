# desr

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17671203.svg)](https://doi.org/10.5281/zenodo.17671203)

*desr* is a package used for Differential Equation Symmetry Reduction and is particularly useful for reducing the number of parameters in dynamical systems. It implements algorithms outlined by Evelyne Hubert and George Labahn <sup>[1](#myfootnote1)</sup>. See the full documentation on [readthedocs](https://desr.readthedocs.io/en/latest/).

This package accompanies a paper "Nondimensionalisation is more science than art".  The Masters dissertation [Differential Algebra and Applications](http://tanbur.github.io/desr/dissertation/differential_algebra_and_applications.pdf) that this code originally accompanied places the algorithms into the theoretical framework of <em>differential algebraic geometry</em> and shows how to extend them to parameter reduction of arbitrary systems of partial differential equations, though this is not yet implemented.

<a name="myfootnote1">1</a>: Hubert, E., & Labahn, G. (2013). Scaling Invariants and Symmetry Reduction of Dynamical Systems. Foundations of Computational Mathematics, 13(4), 479–516. http://doi.org/10.1007/s10208-013-9165-9


## Getting Started

### Installing

From PyPI:

`$ pip install desr`

The figures also need matplotlib, which desr does not install:

`$ pip install matplotlib`

From a clone of the repository:

`$ pip install .`

## Running the tests

From a clone of the repository, with matplotlib and pytest installed:

`$ python -m pytest tests/`

The worked examples in the documentation are doctests too: `python -m sphinx -b doctest docs/source docs/build/doctest`.

## Built With

* [Sphinx](http://www.sphinx-doc.org/en/stable/) - Used to generate the docs.

## Contributing

Submissions for contribution are always welcome.

## Authors

* **Richard Tanburn** - *Initial work* - [tanbur](https://github.com/tanbur)
* **Silviana Amethyst** - *Revisions, bugfixes, maintenance, 2025-* - [ofloveandhate](https://github.com/ofloveandhate) 

## License

This project is licensed under the Apache 2.0 License - see the [LICENSE](https://github.com/ofloveandhate/desr/blob/main/LICENSE) file for details

## Acknowledgments

* Dr. Heather Harrington and Dr. Emilie Dufresne for their supervision of the dissertation.
