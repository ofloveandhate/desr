Introduction
============

*desr* is a package used for Differental Equation Symmetry Reduction and is particularly useful for reducing the number of parameters in dynamical systems.


*desr* implements and extends algorithms originally outlined by Evelyne Hubert and George Labahn :cite:`Hubert2013c`.  It further implements reduction on systems with equality constraints and initial conditions, as described in :cite:`ArxivLinkHere`.


The Masters dissertation `Differential Algebra and Applications <http://tanbur.github.io/desr/dissertation/differential_algebra_and_applications.pdf>`_ that inspired this project places the algorithms into the theoretical framework of *differential algebraic geometry*.  The dissertation further shows how to extend them to parameter reduction of arbitrary systems of partial differential equations, though this is not yet implemented in this software.


Prerequisites
-------------

This package requires the Sympy package.

Installing
----------

To install, download the package and run:

.. code-block:: shell

    pip install .


Running the tests
-----------------

Doctests are included in most files. To run them, from top level run ``sphinx-build -b doctest docs/source/ build/``, or run the individual submodule of desr as a script. e.g. ``python -m doctest -v module.py``

Building the documentation
----------

- The documentation is hosted online at `readthedocs <https://desr.readthedocs.io/latest/>`_.
- `Sphinx <http://www.sphinx-doc.org/en/stable/>`_ - Used to generate the docs.


Contributing
------------

Pull requests and issues on `the GitHub repo <https://github.com/ofloveandhate/desr/>`_ are most certainly welcome.

Authors
-------

- **Richard Tanburn** - *Initial work*
- **Silviana Amethyst** - *Additional implementation, bugfixing, maintenance*

License
-------

This project is licensed under the Apache 2.0 License.

Acknowledgments
---------------

- Dr. Heather Harrington and Dr. Emilie Dufresne for their supervision of the dissertation.
- Thomas Close for writing his diophantine module, which is included in this package.