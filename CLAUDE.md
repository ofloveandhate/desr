# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

`desr` (Differential Equation Symmetry Reduction) is a Python library for reducing the number of parameters in dynamical systems by finding scaling symmetries. It implements algorithms from Hubert & Labahn (2013). The primary use case is nondimensionalization of ODE systems.

## Commands

**Install:**
```bash
pip install .
```

**Run doctests for a single module:**
```bash
python -m doctest -v desr/ode_system.py
```

**Run all doctests via Sphinx:**
```bash
cd docs && make doctest
```

**Run unit tests:**
```bash
python -m pytest desr/unittests.py
# or
python desr/unittests.py
```

**Run an example:**
```bash
python examples/michaelis_menten.py
```

## Architecture

The reduction pipeline flows: `ODESystem` → `ODETranslation` → reduced `ODESystem`.

### Core modules

**`desr/ode_system.py` — `ODESystem`**
Represents a system of ODEs as parallel tuples of `variables` and `derivatives`. The independent variable (default `t`) is always included with derivative `1`. Key entry points:
- `ODESystem.from_tex(tex_str)` — parse LaTeX ODE system
- `ODESystem.from_dict(deriv_dict)` — construct from `{symbol: derivative_expr}` dict
- `exponent_matrix()` — build the exponent (power) matrix used by the scaling algorithm
- `maximal_scaling_matrix()` — compute the scaling matrix `A` from the system
- `reorder_variables(order)` — variable order affects which parameters get normalized out
- `update_initial_conditions(...)`, `add_constraint(lhs, rhs)` — enrich the system before reduction

**`desr/ode_translation.py` — `ODETranslation`**
The reduction engine. Holds a `scaling_matrix` (rows = scaling dimensions, columns = variables) and its column Hermite normal form multiplier `V`. Key methods:
- `ODETranslation.from_ode_system(system)` — auto-construct from a system
- `translate(system)` — produce the reduced `ODESystem`
- `invariants()` — the scaling-invariant monomials (the new variables)
- `auxiliaries()` — the auxiliary variables that get normalized to 1
- `multiplier_add_columns(i, j, alpha)` — manually adjust which invariants are chosen
- `translate_parameter_substitutions(system)` — get substitution rules for original parameters

**`desr/matrix_normal_forms.py`**
Integer matrix algorithms: Hermite Normal Form (`hnf_row`, `hnf_col`, `normal_hnf_col`, `normal_hnf_row`, `hnf_row_lll`) and Smith Normal Form (`smf`). These are the computational backbone.

**`desr/chemical_reaction_network.py` — `ChemicalReactionNetwork`**
Builds `ODESystem` objects from chemical reaction networks specified via `ChemicalSpecies`, `Complex`, and `Reaction` objects.

**`desr/diophantine.py`**
Third-party module (Thomas Close) for solving Diophantine equations; used internally by `matrix_normal_forms`.

**`desr/sympy_helper.py`, `desr/tex_tools.py`**
Utilities: SymPy expression manipulation, LaTeX parsing/rendering.

### Key design notes

- Variables are ordered; the order matters because the last variables tend to be "normalized out" (set to 1) during reduction. Use `reorder_variables()` to control this.
- Doctests are the primary test mechanism — most functions have embedded `>>>` examples. The `unittests.py` file has `unittest.TestCase` classes for matrix normal form routines.
- The naming scheme `('tau', 'nu', 'kappa')` controls output variable names: index-0 = new independent variable, index-1 = pattern for dependent variables, index-2 = pattern for invariant constants.
- `exponent_matrix` was previously called `power_matrix` — watch for stale references in docs/comments.
