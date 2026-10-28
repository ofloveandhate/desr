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
python -m pytest tests/
```
Tests live in the top-level `tests/` directory, split by subject, and are not shipped in
the wheel (`setup.py` packages `desr` only). `tests/test_edge_cases.py` holds adversarial
systems; anything marked `unittest.expectedFailure` there is a known, unfixed bug and will
turn into a failure when it is fixed.

numpy and scipy are optional for *using* desr but required for *testing* it: the suite
imports them unconditionally rather than skipping, since a suite that skips a module when
its dependency is missing reports a success it has not earned.

**Run an example:**
```bash
python examples/michaelis_menten.py
```

## Architecture

The reduction pipeline flows: `ODESystem` → `ODETranslator` → reduced `ODESystem`.

The public names are exported from the package, so a whole workflow is one import:
`from desr import ODESystem, ODETranslator`. `desr/__init__.py` holds the list in
`__all__`; submodules stay reachable as attributes (`desr.matrix_normal_forms`) for
everything finer-grained. The numeric names (`NumericTranslation`, `scheme_of`,
`InsufficientKnownValues`, `ConflictingKnownValues`) resolve through a module-level
`__getattr__` rather than an eager import, because `desr.numerics` needs numpy and
`import desr` must not. `tests/test_namespace.py` pins both halves of that.

### Core modules

**`desr/ode_system.py` — `ODESystem`**
Represents a system of ODEs as parallel tuples of `variables` and `derivatives`. The independent variable (default `t`) is always included with derivative `1`. Key entry points:
- `ODESystem.from_tex(tex_str)` — parse LaTeX ODE system
- `ODESystem.from_dict(deriv_dict)` — construct from `{symbol: derivative_expr}` dict
- `exponent_matrix()` — build the exponent (power) matrix used by the scaling algorithm
- `maximal_scaling_matrix()` — compute the scaling matrix `A` from the system
- `reorder_variables(order)` — variable order affects which parameters get normalized out
- `update_initial_conditions(...)`, `add_constraint(lhs, rhs)` — enrich the system before reduction

**`desr/ode_translation.py` — `ODETranslator`**
The reduction engine. Holds a `scaling_matrix` (rows = scaling dimensions, columns = variables) and its column Hermite normal form multiplier `V`. Key methods:
- `ODETranslator.from_ode_system(system)` — auto-construct from a system
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
- **Parameters versus dependent variables.** A variable with no equation is a *parameter*
  (`system.parameters`); it is stored as `None` in `system.derivatives`. A variable with an
  equation is a *dependent variable* (`system.dependent_variables`), and that includes a
  conserved variable whose equation says `0` — writing `dw/dt = 0` is a deliberate
  statement, not a way of declaring a constant. `parameters=[...]` at construction, or
  `system.make_parameter(name)` afterwards, discards such an equation and makes the
  variable a parameter. Say "parameter", never "constant"; `constant_variables`,
  `num_constants` and `non_constant_variables` were removed in 2.0 without aliases.
- `repr` lists only variables that have equations, then names the parameters on a
  `parameters:` line. `to_tex` omits parameters entirely.
- **pytest is the primary test mechanism**, in the top-level `tests/`. Doctests are
  documentation and a backstop: they are string comparison, which is how `rational_section()`
  once printed its result and returned `None` with its doctest passing. Behaviour shown in
  `docs/source/examples/*.rst` should also be asserted in `tests/`, so a regression surfaces
  before anyone builds the docs.
- The scripts in `examples/` open plot windows and block; they cannot be run unattended.
- The naming scheme controls output variable names. It is a **dict** over the keys
  `parameters`, `independent`, `dependents` and `auxiliaries`; all are optional, and one left
  out takes the default for whichever reduction scheme runs — `desr`'s own `tau`/`nu`/`kappa`
  for the parameter scheme, Hubert & Labahn's `x` and `y` for the general and
  dependent-variable schemes, which also keep the system's own clock unless `independent` says
  otherwise. Whatever is set is honoured by all of them. The positional 3-tuple was removed in
  2.0: its order followed the variable order, which changed, so re-reading it would have
  silently re-pointed every entry.
- `exponent_matrix` was previously called `power_matrix` — watch for stale references in docs/comments.
