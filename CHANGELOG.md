# Changelog

Notes for releases before 3.0.0 are on the
[releases page](https://github.com/ofloveandhate/desr/releases).

## 3.0.0

Scaling symmetry reduction of ODE systems, after Hubert & Labahn (2013).

Not source-compatible with 2.x.  A general-scheme reduction now produces a
different system: two steps that 2.x left undone are done by default, and
each can be turned off.

### Reductions

  - `translate_general` normalises the parameter auxiliaries and identifies a
    clock invariant with the independent variable, both by default, both
    controlled by `normalise_parameters` and `identify_indep_var`.  Passing
    `False` to each reproduces Hubert & Labahn's reduction exactly.  The
    options reach `translate` and `translate_dep_var` too.
  - An auxiliary that is a monomial in the original parameters alone is
    constant, so it can be set to 1.  That is a choice of units, not an
    assumption, and the choice is recorded in the caller's own symbols on
    the reduced system's `normalisations`.
  - A reduced system reports as parameters only the ones its equations use.
  - The reduced clock is called `t` only when it denotes `t`.  Where an
    invariant has been identified with it, the clock is a rescaled time and
    is named `tau`, as the parameter scheme already named the same quantity.
    A reduction that identifies nothing keeps `t` — including every general
    reduction that normalises nothing, and Hubert & Labahn 6.6, whose clock
    invariants all involve a dependent variable.
  - `substitutions()` writes each original symbol in the reduced system's own
    symbols, under every scheme.  It applies both of a reduction's steps, so
    no image names a column the reduction spent.
  - A system with no scaling symmetry has a scaling matrix with no rows, as
    Hubert & Labahn's Proposition 5.1 gives it, and its reduction is the
    identity.  It used to be given a row of zeros, counted as one symmetry,
    and reduced to equations that were not true.  Reducing a reduced system
    therefore returns it unchanged.
  - A scaling matrix passed to `ODETranslator` must have full row rank, and
    one that does not is refused.
  - A reduced symbol that denotes the original symbol of the same name is
    not a name clash.
  - An invariant that is one of the original symbols keeps that symbol's
    name: a dimensionless parameter, a clock the reduction does not rescale,
    and every symbol of a system with no scaling symmetry, which therefore
    reduces to itself.  An auxiliary is always given a new name.  Setting
    `rename_unchanged_symbols` to `True` in the naming scheme gives every
    invariant a generated name instead.
  - `translate_dep_var` accepts a scaling matrix written over the symbols
    other than the independent variable, as the same scaling with a zero
    column there.  This raised `IndexError` before.  A matrix of any other
    width is refused with a message that says so.
  - `expected_variable_count` is right for the dependent-variable scheme
    when an auxiliary is a monomial in the parameters alone and is
    normalised away.  It returned `n` whatever the reduction did.
  - `ODESystem.from_dict` and `from_equations` accept `initial_conditions`.
    They raised `TypeError` for any value of it.  `from_dict` no longer
    adds the independent variable to the dict it is given.
  - `smith_normal_form` is right for every integer matrix.  It failed an
    internal assertion on about one small matrix in ten.
  - The general scheme takes an invariant `k*t`, with `k` a monomial in the
    parameters, as its clock `tau` whenever there is one, and the reduced
    system then has no `t`.  This is a change of independent variable and
    no longer waits on the normalisation making `k` equal to 1.  With it
    the general scheme gives the parameter scheme's system wherever that
    scheme applies.  `identify_indep_var=False` keeps `t`.
  - An invariant built from the parameters alone is named from the naming
    scheme's `parameters` entry, `kappa` by default, under every scheme.
    The general and dependent-variable schemes used to call it `y`.
  - `translate(..., isolate_clock=True)` and `with_isolated_clock()`
    choose a basis of the invariants in which the independent variable
    appears in only one.  Where that one is `k*t`, an autonomous system
    reduces to an autonomous one.  Off by default.

### A translator remembers its system

  - `ODETranslator.from_ode_system` keeps the system, and the methods that
    report on a reduction no longer take one: `scheme`, `invariant_names`,
    `auxiliary_names`, `column_names`, `column_of`, `named_symbols`,
    `symbol_kinds`, `variable_map`, `parameter_invariants` and
    `dependent_invariants` lost their `system` argument.  A translator built
    from a bare matrix says so rather than guessing.
  - `translate` no longer takes `naming_scheme`.  Set it at construction or
    with `set_naming_scheme`, so one reduction cannot be named two ways.
  - A reduced system carries what it needs to describe itself:
    `named_symbols`, `symbol_kinds`, `normalisations`, `auxiliary_names`,
    `forward_monomials` and `original_system`.

### Dimensions

  - `ODESystem.dimensions_of()` gives each symbol's structural dimension,
    derived from the scaling matrix with nothing declared.  They are written
    `[t]`, after the symbol that names the direction.
  - `set_dimension`, `set_dimensions` and `physical_dimensions_of()` are the
    overlay: declare `[t] = time` and the rest follows.  Naming changes no
    matrix.  A dimension no declaration reaches is `None`, not 1.
  - `scaling_symmetries()` reports a symmetry as an action —
    `t -> lamda*t, k_2 -> k_2/lamda` — checked against the equations.
  - `unit_table()` is the reduction's symbols with what each denotes, what it
    measures and what it is measured in.  `num_units_and_dimensions()` counts
    them: a reduction coins `r` units measuring however many dimensions.
  - Also `dimension_of`, `dimension_matrix`, `symbols_by_dimension`,
    `dimensionless_expressions`, `exponent_matrix_hnf` and
    `scaling_symmetry_rows`.

### Figures

  - `desr.figures.flowchart_figure` draws a reduction as nine numbered boxes
    of matrices, page-text wide and as tall as it needs.
  - `desr.figures.symbol_map_figure` draws it as a chain of five systems and
    the four maps between them.
  - `desr.symbol_graph` is the exponent graph both are built on, with
    `ODETranslator.symbol_graph` and `draw_symbol_graph` as the way in.
  - Neither figure computes anything.  Every matrix, monomial, dimension and
    substitution on either page comes from a getter, so a figure drawn for a
    model desr has never seen is as right as the library is.

### Renames

  - `NumericTranslation` is `NumericTranslator`, and is built from a reduced
    system rather than from a system and a translator.
    `ODESystem.numeric_translator()` is the short way to one.
  - `symbols_that_must_not_be_zero` and
    `reduced_symbols_that_must_not_be_zero` are
    `forward_symbols_that_must_not_be_zero` and
    `reverse_symbols_that_must_not_be_zero`, named for the direction rather
    than for the symbol set.  The old `reduced_` name claimed the reduced
    system's symbols and returned the column names, which include the columns
    a reduction spends.
  - `column_monomials` is `forward_monomials`, beside `reverse_monomials`.
  - `indep_var_invariant` is `indep_var_invariants`, and returns all of them.
  - `units_and_dimensions` is `num_units_and_dimensions`, because it counts.
  - `r` keeps working, and `num_scaling_symmetries` says what it is.

### Packaging

  - Two optional extras.  `pip install desr[numerics]` adds numpy and scipy
    for `desr.numerics`; `pip install desr[graphs]` adds networkx and
    matplotlib for the symbol graph and the figures.  `import desr` needs
    neither, and CI proves it on a bare install.
  - `desr.figures` ships with the package.

Tested on 693 unit tests over 1016 subtests, 131 module doctests, 1369
documentation doctests, an HTML documentation build with warnings as errors,
the example scripts, and both figures rendered for four systems each.

See pull requests #40, #41 and #44 for what this is made of.
