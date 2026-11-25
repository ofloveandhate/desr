# Scripts for the paper

Each script here reproduces one example of the paper that describes this library, *Nondimensionalization is more science than art* (Tanburn, Hendron, Maini, Amethyst, Dufresne and Harrington; [arXiv:2512.13455](https://arxiv.org/abs/2512.13455)), and prints the matrices, invariants, substitutions and reduced systems that the paper shows. A revision of the paper is in preparation, and these scripts follow the revision.

| script | example | source of the model |
|---|---|---|
| `michaelis_menten_two_variable.py` | Michaelis–Menten in two variables | Michaelis & Menten (1913) |
| `michaelis_menten_initial_conditions.py` | the same, with initial conditions and the Michaelis constant, brought to Murray's form | Murray, *Mathematical Biology I*, (6.20) |
| `michaelis_menten_four_variable.py` | Michaelis–Menten in four variables, with a choice of invariants imposed | Michaelis & Menten (1913) |
| `compartment_model.py` | a two-compartment model of a vaccine injection | DiStefano (2015), Example 13.6; Weissburg et al. (1995) |
| `segel_slemrod.py` | Michaelis–Menten on two timescales, and alternative nondimensionalisations | Segel & Slemrod (1989) |
| `cell_cycle.py` | a cell cycle control network | Sible & Tyson (2007), Figure 2 |

## Running them

From the root of the repository, with desr installed:

```bash
python paper/cell_cycle.py
```

Two of the scripts, `michaelis_menten_two_variable.py` and `segel_slemrod.py`, have no code of their own. Those examples are worked out on pages of the documentation, `docs/source/examples/mm_basic.rst` and `segel_slemrod.rst`, and the scripts run the computations on those pages in order and print each with its result. `_from_page.py` is what does that.

## Figures

A script that has figures in the paper draws them when it is run, into `figures/`, each as a PDF for the paper and a PNG for looking at. The PDF is written without a creation date, so the same figure gives the same file. `figures/` is not kept in the repository; run the script, or `regenerate_expected.py`, to make them.

| script | figures |
|---|---|
| `michaelis_menten_initial_conditions.py` | `michaelis_menten_initial_conditions_flowchart`, `michaelis_menten_initial_conditions_symbol_map` |

The figures need matplotlib, which desr does not install: `pip install matplotlib`.

## The expected output

`expected/` holds what each script printed when it was last looked at. `tests/test_paper_scripts.py` runs every script and compares. If a change to the library alters something a script prints, the test fails and shows which lines moved, so a result that is in the paper cannot change without anyone noticing.

When a change is intended:

```bash
python paper/regenerate_expected.py
```

Then read the difference it made to `expected/` as you would any other change, check it against the paper, and commit both.

## Which version

The output in `expected/` is that of the library at the commit it was written with. To reproduce the paper as published, check out the version of desr the paper cites.
