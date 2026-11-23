# Releasing desr

The steps for cutting a release, in order. Nothing here happens on a merge: a release is made by pushing a tag.

## 1. Before the tag

1. **The version.** It is written in one place, `__version__` in `desr/__init__.py`. `setup.py` and `docs/source/conf.py` read it from there.
2. **The changelog.** `CHANGELOG.md` has a section for the version. Its last paragraph gives test counts; update them to what `main` runs now. The section is what goes into the GitHub release.
3. **The milestone.** No open issues on the milestone for this version.
4. **CI is green on `main`** for the commit to be tagged. The `tests` workflow runs the unit tests on Python 3.10 to 3.14, the module doctests, a bare install without the optional dependencies, the Sphinx doctests, the HTML build with warnings as errors, the example scripts, and coverage.

To run the same checks locally:

```bash
python -m pytest tests/
python -m pytest --doctest-modules desr/
python -m sphinx -b doctest docs/source docs/build/doctest
python -m sphinx -b html -W docs/source docs/build/html
python .github/scripts/run_examples.py
```

## 2. Tag and release on GitHub

Tags are `v` followed by the version: `v3.0.0`.

```bash
git switch main
git pull
git tag -a v3.0.0 -m "desr 3.0.0"
git push origin v3.0.0
```

Then create the GitHub release from that tag. Title: `desr 3.0.0`. Body: the version's section of `CHANGELOG.md`.

```bash
gh release create v3.0.0 --title "desr 3.0.0" --notes-file <file holding the changelog section>
```

A pushed tag and a published release are public. Check the commit before pushing the tag.

## 3. Documentation

Read the Docs builds from `.readthedocs.yaml`, and builds a version for a new tag on its own. After pushing the tag, check there that the version appears and that its build passed.

## 4. Zenodo

The DOI on the README, `10.5281/zenodo.17671203`, is the concept record: it means desr rather than one version of it. Each release adds a version to that record.

The deposit is made by hand. The automated job in `.gitlab-ci.yml` is not working at this time.

From the repository root, with a zenodo.org token:

```bash
pip install gitlab2zenodo==0.0b8
git archive --format zip --output 3.0.0.zip v3.0.0
g2z-send -i 17671203 -t <zenodo.org token> -m .zenodo.json 3.0.0.zip
```

This leaves a **draft** on the record. Finish on zenodo.org: open the draft, read the metadata, and publish. **Publishing a Zenodo record cannot be undone.**

Two flags are deliberately absent from `g2z-send`: `-s` would send to the sandbox, a separate service with its own tokens and records, and `-p` would publish without anyone reading the draft.

`.zenodo.json` is known to be incomplete; see issue #32.

## 5. PyPI

desr is not on PyPI; see issue #35. There is no publishing step.

## 6. Afterwards

1. Close the milestone.
2. Delete merged branches.
3. If work on the next version starts at once, bump `__version__` and open a new section in `CHANGELOG.md`.
