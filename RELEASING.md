# Releasing desr

The steps for cutting a release, in order. Nothing here happens on a merge: a release is made by pushing a tag.

## 1. Before the tag

1. **The version.** It is written in one place, `__version__` in `desr/__init__.py`. `setup.py` and `docs/source/conf.py` read it from there.
2. **The changelog.** `CHANGELOG.md` has a section for the version. Its last paragraph gives test counts; update them to what `main` runs now. The section is what goes into the GitHub release.
3. **The milestone.** No open issues on the milestone for this version.
4. **CI is green on `main`** for the commit to be tagged. The `tests` workflow runs the unit tests on Python 3.10 to 3.14, the module doctests, a bare install without the optional dependencies, the Sphinx doctests, the HTML build with warnings as errors, and coverage.

To run the same checks locally:

```bash
python -m pytest tests/
python -m pytest --doctest-modules desr/
python -m sphinx -b doctest docs/source docs/build/doctest
python -m sphinx -b html -W docs/source docs/build/html
```

## 2. Tag

Tags are `v` followed by the version: `v3.0.0`. Pushing the tag starts the release: section 5 says what it does.

```bash
git switch main
git pull
git tag -a v3.0.0 -m "desr 3.0.0"
git push origin v3.0.0
```

A pushed tag is public, and for a final version it ends in a PyPI release and a GitHub release that cannot be taken back. Check the commit before pushing the tag.

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

**A pushed tag publishes.** `.github/workflows/publish.yml` runs every check in `tests.yml`, builds both distributions, checks them with `twine`, checks that the tag is `v` followed by exactly the version in `desr/__init__.py`, and uploads:

| version | example tag | goes to |
|---|---|---|
| development or pre-release | `v3.0.0.dev1`, `v3.0.0rc1` | TestPyPI |
| final | `v3.0.0` | PyPI, after the required reviewer approves the `pypi` environment; then a GitHub release |

The GitHub release is titled `desr 3.0.0`, has the version's section of `CHANGELOG.md` as its notes, and has the wheel and the source distribution attached, the same files PyPI has. A tag whose final version has no section in `CHANGELOG.md` still reaches PyPI but gets no release, and the run fails at that step; make the release by hand with `gh release create`. A development version gets no GitHub release.

A tag that disagrees with the version publishes nothing. **A version number on PyPI is single-use**: a release can be deleted, but the same version can never be uploaded again.

### Trying a version on TestPyPI first

On the branch being released:

1. Set `__version__` to a development version, `3.0.0.dev1`, and commit.
2. Tag and push the tag: `git tag -a v3.0.0.dev1 -m "desr 3.0.0.dev1"` and `git push origin v3.0.0.dev1`.
3. When the workflow has run, read https://test.pypi.org/project/desr/ and install from there into a fresh environment. The dependencies are not on TestPyPI, so the real index is given as well:

   ```bash
   python -m venv try-desr
   try-desr/bin/pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ desr==3.0.0.dev1
   ```

4. If something needs fixing, fix it and repeat with `3.0.0.dev2`.
5. Set `__version__` back to `3.0.0`, merge, and tag `v3.0.0` on `main` as in section 2. That tag goes to PyPI.

The development tags can be deleted afterwards, locally and on GitHub; what was uploaded to TestPyPI stays there.

### One-time setup: trusted publishing

No token is stored. PyPI and TestPyPI each accept uploads from this repository's `publish.yml`, in a named GitHub environment.

1. On https://test.pypi.org, under Account settings, Publishing, add a pending publisher: project `desr`, owner `ofloveandhate`, repository `desr`, workflow `publish.yml`, environment `testpypi`.
2. On https://pypi.org, the same with environment `pypi`.
3. In the repository's Settings, Environments, create `testpypi` and `pypi`. Giving `pypi` a required reviewer makes every real upload wait for a click.

### By hand

If the workflow cannot be used, the same steps by hand, from the tagged commit:

```bash
pip install build twine
python -m build --outdir dist-3.0.0
python -m twine check dist-3.0.0/*
python -m twine upload --repository testpypi dist-3.0.0/*   # or, for real: python -m twine upload dist-3.0.0/*
```

## 6. Afterwards

1. Close the milestone.
2. Delete merged branches.
3. If work on the next version starts at once, bump `__version__` and open a new section in `CHANGELOG.md`.
