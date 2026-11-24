'''
The repository's own tooling, which the package does not ship: where the version is
written, and the two scripts `.github/workflows/publish.yml` runs.

These read files outside `desr/` and `tests/` -- `setup.py`, `docs/source/conf.py`,
`.github/scripts/`, `CHANGELOG.md` -- so this file is left out of the source
distribution (see `MANIFEST.in`), whose tests must run from the unpacked archive alone.

`release_index.py`: a pushed tag publishes only if it is the version, and a pre-release
goes to TestPyPI, never to PyPI.  `changelog_section.py`: a GitHub release's notes are its
version's section of CHANGELOG.md, and there is one.
'''
import pathlib
import sys
from unittest import TestCase

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / '.github' / 'scripts'))

from changelog_section import changelog_section              # noqa: E402
from release_index import package_version, release_index    # noqa: E402

import desr    # noqa: E402


class TestTheVersion(TestCase):

    def test_the_version_is_written_in_exactly_one_place(self):
        '''
        `desr.__version__` is the source; `setup.py` parses it and `docs/source/conf.py`
        imports it.  They had drifted before this -- setup.py said 1.1.0 while the docs
        said version 1, release 0 -- so the agreement is asserted rather than trusted.
        '''
        import ast
        import pathlib
        import re

        root = pathlib.Path(__file__).parent.parent
        self.assertRegex(desr.__version__, r'^\d+\.\d+')

        # setup.py must read it from the package, not carry its own literal.
        setup_source = (root / 'setup.py').read_text()
        self.assertIn('read_version()', setup_source)
        self.assertNotRegex(setup_source, r"version\s*=\s*['\"]\d")

        # And the parse setup.py performs must agree with the imported value.
        tree = ast.parse((root / 'desr' / '__init__.py').read_text())
        parsed = next(ast.literal_eval(node.value)
                      for node in tree.body if isinstance(node, ast.Assign)
                      for target in node.targets
                      if getattr(target, 'id', None) == '__version__')
        self.assertEqual(desr.__version__, parsed)

        # conf.py must derive both of its strings from the package.
        conf_source = (root / 'docs' / 'source' / 'conf.py').read_text()
        self.assertIn('_desr.__version__', conf_source)
        self.assertNotRegex(conf_source, r"^release\s*=\s*u?['\"]\d", )
        self.assertFalse(re.search(r"^version\s*=\s*u?['\"]\d", conf_source,
                                   flags=re.MULTILINE))


class TestReleaseIndex(TestCase):

    def test_a_final_version_goes_to_pypi(self):
        self.assertEqual('pypi', release_index('v3.0.0', '3.0.0'))
        self.assertEqual('pypi', release_index('v3.0.1.post1', '3.0.1.post1'))

    def test_a_pre_release_or_development_version_goes_to_testpypi(self):
        for version in ('3.0.0.dev1', '3.0.0a1', '3.0.0b2', '3.0.0rc1', '3.1.0.dev0'):
            with self.subTest(version=version):
                self.assertEqual('testpypi', release_index('v' + version, version))

    def test_a_tag_that_is_not_the_version_publishes_nothing(self):
        for tag, version in (('v3.0.1', '3.0.0'), ('v3.0.0', '3.0.0.dev1'),
                             ('3.0.0', '3.0.0'), ('v3.0', '3.0.0'),
                             ('v3.0.0.dev01', '3.0.0.dev1'), ('vnext', '3.0.0')):
            with self.subTest(tag=tag, version=version):
                with self.assertRaises(ValueError):
                    release_index(tag, version)

    def test_the_version_it_reads_is_the_packages(self):
        self.assertEqual(desr.__version__, package_version())


class TestChangelogSection(TestCase):
    '''The notes of a GitHub release are the version's section of CHANGELOG.md.'''

    TEXT = '''# Changelog

Intro.

## 3.1.0

### Added

  - a thing

## 3.0.0

Scaling symmetry reduction.

### Reductions

  - another thing
'''

    def test_the_section_runs_to_the_next_version(self):
        self.assertEqual('### Added\n\n  - a thing\n',
                         changelog_section('3.1.0', self.TEXT))

    def test_the_last_section_runs_to_the_end(self):
        self.assertEqual('Scaling symmetry reduction.\n\n### Reductions\n\n  - another thing\n',
                         changelog_section('3.0.0', self.TEXT))

    def test_a_missing_or_empty_section_is_refused(self):
        with self.assertRaises(ValueError):
            changelog_section('2.9.9', self.TEXT)
        with self.assertRaises(ValueError):
            changelog_section('3.0.0', '## 3.0.0\n\n## 2.0.0\n\nold\n')
        # A heading that merely starts with the version is a different version.
        with self.assertRaises(ValueError):
            changelog_section('3.0', self.TEXT)

    def test_the_version_being_released_has_its_section(self):
        # Whatever __version__ is, its final version has notes waiting in CHANGELOG.md.
        from packaging.version import Version
        changelog = (pathlib.Path(__file__).resolve().parents[1] / 'CHANGELOG.md').read_text(
            encoding='utf-8')
        self.assertTrue(changelog_section(Version(desr.__version__).base_version, changelog))
