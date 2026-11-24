'''
Which index a pushed tag publishes to, after checking that the tag is the version.

Usage: python .github/scripts/release_index.py <tag>

Prints `index=testpypi` or `index=pypi`, for `$GITHUB_OUTPUT`.  A pre-release or
development version (`3.0.0.dev1`, `3.0.0rc1`) goes to TestPyPI, where it can be looked
at and installed before anything is published for real.  A final version goes to PyPI.

Exits with an error, and publishes nothing, if the tag is not `v` followed by exactly the
version written in `desr/__init__.py`.  A version number on PyPI can be used once, so a
tag that disagrees with the package would spend a number on the wrong build.
'''
import ast
import pathlib
import sys

from packaging.version import InvalidVersion, Version

ROOT = pathlib.Path(__file__).resolve().parents[2]


def package_version(root=ROOT):
    ''' `__version__` from `desr/__init__.py`, parsed as `setup.py` parses it. '''
    source = root / 'desr' / '__init__.py'
    for node in ast.parse(source.read_text()).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if getattr(target, 'id', None) == '__version__':
                    return ast.literal_eval(node.value)
    raise RuntimeError('no __version__ in {}'.format(source))


def release_index(tag, version):
    '''
    `'testpypi'` or `'pypi'` for a tag and the package's version.

    Raises:
        ValueError: If the tag is not `v` followed by exactly the version.
    '''
    if not tag.startswith('v'):
        raise ValueError('the tag {!r} does not start with v'.format(tag))
    try:
        tagged = Version(tag[1:])
    except InvalidVersion:
        raise ValueError('the tag {!r} is not v and a version'.format(tag))
    if tag[1:] != version or tagged != Version(version):
        raise ValueError('the tag {!r} is not the package version {!r}'.format(tag, version))
    return 'testpypi' if tagged.is_prerelease else 'pypi'


if __name__ == '__main__':
    try:
        print('index={}'.format(release_index(sys.argv[1], package_version())))
    except (IndexError, ValueError) as error:
        sys.exit('not publishing: {}'.format(error))
