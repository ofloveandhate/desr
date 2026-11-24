'''
A version's section of CHANGELOG.md, for the notes of its GitHub release.

Usage: python .github/scripts/changelog_section.py <version> > notes.md

The section is everything under the heading `## <version>`, up to the next `## `
heading, without the heading itself (the release has its own title).  Exits with an
error if there is no such section, so a release is never made with empty notes.
'''
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]


def changelog_section(version, text):
    '''
    The body of the section headed `## <version>` in `text`.

    Raises:
        ValueError: If there is no such section, or it is empty.
    '''
    lines = text.splitlines()
    heading = '## {}'.format(version)
    try:
        start = next(index for index, line in enumerate(lines) if line.strip() == heading)
    except StopIteration:
        raise ValueError('CHANGELOG.md has no section {!r}'.format(heading))
    end = next((index for index in range(start + 1, len(lines))
                if lines[index].startswith('## ')), len(lines))
    body = '\n'.join(lines[start + 1:end]).strip()
    if not body:
        raise ValueError('the section {!r} of CHANGELOG.md is empty'.format(heading))
    return body + '\n'


if __name__ == '__main__':
    try:
        sys.stdout.write(changelog_section(sys.argv[1],
                                           (ROOT / 'CHANGELOG.md').read_text(encoding='utf-8')))
    except (IndexError, ValueError) as error:
        sys.exit('no release notes: {}'.format(error))
