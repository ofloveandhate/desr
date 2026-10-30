#!/usr/bin/env python
"""Run every script in examples/ and report which ones fail.

Nothing else runs most of these.  Two are referenced by `.. plot::` directives and so are
exercised by the HTML docs build; the other six are exercised by nothing at all, which is
how a parameter-ordering bug survived until examples/linear_compartment.py tripped over it.

They can be run unattended: with the Agg backend `plt.show()` is a no-op.

Every script is run even after one fails, so a single breakage does not hide the rest, and
the exit status is non-zero if any failed.  `SKIP` names scripts that are known to fail and
are tracked by an issue -- a skipped script is reported, never silently passed over.
"""
import io
import os
import pathlib
import runpy
import sys
import traceback

os.environ.setdefault('MPLBACKEND', 'Agg')

try:
    import matplotlib
    matplotlib.use('Agg')
except ImportError:  # pragma: no cover - matplotlib is installed in CI
    pass

EXAMPLES = pathlib.Path(__file__).resolve().parents[2] / 'examples'

#: Scripts known to fail, each with the issue that tracks it.  Remove an entry when its
#: issue closes; a script that starts passing while still listed here is reported as a
#: surprise rather than quietly accepted.
SKIP = {
    'michaelis_menten.py': 'https://github.com/ofloveandhate/desr/issues/19',
}


def run(script):
    """Run one script, returning None on success or a one-line description of the failure."""
    captured = io.StringIO()
    stdout, sys.stdout = sys.stdout, captured
    try:
        runpy.run_path(str(script), run_name='__main__')
        return None
    except Exception as exc:
        first = str(exc).strip().splitlines()
        detail = first[0] if first else ''
        return '{}: {}'.format(type(exc).__name__, detail)
    finally:
        sys.stdout = stdout


def main():
    scripts = sorted(p for p in EXAMPLES.glob('*.py') if p.name != '__init__.py')
    if not scripts:
        print('no example scripts found in {}'.format(EXAMPLES), file=sys.stderr)
        return 1

    failed, skipped, unexpectedly_fine = [], [], []

    for script in scripts:
        failure = run(script)
        expected = SKIP.get(script.name)

        if failure and expected:
            skipped.append((script.name, failure, expected))
            print('SKIP {:42s} {}'.format(script.name, failure))
        elif failure:
            failed.append((script.name, failure))
            print('FAIL {:42s} {}'.format(script.name, failure))
            traceback.print_exc()
        elif expected:
            unexpectedly_fine.append(script.name)
            print('PASS {:42s} (listed in SKIP -- remove it and close {})'.format(
                script.name, expected))
        else:
            print('ok   {}'.format(script.name))

    print('\n{} ran, {} ok, {} failed, {} skipped'.format(
        len(scripts), len(scripts) - len(failed) - len(skipped), len(failed), len(skipped)))

    for name, _failure, issue in skipped:
        print('  skipped {} -- tracked by {}'.format(name, issue))

    if unexpectedly_fine:
        print('\n{} now pass(es) but are still listed in SKIP; remove the entry.'.format(
            ', '.join(unexpectedly_fine)), file=sys.stderr)
        return 1

    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
