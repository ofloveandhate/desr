'''
The scripts in `paper/` print what their files in `paper/expected/` say they print.

Each script reproduces one example of the paper, and its printed output is what is
copied into the paper.  So a change to the library that alters that output alters
something in print, and this is where it is noticed: the test fails, naming the script
and showing the lines that moved.

When the change is intended, run `python paper/regenerate_expected.py`, read what it
did to `paper/expected/`, and commit both.
'''
import difflib
import pathlib
import sys
from unittest import TestCase

PAPER = pathlib.Path(__file__).resolve().parent.parent / 'paper'
sys.path.insert(0, str(PAPER))

from regenerate_expected import EXPECTED, output_of, scripts    # noqa: E402

#: The examples the paper has.  A script added or removed has to be added or removed
#: here too, so that the list is a decision and not an accident of what is in the folder.
EXAMPLES = ['cell_cycle', 'compartment_model', 'michaelis_menten_four_variable',
            'michaelis_menten_initial_conditions', 'michaelis_menten_two_variable',
            'segel_slemrod']


class TestThePaperScripts(TestCase):

    def test_there_is_a_script_and_an_expected_output_for_each_example(self):
        self.assertEqual(EXAMPLES, [script.stem for script in scripts()])
        self.assertEqual(EXAMPLES, sorted(path.stem for path in EXPECTED.glob('*.txt')))

    def test_each_script_prints_what_is_expected(self):
        for script in scripts():
            with self.subTest(script=script.name):
                expected = (EXPECTED / (script.stem + '.txt')).read_text()
                # Another hash seed than the one the expected output was written with,
                # so that output depending on the order a set iterates in fails here.
                printed = output_of(script, hash_seed='1')
                if printed != expected:
                    difference = '\n'.join(list(difflib.unified_diff(
                        expected.splitlines(), printed.splitlines(),
                        'expected', 'printed', lineterm='', n=1))[:60])
                    self.fail('{} no longer prints what paper/expected/{}.txt holds:\n{}'
                              .format(script.name, script.stem, difference))

    def test_an_expected_output_is_not_empty(self):
        for name in EXAMPLES:
            with self.subTest(example=name):
                self.assertGreater(len((EXPECTED / (name + '.txt')).read_text().splitlines()), 20)
