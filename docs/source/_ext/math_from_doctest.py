'''
Typeset TeX that a page's own doctests compute.

    .. math-from-doctest:: reduced.to_tex()

runs the ``>>>`` examples of the page above the directive, in a namespace prepared by
``doctest_global_setup`` as the doctest builder prepares it, then evaluates the argument
there and sets the string it gives as an ``align`` block.  The page shows the system
the doctests check, so it cannot go stale the way a block copied from their output can.

The examples run during every build that reads the page, the HTML build included, so
the page costs its doctests' time.  Their output is not compared to anything here: that
is the doctest builder's job.  A failing example stops the build, naming the page and
the line.
'''
import contextlib
import doctest
import io

from docutils import nodes
from docutils.statemachine import StringList
from sphinx.util.docutils import SphinxDirective


class MathFromDoctest(SphinxDirective):
    required_arguments = 1
    final_argument_whitespace = True
    has_content = False

    def run(self):
        path = self.env.doc2path(self.env.docname)
        with open(path, encoding='utf-8') as handle:
            lines = handle.read().splitlines()
        above = '\n'.join(lines[:self.lineno - 1])

        namespace = {}
        exec(self.config.doctest_global_setup, namespace)
        with contextlib.redirect_stdout(io.StringIO()):
            for example in doctest.DocTestParser().get_examples(above):
                where = '{}, line {}'.format(path, example.lineno + 1)
                try:
                    exec(compile(example.source, where, 'exec'), namespace)
                except Exception as error:
                    raise self.error('The example at {} failed: {!r}'.format(where, error))
            tex = eval(self.arguments[0], namespace)

        if not isinstance(tex, str):
            raise self.error('{} gave a {}, not TeX.'.format(
                self.arguments[0], type(tex).__name__))

        block = ['.. math::', '    :nowrap:', '', '    \\begin{align}']
        block += ['    ' + line for line in tex.splitlines()]
        block += ['    \\end{align}']
        container = nodes.container()
        self.state.nested_parse(StringList(block, source=path), self.content_offset, container)
        return container.children


def setup(app):
    app.add_directive('math-from-doctest', MathFromDoctest)
    return {'parallel_read_safe': True, 'parallel_write_safe': True}
