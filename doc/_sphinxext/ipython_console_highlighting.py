"""reST directive for syntax-highlighting ipython interactive sessions.

XXX - See what improvements can be made based on the new (as of Sept 2009)
'pycon' lexer for the python console.  At the very least it will give better
highlighted tracebacks.
"""

#-----------------------------------------------------------------------------
# Needed modules

# Standard library
import re

# Third party
from pygments.lexer import Lexer, do_insertions
from pygments.lexers.agile import (PythonConsoleLexer, PythonLexer,
                                   PythonTracebackLexer)
from pygments.token import Comment, Generic

from sphinx import highlighting

#-----------------------------------------------------------------------------
# Global constants
line_re = re.compile('.*?\n')

#-----------------------------------------------------------------------------
# Code begins - classes and functions

class IPythonConsoleLexer(Lexer):
    """
    For IPython console output or doctests, such as:

    .. sourcecode:: ipython

      In [1]: a = 'foo'

      In [2]: a
      Out[2]: 'foo'

      In [3]: print a
      foo

      In [4]: 1 / 0

    Notes:

      - Tracebacks are not currently supported.

      - It assumes the default IPython prompts, not customized ones.
    """

    name = 'IPython console session'
    aliases = ['ipython']
    mimetypes = ['text/x-ipython-console']
    input_prompt = re.compile("(In \[[0-9]+\]: )|(   \.\.\.+: )")
    output_prompt = re.compile("Out\[[0-9]+\]: ")
    tb_start = re.compile("\-+")

    def get_tokens_unprocessed(self, text):
        pylexer = PythonLexer(**self.options)
        tblexer = PythonTracebackLexer(**self.options)

        # print '\nTEXT > \n', text, '\n TEXT'
        for line in text.splitlines():
            lstrip = line.lstrip()
            if lstrip.startswith('Out'):
                line = lstrip + '\n'
            elif lstrip.startswith('...'):
                line = line + '\n'
            else:
                line = line + '\n'
            input_prompt = self.input_prompt.match(line)
            output_prompt = self.output_prompt.match(line)

            if input_prompt is not None:
                yield (0, Generic.Prompt, input_prompt.group())
                code = line[input_prompt.end():]
                for item in pylexer.get_tokens_unprocessed(code):
                  yield item
            elif output_prompt is not None:
                # Use the 'error' token for output.  We should probably make
                # our own token, but error is typicaly in a bright color like
                # red, so it works fine for our output prompts.
                yield (0, Generic.Error, output_prompt.group())
                index = output_prompt.end()
                yield index, Generic.Output, line[index:]
            else:
                yield 0, Generic.Output, line


def setup(app):
    """Setup as a sphinx extension."""

    # This is only a lexer, so adding it below to pygments appears sufficient.
    # But if somebody knows that the right API usage should be to do that via
    # sphinx, by all means fix it here.  At least having this setup.py
    # suppresses the sphinx warning we'd get without it.
    pass

#-----------------------------------------------------------------------------
# Register the extension as a valid pygments lexer
highlighting.lexers['ipython'] = IPythonConsoleLexer()
