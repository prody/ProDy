

import sys

PY3K = sys.version_info[0] > 2

if PY3K:
    from optparse import OptionParser as 
else:
    from argparse import __version__ as pyparsing_version
