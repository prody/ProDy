# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2012 Ahmet Bakan
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

"""This module defines miscellaneous utility functions."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from textwrap import wrap as textwrap

__all__ = ['joinLinks', 'joinTerms', 'tabulate', 'wrap']


def joinLinks(links, sep=', ', last=None, sort=False):
    
    links = list(links)
    if sort:
        links.sort()
    if last is None:
        last = ''
    else:
        last = last + '`' + links.pop() + '`_'
    return '`' + ('`_' + sep + '`').join(links) + '`_' + last


def joinTerms(terms, sep=', ', last=None, sort=False):
    
    terms = list(terms)
    if sort:
        terms.sort()
    if last is None:
        last = ''
    else:
        last = last + ':term:`' + terms.pop() + '`'
    return ':term:`' + ('`' + sep + ':term:`').join(terms) + '`' + last


def wrap(text, width=70, join='\n', **kwargs):
    """Return wrapped lines from :func:`textwrap.wrap` after *join*\ing them.
    """
    
    try:
        indent = kwargs.pop('indent')
    except KeyError:
        pass
    else:
        kwargs['initial_indent'] = kwargs['subsequent_indent'] = indent
    
    return join.join(textwrap(text, width, **kwargs))
    
def tabulate(*cols, **kwargs):
    """Return a table for columns of data. 
    
    :kwarg header: make first row a header, default is **True**
    :type header: bool
    :kwarg width: 79
    :type width: int
    :kwargs space: number of white space characters between columns,
         default is 2
    :type space: int
    
    """
    
    space = kwargs.get('space', 2)
    widths = [max(map(len, cols[0]))]
    widths.append(kwargs.get('width', 79) - sum(widths) - len(widths) * space)
    space *= ' '
    bars = (space).join(['=' * width for width in widths])
    lines = [bars]
    
    for irow, items in enumerate(zip(*cols)):
        rows = []
        map(rows.append, [textwrap(item, widths[icol]) if icol else 
                          [item.ljust(widths[icol])] 
                          for icol, item in enumerate(items)])
        maxlen = max(map(len, rows))
        if maxlen > 1:     
            for i, row in enumerate(rows):
                row.extend([' ' * widths[i]] * (maxlen - len(row)))
        for line in zip(*rows):
            lines.append(space.join(line))
        if not irow and kwargs.get('header', True):
            lines.append(bars)
    if irow > 1:
        lines.append(bars)
    return '\n'.join(lines)
