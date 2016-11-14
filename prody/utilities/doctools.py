"""This module defines miscellaneous utility functions."""

from textwrap import wrap

__all__ = ['joinLinks', 'joinRepr', 'joinTerms', 'tabulate', 'wrapText']


def joinLinks(links, sep=', ', last=None, sort=False):
    """Returns a string joining *links* as reStructuredText."""

    links = list(links)
    if sort:
        links.sort()
    if last is None:
        last = ''
    else:
        last = last + '`' + links.pop() + '`_'
    return '`' + ('`_' + sep + '`').join(links) + '`_' + last


def joinRepr(items, sep=', ', last=None, sort=False):
    """Returns a string joining representations of *items*."""

    items = [repr(item) for item in items]
    if sort:
        items.sort()
    if last is None:
        last = ''
    else:
        last = last + items.pop()
    return sep.join(items) + last


def joinTerms(terms, sep=', ', last=None, sort=False):
    """Returns a string joining *terms* as reStructuredText."""

    terms = list(terms)
    if sort:
        terms.sort()
    if last is None:
        last = ''
    else:
        last = last + ':term:`' + terms.pop() + '`'
    return ':term:`' + ('`' + sep + ':term:`').join(terms) + '`' + last


def wrapText(text, width=70, join='\n', **kwargs):
    """Returns wrapped lines from :func:`textwrap.wrap` after *join*\ing them.
    """

    try:
        indent = kwargs.pop('indent')
    except KeyError:
        pass
    else:
        kwargs['initial_indent'] = kwargs['subsequent_indent'] = ' ' * indent
    if join:
        return join.join(wrap(text, width, **kwargs))
    else:
        return wrap(text, width, **kwargs)

def tabulate(*cols, **kwargs):
    """Returns a table for columns of data.

    :kwarg header: make first row a header, default is **True**
    :type header: bool
    :kwarg width: 79
    :type width: int
    :kwargs space: number of white space characters between columns,
         default is 2
    :type space: int

    """

    indent = kwargs.get('indent', 0)
    space = kwargs.get('space', 2)
    widths = [max(map(len, cols[0]))] # PY3K: OK
    widths.append(kwargs.get('width', 79) - sum(widths) -
                  len(widths) * space)
    space *= ' '
    bars = (space).join(['=' * width for width in widths])
    lines = [bars]

    for irow, items in enumerate(zip(*cols)): # PY3K: OK
        rows = []
        rows.extend([wrap(item, widths[icol]) if icol else
                      [item.ljust(widths[icol])]
                         for icol, item in enumerate(items)])
        maxlen = max(map(len, rows)) # PY3K: OK
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
