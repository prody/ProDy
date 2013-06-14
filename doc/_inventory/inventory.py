"""A module for managing ProDy and external inventories."""

import zlib

__all__ = ['trim_labels', 'remove_api']


def parse_inventory(inv):

    inp = open(inv, 'rb')
    header = [inp.readline().decode('utf-8') for i in range(4)]
    if 'zlib' not in header[-1]:
        raise ValueError
    lines = zlib.decompress(inp.read()).decode('utf-8').splitlines()
    inp.close()
    return header, lines


def trim_labels(inv, out):

    header, lines = parse_inventory(inv)
    trimmed = []
    for line in lines:
        items = line.split(' py:')

        if len(items) > 1:
            if (items[1].startswith('method') or
                    items[1].startswith('attribute')):
                items[0] = '.'.join(items[0].split('.')[-2:])
            elif items[1].startswith('module'):
                _ = items[0].split('.')
                if len(_) > 1:
                    if _[-1] == _[-2]:
                        items[0] = '.'.join(_[-2:] or _)
                    else:
                        items[0] = _[-1]
            else:
                items[0] = items[0].split('.')[-1]
            trimmed.append(items[0] + ' py:' + items[1])
        else:
            trimmed.append(line)

    out = open(out, 'wb')
    for line in header:
        out.write(line.encode('utf-8'))
    compressor = zlib.compressobj(9)

    for line in trimmed:
        out.write(compressor.compress((line + '\n').encode('utf-8')))
    out.write(compressor.flush())
    out.close()


def remove_api(inv, out):

    header, lines = parse_inventory(inv)
    trimmed = []
    for line in lines:
        items = line.split(' py:')
        if len(items) == 1:
            trimmed.append(line)

    out = open(out, 'wb')
    for line in header:
        out.write(line.encode('utf-8'))
    compressor = zlib.compressobj(9)

    for line in trimmed:
        out.write(compressor.compress((line + '\n').encode('utf-8')))
    out.write(compressor.flush())
    out.close()
