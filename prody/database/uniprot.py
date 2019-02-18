import os.path
import numpy as np
from prody import LOGGER, PY3K
from prody.utilities import dictElement, dictElementLoop, openURL, which

import platform, os, re, sys, time, urllib

if PY3K:
    import urllib.parse as urllib
    import urllib.request as urllib2
else:
    import urllib
    import urllib2

from xml.etree.cElementTree import XML

__all__ = ['queryUniprot', ]

def queryUniprot(id, expand=[], regex=True):
    """Query Uniprot with *id* and return a `dictionary` containing the results
    
    :arg expand: entries through which you want to loop dictElements
        until there aren't any elements left
    :type expand: list
    """

    if not isinstance(id, str):
        raise TypeError('id should be a string')

    try:
        record_file = openURL('http://www.uniprot.org/uniprot/{0}.xml'.format(id))
    except:
        raise ValueError('No Uniprot record found with that id')
    
    data = record_file.read()
    record_file.close()
    data = XML(data)

    data = dictElement(data.getchildren()[0], '{http://uniprot.org/uniprot}', number_multiples=True)

    for key in data:
        value = data[key]
        if not key.startswith('dbReference'):
            continue
        
        try:
            if value.get('type') != 'PDB':
                continue
        except AttributeError:
            continue

        pdbid = value.get('id')
        refdata = {'PDB': pdbid}
        for prop in value:
            prop_key = prop.get('type')
            prop_val = prop.get('value')
            refdata[prop_key] = prop_val
        data[key] = refdata
            
    if expand:
        keys = []
        if regex:
            for lt in expand:
                lt_re = re.compile(lt)
                for key in data.keys():
                    if lt_re.match(key):
                        keys.append(key)
        else:
            keys = expand
        data = dictElementLoop(data, keys, '{http://uniprot.org/uniprot}')
    
    return data