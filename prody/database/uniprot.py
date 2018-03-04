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

import xml.etree.cElementTree as ET

from prody.sequence import Sequence

__all__ = ['queryUniprot', ]

def queryUniprot(id, loop_through=[]):
    """Query Uniprot with *id* and return a `dictionary` containing the results
    
    :arg loop_through: entries through which you want to loop dictElements
        until there aren't any elements left
    :type loop_through: list
    """

    if not isinstance(id, str):
        raise TypeError('id should be a string')

    try:
        record_file = urllib2.urlopen('http://www.uniprot.org/uniprot/{0}.xml'.format(id))
    except:
        raise ValueError('No Uniprot record found with that id')
    
    data = record_file.read()
    record_file.close()
    data = ET.XML(data)

    data = dictElement(data.getchildren()[0], '{http://uniprot.org/uniprot}')

    if loop_through != []:
        data = dictElementLoop(data, loop_through, '{http://uniprot.org/uniprot}')
    
    return data