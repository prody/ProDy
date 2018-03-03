import os.path
import numpy as np
from prody import LOGGER, PY3K
from prody.utilities import dictElement, openURL, which

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

def queryUniprot(id):
    """Query Uniprot with *id* and return a `dictionary` containing the results"""

    if not isinstance(id, str):
        raise TypeError('id should be a string')

    try:
        record_file = urllib2.urlopen('http://www.uniprot.org/uniprot/{0}.xml'.format(id))
    except:
        raise ValueError('No Uniprot record found with that id')
    
    data = record_file.read()
    record_file.close()
    data = ET.XML(data)

    return data