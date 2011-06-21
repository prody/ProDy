#!/usr/bin/python
"""Plot monthly download statistics obtained from PyPI.

Information: 

  * http://pypi.python.org/stats/months/
  * http://www.python.org/dev/peps/pep-0381/

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'


import bz2
import datetime
import urllib

def incrMonth(dt):
    dt += datetime.timedelta(31)
    return datetime.date(dt.year, dt.month, 1)

# Package info
release_date = datetime.date(2010, 11, 1)
package_name = 'ProDy'

# Parse statistics
labels = []
counts = []
month = release_date
while month <= datetime.date.today(): 
    counts.append(0)
    stats = bz2.decompress(urllib.urlopen(
                           'http://pypi.python.org/stats/months/{0:s}.bz2'
                           .format(month.isoformat()[:7])).read()).split('\n')
    for line in stats:
        if line.startswith(package_name):
            counts[-1] += int(line.split(',')[-1])
    labels.append(month.strftime('%b %y'))
    month = incrMonth(month)

# Make figure
import numpy as np
import matplotlib.pyplot as plt
plt.figure(figsize=(7.5,4))
plt.bar(range(len(counts)), counts, color='black')
plt.xticks(np.arange(len(labels))+0.5, labels, rotation=15, fontsize=10)
plt.yticks(plt.yticks()[0], fontsize=10)
plt.grid()
plt.title('{0:s} monthly download statistics'.format(package_name), fontsize=12)
plt.ylabel('Number of downloads', fontsize=11)
plt.savefig('_static/pypi_monthly.png', dpi=72)
