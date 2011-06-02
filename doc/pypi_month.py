#!/usr/bin/python
"""Plots monthly download statistics parsed from PyPI.

Information: 

  * http://pypi.python.org/stats/months/
  * http://www.python.org/dev/peps/pep-0381/

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'


from collections import defaultdict
import bz2
import datetime
import urllib
import xmlrpclib
import matplotlib.pyplot as plt

release_date = datetime.date(2010, 11, 1)
package_name = 'ProDy'
show_hidden = True

def incrMonth(dt):
    dt += datetime.timedelta(31)
    return datetime.date(dt.year, dt.month, 1)

pypi = xmlrpclib.Server('http://pypi.python.org/pypi')
releases = pypi.package_releases(package_name, show_hidden)
month = release_date
today = datetime.date.today()
month_labels = []
counts_total = []
counts_release = defaultdict(list)
while month <= today: 
    stats = bz2.decompress(urllib.urlopen('http://pypi.python.org/stats/months/{0:s}.bz2'
                           .format(month.isoformat()[:7])).read()).split('\n')
    counts = defaultdict(int)
    total = 0 
    for l in stats:
        for rel in releases:
            if l.startswith(package_name) and rel in l:
                c = int(l.split(',')[-1])
                counts[rel] += c
                total += c
    for rel in releases:
        counts_release[rel].append(counts[rel])
    month_labels.append(month.strftime('%b %y'))
    counts_total.append(total)
    month = incrMonth(month)

plt.figure(figsize=(7.5,4))
plt.bar(range(len(total_counts)), total_counts, color='black')
plt.xticks(arange(len(month_labels))+0.5, month_labels, rotation=15, fontsize=10)
plt.yticks(plt.yticks()[0], fontsize=10)
plt.grid()
plt.title('{0:s} monthly download statistics'.format(package_name), fontsize=12)
plt.ylabel('Number of downloads', fontsize=11)
plt.savefig('_static/pypi_monthly.png', dpi=72)

"""
plt.figure(figsize=(8,4))
for rel in releases:
    plt.plot(range(len(total_counts))[:-1], counts_release[rel][:-1], '-', marker='o', label=rel)

plt.xticks(arange(len(month_labels))[:-1]+0.5, month_labels[:-1], rotation=15, fontsize=10)
plt.yticks(plt.yticks()[0], fontsize=10)
plt.grid()
plt.legend()
plt.savefig('_static/pypi_monthly_release.png', dpi=72)
"""
