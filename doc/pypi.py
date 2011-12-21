#!/usr/bin/python
"""Update Python package release and monthly download statistics.

Information: 

  * http://pypi.python.org/stats/months/
  * http://www.python.org/dev/peps/pep-0381/

"""
__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2011  Ahmet Bakan'

import time
import datetime
import xmlrpclib
import cPickle
import sys
import os.path

PACKAGE_NAME = 'ProDy'
RELEASE_DATE = datetime.date(2010, 11, 1)

def loadStats():
    # Gather download statistics
    if os.path.isfile(PACKAGE_NAME.lower() + '_pypi_stats.pkl'):
        pkl = open(PACKAGE_NAME.lower() + '_pypi_stats.pkl')
        stats = cPickle.load(pkl)
        pkl.close()
        return stats
    else:
        from collections import defaultdict
        return defaultdict(int)

def saveStats(stats):
    pkl = open(PACKAGE_NAME.lower() + '_pypi_stats.pkl', 'wb')
    cPickle.dump(stats, pkl)
    pkl.close()
    
def updateReleaseStats():

    pypi = xmlrpclib.Server('http://pypi.python.org/pypi')

    STATS = loadStats()    

    show_hidden = True
    releases = pypi.package_releases(PACKAGE_NAME, show_hidden)
    current_urls = None
    for release in releases:
        urls = pypi.release_urls(PACKAGE_NAME, release)
        if current_urls is None:
            current_urls = urls
        downloads = 0
        for url in urls:
            downloads += url['downloads']
        date = datetime.datetime.strptime(url['upload_time'].value, 
                                          "%Y%m%dT%H:%M:%S")
        STATS[release] = (date.strftime('%B %e, %Y'), downloads)
        
    saveStats(STATS)
    
    # Write statistics in a CSV file
    releases = STATS.keys()
    releases.sort(reverse=True)
    downloads_total = 0
    stats = open('_static/pypi_statistics.csv', 'w')
    stats.write('Release;Date;Downloads\n')
    for release in releases:
        item = STATS[release]
        if isinstance(item, tuple) and len(item) == 2:
            date, downloads = item
            stats.write('{0:s};{1:s};{2:d}\n'.format(release, date, downloads))
            downloads_total += downloads
    stats.write('Total;{0:s};{1:d}\n'.format(time.strftime('%B %e, %Y'), 
                                             downloads_total))
    stats.close()

    # Write a CSV file with info on and links to the current downloads 
    packagetype_map = {'sdist': 'Source', 
                       'bdist_wininst': 'MS Windows installer'}
    python_version_map = {'source': ''} 
    files = open('_static/pypi_downloads.csv', 'w')
    files.write('File;Type;Py Version;Size;Downloads\n')
    for url in current_urls:
        url['packagetype'] = packagetype_map.get(url['packagetype'], 
                                                 url['packagetype'])
        url['python_version'] = python_version_map.get(url['python_version'], 
                                                       url['python_version'])
        files.write('`{0[filename]:s} <{0[url]:s}>`_ '.format(url))
        files.write('(`md5 <http://pypi.python.org/pypi?:action=show_md5&'
                    'digest={0[md5_digest]:s}>`_);'.format(url))
        files.write('{0[packagetype]:s};'.format(url))
        files.write('{0[python_version]:s};'.format(url))
        files.write('{0:d}KB;'.format(int(url['size']/1024)))
        files.write('{0[downloads]:d}'.format(url))
        files.write('\n')
    files.close()

    # Write an HTML sidebar to show the total number of downloads in index
    html = open('_templates/getprody.html', 'w')
    html.write('''
    <h3>Download</h3>
    <p>Current version: <b>{{ version }}</b></p>
    <p><b>{{ docrelease }}</b></p>
    <p>See <a href="{{ pathto("getprody") }}">Getting ProDy</a> for installing 
       instructions, or install it with: 
       <pre>easy_install -U ProDy</pre>   
    ''')
    html.write('Total of {0:d} downloads '.format(downloads_total)) 
    html.write('(see <a href="{{ pathto("reports/pypi_statistics") }}">'
               'details</a>).</p>')

def updateMonthlyStats():
    import bz2
    import datetime
    import urllib

    def incrMonth(dt):
        dt += datetime.timedelta(31)
        return datetime.date(dt.year, dt.month, 1)

    package_name = PACKAGE_NAME

    STATS = loadStats()
    # Parse statistics
    #labels = []
    #counts = []
    month = RELEASE_DATE
    today = datetime.date.today()
    while month <= today: 
        key = month.isoformat()[:7]
        if not key in STATS:
            count = 0
            try:
                stats = bz2.decompress(urllib.urlopen(
                               'http://pypi.python.org/stats/months/{0:s}.bz2'
                               .format(key)).read()).split('\n')
            except:
                break
            current = defaultdict(int)
            for line in stats:
                if line.startswith(package_name):
                    items = line.split(',')
                    incr = int(items[-1])
                    count += incr
                    release = ''
                    for char in items[1][len(package_name)+1:]:
                        if char.isalpha():
                            break
                        release += char
                    release = release[:-1]
                    if release[-1] == '.':
                        release = release[:-1]
                    if release.endswith('.0'):
                        release = release[:-2]
                    current[release] += incr            
                    
            if month.month == today.month and month.year == today.year: 
                STATS['THIS-MO'] = current
            else:
                STATS[key] = current
        #labels.append(month.strftime('%b %y'))
        month = incrMonth(month)
    saveStats(STATS)
    return None
    # Make figure
    import numpy as np
    import matplotlib.pyplot as plt
    plt.figure(figsize=(7.5,4))
    plt.bar(range(len(counts)), counts, color='black')
    plt.xticks(np.arange(len(labels))+0.5, labels, rotation=15, fontsize=10)
    plt.yticks(plt.yticks()[0], fontsize=10)
    plt.grid()
    plt.title('{0:s} monthly download statistics'.format(package_name), 
              fontsize=12)
    plt.ylabel('Number of downloads', fontsize=11)
    plt.savefig('_static/pypi_monthly.png', dpi=72)

if __name__ == '__main__':
    
    if len(sys.argv) > 1:
        if sys.argv[1].lower().startswith('m'):
            updateMonthlyStats()
        elif sys.argv[1].lower().startswith('r'):
            updateReleaseStats()
    else:
        stats = loadStats()
        keys = stats.keys()
        keys.sort(reverse=True)
        month = 0
        release = 0
        for key in keys: 
            value = stats[key]
            if isinstance(value, int):
                month += value
            else:
                release += value[1]
            print key, value
        print 'Total (month):', month
        print 'Total (release):', release
