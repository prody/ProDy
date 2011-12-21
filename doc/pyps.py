#!/usr/bin/python
"""Python package download statistics.

Information: 

  * http://pypi.python.org/stats/months/
  * http://www.python.org/dev/peps/pep-0381/

"""
__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2011  Ahmet Bakan'

import sys
import csv
import bz2
import cPickle
import logging
import os.path
import urllib2
import datetime
from HTMLParser import HTMLParser
from collections import defaultdict

PACKAGE = 'ProDy'
LOGGER = logging.getLogger('.pyps')
LOGGER.setLevel(logging.INFO)
LOGGER.addHandler(logging.StreamHandler())

class PyPIStatsFile(object):
    
    def __init__(self, url, month, modified=None, size=None):
        
        self.url = url
        self.month = month
        self.modified = modified
        self.size = size
    
    def read(self):
        url = urllib2.urlopen(self.url)
        stats = bz2.decompress(url.read())
        url.close()
        return stats
                

class PyPIStatsURLParser(HTMLParser):

    def __init__(self, url):
        
        HTMLParser.__init__(self)
        self.url = url
        self.files = []
        self._temp = None
    
    def handle_starttag(self, tag, attrs):
        
        if tag == 'a':
            self._temp = None

    def handle_endtag(self, tag):
        
        pass

    def handle_data(self, data):
        
        data = data.strip()
        if data:
            if self._temp:
                items = data.split()
                if len(items) == 3:
                    _ = ' '.join(items[:2])
                    date = datetime.datetime.strptime(_, '%d-%b-%Y %H:%M')
                    self.files.append(
                        PyPIStatsFile(self.url + self._temp,
                        os.path.splitext(self._temp)[0],
                        date, int(items[2]))
                    )               
            elif data.startswith('2') and data.endswith('.bz2'):
                self._temp = data
        
def fetchURLs(url='http://pypi.python.org/stats/months/'):

    LOGGER.info("Fetching content from '{0:s}'.".format(url))
    import urllib2
    stats_months = urllib2.urlopen(url)
    feed = stats_months.read()
    stats_months.close()
    
    LOGGER.info("Parsing monthly statistics file details.")
    parser = PyPIStatsURLParser(url)
    parser.feed(feed)
    return parser


def package_filename(package):
    
    return package + '_stats.pkl'

def load_stats(filename):

    # Gather download statistics
    if filename and os.path.isfile(filename):
        LOGGER.info("Loading statistics from '{0:s}'".format(filename))
        pkl = open(filename)
        stats = cPickle.load(pkl)
        pkl.close()
        return stats
    else:
        from collections import defaultdict
        return defaultdict(int)

def save_stats(filename, stats):

    pkl = open(filename, 'wb')
    cPickle.dump(stats, pkl)
    pkl.close()
    return filename

def get_version(package, filename):
    
    ndot = 0
    version = ''
    i = len(package) + 1
    while i < len(filename):
        char = filename[i]
        if char == '.':
            if not filename[i+1].isdigit():
                break
        version += char
        i += 1
    return version

def update_stats(package=PACKAGE, filename=None):
    
    if filename is None:
        filename = package_filename(package)
    stats = load_stats(filename)
    p = fetchURLs()
    noupdates = True
    for f in p.files:
        if f.month in stats and stats[f.month]['modified'] == f.modified: 
            continue
        noupdates = False
        LOGGER.info("Updating statistics for " + f.month + '.')
        stats[f.month] = defaultdict(int)
        month = stats[f.month]
        month['modified'] = f.modified
        for line in f.read().split('\n'):
            if not line.startswith(package):
                continue
            items = line.split(',')
            version = get_version(package, items[1])
            count = int(items[-1])
            month[version] += count
    if noupdates:
        LOGGER.info("Nothing to update.")
    else:
        save_stats(filename, stats)
        LOGGER.info("Package statistics are updated ({0:s}).".format(filename))
        return filename
    
def release_stats(stats):
    
    stats = load_stats(stats)
    releases = defaultdict(int)
    for month in stats.itervalues(): 
        for key, value in month.iteritems():
            if key == 'modified':
                continue
            releases[key] += value
    releases = releases.items()
    releases.sort()
    return releases

def release_stats_csv(stats, outname, delimiter=','):
    
    ostream = open(outname, 'wb')
    out = csv.writer(ostream, delimiter=delimiter)
    
    out.writerow(['Release', 'Downloads'])
    for month in release_stats(stats):
        out.writerow(month)
    ostream.close()
    LOGGER.info("Monthly statistics are written in '{0:s}'".format(outname))
    return outname

def total_downloads(stats):

    stats = load_stats(stats)
    total = 0
    for month in stats.itervalues(): 
        for key, value in month.iteritems():
            if key == 'modified':
                continue
            total += value
    return total

def monthly_stats(stats):

    stats = load_stats(stats)
    months = []
    for month, stats in stats.iteritems(): 
        counts = 0
        for key, value in stats.iteritems():
            if key == 'modified':
                continue
            counts += value
        if counts > 0:
            months.append((month, counts))
    months.sort()
    return months

def monthly_stats_csv(stats, outname, delimiter=','):
    
    ostream = open(outname, 'wb')
    out = csv.writer(ostream, delimiter=delimiter)
    
    out.writerow(['Month', 'Downloads'])
    for month in monthly_stats(stats):
        out.writerow(month)
    ostream.close()
    LOGGER.info("Monthly statistics are written in '{0:s}'".format(outname))
    return outname

def monthly_stats_plot(stats, outname, dpi=72, xlabelstep=2):
    labels = []
    counts = []
    for m, c in monthly_stats(stats):
        labels.append(m)
        counts.append(c)
    import numpy as np
    import matplotlib.pyplot as plt
    plt.figure(figsize=(7.5,4))
    plt.bar(range(len(counts)), counts, color='black')
    plt.xticks(np.arange(len(labels))[::xlabelstep]+0.5, 
               labels[::xlabelstep], rotation=15, fontsize=10)
    plt.yticks(plt.yticks()[0], fontsize=10)
    plt.grid()
    plt.title('Monthly downloads', fontsize=12)
    plt.ylabel('Number of downloads', fontsize=11)
    plt.savefig(outname, dpi=dpi)
    LOGGER.info("Monthly downloads plot is saved as '{0:s}'.".format(outname))
    return outname

def current_release(package=PACKAGE):
    
    url = 'http://pypi.python.org/pypi'
    LOGGER.info("Connecting to '{0:s}'.".format(url))
    import xmlrpclib
    pypi = xmlrpclib.Server(url)
    
    show_hidden = False
    releases = pypi.package_releases(PACKAGE, show_hidden)
    
    for release in releases:
        urls = pypi.release_urls(PACKAGE, release)
        break
    return urls

def current_release_csv(package=PACKAGE, outname=None, delimiter=';', 
                        rst=False):
    
    if outname is None:
        outname = package + '_current_release.csv'
    elif not outname.endswith('.csv'):
        outname += '.csv'
    ostream = open(outname, 'wb')
    out = csv.writer(ostream, delimiter=delimiter)
    
    # Write a CSV file with info on and links to the current downloads 
    packagetype_map = {'sdist': 'Source', 
                       'bdist_wininst': 'MS Windows installer'}
    python_version_map = {'source': ''} 
    if rst:
        out.writerow(['File', 'Type', 'Py Version', 'Size', 'Downloads'])
    else:
        out.writerow(['File', 'URL', 'md5', 'Type', 
                      'Py Version', 'Size', 'Downloads'])
    for url in current_release(package):
        url['packagetype'] = packagetype_map.get(url['packagetype'], 
                                                 url['packagetype'])
        url['python_version'] = python_version_map.get(url['python_version'], 
                                                       url['python_version'])
        if rst:
            out.writerow(
                ['`{0[filename]:s} <{0[url]:s}>`_ '
                 '(`md5 <http://pypi.python.org/pypi?:action=show_md5&'
                 'digest={0[md5_digest]:s}>`_)'.format(url),
                 '{0[packagetype]:s}'.format(url),
                 '{0[python_version]:s}'.format(url),
                 '{0:d}KB'.format(int(url['size']/1024)),
                 '{0[downloads]:d}'.format(url)]            
            )
        else:
            out.writerow(
                ['{0[filename]:s}'.format(url),
                 '{0[url]:s}'.format(url),
                 '{0[md5_digest]:s}'.format(url),
                 '{0[packagetype]:s}'.format(url),
                 '{0[python_version]:s}'.format(url),
                 '{0:d}KB'.format(int(url['size']/1024)),
                 '{0[downloads]:d}'.format(url)]            
            )

    ostream.close()
    LOGGER.info("Current release details are written to '{0:s}'."
                .format(outname))
    return outname


import argparse
parser = argparse.ArgumentParser(
    description="Fetch package download statistics from Python Package "
                "Index (PyPI). Package needs to be distributed via PyPI.",
    )#usage='%(prog)s [--help] <command> [options] pkg')
    
parser.add_argument('-q', '--quiet', help="don't print log messages to stderr",
    action='store_true', default=False)

subparsers = parser.add_subparsers(
    title='subcommands')
        
subparser = subparsers.add_parser('update', 
    help='update download statistics')

subparser.add_argument('-s', '--stats', default=None, metavar='FILENAME',
        help="filename for storing package statistics (default: "
             "'package_stats.pkl')")

subparser.add_argument('package', help='Python package name')

subparser.set_defaults(
    func=lambda args: update_stats(args.package, args.stats))


subparser = subparsers.add_parser('current', 
    help='write current release details in a CSV file')

subparser.add_argument('-o', default=None, metavar='FILENAME',
        help="CSV filename for writing current (default: "
             "'package_current_release.csv')")
             
subparser.add_argument('-d', default=',', metavar='DELIMITER',
        help="column delimiter (default: '%(default)s')")

subparser.add_argument('--rst', default=False, action='store_true',
        help="write reStructured text")

subparser.add_argument('package', help='Python package name')

subparser.set_defaults(
    func=lambda args: 
        current_release_csv(args.package, outname=args.o, 
                            delimiter=args.d, rst=args.rst))

subparser = subparsers.add_parser('monthly', 
    help='write monthly download statistics in a CSV file')

subparser.add_argument('-o', default='monthly_statistics.pkl',
    metavar='FILENAME', 
    help="CSV filename (default: '%(default)s')")
             
subparser.add_argument('-d', default=',', metavar='DELIMITER',
        help="column delimiter (default: '%(default)s')")

subparser.add_argument('stats', help='package statistics filename')

subparser.set_defaults(
    func=lambda args: 
        monthly_stats_csv(args.stats, args.o, args.d))

subparser = subparsers.add_parser('msplot', 
    help='plot monthly statistics, requires Matplotlib')

subparser.add_argument('-o', default='monthly_statistics.png',
    metavar='FILENAME', help="figure filename (default: '%(default)s')")
             
subparser.add_argument('-d', default=',', metavar='DELIMITER',
        help="column delimiter (default: '%(default)s')")

subparser.add_argument('stats', help='package statistics filename')

subparser.set_defaults(
    func=lambda args: monthly_stats_plot(args.stats, args.o))


subparser = subparsers.add_parser('release', 
    help='write release download statistics in a CSV file')

subparser.add_argument('-o', default='release_statistics.pkl',
    metavar='FILENAME', 
    help="CSV filename (default: '%(default)s')")
             
subparser.add_argument('-d', default=',', metavar='DELIMITER',
        help="column delimiter (default: '%(default)s')")

subparser.add_argument('stats', help='package statistics filename')

subparser.set_defaults(
    func=lambda args: release_stats_csv(args.stats, args.o, args.d))


if __name__ == '__main__':

    if len(sys.argv) == 1:    
        parser.print_help()
    else:
        args = parser.parse_args()
        if args.quiet: 
            LOGGER.setLevel(logging.WARNING)
        args.func(args)
