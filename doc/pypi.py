#!/usr/bin/python
__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

from datetime import datetime
import xmlrpclib

package_name = 'ProDy'

pypi = xmlrpclib.Server('http://pypi.python.org/pypi')

# Write statistics in a CSV file
stats = open('_static/pypi_statistics.csv', 'w')
stats.write('Release;Date;Downloads\n')
show_hidden = True
releases = pypi.package_releases(package_name, show_hidden)
downloads_total = 0
current_urls = None
for release in releases:
    urls = pypi.release_urls(package_name, release)
    if current_urls is None:
        current_urls = urls
    downloads = 0
    for url in urls:
        downloads += url['downloads']
    downloads_total += downloads
    date = datetime.strptime(url['upload_time'].value, "%Y%m%dT%H:%M:%S")
    stats.write('{0:s};{1:s};{2:d}\n'.format(
        release, 
        date.strftime('%B %e, %Y'), 
        downloads))
    
stats.write('Total; ;{0:d}\n'.format(downloads_total))
stats.close()


# Write a CSV file with info on and links to the current downloads 
packagetype_map = {'sdist': 'Source', 'bdist_wininst': 'MS Windows installer'}
python_version_map = {'source': ''} 
files = open('_static/pypi_downloads.csv', 'w')
files.write('File;Type;Py Version;Size;Downloads\n')
for url in current_urls:
    url['packagetype'] = packagetype_map.get(url['packagetype'], url['packagetype'])
    url['python_version'] = python_version_map.get(url['python_version'], url['python_version'])
    files.write('`{0[filename]:s} <{0[url]:s}>`_ '.format(url))
    files.write('(`md5 <http://pypi.python.org/pypi?:action=show_md5&digest={0[md5_digest]:s}>`_);'.format(url))
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
html.write('''(see <a href="{{ pathto("reports/pypi_statistics") }}">details</a>).
</p>''')
