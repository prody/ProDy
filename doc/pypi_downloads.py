#!/usr/bin/python
__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

from datetime import datetime
import xmlrpclib

package_name = 'ProDy'

pypi = xmlrpclib.Server('http://pypi.python.org/pypi')
# Write statistics in a CSV file
rst = open('reports/pypi_downloads.csv', 'w')
rst.write('Release;Date;Downloads;\n')
show_hidden = True
releases = pypi.package_releases(package_name, show_hidden)
downloads_total = 0
for release in releases:
    urls = pypi.release_urls(package_name, release)
    downloads = 0
    for url in urls:
        downloads += url['downloads']
    downloads_total += downloads
    date = datetime.strptime(url['upload_time'].value, "%Y%m%dT%H:%M:%S")
    rst.write('{0:s};{1:s};{2:d};{3:s}\n'.format(
        release, 
        date.strftime('%B %e, %Y'), 
        downloads,
        '#' + '=' * int(round(downloads / 10.)) + '>'))
    
downloads_average = int(round((downloads_total * 1.0 / len(releases))))
rst.write('Average; ;{0:d};{1:s}\n'.format(
    downloads_average, '#' + '=' * int(round(downloads_average / 10.)) + '>'))
rst.write('Total; ;{0:d};\n'.format(downloads_total))
rst.close()


# Now write an HTML sidebar to show the total number of downloads in index 

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
html.write('''(see <a href="{{ pathto("reports/pypi_downloads") }}">details</a>).
</p>''')
