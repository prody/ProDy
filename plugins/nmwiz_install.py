# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2011 Ahmet Bakan
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

__version__ = '1.0'

import os
import os.path
import sys
import shutil
import glob
from types import StringType, UnicodeType

PY3K = sys.version_info[0] > 2

def getVMDpaths():
    """Return VMDDIR, if bin=True, return path to the executable."""
    vmdbin = None
    vmddir = None
    if sys.platform == 'win32': 
        if PY3K:
            import winreg as _winreg
        else:
            import _winreg
        for vmdversion in ('1.8.7', '1.9'): 
            try:
                key = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE, 
                    'Software\\University of Illinois\\VMD\\' + vmdversion)
                vmddir = _winreg.QueryValueEx(key, 'VMDDIR')[0]
                vmdbin = os.path.join(vmddir, 'vmd.exe') 
            except:    
                pass
            try:
                key = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE, 
                    'Software\\WOW6432node\\University of Illinois\\VMD\\' + 
                    vmdversion)
                vmddir = _winreg.QueryValueEx(key, 'VMDDIR')[0]
                vmdbin = os.path.join(vmddir, 'vmd.exe') 
            except:    
                pass
    else:
        try:
            pipe = os.popen('which vmd')
            vmdbin = pipe.next().strip()
            vmdfile = open(vmdbin)
            for line in vmdfile:
                if 'defaultvmddir' in line:
                    exec(line.strip())
                    vmddir = defaultvmddir
                    break
            vmdfile.close()
        except:
            pass
    if isinstance(vmdbin, (StringType, UnicodeType)) and \
       isinstance(vmddir, (StringType, UnicodeType)) and \
       os.path.isfile(vmdbin) and os.path.isdir(vmddir):  
        return vmdbin, vmddir
    return None, None

def installNMWiz(vmddir):
    """Copy NMWiz plugin files to $VMDDIR/plugins/noarch/tcl folder."""
    plugindir = os.path.join(vmddir, 'plugins', 'noarch', 'tcl')
    nmwiz = 'nmwiz' + __version__[:3]
    nmwizdir = os.path.join(plugindir, nmwiz)
    if not os.path.isdir(nmwizdir):
        os.mkdir(nmwizdir)
    print('installing NMWiz into ' + plugindir)
    for fn in ('nmwiz.tcl', 'pkgIndex.tcl'):
        print('copying ' + os.path.join(nmwiz, fn) + ' -> ' + os.path.join(nmwizdir, fn))
        shutil.copy(os.path.join(nmwiz, fn), os.path.join(nmwizdir, fn))
    loadplugins = os.path.join(vmddir, 'scripts', 'vmd', 'loadplugins.tcl') 
    tcl = open(loadplugins)
    oldlines = tcl.readlines()
    newlines = []
    update = True
    for line in oldlines:
        newlines.append(line)
        if 'nmwiz_tk' in line:
            update = False
            break
        if 'namdplot_tk' in line:
            newlines.append('  vmd_install_extension nmwiz   nmwiz_tk   "Analysis/Normal Mode Wizard"\n')
    tcl.close()
    if update:
        print('updating ' + loadplugins)
        tcl = open(loadplugins, 'w')
        for line in newlines:        
            tcl.write(line)
        tcl.close()
    else:
        print('skipping update of ' + loadplugins)
    
def removeNMWiz(vmddir):
    """Remove older versions of NMWiz from $VMDDIR/plugins/noarch/tcl folder."""
    plugindir = os.path.join(vmddir, 'plugins', 'noarch', 'tcl')
    nmwiz = 'nmwiz' + __version__[:3]
    for nmwizdir in glob.glob(os.path.join(plugindir, 'nmwiz*')): 
        if nmwiz in nmwizdir: 
            continue
        print('removing previous NMWiz release from ' + nmwizdir)
        for nmwizfile in glob.glob(os.path.join(nmwizdir, '*')):
            print('removing ' + nmwizfile) 
            os.remove(nmwizfile)
        print('removing ' + nmwizdir)
        os.rmdir(nmwizdir)
    
if __name__ == '__main__':
    vmdbin, vmddir = getVMDpaths()
    if vmddir is not None:
        try:
            installNMWiz(vmddir)
        except:
            print('NMWiz could not be installed. User might not have '
                  'rights to write in the VMD path {0:s}.'
                  .format(vmddir))
        else:
            removeNMWiz(vmddir)
    else:
        print('NMWiz could not be installed, VMD could not be located.')
    raw_input('Press Enter to exit.')
