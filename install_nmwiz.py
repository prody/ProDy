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

import os
import os.path
import sys
import shutil

def getVMDDIR():
    if sys.platform == 'win32': 
        import _winreg
        for vmdversion in ('1.8.7',): 
            try:
                key = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE, 
                    'Software\\University of Illinois\\VMD\\' + vmdversion)
                vmddir = _winreg.QueryValueEx(key, 'VMDDIR')[0]
                if os.path.isdir(vmddir):
                    return vmddir
            except:    
                pass
            try:
                key = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE, 
                    'Software\\WOW6432node\\University of Illinois\\VMD\\' + 
                    vmdversion)
                vmddir = _winreg.QueryValueEx(key, 'VMDDIR')[0]
                if os.path.isdir(vmddir):
                    return vmddir
            except:    
                pass
    else:
        try:
            pipe = os.popen('which vmd')
            vmdbin = pipe.next().strip()
            vmdbin = open(vmdbin)
            for line in vmdbin:
                if 'defaultvmddir' in line:
                    exec(line.strip())
                    vmddir = defaultvmddir
                    break
            vmdbin.close()
        except:
            pass
        else:
            if os.path.isdir(vmddir):
                return vmddir
    return None

def installNMWiz(vmddir):
    """Copy NMWiz plug-in files to $VMDDIR/plugins/noarch/tcl folder."""
    plugindir = os.path.join(vmddir, 'plugins', 'noarch', 'tcl')
    nmwiz = 'nmwiz' + __version__[:3]
    nmwizdir = os.path.join(plugindir, nmwiz)
    if not os.path.isdir(nmwizdir):
        os.mkdir(nmwizdir)
    print('Copying NMWiz files to ' + nmwizdir)
    for fn in ('nmwiz.tcl', 'pkgIndex.tcl'):
        print(os.path.join(nmwiz, fn) + ' -> ' + os.path.join(nmwizdir, fn))
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
        tcl = open(loadplugins, 'w')
        for line in newlines:        
            tcl.write(line)
        tcl.close()

if __name__ == '__main__':
    vmddir = getVMDDIR()
    if vmddir is not None:
        installNMWiz(vmddir)
