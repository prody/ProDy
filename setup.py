import glob
import os
import os.path
import sys
import shutil
import cPickle

from distutils.core import setup
from distutils.extension import Extension
from distutils.command.install import install

readme = open('README.txt')
long_description = ''.join(readme.readlines())
readme.close()

__version__ = ''
for line in open('prody/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())
        break
    
def isInstalled(module_name):
    """Check if a required package is installed, by trying to import it."""
    try:
        return __import__(module_name)
    except ImportError:
        return False
    else:
        return True

if not isInstalled('numpy'):    
    print("""NumPy is not installed. This package is required for main ProDy 
features and needs to be installed before you can use ProDy.  
You can find NumPy at: http://numpy.scipy.org""")
    
PACKAGES = ['prody']

EXTENSIONS = []

if os.name != 'java' and sys.version_info[0] == 2:
    pairwise2 = ['cpairwise2module.c', 'pairwise2.py']
    if all([os.path.isfile(os.path.join('prody', fn)) for fn in pairwise2]):  
        EXTENSIONS.append(
            Extension('prody.cpairwise2',
                      ['prody/cpairwise2module.c'],
                      include_dirs=["prody"]
                      ))
    if isInstalled('numpy'):
        import numpy
        kdtree = ['__init__.py', 'KDTree.c', 'KDTree.h', 'KDTree.py', 
                  'KDTreemodule.c', 'Neighbor.h']
        if all([os.path.isfile(os.path.join('prody/KDTree', fn)) for fn in kdtree]):
            EXTENSIONS.append(
                Extension('prody.KDTree._CKDTree',
                          ['prody/KDTree/KDTree.c',
                           'prody/KDTree/KDTreemodule.c'],
                          include_dirs=[numpy.get_include()],
                          ))
        PACKAGES.append('prody.KDTree')

SCRIPTS = glob.glob('scripts/*py')

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

def updateVMDpath(vmdpath):
    """Update VMD path in ProDy settings file (.prody)."""
    
    fn = os.path.join(os.path.expanduser('~'), '.prody')
    settings = None
    if os.path.isfile(fn):
        inf = open(fn)
        settings = cPickle.load(inf)
        inf.close()
        if isinstance(settings, dict):
            settings['vmd'] = vmdpath
        else:
            settings = None
        
    if settings is None:
        settings = {
            'loglevel': 'debug',
            'local_pdb_folder': None,
            'vmd': vmdpath}
    out = open(fn, 'w')
    cPickle.dump(settings, out)
    out.close()
        
class installProDy(install):
    """Override the standard install to install VMD plug-in NMWiz."""
    
    def run(self):
        """Try installing NMWiz and then continue normal installation."""
        vmddir = getVMDDIR()
        if vmddir is not None:
            installNMWiz(vmddir)
            updateVMDpath(os.path.join(vmddir, 'vmd.exe'))
        else:
            print('NMWiz could not be installed, VMD could not be located.')

        
        install.run(self)

setup(
    name='ProDy',
    version=__version__,
    author='Ahmet Bakan',
    author_email='ahb12 at pitt dot edu',
    description='A Python Package for Protein Dynamics Analysis',
    long_description=long_description,
    url='http://www.csb.pitt.edu/ProDy',
    cmdclass={'install' : installProDy,},
    packages=PACKAGES,
    ext_modules=EXTENSIONS,
    license='GPLv3',
    keywords=('protein, dynamics, elastic network model, Gaussian network model, '
              'anisotropic network model, essential dynamics analysis, '
              'principal component analysis, Protein Data Bank, PDB, GNM, ANM, PCA'),
    classifiers=[
                 'Development Status :: 4 - Beta',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: GNU General Public License (GPL)',
                 'Operating System :: MacOS',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 2',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Topic :: Scientific/Engineering :: Chemistry',
                ],
    scripts=SCRIPTS,
    requires=['NumPy', ],
    provides=['ProDy({0:s})'.format(__version__)]
    )
