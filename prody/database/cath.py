
from os import getcwd, remove
# from glob import glob
from os.path import sep as pathsep
from os.path import isdir, isfile, join, split, splitext, normpath

from prody import LOGGER, SETTINGS
from prody.utilities import makePath, gunzip, relpath, copyFile, openURL
from prody.utilities import sympath

__all__ = ['fetchCATH', 'buildCATHNameDict', 'buildPDBChainCATHDict']

def fetchCATH(filename, ftp_host=None, ftp_path=None, **kwargs):
    """Downloads CATH file via FTP."""
    if ftp_host == None:
        ftp_host = 'orengoftp.biochem.ucl.ac.uk'
    if ftp_path == None:
        ftp_path = '/cath/releases/daily-release/newest/'
    from ftplib import FTP
    report = kwargs.get('report', True)
    output_folder = kwargs.pop('folder', None)
    ftp_fn = filename
    try:
        ftp = FTP(ftp_host)
    except Exception as error:
        raise type(error)('FTP connection problem, potential reason: '
                          'no internet connectivity')
    else:
        success = 0
        failure = 0
        filenames = []
        ftp.login('')
        
        data = []
        try:
            ftp.cwd(ftp_path)
            ftp.retrbinary('RETR ' + ftp_fn, data.append)
        except Exception as error:
            if ftp_fn in ftp.nlst():
                LOGGER.warn('{0} download failed ({1}). It is '
                            'possible that you do not have rights to '
                            'download .gz files in the current network.'
                            .format(ftp_fn, str(error)))
            else:
                LOGGER.warn('{0} download failed. {1} does not exist '
                            'on {2}.'.format(ftp_fn, ftp_fn, ftp_host))
            failure += 1
            filenames.append(None)
        else:
            if len(data):
                if output_folder is None:
                    output_folder = getcwd()
                    filename_full = join(output_folder, ftp_fn)

                    with open(filename_full, 'w+b') as pdbfile:
                        write = pdbfile.write
                        [write(block) for block in data]

                    filename_full = normpath(relpath(filename_full))
                    if report: 
                        LOGGER.debug('{0} downloaded ({1})'
                                     .format(ftp_fn, sympath(filename_full)))
                    success += 1
                    filenames.append(filename_full)
                else:
                    LOGGER.warn('{0} download failed, reason unknown.'
                                .format(ftp_fn))
                    failure += 1
                    filenames.append(None)
        ftp.quit()


def buildCATHNameDict(cath_file, iscommpressed=True):
    """Returns a dictionary for CATH names with key of CATH ID."""
    if iscommpressed:
        gunzip(cath_file, 'cath_b.names.temp')
        cath_file = 'cath_b.names.temp'
        
    cath_id2name = dict()
    with open(cath_file, 'r') as file_temp:
        for line in file_temp:
            ind_temp = line.find(' ')
            cath_id2name[line[:ind_temp]] = line[ind_temp:].strip()
    if iscommpressed:
        remove(cath_file) 
    return cath_id2name
    
def buildPDBChainCATHDict(cath_file, iscommpressed=True):
    """Returns a dictionary for CATH info (ID and version) with key of PDB chain."""
    if iscommpressed:
        gunzip(cath_file, 'cath_b.all.temp')
        cath_file = 'cath_b.all.temp'
    
    cath_dict_temp = dict()
    cath_i_dict = dict()
    with open(cath_file, 'r') as file_temp:
        for line in file_temp:
            line = line.strip()
            if line != '':
                line_list = line.split(' ')
                cath_dict_temp[line_list[0]] = line_list[1:]
                key, value = line[0:5], line[5:7]
                if key in cath_i_dict:
                    cath_i_dict[key].append(value)
                else:
                    cath_i_dict[key] = [value]
    pdbChain2CATH = dict()
    for key, values in cath_i_dict.items():
        pdbChain2CATH[key] = []
        for v in values:
            pdbChain2CATH[key].append(cath_dict_temp[key+v])
    if iscommpressed:
        remove(cath_file) 
    return pdbChain2CATH
    
    
# ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/daily-release/newest/

# fetchCATH('cath-b-newest-names.gz')
# cath_id2name = buildCATHNameDict('cath-b-newest-names.gz')

# fetchCATH('cath-b-newest-all.gz')
# pdbChain2CATH = buildPDBChainCATHDict('cath-b-newest-all.gz')
