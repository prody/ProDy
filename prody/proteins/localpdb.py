# -*- coding: utf-8 -*-
"""This module defines functions for handling local PDB folders."""

from glob import glob, iglob
from os.path import sep as pathsep
from os.path import abspath, isdir, isfile, join, split, splitext, normpath

from prody import LOGGER, SETTINGS
from prody.utilities import makePath, gunzip, relpath, copyFile, isWritable
from prody.utilities import sympath, isListLike

from . import wwpdb
from .wwpdb import checkIdentifiers, fetchPDBviaFTP, fetchPDBviaHTTP


__all__ = ['pathPDBFolder', 'pathPDBMirror',
           'fetchPDB', 'fetchPDBfromMirror',
           'iterPDBFilenames', 'findPDBFiles',
           'fetchPDBs']

def pathPDBFolder(folder=None, divided=False):
    """Returns or specify local PDB folder for storing PDB files downloaded from
    `wwPDB <http://www.wwpdb.org/>`_ servers.  Files stored in this folder can
    be accessed via :func:`.fetchPDB` from any working directory.  To release
    the current folder, pass an invalid path, e.g. ``folder=''``.

    If *divided* is **True**, the divided folder structure of wwPDB servers
    will be assumed when reading from and writing to the local folder.  For
    example, a structure with identifier **1XYZ** will be present as
    :file:`pdblocalfolder/yz/pdb1xyz.pdb.gz`.

    If *divided* is **False**, a plain folder structure will be expected and
    adopted when saving files.  For example, the same structure will be
    present as :file:`pdblocalfolder/1xyz.pdb.gz`.

    Finally, in either case, lower case letters will be used and compressed
    files will be stored."""

    if folder is None:
        folder = SETTINGS.get('pdb_local_folder')
        if folder:
            if isdir(folder):
                return folder, SETTINGS.get('pdb_local_divided', True)
            else:
                LOGGER.warn('PDB local folder {0} is not a accessible.'
                            .format(repr(folder)))
    else:
        if isdir(folder):
            folder = abspath(folder)
            LOGGER.info('Local PDB folder is set: {0}'.format(repr(folder)))
            if divided:
                LOGGER.info('wwPDB divided folder structure will be assumed.')
            else:
                LOGGER.info('A plain folder structure will be assumed.')
            SETTINGS['pdb_local_folder'] = folder
            SETTINGS['pdb_local_divided'] = bool(divided)
            SETTINGS.save()
        else:
            current = SETTINGS.pop('pdb_local_folder')
            if current:
                LOGGER.info('PDB folder {0} is released.'
                            .format(repr(current)))
                SETTINGS.pop('pdb_local_divided')
                SETTINGS.save()
            else:
                raise IOError('{0} is not a valid path.'.format(repr(folder)))

wwpdb.pathPDBFolder = pathPDBFolder

def pathPDBMirror(path=None, format=None):
    """Returns or specify PDB mirror path to be used by :func:`.fetchPDB`.
    To release the current mirror, pass an invalid path, e.g. ``path=''``.
    If you are keeping a partial mirror, such as PDB files in
    :file:`/data/structures/divided/pdb/` folder, specify *format*, which is
    ``'pdb'`` in this case."""

    if path is None:
        path = SETTINGS.get('pdb_mirror_path')
        format = SETTINGS.get('pdb_mirror_format', None)
        if path:
            if isdir(path):
                if format is None:
                    return path
                else:
                    return path, format
            else:
                LOGGER.warning('PDB mirror path {0} is not a accessible.'
                               .format(repr(path)))
    else:
        if isdir(path):
            path = abspath(path)
            LOGGER.info('Local PDB mirror path is set: {0}'
                        .format(repr(path)))
            SETTINGS['pdb_mirror_path'] = path
            SETTINGS['pdb_mirror_format'] = format
            SETTINGS.save()
        else:
            current = SETTINGS.pop('pdb_mirror_path')
            if current:
                LOGGER.info('PDB mirror {0} is released.'
                            .format(repr(current)))
                SETTINGS.save()
            else:
                raise IOError('{0} is not a valid path.'.format(repr(path)))


def fetchPDBfromMirror(*pdb, **kwargs):
    """Returns path(s) to PDB (default), PDBML, or mmCIF file(s) for specified
    *pdb* identifier(s).  If a *folder* is specified, files will be copied
    into this folder.  If *compressed* is **False**, files will decompressed.
    *format* argument can be used to get `PDBML <http://pdbml.pdb.org/>`_ and
    `mmCIF <http://mmcif.pdb.org/>`_ files: ``format='cif'`` will fetch an
    mmCIF file, and ``format='xml'`` will fetch a PDBML file.  If PDBML header
    file is desired, ``noatom=True`` argument will do the job."""

    mirror = pathPDBMirror()
    if mirror is None:
        raise IOError('no mirror path is set')


    try:
        mirror, mirror_format = mirror
    except ValueError:
        mirror_format = None

    format = str(kwargs.pop('format', 'pdb')).lower()

    if kwargs.get('check', True):
        identifiers = checkIdentifiers(*pdb)
    else:
        identifiers = list(pdb)

    if format == 'pdb':
        ftp_divided = 'data/structures/divided/pdb'
        ftp_pdbext = '.ent.gz'
        ftp_prefix = 'pdb'
        extension = '.pdb'
    elif format == 'xml':
        if bool(kwargs.pop('noatom', False)):
            ftp_divided = 'data/structures/divided/XML-noatom'
            ftp_pdbext = '-noatom.xml.gz'
            extension = '-noatom.xml'
        else:
            ftp_divided = 'data/structures/divided/XML'
            ftp_pdbext = '.xml.gz'
            extension = '.xml'
        ftp_prefix = ''
    elif format == 'cif':
        ftp_divided = 'data/structures/divided/mmCIF'
        ftp_pdbext = '.cif.gz'
        ftp_prefix = ''
        extension = '.cif'
    else:
        if format:
            raise ValueError('{0} is not a recognized format'
                             .format(repr(format)))
        else:
            raise ValueError('please specify a valid format')

    if mirror_format:
        if mirror_format.lower() != format:
            raise IOError('mirror contains only ' + mirror_format + ' files')
        ftp_divided = ''
    else:
        ftp_divided = join(*ftp_divided.split('/'))
    folder = kwargs.get('folder')
    compressed = kwargs.get('compressed', True)
    filenames = []
    append = filenames.append
    success = 0
    failure = 0
    for pdb in identifiers:
        if pdb is None:
            append(None)
            continue
        fn = join(mirror, ftp_divided, pdb[1:3],
                  ftp_prefix + pdb + ftp_pdbext)
        if isfile(fn):
            if folder or not compressed:
                if compressed:
                    fn = copyFile(fn, join(folder or '.',
                                             pdb + extension + '.gz'))
                else:
                    fn = gunzip(fn, join(folder or '.', pdb + extension))
            append(normpath(fn))
            success += 1
        else:
            append(None)
            failure += 1

    if len(identifiers) == 1:
        fn = filenames[0]
        if success:
            LOGGER.debug('PDB file is found in the local mirror ({0}).'
                            .format(sympath(fn)))
        return fn
    else:
        LOGGER.debug('PDB files found in the local mirror ({0} found, '
                        '{1} missed).'.format(success, failure))
        return filenames


def fetchPDB(*pdb, **kwargs):
    """Returns path(s) to PDB file(s) for specified *pdb* identifier(s).  Files
    will be sought in user specified *folder* or current working directory, and
    then in local PDB folder and mirror, if they are available.  If *copy*
    is set **True**, files will be copied into *folder*.  If *compressed* is
    **False**, all files will be decompressed into *folder*.  See :func:`pathPDBFolder` 
    and :func:`pathPDBMirror` for managing local resources, :func:`.fetchPDBviaFTP`
    and :func:`.fetchPDBviaHTTP` for downloading files from PDB servers."""

    if len(pdb) == 1 and isinstance(pdb[0], list):
        pdb = pdb[0]

    folder = kwargs.get('folder', '.')
    compressed = kwargs.get('compressed')
    format_ = kwargs.get('format', 'pdb')

    if format_ != 'emd':
        identifiers = checkIdentifiers(*pdb)
    else:
        identifiers = pdb

    # check *folder* specified by the user, usually pwd ('.')
    filedict = findPDBFiles(folder, compressed=compressed, 
                            format=format_)

    filenames = []
    not_found = []
    exists = 0
    for i, pdb in enumerate(identifiers):
        if pdb is None:
            filenames.append(None)
        elif pdb in filedict:
            filenames.append(filedict[pdb])
            exists += 1
        else:
            filenames.append(None)
            not_found.append((i, pdb))

    if not not_found:
        if len(filenames) == 1:
            filenames = filenames[0]
            if exists:
                LOGGER.debug('{0} file is found in working directory ({1}).'
                             .format(format_.upper(), sympath(filedict[pdb])))
        return filenames

    if not isWritable(folder):
        raise IOError('permission to write in {0} is denied, please '
                      'specify another folder'.format(folder))

    if compressed is not None and not compressed:
        filedict = findPDBFiles(folder, compressed=True)
        not_found, decompress = [], not_found
        for i, pdb in decompress:
            if pdb in filedict:
                fn = filedict[pdb]
                filenames[i] = gunzip(fn, splitext(fn)[0])
            else:
                not_found.append((i, pdb))

    if not not_found:
        return filenames[0] if len(identifiers) == 1 else filenames

    local_folder = pathPDBFolder()
    copy = kwargs.setdefault('copy', False)
    if local_folder:
        local_folder, is_divided = local_folder
        temp, not_found = not_found, []
        for i, pdb in temp:
            if is_divided:
                fn = join(local_folder, pdb[1:3], 'pdb' + pdb + '.pdb.gz')
            else:
                fn = join(local_folder, pdb + '.pdb.gz')
            if isfile(fn):
                if copy or not compressed and compressed is not None:
                    if compressed:
                        fn = copyFile(fn, join(folder, pdb + 'pdb.gz'))
                    else:
                        fn = gunzip(fn, join(folder, pdb + '.pdb'))
                filenames[i] = normpath(fn)
            else:
                not_found.append((i, pdb))

    if not not_found:
        if len(identifiers) == 1:
            fn = filenames[0]
            items = fn.split(pathsep)
            if len(items) > 5:
                fndisp = pathsep.join(items[:3] + ['...'] + items[-1:])
            else:
                fndisp = relpath(fn)
            LOGGER.debug('PDB file is found in the local folder ({0}).'
                            .format(fndisp))
            return fn
        else:
            return filenames

    if kwargs['copy'] or (compressed is not None and not compressed):
        kwargs['folder'] = folder

    downloads = [pdb for i, pdb in not_found]
    fns = None

    try:
        fns = fetchPDBfromMirror(*downloads, **kwargs)
    except IOError:
        pass
    else:
        if len(downloads) == 1: fns = [fns]
        temp, not_found = not_found, []
        for i, fn in enumerate(fns):
            if fn is None:
                not_found.append(temp[i])
            else:
                i, _ = temp[i]
                filenames[i] = fn

    if not not_found:
        return filenames[0] if len(identifiers) == 1 else filenames

    if fns:
        downloads = [pdb for i, pdb in not_found]

    fns = None

    tp = kwargs.pop('tp', None)
    if tp is not None:
        tp = tp.lower()
    
    if tp == 'http':
        try:
            fns = fetchPDBviaHTTP(*downloads, check=False, **kwargs)
        except Exception as err:
            LOGGER.warn('Downloading PDB files via HTTP failed '
                        '({0}).'.format(str(err)))
    elif tp == 'ftp':
        try:
            fns = fetchPDBviaFTP(*downloads, check=False, **kwargs)
        except Exception as err:
            LOGGER.warn('Downloading PDB files via FTP failed '
                        '({0}).'.format(str(err)))
    else:
        tryHTTP = False
        try:
            fns = fetchPDBviaFTP(*downloads, check=False, **kwargs)
        except Exception as err:
            tryHTTP = True
   
        if fns is None or isinstance(fns, list) and None in fns:
            tryHTTP = True
        elif isinstance(fns, list): 
            downloads = [not_found[i][1] for i in range(len(fns)) if fns[i] is None]
            if len(downloads) > 0: 
                tryHTTP = True
        if tryHTTP:
            LOGGER.info('Downloading PDB files via FTP failed, '
                        'trying HTTP.')
            try:
                fns = fetchPDBviaHTTP(*downloads, check=False, **kwargs)
            except Exception as err:
                LOGGER.warn('Downloading PDB files via HTTP also failed '
                            '({0}).'.format(str(err)))
    
    if len(downloads) == 1: fns = [fns]
    if fns:
        for i, fn in zip([i for i, pdb in not_found], fns):
            filenames[i] = fn

    return filenames[0] if len(identifiers) == 1 else filenames

def fetchPDBs(*pdb, **kwargs):
    """"Wrapper function to fetch multiple files from the PDB. 
    If no format is given, it tries PDB then mmCIF then EMD.
    
    :arg pdb: one PDB identifier or filename, or a list of them.
        If needed, PDB files are downloaded using :func:`.fetchPDB()` function.
    """

    n_pdb = len(pdb)
    if n_pdb == 0:
        raise ValueError('Please provide a PDB ID or filename')

    if n_pdb == 1:
        if isListLike(pdb[0]):
            pdb = pdb[0]
            n_pdb = len(pdb)

    fnames = []
    for p in pdb:
        format = kwargs.pop('format', None)
        
        if format is not None:
            filename = fetchPDB(p, format=format, **kwargs)

        else:
            filename = fetchPDB(p, **kwargs)

            if filename is None:
                filename = fetchPDB(p, format='cif', **kwargs)

            if filename is None:
                filename = fetchPDB(p, format='emd', **kwargs)

        fnames.append(filename)

    return fnames

def iterPDBFilenames(path=None, sort=False, unique=True, **kwargs):
    """Yield PDB filenames in *path* specified by the user or in local PDB
    mirror (see :func:`.pathPDBMirror`).  When *unique* is **True**, files
    one of potentially identical files will be yielded (e.g. :file:`1mkp.pdb`
    and :file:`pdb1mkp.ent.gz1`).  :file:`.pdb` and :file:`.ent` extensions,
    and compressed files are considered."""

    from re import compile, IGNORECASE

    format = kwargs.get('format')

    if path is None or kwargs.get('mirror') is True:
        if path is None:
            path = pathPDBMirror()
        if path is None:
            raise ValueError('path must be specified or PDB mirror path '
                             'must be set')
        if sort:
            pdbs = glob(join(path, 'data/structures/divided/pdb/',
                        '*/*.ent.gz'))
            pdbs.sort(reverse=kwargs.get('reverse'))
        else:
            pdbs = iglob(join(path, 'data/structures/divided/pdb/',
                        '*/*.ent.gz'))
        for fn in pdbs:
            yield fn
    else:
        unique=bool(unique)
        if unique:
            yielded = set()
        compressed = kwargs.get('compressed')
        if compressed is None:
            pdbext = compile('[.](pdb|ent)([.]gz)?$', IGNORECASE)
            cifext = compile('[.](cif)([.]gz)?$', IGNORECASE)
            emdext = compile('[.](emd|map|mrc)([.]gz)?$', IGNORECASE)
        elif compressed:
            pdbext = compile('[.](pdb|ent)[.]gz$', IGNORECASE)
            cifext = compile('[.](cif)[.]gz$', IGNORECASE)
            emdext = compile('[.](emd|map|mrc)[.]gz$', IGNORECASE)
        else:
            pdbext = compile('[.](pdb|ent)$', IGNORECASE)
            cifext = compile('[.](cif)$', IGNORECASE)
            emdext = compile('[.](emd|map|mrc)$', IGNORECASE)
        if format == 'cif':
            pdbext = cifext
        if format == 'emd':
            pdbext = emdext
        pdbs = [pdb for pdb in iglob(join(path, '*')) if pdbext.search(pdb)]
        if sort:
            pdbs.sort(reverse=kwargs.get('reverse'))
        for fn in pdbs:
            if unique:
                pdb = splitext(splitext(split(fn)[1])[0])[0]
                if len(pdb) == 7 and pdb.startswith('pdb'):
                    pdb = pdb[3:]
                if pdb in yielded:
                    continue
                else:
                    yielded.add(pdb)
            yield fn


def findPDBFiles(path, case=None, **kwargs):
    """Returns a dictionary that maps PDB filenames to file paths.  If *case*
    is specified (``'u[pper]'`` or ``'l[ower]'``), dictionary keys (filenames)
    will be modified accordingly.  If a PDB filename has :file:`pdb` prefix,
    it will be trimmed, for example ``'1mkp'`` will be mapped to file path
    :file:`./pdb1mkp.pdb.gz`).  If a file is present with multiple extensions,
    only one of them will be returned. See also :func:`.iterPDBFilenames`."""

    case = str(case).lower()
    upper = lower = False
    if case.startswith('u'):
        upper = True
    elif case.startswith('l'):
        lower = True

    pdbs = {}
    for fn in iterPDBFilenames(path, sort=True, reverse=True, **kwargs):
        fn = normpath(fn)
        pdb = splitext(split(fn)[1])[0]
        ending = splitext(split(fn)[1])[1]
        if ending == '.gz':
            pdb = splitext(pdb)[0]
        if len(pdb) == 7 and pdb.startswith('pdb'):
            pdb = pdb[3:]
        if upper:
            pdbs[pdb.upper()] = fn
        elif lower:
            pdbs[pdb.lower()] = fn
        else:
            pdbs[pdb] = fn

    return pdbs

