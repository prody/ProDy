import numpy as np

def msaeye(msa, unique, turbo):
    tic1 = timeit.default_timer()
    length = msa.shape[1]
    number = msa.shape[0]
    # number = 5
    array = np.eye(int(number))

    seqs = []
    for i in range(number):
        seqs.append(msa[i,:])
    iseq = np.zeros((number, length), dtype=int)

    for i in range(0,number-1):
        if i == 0:
            for k in range(length):
                if ord(seqs[i][k])>90:
                    iseq[i,k]=ord(seqs[i][k])-96 if ord(seqs[i][k])-96 > 0 \
                              and ord(seqs[i][k])-96 < 26 else 0
                else:
                    iseq[i,k]=ord(seqs[i][k])-64 if ord(seqs[i][k])-64 > 0 \
                              and ord(seqs[i][k])-64 < 26 else 0
            for j in range(i+1,number):
                score=0.
                ncols=0.
                for k in range(length):
                    if ord(seqs[j][k])>90:
                        iseq[j,k]=ord(seqs[j][k])-96 if ord(seqs[j][k])-96 > 0 \
                                  and ord(seqs[j][k])-96 < 26 else 0
                    else:
                        iseq[j,k]=ord(seqs[j][k])-64 if ord(seqs[j][k])-64 > 0 and ord(seqs[j][k])-64 < 26 else 0
                    if iseq[i,k] or iseq[j,k]:
                        ncols += 1
                        if iseq[i,k]==iseq[j,k]:
                            score+=1
                array[i,j]=float(score)/ncols
                array[j,i]=array[i,j]
            # print iseq[0]
            # print seqs[0]
            # raw_input()
        else:
            for j in range(i+1,number):
                score=0.
                ncols=0.
                for k in range(length):
                    if iseq[i,k] or iseq[j,k]:
                        ncols += 1
                        if iseq[i,k]==iseq[j,k]:
                            score+=1
                array[i,j]= float(score)/ncols#float(sum((iseq[i] == iseq[j])*(iseq[i]*iseq[j]!=0))) / sum(iseq[i]*iseq[j]!=0)
                array[j,i]=array[i,j]

    toc1 = timeit.default_timer()
    elapsed1 = toc1 - tic1
    LOGGER.debug('Elapsed: %4.2fs'%elapsed1)

def showPerturbResponse(**kwargs):
    """ Plot the PRS matrix with the profiles along the right and bottom.

    If no PRS matrix or profiles are provided, these will be calculated first
    using the provided options with a provided model (e.g. ANM, GNM or EDA).
    So as to obtain different sensors and effectors, normMatrix=True by default.

    If atoms are provided then residue numbers can be used from there.
    *model* and *atoms* must have the same number of atoms. *atoms* must be an
    :class:`.AtomGroup` instance.

    :arg prs_matrix: a perturbation response matrix
    :type prs_matrix: :class:`~numpy.array`

    :arg effectiveness: an effectiveness profile from a PRS matrix
    :type effectiveness: :class:`~numpy.array`

    :arg sensitivity: a sensitivity profile from a PRS matrix
    :type sensitivity: :class:`~numpy.array`

    :arg model: any object with a calcCovariance method
        e.g. :class:`.ANM` instance
    :type model: NMA

    :arg atoms: a :class: `AtomGroup` instance
    :type atoms: AtomGroup

    :arg returnData: whether to return data for further analysis
        default is False
    :type returnData: bool
    
    :arg percentile: percentile argument for showAtomicMatrix
    :type percentile: float

    Return values are prs_matrix, effectiveness, sensitivity, ax1, ax2, im, ax3, ax4
    The PRS matrix, effectiveness and sensitivity will not be returned if provided. 
    If returnData is False then only the last five objects are returned.
    """

    import matplotlib.pyplot as plt
    import matplotlib

    prs_matrix = kwargs.get('prs_matrix')
    effectiveness = kwargs.get('effectiveness')
    sensitivity = kwargs.get('sensitivity')
    model = kwargs.pop('model')
    atoms = kwargs.get('atoms')
    returnData = kwargs.pop('returnData', False)

    if atoms is None:

        if prs_matrix is None:
            if model is None:
                raise ValueError('Please provide a PRS matrix or model.')
            else:
                prs_matrix = calcPerturbResponse(model)

        if effectiveness is None or sensitivity is None:
            effectiveness, sensitivity = calcPerturbResponseProfiles(prs_matrix)
        else:
            returnData = False

        showMatrix_returns = showAtomicMatrix(prs_matrix, effectiveness, sensitivity, **kwargs)

    else:
        if not isinstance(atoms, AtomGroup) and not isinstance(atoms, Selection):
            raise TypeError('atoms must be an AtomGroup instance')
        elif model is not None and atoms.numAtoms() != model.numAtoms():
            raise ValueError('model and atoms must have the same number atoms')

        if prs_matrix is None: 
            if model is None:
                raise ValueError('Please provide a PRS matrix or model.')
            atoms, prs_matrix = calcPerturbResponse(model=model,atoms=atoms)

        if effectiveness is None or sensitivity is None:
            effectiveness, sensitivity = calcPerturbResponseProfiles(prs_matrix, atoms=atoms)

        showMatrix_returns = showAtomicMatrix(prs_matrix, effectiveness, sensitivity, **kwargs)

    if not returnData:
        return showMatrix_returns
    elif kwargs.get('prs_matrix') is not None:
       if atoms is not None:
           return atoms, effectiveness, sensitivity, showMatrix_returns
       else:
           return effectiveness, sensitivity, showMatrix_returns
    else:
        return prs_matrix, effectiveness, sensitivity, showMatrix_returns

def calcPerturbResponseProfiles(prs_matrix, atoms=None):
    """ Calculate the effectiveness and sensitivity
    profiles, which are the averages over the rows
    and columns of the PRS matrix.

    :arg prs_matrix: a perturbation response matrix
    :type prs_matrix: ndarray 

    When an *atoms* instance is given, the profiles will be added as data, 
    which can be retrieved with ``atoms.getData('effectiveness')`` and 
    ``atoms.getData('sensitivity')``. 
    """

    effectiveness = np.mean(prs_matrix, axis=1)
    sensitivity = np.mean(prs_matrix, axis=0)

    if atoms is not None:
        try:
            ag = atoms.getAtomGroup()
            defdata = np.zeros(ag.numAtoms(), dtype=float)
            ag.setData('effectiveness', defdata.copy())
            ag.setData('sensitivity', defdata.copy())
        except AttributeError:
            pass
        atoms.setData('effectiveness', effectiveness)
        atoms.setData('sensitivity', sensitivity)
        
    return effectiveness, sensitivity

def writePerturbResponsePDB(prs_matrix, pdbIn=None, **kwargs):
    """ Write the average response to perturbation of
    a particular residue (a row of a perturbation response matrix)
    or the average effect of perturbation of a particular residue
    (a column of a normalized perturbation response matrix)
    into the b-factor field of a PDB file for visualisation in a
    molecular graphics program.
    If no chain is given this will be done for that residue in all chains.

    If no residue number is given then the effectiveness and sensitivity
    profiles will be written out instead. These two profiles are also returned
    as arrays for further analysis if they aren't already provided.

    :arg prs_matrix: a perturbation response matrix 
        or a :class:`.AtomGroup` object with a PRS matrix associated as data
    :type prs_matrix: array or :class:`.AtomGroup`

    :arg pdbIn: file name for the input PDB file where you would like the PRS
        data mapped
    :type pdbIn: str

    :arg pdbOut: a list of file names (enclosed in square
        brackets) for the output PDB file, default is to append
        the chain and residue info (name and number) onto the pdbIn stem.
        The input for pdbOut can also be used as a stem if you enter a 
        single string enclosed in quotes.
        If no residue number is supplied, chain is ignored and the default 
        is to append '_effectiveness' and '_sensitivity' onto the stem.
    :type pdbOut: list

    :arg chain: chain identifier for the residue of interest, default is all chains
        If you want to analyse residues in a subset of chains, concatentate them
        together e.g. 'AC'
    :type chain: str

    :arg resnum: residue number for the residue of interest
    :type resnum: int

    :arg direction: the direction you want to use to read data out
        of the PRS matrix for plotting: the options are 'effect' or 'response'.
        Default is 'effect'.
        A row gives the effect on each residue of peturbing the specified 
        residue.
        A column gives the response of the specified residue to perturbing 
        each residue.
        If no residue number is provided then this option will be ignored
    :type direction: str

    :arg returnData: whether to return effectiveness and sensitivity for analysis
        default is False
    :type returnProfiles: bool

    :arg effectiveness: effectiveness profile
    :type array

    :arg sensitivity: sensitivity profile
    :type array
    """

    if not isinstance(prs_matrix, np.ndarray):
        try:
            prs_matrix = prs_matrix.getData('prs_matrix')
        except:
            raise TypeError('Please provide a valid PRS matrix in numpy ndarray format.')

    try:
        fi = open(pdbIn,'r')
        lines = fi.readlines()
        fi.close()
    except:
        raise PRSMatrixParseError('Please provide a valid file name for the input PDB.')
 
    chain = kwargs.get('chain', None)

    structure = parsePDB(pdbIn,subset='ca')
    structure.setData('prs_matrix',prs_matrix)

    hv = structure.getHierView()
    chains = []
    for i in range(len(list(hv))):
        chainAg = list(hv)[i]
        chains.append(chainAg.getChids()[0])

    chains = np.array(chains)
    if chain is None:
        chain = ''.join(chains)

    resnum = kwargs.get('resnum', None)
    pdbOut = kwargs.get('pdbOut', None)
    if pdbOut is None:
        out_stem = pdbIn.split('.')[0]
    elif type(pdbOut) is str:
        out_stem = pdbOut.split('.')[0]
        pdbOut = None

    if resnum is None:
        effectiveness = kwargs.get('effectiveness', None)
        sensitivity = kwargs.get('sensitivity', None)
        if effectiveness is None or sensitivity is None:
            effectiveness = np.mean(prs_matrix, axis=1)
            sensitivity = np.mean(prs_matrix, axis=0)

        structure.setData('effectiveness', effectiveness)
        structure.setData('sensitivity', sensitivity)

        file_effs_name = '{0}_effectiveness.pdb'.format(out_stem)
        file_sens_name = '{0}_sensitivity.pdb'.format(out_stem)
        fileEffs = open(file_effs_name,'w')
        fileSens = open(file_sens_name,'w')

        for line in lines:            
            if line.find('ATOM') != 0 and line.find('HETATM') != 0 and line.find('ANISOU') != 0:
                fileEffs.write(line)                    
                fileSens.write(line)
            elif line.find('ATOM') == 0:
                fileEffs.write(line[:60] + '{:6.2f}'.format(float(structure.select( \
                               'chain {0} and resnum {1}'.format(line[21],line[22:26])) \
                               .getData('effectiveness')) * 100/np.max( \
                               structure.getData('effectiveness'))) + line[66:])
                fileSens.write(line[:60] + '{:6.2f}'.format(float(structure.select( \
                               'chain {0} and resnum {1}'.format(line[21],line[22:26])) \
                               .getData('sensitivity')) * 100/np.max( \
                               structure.getData('sensitivity'))) + line[66:])
            elif line.find('HETATM') == 0:
                fileEffs.write(line[:60] + '  0.00' + line[66:])
                fileSens.write(line[:60] + '  0.00' + line[66:])
                      
        fileEffs.close()
        fileSens.close()
        LOGGER.info('The effectiveness and sensitivity profiles were written' \
                    ' to {0} and {1}.'.format(file_effs_name,file_sens_name))

        returnData = kwargs.get('returnData',False)
        if returnData:
            return structure, effectiveness, sensitivity
        else:
            return
 
    direction = kwargs.get('direction','effect')
    for n in range(len(chain)):
        if not chain[n] in chains:
            raise PRSMatrixParseError('Chain {0} was not found in {1}'.format(chain[n], pdbIn))

    if pdbOut is None:
        pdbOut = []
        for c in chain:
            pdbOut.append('{0}_{1}_{2}{3}_{4}.pdb' \
                          .format(out_stem, c, \
                                  str(structure.select('chain {0} and resnum {1}' \
                                      .format(c, resnum)).getResnames()), \
                                  resnum, direction))

    for c in chain:
        fo = open(pdbOut[n],'w')
        for line in lines:
            if line.find('ATOM') != 0 and line.find('HETATM') != 0 and line.find('ANISOU') != 0:
                fo.write(line)
            elif line.find('ATOM') == 0:
                if direction == 'effect':
                    fo.write(line[:60] + '{:6.2f}'.format(float(structure.getData('prs_matrix') \
                                         [structure.select('chain {0} and resnum {1}' \
                                          .format(c, resnum)).getResindices(), \
                                          structure.select('chain {0} and resnum {1}' \
                                          .format(line[21], line[22:26])).getResindices()])*100) \
                             + line[66:])
                else:
                    fo.write(line[:60] + '{:6.2f}'.format(float(structure.getData('prs_matrix') \
                                         [structure.select('chain {0} and resnum {1}' \
                                          .format(line[21], line[22:26])).getResindices(), \
                                          structure.select('chain {0} and resnum {1}' \
                                          .format(c, resnum)).getResindices()])*100) \
                             + line[66:])
            elif line.find('HETATM') == 0:
                fo.write(line[:60] + '  0.00' + line[66:])

        LOGGER.info('Perturbation responses for specific residues were written' \
                    ' to {0}.'.format(', '.join(pdbOut)))

def parsePerturbResponseMatrix(prs_matrix_file, norm=False):
    """Parses a perturbation response matrix from a file into a numpy ndarray.

    :arg prs_matrix_file: name of the file containing a PRS matrix
    :type prs_matrix_file: str

    :arg norm: whether to normalize the PRS matrix after parsing it.
        Default is False. If you used an old version of the script 
        and didn't normalize before saving, set this to True.
    :type norm: bool

    """
    fmat = open(prs_matrix_file, 'rb')
    matlines = fmat.readlines()
    fmat.close()

    prs_matrix = []
    for line in matlines:
        prs_matrix.append([float(entry) for entry in line.split()])

    prs_matrix = np.array(prs_matrix)

    if norm:
       # normalize the PRS matrix
       self_dp = np.diag(prs_matrix)  # using self displacement (diagonal of
                              # the original matrix) as a
                              # normalization factor
       self_dp = self_dp.reshape(len(prs_matrix), 1)
       prs_matrix = prs_matrix / np.repeat(self_dp, len(prs_matrix), axis=1)

    return prs_matrix

class PRSMatrixParseError(Exception):
    pass

def buildDaliEnsemble(PDBs, record):
    daliInfo = record._alignPDB

    n_confs = len(PDBs)
    
    ref_pdb_ca = PDBs[0]
    ref_chain = list(ref_pdb_ca.getHierView().iterChains())[0]
    ref_indices_set = set(range(len(ref_chain)))
    ensemble = PDBEnsemble('Dali ensemble - ' + record.getTitle())
    ensemble.setAtoms(ref_chain)
    ensemble.setCoords(ref_chain)
    
    LOGGER.progress('Building PDB ensemble for {0} conformations from Dali...'
                    .format(n_confs), n_confs, '_prody_buildDaliEnsemble')

    for i, pdb in enumerate(PDBs):
        pdb_chain = pdb.getTitle()[:5]
        temp_dict = daliInfo[pdb_chain]
        
        sel_pdb_ca = PDBs[i]
        map_ref = temp_dict['map_ref']
        map_sel = temp_dict['map_sel']
        dum_sel = list(ref_indices_set - set(map_ref))
        atommap = AtomMap(sel_pdb_ca, indices=map_sel, mapping=map_ref, dummies=dum_sel)
        ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'), degeneracy=True)

        LOGGER.update(i, label='_prody_buildDaliEnsemble')
    LOGGER.finish()

    try:
        ensemble.iterpose()
    except:
        LOGGER.warn('failed to iterpose the ensemble.')
        
    return ensemble

def fetchCATH(filename, ftp_host=None, ftp_path=None, **kwargs):
    """Downloads CATH file via FTP."""
    if ftp_host == None:
        ftp_host = 'orengoftp.biochem.ucl.ac.uk'
    if ftp_path == None:
        ftp_path = '/cath/releases/daily-release/newest/'
    from ftplib import FTP
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


def fetchCATH(filename, ftp_host=None, ftp_path=None, **kwargs):
    """Downloads CATH file via FTP."""
    if ftp_host == None:
        ftp_host = 'orengoftp.biochem.ucl.ac.uk'
    if ftp_path == None:
        ftp_path = '/cath/releases/daily-release/newest/'
    from ftplib import FTP
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

# ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/daily-release/newest/

# fetchCATH('cath-b-newest-names.gz')
# cath_id2name = buildCATHNameDict('cath-b-newest-names.gz')

# fetchCATH('cath-b-newest-all.gz')
# pdbChain2CATH = buildPDBChainCATHDict('cath-b-newest-all.gz')

def extend(model, nodes, atoms):
    """Returns mapping indices and an :class:`.AtomMap`."""

    try:
        n_atoms = model.numAtoms()
        is3d = model.is3d()
    except AttributeError:
        raise ValueError('model must be an NMA instance')

    try:
        n_nodes = nodes.numAtoms()
        i_nodes = nodes.iterAtoms()
    except AttributeError:
        raise ValueError('nodes must be an Atomic instance')

    if n_atoms != n_nodes:
        raise ValueError('atom numbers must be the same')

    if not nodes in atoms:
        raise ValueError('nodes must be a subset of atoms')

    atom_indices = []
    indices = []
    get = HierView(atoms).getResidue

    for i, node in enumerate(i_nodes):
        res = get(node.getChid() or None, node.getResnum(),
                  node.getIcode() or None, node.getSegname() or None)
        if res is None:
            raise ValueError('atoms must contain a residue for all atoms')
        atom_indices.append(res._getIndices())
        if is3d:
            indices.append(list(range(i*3, (i+1)*3)) * len(res))
        else:
            indices.append([i] * len(res))
    atom_indices = np.concatenate(atom_indices)
    indices = np.concatenate(indices)

    try:
        ag = atoms.getAtomGroup()
    except AttributeError:
        ag = atoms
    atommap = AtomMap(ag, atom_indices, atoms.getACSIndex(),
                      title=str(atoms), intarrays=True)
    return indices, atommap

def extendAtomicData(data, nodes, atoms):
    """Extend a coarse grained data obtained for *nodes* to *atoms*.

    :arg data: any data array
    :type data: :class:`~numpy.ndarray`

    :arg nodes: a set of atoms that has been used
        as nodes in data generation
    :type nodes: :class:`.Atomic`

    :arg atoms: atoms to be selected from
    :type atoms: :class:`.Atomic`

    """
    from collections import Counter

    try:
        data = np.asarray(data)
    except:
        raise TypeError('The data must be array-like.')

    if not isinstance(nodes, Atomic):
        raise TypeError('nodes must be an Atomic instance')

    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance')

    nnodes = nodes.numAtoms()

    is3d = False
    if len(data) != nnodes:
        if data.shape[0] == nnodes * 3:
            is3d = True
        else:
            raise ValueError('data and atoms must have the same size')

    indices = nodes.getResindices()
    if is3d:
        indices = np.array([[i*3, i*3+1, i*3+2] 
                        for i in indices]
                        ).reshape(3*len(indices))

    data_ext = []
    resid_counter = Counter(atoms.getResindices())
    for i in indices:
        data_ext.extend(resid_counter.values()[i]*[data[i]])

    resid_selstr = ' '.join([str(resid) for resid in nodes.getResindices()])
    rest = atoms.select('not resid {0}'.format(resid_selstr))
    data_ext.extend(np.zeros(rest.numAtoms()))
        
    return data_ext


def refineEnsemble(ens, lower=.5, upper=10.):
    """Refine a PDB ensemble based on RMSD criterions.""" 

    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    from collections import Counter

    ### calculate pairwise RMSDs ###
    RMSD = ens.getRMSDs(pairwise=True)

    # convert the RMSD table to the compressed form
    v = squareform(RMSD)

    ### apply upper threshold ###
    Z_upper = linkage(v, method='complete')
    labels = fcluster(Z_upper, upper, criterion='distance')
    most_common_label = Counter(labels).most_common(1)[0][0]
    I = np.where(labels==most_common_label)[0]

    ### apply lower threshold ###
    Z_lower = linkage(v, method='single')
    labels = fcluster(Z_lower, lower, criterion='distance')
    uniq_labels = np.unique(labels)

    clusters = []
    for label in uniq_labels:
        indices = np.where(labels==label)[0]
        clusters.append(indices)

    J = np.ones(len(clusters), dtype=int) * -1
    rmsd = None
    for i, cluster in enumerate(clusters):
        if len(cluster) > 0:
            # find the conformations with the largest coverage 
            # (the weight of the ref should be 1)
            weights = [ens[j].getWeights().sum() for j in cluster]
            js = np.where(weights==np.max(weights))[0]

            # in the case where there are multiple structures with the same weight,
            # the one with the smallest rmsd wrt the ens._coords is selected. 
            if len(js) > 1:
                # rmsd is not calulated unless necessary for the sake of efficiency
                rmsd = ens.getRMSDs() if rmsd is None else rmsd
                j = js[np.argmin(rmsd[js])]
            else:
                j = js[0]
            J[i] = cluster[j]
        else:
            J[i] = cluster[0]

    ### refine ensemble ###
    K = np.intersect1d(I, J)

    reens = ens[K]

    return reens

def showVarianceBar(mode_ensemble, highlights=None, **kwargs):

    from matplotlib.pyplot import figure, gca, annotate, subplots_adjust, plot
    from matplotlib.figure import Figure
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize, NoNorm
    from matplotlib import cm, colors
    
    fig = kwargs.pop('figure', None)

    if isinstance(fig, Figure):
        fig_num = fig.number
    elif fig is None or isinstance(fig, (int, str)):
        fig_num = fig
    else:
        raise TypeError('figure can be either an instance of matplotlib.figure.Figure '
                        'or a figure number.')
    if SETTINGS['auto_show']:
        if fig_num is None:
            figure(figsize=(6, 2))
        else:
            figure(fig_num)
    elif fig_num is not None:
        figure(fig_num)
    ax = gca()

    # adjust layouts
    box = ax.get_position()
    _, _, _, height = box.bounds
    ratio = 2.5
    box.y1 = box.y0 + height/ratio
    #box.y0 += height/7.
    ax.set_position(box)

    fract = kwargs.pop('fraction', True)

    #defarrow = {'width':1, 'headwidth':2, 
    #            'facecolor':'black',
    #            'headlength': 4}
    defarrow = {'arrowstyle': '->'}
    arrowprops = kwargs.pop('arrowprops', defarrow)

    if fract:
        sig = calcSignatureFractVariance(mode_ensemble)
    else:
        sig = mode_ensemble.getVariances() 

    variances = sig.getArray().sum(axis=1)
    #meanVar = variances.mean()
    #stdVar = variances.std()
    
    #variances = (variances - meanVar)/stdVar

    maxVar = variances.max()
    minVar = variances.min()

    cmap = kwargs.pop('cmap', 'jet')
    norm = Normalize(vmin=minVar, vmax=maxVar)
    cb = ColorbarBase(ax, cmap=cmap, norm=norm,
                      orientation='horizontal')

    if not highlights:
        highlights = []

    indices = []; labels = []
    ens_labels = mode_ensemble.getLabels()
    for hl in highlights:
        if isinstance(hl, str):
            if not ens_labels:
                raise TypeError('highlights should be a list of integers because '
                                    'mode_ensemble has no label')
            indices.append(ens_labels.index(hl))
            labels.append(hl)
        else:
            try:
                index = int(hl)
            except:
                raise TypeError('highlights should be a list of integers or strings') 
            indices.append(index)
            if ens_labels:
                labels.append(ens_labels[index])
            else:
                labels.append(str(index))

    annotations = []
    for i, label in zip(indices, labels):
        x = norm(variances[i])
        an = annotate(label, xy=(x, 1), xytext=(x, ratio), arrowprops=arrowprops)
        annotations.append(an)

    for i in range(len(variances)):
        x = norm(variances[i])
        plot([x, x], [0, 1], 'w')

    cb.set_label('Variances')

    if SETTINGS['auto_show']:
        showFigure()
    return cb, annotations
    
def mapChainByChain(atoms, target, **kwargs):
    """This function is similar to :func:`.mapOntoChain` but correspondence 
    of chains is found by their chain identifiers. 
    
    :arg atoms: atoms to be mapped onto *target*
    :type atoms: :class:`.Atomic`
    
    :arg target: reference structure for mapping
    :type target: :class:`.Atomic`
    
    :arg return_all: whether to return all mappings.
        If False, only mappings for the first chain will be returned. 
        Default is **True**
    :arg return_all: bool

    :arg correspondence: chain IDs in atoms corresponding to those in ref
        Default is to use the same chain IDs as in ref.
    :type correspondence: str, list, dict
    """

    mappings = []

    if isinstance(target, AtomGroup):
        chs_ref_ag = target.iterChains()
    else:
        chs_ref_ag = target.getAtomGroup().iterChains()

    id_atm = atoms.getTitle()
    id_ref = target.getTitle()
    
    chs_atm = [chain for chain in atoms.getHierView().iterChains()]
    chs_ref = [chain for chain in target.getHierView().iterChains()]

    corr_input = kwargs.get('correspondence', None)

    if isinstance(corr_input, dict):
        correspondence = corr_input
    elif corr_input is None:
        correspondence = {}
    elif isinstance(corr_input, str):
        correspondence = {}
        correspondence[atoms.getTitle()] = corr_input
    else:
        correspondence = {}
        try:
            correspondence[id_atm] = corr_input[0]
            correspondence[id_ref] = corr_input[1]
        except (IndexError, TypeError):
            raise TypeError('correspondence should be a dict with keys being titles of atoms and ref, '
                            'and values are str indicating chID correspondences')

    if not id_atm in correspondence:
        correspondence[id_atm] = ''.join([chain.getChid() for chain in chs_atm])

    if not id_ref in correspondence:
        correspondence[id_ref] = ''.join([chain.getChid() for chain in chs_ref_ag])

    corr_tar = correspondence[id_atm]
    corr_ref = correspondence[id_ref]
    for chain in chs_ref:
        try:
            i = corr_ref.index(chain.getChid())
            chid = corr_tar[i]
        except ValueError:
            continue

        for target_chain in chs_atm:
            if target_chain.getChid() == chid:
                mappings_ = mapOntoChainByAlignment(target_chain, chain, **kwargs)
                if len(mappings_):
                    mappings.append(mappings_[0])
        
    return mappings
    