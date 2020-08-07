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

def _extend(self, arr, defval=0):
    mask = self.mask#.copy()
    if self.is3d():
        mask = np.repeat(mask, 3)

    n_true = np.sum(mask)
    N = len(mask)

    if arr.ndim == 1:
        whole_array = np.empty(N, dtype=arr.dtype)
        whole_array.fill(defval)
        whole_array[mask] = arr[:n_true]
    elif arr.ndim == 2:
        n, m = arr.shape
        whole_array = np.empty((N, m), dtype=arr.dtype)
        whole_array.fill(defval)
        #mask = np.expand_dims(mask, axis=1)
        #mask = mask.repeat(m, axis=1)
        whole_array[mask] = arr[:n_true, :]
    else: # only developers can trigger this case
        raise ValueError('arr can only be either 1D or 2D')
    return whole_array

def alignPDBEnsemble(ensemble, suffix='_aligned', outdir='.', gzip=False):
    """Align PDB files using transformations from *ensemble*, which may be
    a :class:`.PDBEnsemble` or a :class:`.PDBConformation` instance. Label of
    the conformation (see :meth:`~.PDBConformation.getLabel`) will be used to
    determine the PDB structure and model number.  First four characters of
    the label is expected to be the PDB identifier and ending numbers to be the
    model number.  For example, the :class:`.Transformation` from conformation
    with label *2k39_ca_selection_'resnum_<_71'_m116* will be applied to 116th
    model of structure **2k39**.  After applicable transformations are made,
    structure will be written into *outputdir* as :file:`2k39_aligned.pdb`.
    If ``gzip=True``, output files will be compressed.  Return value is
    the output filename or list of filenames, in the order files are processed.
    Note that if multiple models from a file are aligned, that filename will
    appear in the list multiple times."""

    if not isinstance(ensemble, (PDBEnsemble, PDBConformation)):
        raise TypeError('ensemble must be a PDBEnsemble or PDBConformation')
    if isinstance(ensemble, PDBConformation):
        ensemble = [ensemble]
    if gzip:
        gzip = '.gz'
    else:
        gzip = ''
    output = []
    pdbdict = {}
    for conf in ensemble:
        trans = conf.getTransformation()
        if trans is None:
            raise ValueError('transformations are not calculated, call '
                             '`superpose` or `iterpose`')
        label = conf.getLabel()

        pdb = label[:4]
        filename = pdbdict.get(pdb, fetchPDB(pdb))
        if filename is None:
            LOGGER.warning('PDB file for conformation {0} is not found.'
                           .format(label))
            output.append(None)
            continue
        LOGGER.info('Parsing PDB file {0} for conformation {1}.'
                    .format(pdb, label))

        acsi = None
        model = label.rfind('m')
        if model > 3:
            model = label[model+1:]
            if model.isdigit():
                acsi = int(model) - 1
            LOGGER.info('Applying transformation to model {0}.'
                        .format(model))

        if isinstance(filename, str):
            ag = parsePDB(filename)
        else:
            ag = filename

        if acsi is not None:
            if acsi >= ag.numCoordsets():
                LOGGER.warn('Model number {0} for {1} is out of range.'
                            .format(model, pdb))
                output.append(None)
                continue
            ag.setACSIndex(acsi)
        trans.apply(ag)
        outfn = os.path.join(outdir, pdb + suffix + '.pdb' + gzip)
        if ag.numCoordsets() > 1:
            pdbdict[pdb] = ag
        else:
            writePDB(outfn, ag)
        output.append(os.path.normpath(outfn))

    for pdb, ag in pdbdict.items():  # PY3K: OK
        writePDB(os.path.join(outdir, pdb + suffix + '.pdb' + gzip), ag)
    if len(output) == 1:
        return output[0]
    else:
        return output