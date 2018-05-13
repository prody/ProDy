def showAtomicMatrix(matrix=None, x_array=None, y_array=None, **kwargs):
    """Show a matrix using :meth:`~matplotlib.axes.Axes.imshow`. Curves on x- and y-axis can be added.
    The first return value is the :class:`~matplotlib.axes.Axes` object for the upper plot, and the second
    return value is equivalent object for the left plot. The third return value is 
    the :class:`~matplotlib.image.AxesImage` object for the matrix plot. The last return value is the 
    :class:`~matplotlib.axes.Axes` object for the color bar.

    :arg matrix: Matrix to be displayed.
    :type matrix: :class:`~numpy.ndarray`

    :arg x_array: Data to be plotted above the matrix.
    :type x_array: :class:`~numpy.ndarray`

    :arg y_array: Data to be plotted on the left side of the matrix.
    :type y_array: :class:`~numpy.ndarray`

    :arg percentile: A percentile threshold to remove outliers, i.e. only showing data within *p*-th 
                     to *100-p*-th percentile.
    :type percentile: float

    :arg vmin: Minimum value that can be used together with vmax 
               as an alternative way to remove outliers
    :type vmin: float

    :arg vmax: Maximum value that can be used together with vmin 
               as alternative way to remove outliers
    :type vmax: float

    :arg atoms: a :class: `AtomGroup` instance for matching 
        residue numbers and chain IDs. 
    :type atoms: :class: `AtomGroup`

    :arg num_div: the number of divisions for each chain
        default 2
    :type num_div: int

    :arg resnum_tick_labels: residue number labels in place of num_div.
         A list can be used to set the same labels on all chains or 
         a dictionary of lists to set different labels for each chain
    :type resnum_tick_labels: list or dictionary

    :arg add_last_resi: whether to add a label for the last residue
        default False
    :type add_last_resi: bool

    :arg label_size: size for resnum labels
        default is 6, which works well for 4 residues on 4 chains
    :type label_size: int
    """ 

    num_div = kwargs.pop('num_div',2)
    resnum_tick_labels = kwargs.pop('resnum_tick_labels',None)
    add_last_resi = kwargs.pop('add_last_resi',False)

    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from matplotlib.collections import LineCollection
    from matplotlib.pyplot import figure, imshow

    if matrix is None:
        raise TypeError('You need to provide a matrix.') 
    elif len(np.shape(matrix)) != 2:
        raise ValueError('The matrix must be a 2D array.')

    p = kwargs.pop('percentile', None)
    vmin = kwargs.pop('vmin', None)
    vmax = kwargs.pop('vmax', None)

    if vmin is None and vmax is None and p is not None:
        vmin = np.percentile(matrix, p)
        vmax = np.percentile(matrix, 100-p)
    
    W = 10
    H = 10
    aspect = 'auto'

    if x_array is not None and y_array is not None:
        nrow = 3; ncol = 5
        i = 1; j = 1
        width_ratios = [1, W, 0.2]
        height_ratios = [1, H, 0.2]
    elif x_array is not None and y_array is None:
        nrow = 3; ncol = 4
        i = 1; j = 0
        width_ratios = [W, 0.2]
        height_ratios = [1, H, 0.2]
    elif x_array is None and y_array is not None:
        nrow = 2; ncol = 5
        i = 0; j = 1
        width_ratios = [1, W, 0.2]
        height_ratios = [H, 0.2]
    else:
        nrow = 2; ncol = 4
        i = 0; j = 0
        width_ratios = [W, 0.2]
        height_ratios = [H, 0.2]

    main_index = (i,j)
    upper_index = (i-1,j)
    lower_index = (i+1,j)
    left_index = (i,j-1)
    right_index = (i,j+1)

    outer = GridSpec(1, 3, width_ratios = [sum(width_ratios), 1, 4], hspace=0., wspace=0.2) 

    gs = GridSpecFromSubplotSpec(nrow, ncol-2, subplot_spec = outer[0], width_ratios=width_ratios,
                                 height_ratios=height_ratios, hspace=0., wspace=0.)

    gs_bar = GridSpecFromSubplotSpec(nrow-1, 1, subplot_spec = outer[1], height_ratios=height_ratios[:-1], hspace=0., wspace=0.)
    gs_legend = GridSpecFromSubplotSpec(nrow-1, 1, subplot_spec = outer[2], height_ratios=height_ratios[:-1], hspace=0., wspace=0.)

    if SETTINGS['auto_show']:
        fig = plt.figure(figsize=[9.5,6]) 
    axes = []

    atoms = kwargs.pop('atoms', None)
    if atoms is not None:
        if not isinstance(atoms, AtomGroup) and not isinstance(atoms, Selection):
            raise TypeError('atoms must be an AtomGroup instance')
    
    cmap = kwargs.pop('cmap', 'jet')
    label_size = kwargs.pop('label_size', 6)
 
    ax1 = ax2 = ax3 = ax4 = ax5 = ax6 = ax7 = None
    if nrow > 2:
        y1 = x_array
        x1 = np.arange(len(y1))
        ax1 = plt.subplot(gs[upper_index])
        points1 = np.array([x1, y1]).T.reshape(-1, 1, 2)
        segments1 = np.concatenate([points1[:-1], points1[1:]], axis=1)
        lc1 = LineCollection(segments1, array=y1, linewidths=1, cmap=cmap)
        ax1.add_collection(lc1)

        ax1.set_xlim(x1.min(), x1.max())
        ax1.set_ylim(y1.min(), y1.max())
        ax1.axis('off')

    if ncol > 4:
        x2 = y_array
        y2 = np.arange(len(x2))
        ax2 = plt.subplot(gs[left_index])
        points2 = np.array([x2, y2]).T.reshape(-1, 1, 2)
        segments2 = np.concatenate([points2[:-1], points2[1:]], axis=1)
        lc2 = LineCollection(segments2, array=x2, linewidths=1, cmap=cmap)
        ax2.add_collection(lc2)

        ax2.set_xlim(x2.min(), x2.max())
        ax2.set_ylim(y2.min(), y2.max())
        ax2.axis('off')
        ax2.invert_xaxis()

    ax3 = plt.subplot(gs[main_index])
    im = imshow(matrix, aspect=aspect, vmin=vmin, vmax=vmax, cmap=cmap, **kwargs)
    ax3.yaxis.tick_right()

    ax4 = plt.subplot(gs_bar[-1])
    plt.colorbar(cax=ax4)
 
    if atoms is not None:

        # Add bars along the bottom and right that are colored by chain and numbered with residue

        ax5 = plt.subplot(gs[lower_index])
        ax6 = plt.subplot(gs[right_index])

        n = 0
        resnum_tick_locs = []
        resnum_tick_labels_list = []

        if resnum_tick_labels is None:
            resnum_tick_labels = []
            user_set_labels = False
        elif type(resnum_tick_labels) is list:
            user_set_labels = list
        elif type(resnum_tick_labels) is dict:
            user_set_labels = dict
        else:
            raise TypeError('The resnum tick labels should be a list or dictionary of lists')

        chain_colors = 'gcmyrwbk'
        chain_handles = []
        for i in atoms.getHierView().iterChains():
            
            chain_handle, = ax5.plot([i.getResindices()[0], i.getResindices()[-1]], [0, 0], \
                                     '-', linewidth=3, color=chain_colors[n], label=str(i))
            chain_handles.append(chain_handle)

            ax6.plot([0,0], [np.flip(i.getResindices(),0)[0], np.flip(i.getResindices(),0)[-1]], \
                     '-', linewidth=3, color=chain_colors[n], label=str(i))

            if not user_set_labels:
                for j in range(num_div):
                    resnum_tick_locs.append(i.getResindices()[i.numAtoms()/num_div*j])
                    resnum_tick_labels.append(i.getResnums()[i.numAtoms()/num_div*j])
            elif user_set_labels is list:
                for j in resnum_tick_labels:
                    resnum_tick_locs.append(i.getResindices()[np.where(i.getResnums() == j)[0][0]])
                    resnum_tick_labels_list.append(j)
            else:
                for k in resnum_tick_labels.keys():
                    if i.getChids()[0] == k:
                       for j in resnum_tick_labels[k]:
                           resnum_tick_locs.append(i.getResindices()[np.where(i.getResnums() == j)[0][0]])
                           resnum_tick_labels_list.append(j)

            n += 1

        ax7 = plt.subplot(gs_legend[-1])
        plt.legend(handles=chain_handles, loc=2, bbox_to_anchor=(0.25, 1))
        ax7.axis('off')

        if add_last_resi:
            resnum_tick_locs.append(atoms.getResindices()[-1])
            resnum_tick_labels_list.append(atoms.getResnums()[-1])

        resnum_tick_locs = np.array(resnum_tick_locs)
        resnum_tick_labels = np.array(resnum_tick_labels_list)

        ax3.axis('off')

        ax5.set_xticks(resnum_tick_locs)
        ax5.set_xticklabels(resnum_tick_labels)
        ax5.tick_params(labelsize=label_size)
        ax5.set_yticks([])

        ax6.set_xticks([])
        ax6.yaxis.tick_right()
        ax6.set_yticks(resnum_tick_locs)
        ax6.set_yticklabels(resnum_tick_labels)
        ax6.tick_params(labelsize=label_size)

        ax5.set_xlim([-0.5, len(matrix)+0.5])
        ax6.set_ylim([-0.5, len(matrix.T)+0.5])

    if SETTINGS['auto_show']:
        showFigure()
 
    plt.sca(ax3)

    return ax1, ax2, im, ax3, ax4, ax5, ax6, ax7

def showPlot(y, **kwargs):

    """
    Show a plot with the option to include chain color bars using provided atoms.
    
    :arg atoms: a :class: `AtomGroup` instance for matching 
        residue numbers and chain IDs. 
    :type atoms: :class: `AtomGroup`
    
    :arg num_div: the number of divisions for each chain
        default 2
    :type num_div: int

    :arg resnum_tick_labels: residue number labels in place of num_div.
         A list can be used to set the same labels on all chains or 
         a dictionary of lists to set different labels for each chain
    :type resnum_tick_labels: list or dictionary

    :arg add_last_resi: whether to add a label for the last residue
        default False
    :type add_last_resi: bool

    :arg label_size: size for resnum labels
        default is 6, which works well for 4 residues on 4 chains
    :type label_size: int

    :arg overlay_chains: overlay the chains rather than having them one after another
        default False
    :type overlay_chains: bool

    :arg domains: color the bar at the bottom by domains rather than chains
        default False
    :type domains: bool
    """
    atoms = kwargs.pop('atoms',None)
    overlay_chains = kwargs.pop('overlay_chains',False)
    domains = kwargs.pop('domains',False)

    num_div = kwargs.pop('num_div',2)
    resnum_tick_labels = kwargs.pop('resnum_tick_labels',None)
    add_last_resi = kwargs.pop('add_last_resi',False)
    label_size = kwargs.pop('label_size',6)

    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from matplotlib.collections import LineCollection
    from matplotlib.pyplot import figure, imshow

    if y is None:
        raise TypeError('You need to provide data for the y-axis.')
    elif len(np.shape(y)) != 1:
        raise ValueError('The data must be a 1D array.')

    if SETTINGS['auto_show']:
        fig = plt.figure(figsize=[9.5,6])
    axes = [] 

    if atoms is not None:
        height_ratios = [15,0.2]
        nrows = 2
    else:
        height_ratios = None
        nrows = 1

    outer = GridSpec(1, 2, width_ratios = [16, 4], hspace=0., wspace=0.)

    gs = GridSpecFromSubplotSpec(nrows, 1, subplot_spec = outer[0], \
                                 height_ratios=height_ratios, hspace=0., wspace=0.)

    gs_legend = GridSpecFromSubplotSpec(1, 1, subplot_spec = outer[1], hspace=0., wspace=0.)
    
    ax1 = plt.subplot(gs[0])

    chain_colors = 'gcmyrwbk'
    chain_handles = []

    if overlay_chains:
        n = 0
        for i in atoms.getHierView().iterChains():
            chain_handle, = ax1.plot(y[i.getResindices()[0]:i.getResindices()[-1]], color=chain_colors[n], label=str(i), **kwargs)
            chain_handles.append(chain_handle)
            n += 1
    else:
        ax1.plot(y, **kwargs)

    if nrows > 1:
        ax2 = plt.subplot(gs[1])

        resnum_tick_locs = []
        resnum_tick_labels_list = []

        if resnum_tick_labels is None:
            resnum_tick_labels = []
            user_set_labels = False
        elif type(resnum_tick_labels) is list:
            user_set_labels = list
        elif type(resnum_tick_labels) is dict:
            user_set_labels = dict
        else:
            raise TypeError('The resnum tick labels should be a list or dictionary of lists')

        n = 0
        for i in atoms.getHierView().iterChains():
            if not overlay_chains:
                chain_handle, = ax2.plot([i.getResindices()[0], i.getResindices()[-1]], [0, 0], \
                                         '-', linewidth=3, color=chain_colors[n], label=str(i))
                chain_handles.append(chain_handle)

            if not user_set_labels:
                for j in range(num_div):
                    resnum_tick_locs.append(i.getResindices()[i.numAtoms()/num_div*j])
                    resnum_tick_labels.append(i.getResnums()[i.numAtoms()/num_div*j])
            elif user_set_labels is list:
                for j in resnum_tick_labels:
                    resnum_tick_locs.append(i.getResindices()[np.where(i.getResnums() == j)[0][0]])
                    resnum_tick_labels_list.append(j)
            else:
                for k in resnum_tick_labels.keys():
                    if i.getChids()[0] == k:
                       for j in resnum_tick_labels[k]: 
                           resnum_tick_locs.append(i.getResindices()[np.where(i.getResnums() == j)[0][0]])
                           resnum_tick_labels_list.append(j)

            n += 1

        if domains:
            try:
                atoms.getData('domain')[0]
            except:
                raise ValueError('A domain bar can only be generated if \
                                  there is domain data associated with \
                                  the atoms.')

            borders = {}
            for i in range(atoms.numAtoms()/atoms.getHierView().numChains()):
                if atoms.getData('domain')[i] != atoms.getData('domain')[i-1]:
                    if i != 0:
                        borders[atoms.getData('domain')[i-1]][-1].append(i-1)
                    if not atoms.getData('domain')[i] in borders.keys():
                        borders[atoms.getData('domain')[i]] = []
                    borders[atoms.getData('domain')[i]].append([])
                    borders[atoms.getData('domain')[i]][-1].append(i)

            hsv = plt.get_cmap('hsv')
            colors = hsv(np.linspace(0, 1.0, len(borders.keys())))

            for chain in atoms.getHierView().iterChains():
                domains_found = []
                for i in range(chain.numAtoms()):
                    if not atoms.getData('domain')[i] in domains_found and str(atoms.getData('domain')[i]) is not '':
                        n = 0
                        for j in borders[atoms.getData('domain')[i]]:
                            m = 0
                            if m == 0:
                                domain_handle, = ax2.plot([j[0], j[-1]], [0, 0], '-', linewidth=3, \
                                                          color=colors[n], label=str(atoms.getData('domain')[i]))
                                chain_handles.append(domain_handle)
                            else:
                                ax2.plot([j[0], j[-1]], [0, 0], '-', linewidth=3, color=colors[n])
                            m += 1
                        n += 1
 
        ax3 = plt.subplot(gs_legend[-1])
        plt.legend(handles=chain_handles, loc=2, bbox_to_anchor=(0.25, 1))
        ax3.axis('off')

        if not user_set_labels:
            resnum_tick_labels_list = resnum_tick_labels

        if add_last_resi:
            resnum_tick_locs.append(atoms.getResindices()[-1])
            resnum_tick_labels_list.append(atoms.getResnums()[-1])

        resnum_tick_locs = np.array(resnum_tick_locs)
        resnum_tick_labels = np.array(resnum_tick_labels_list)

        ax1.set_xticks([])

        if overlay_chains:
            ax1.set_xlim(-0.5,atoms.numAtoms()/atoms.getHierView().numChains()+0.5)

        ax2.set_xticks(resnum_tick_locs)
        ax2.set_xticklabels(resnum_tick_labels)
        ax2.tick_params(labelsize=label_size)
        ax2.set_yticks([])

        ax2.set_xlim(ax1.get_xlim())

    if atoms is not None:
        return ax1, ax2, ax3
    else:
        return ax1

def msaeye(msa, unique, turbo):
    tic1 = timeit.default_timer()
    length = msa.shape[1]
    number = msa.shape[0]
    # number = 5
    array = eye(int(number))

    seqs = []
    for i in range(number):
        seqs.append(msa[i,:])
    iseq = zeros((number, length), dtype=int)

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
                if direction is 'effect':
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
    :type data: `~numpy.ndarray`

    :arg nodes: a set of atoms that has been used
        as nodes in data generation
    :type nodes: :class:`

    :arg atoms: atoms to be selected from
    :type atoms: :class:`Atomic`

    """
    from collections import Counter

    try:
        data = asarray(data)
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
        indices = array([[i*3, i*3+1, i*3+2] 
                        for i in indices]
                        ).reshape(3*len(indices))

    data_ext = []
    resid_counter = Counter(atoms.getResindices())
    for i in indices:
        data_ext.extend(resid_counter.values()[i]*[data[i]])

    resid_selstr = ' '.join([str(resid) for resid in nodes.getResindices()])
    rest = atoms.select('not resid {0}'.format(resid_selstr))
    data_ext.extend(zeros(rest.numAtoms()))
        
    return data_ext
