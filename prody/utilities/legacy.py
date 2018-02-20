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

    :arg domain_bar: color the bar at the bottom by domains rather than chains
        default False
    :type domain_bar: bool
    """
    atoms = kwargs.pop('atoms',None)
    overlay_chains = kwargs.pop('overlay_chains',False)
    domain_bar = kwargs.pop('domain_bar',False)

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

        if domain_bar:
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