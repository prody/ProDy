# -*- coding: utf-8 -*-

import time
import math
import os
import numpy as np
from prody import LOGGER
from prody.proteins import parsePDB
from prody.atomic import Atomic
from .anm import ANMBase, ANM
from .nma import NMA
from prody.utilities import checkCoords
#import _calcStiffnessMatrixElement

__all__ = ['calcStiffnessMatrix', 'calcPairDeformationDist','getVMDStiffness']


class SMBase(ANMBase, NMA):

    def __init__(self, name = 'Unknown'):
       
       super(SMBase, self).__init__(name)
       self.VMDname = None
       
    
    def getSMrange(self, stiffness_matrix):
 
        """ Function gives the lower and the higher value of Kij
        from stiffness matrix."""
 
        SMmin = np.amin(stiffness_matrix)
        SMmax = np.amax(stiffness_matrix)
        return SMmin, SMmax
        

    #@autojit
    def StiffnessMatrix(self, numCalphas, n_modes, coords, eigvals, eigvecs):

	"""This function calculates stiffness matrix as described in "Eran Eyal, Ivet Bahar, 
	Toward a Molecular Understanding of the Anisotropic Response of Proteins to External Forces: 
	Insights from Elastic Network Models, Biophysical Journal, 2008 (94) 3424-34355. 
	
	:numCalphas    : Number of C alpha atoms.
	:n_modes       : How many modes you want to include in your computations.
	:coords        : Equilibrium coordinates in a Nx3 matrix. Each element of the matrix coordinates of C alpha atom i: [x_i, y_i, z_i]
	:eigvals       : Eigenvalues 'including' zero modes.
	:eigvecs       : Eigenvetors 'including' zero modes.

	Autor: tekpinar@buffalo.edu"""

	k_B_x_T       = 1.0

	LOGGER.timeit('_sm')
	stiffness_matrix = np.zeros((numCalphas, numCalphas), np.double)
	inv_sqrt_eigvals = np.zeros((3*numCalphas), np.double)

	for k in xrange(6, n_modes):
	    inv_sqrt_eigvals[k] = math.sqrt(k_B_x_T/eigvals[k])

	eigvecs_flat = eigvecs.flatten()
	LOGGER.info('Calculating stiffness matrix.')

        from .smtools import calcStiffnessMatrixElement_calcStiffnessMatrixElement
	for i in xrange (0, numCalphas):
	    pairInd_I=i
	    ind_3i=pairInd_I*3
	    coord_I=coords[pairInd_I]

	    for j in xrange (i+1, numCalphas):

		pairInd_J=j
		coord_J=coords[pairInd_J]
		ind_3j=pairInd_J*3

		R_ij_sup_0=[(coord_J[0]-coord_I[0]), (coord_J[1]-coord_I[1]), (coord_J[2]-coord_I[2])]

		R_ij_sup_0_normalized_vec=R_ij_sup_0/np.linalg.norm(R_ij_sup_0)
		result=calcStiffnessMatrixElement_calcStiffnessMatrixElement(numCalphas, n_modes, ind_3i, ind_3j, R_ij_sup_0_normalized_vec, inv_sqrt_eigvals, eigvals, eigvecs_flat)
		stiffness_matrix[i][j]=result
		stiffness_matrix[j][i]=result

	LOGGER.report('Stiffness matrix calculated in %.2lfs.', label='_sm')

	return stiffness_matrix

    def MeanValueOfStiffness(self, stiffness_matrix):
        meanSiff = stiffness_matrix.sum(axis=0)
        return meanSiff
    
    #@autojit
    def calcPairDeformation(self, coords, selstr1, selstr2, coords_all, numCalphas, eigvals, eigvecs, n_modes):

        """This function calculate distribution of the deformations in distance dij for a pair of amino acids residue number
        and returns list with mode numbers and dij value for selected pair of amino acids.    
        """
        
        SM = SMBase()
        first_aa = SM.FirstLastAA(coords_all)
        selstr1 = selstr1 - int(first_aa[0]) + 1
        selstr2 = selstr2 - int(first_aa[0]) + 1
        LOGGER.timeit('_pairdef')
	inv_sqrt_eigvals = np.zeros((3*numCalphas), np.double)
	R_ij_0_matrix = np.zeros((numCalphas, numCalphas, 3), np.double)
	R_ij_sup_0_norm_matrix = np.zeros((numCalphas, numCalphas, 3), np.double)
	D_ij_k_matrix = np.zeros((numCalphas, numCalphas), np.double)

	k_B_x_T = 1.0
	    
	for k in xrange(6, n_modes):
	    inv_sqrt_eigvals[k] = math.sqrt(k_B_x_T/eigvals[k])

	counter = 0
	for i in xrange(0, numCalphas):
	    for j in xrange(counter+1, numCalphas):
		  R_ij_0_matrix[j][counter] = coords[j] - coords[counter]
		  R_ij_0_matrix[counter][j] = coords[j] - coords[counter]
		  R_ij_sup_0_norm_matrix[j][counter] = (coords[j] - coords[counter])/(np.linalg.norm(coords[j] - coords[counter]))
		  R_ij_sup_0_norm_matrix[counter][j] = (coords[j] - coords[counter])/(np.linalg.norm(coords[j] - coords[counter]))
	    counter = counter + 1

	D_pair_k = []
	mode_nr = []
	for m in xrange(6,n_modes):
	    U_ij_k = [(eigvecs[(selstr1-1)*3][m] - eigvecs[(selstr2-1)*3][m]), (eigvecs[((selstr1-1)*3)+1][m] - eigvecs[((selstr2-1)*3)+1][m]), (eigvecs[((selstr1-1)*3)+2][m] - eigvecs[((selstr2-1)*3)+2][m])]
	    D_ij_k = abs((inv_sqrt_eigvals[m])*(np.vdot(R_ij_sup_0_norm_matrix[selstr1-1][selstr2-1], U_ij_k))) 	
	    D_pair_k.append(D_ij_k)
	    mode_nr.append(m)
        
        LOGGER.report('Deformation was calculated in %.2lfs.', label='_pairdef')
	
	return mode_nr, D_pair_k
    

    def FirstLastAA(self, coords):

        """This function check the first and last residue number from parsed PDB.    
        """

	calphas = coords.select('protein and name CA')
	first = calphas.getResnums()[0]
	last = calphas.getResnums()[-1]
	return first, last


    def VMDdeformationPair(self, AAselstr, AArange, coords, coords_all, stiffness_matrix, selstr):
	
        """This function finding amino acids with specified range (AArange) of effective force constant
        for particular amino acids (AAselstr=[#aa1, #aa2] (from resid #aa1 to #aa2)).
        Its creates VMD file with found amino acids. This file is removed after closing VMD program. 
        
        :AAselstr = [selstr1, selstr2]	: range of amino acid structure using resid from PDB file;
                                            [selstr1] - use when only one resid will be given.
        :AArange = [value, value]	: range of Kij (effective force constant [N/m]);
        :coords				: coordinates of protein (parsePDB(pdf_file) or parse(fetchPDB(PDB_code));
        :coords_all			: this is parsed PDBfile, needed for writing PDB file with all atoms;
        :stiffness_matrix		: using StiffnessMatrix() """
        
        from prody import writePDB
        if len(AAselstr) == 1:
            VMDname = 'defpairVMD_'+str(AAselstr[0])+'aa_r'+str(AArange[0])+'-'+str(AArange[1])
        elif len(AAselstr) == 2:
            VMDname = 'defpairVMD_'+str(AAselstr[0])+'-'+str(AAselstr[1])+'aa_r'+str(AArange[0])+'-'+str(AArange[1])
        
	self.VMDname = VMDname

	writePDB(VMDname+'.pdb', coords_all)
	TextFile = open(VMDname+'.txt', 'w')
	VMDfile = open(VMDname+'.tcl','w')
	
	LOGGER.info('Creating VMD file.')
	
	VMDfile.write('display rendermode GLSL \n')
	VMDfile.write('display projection orthographic\n')
	VMDfile.write('color Display Background white\n')
	VMDfile.write('display shadows on\n')
	VMDfile.write('display depthcue off\n')
	VMDfile.write('axes location off\n')
	VMDfile.write('stage location off\n')
	VMDfile.write('light 0 on\n')
	VMDfile.write('light 1 on\n')
	VMDfile.write('light 2 off\n')
	VMDfile.write('light 3 on\n')
	VMDfile.write('mol addrep 0\n')
	VMDfile.write('display resetview\n')
	VMDfile.write('mol new {./'+str(self.VMDname)+'.pdb} type {pdb} first 0 last -1 step 1 waitfor 1\n')
	VMDfile.write('mol modselect 0 0 protein\n')
	VMDfile.write('mol modstyle 0 0 NewCartoon 0.300000 10.000000 4.100000 0\n')
	VMDfile.write('mol modcolor 0 0 Structure\n')
	VMDfile.write('mol color Structure\n')
	VMDfile.write('mol representation NewCartoon 0.300000 10.000000 4.100000 0\n')
	VMDfile.write('mol selection protein\n')
	VMDfile.write('mol material Opaque\n')

	SM = SMBase()
	first_aa = SM.FirstLastAA(coords_all)
	first_aa = int(first_aa[0])    
	
	colors = ['blue', 'red', 'gray', 'orange','yellow', 'tan','silver', 'green', 'white', 'pink', 'cyan', 'purple', 'lime', 'mauve', 'ochre', 'iceblue', 'black','yellow2','yellow3','green2','green3','cyan2','cyan3','blue2','blue3','violet','violet2','magenta','magenta2','red2','red3','orange2','orange3']*50
	
	if len(AAselstr) == 1:
	    AAselstr1 = int(AAselstr[0])
	    AAselstr2 = int(AAselstr[0])
	
	elif len(AAselstr) == 2:
	    AAselstr1 = int(AAselstr[0])
	    AAselstr2 = int(AAselstr[1])

        color_nr = 1 # starting from red color in VMD
        ResCounter = []
        for r in xrange(AAselstr1, AAselstr2+1):
	    baza_col = [] # Value of Kij is here for each residue
	    nr_baza_col = [] # Resid of aa are here
            VMDfile.write("draw color "+str(colors[color_nr])+"\n")
            
            for nr_i, i in enumerate(stiffness_matrix[r-first_aa]):
                if AArange[0] < float(i) < AArange[1]:
                    baza_col.append(i)
                    nr_baza_col.append(nr_i+first_aa)
                    coords_all_ca = coords_all.select(selstr)
                    resid_r = str(coords_all_ca.getResnames()[r-1])+str(r)
                    resid_r2 = str(coords_all_ca.getResnames()[nr_i])+str(nr_i+first_aa)
                    
                    if len(baza_col) == 0: # if base is empty then it will not change the color
                        color_nr = 0
                    else:
                        VMDfile.write("draw line "+'{'+str(coords[r-first_aa])[1:-1]+'} {'+str(coords[nr_i])[1:-1]+'} width 3 style solid \n')
                        TextFile.write(str(resid_r)+'\t'+resid_r2+'\t'+str(i)+'\n')
                        ResCounter.append(len(baza_col))
                        
                else: pass
	    if len(baza_col) != 0:
                VMDfile.write('mol addrep 0\n')
                VMDfile.write('mol modselect '+str(color_nr+1)+' 0 protein and name CA and resid '+ str(r)+' '+str(nr_baza_col)[1:-1].replace(',','')+'\n')
                VMDfile.write('mol modcolor '+str(color_nr+1)+' 0 ColorID '+str(color_nr)+'\n')
                VMDfile.write('mol modstyle '+str(color_nr+1)+' 0 VDW 0.600000 12.000000\n')
                VMDfile.write('mol color ColorID '+str(color_nr)+'\n')
                VMDfile.write('mol representation VDW 1.000000 12.000000 \n')
                VMDfile.write('mol selection protein and name CA and resid '+ str(r)+' '+str(nr_baza_col)[1:-1].replace(',','')+'\n')
                VMDfile.write('mol material Opaque \n')
                color_nr = color_nr + 1
                
        VMDfile.write('mol addrep 0\n')
        VMDfile.close()
        TextFile.close()
        
        if len(ResCounter) > 0:
            return VMDname
        elif len(ResCounter) == 0:
            LOGGER.info('There is no residue pair in this Kij range.')
            return 'None'   

     
# Definitions:

def calcStiffnessMatrix(coords,  n_modes='None', selstr='calpha', cutoff=13., gamma=1., saveMatrix='True', saveMap='False'):

    """This function calculates stiffness matrix as described in "Eran Eyal, Ivet Bahar, 
    Toward a Molecular Understanding of the Anisotropic Response of Proteins to External Forces: 
    Insights from Elastic Network Models, Biophysical Journal, 2008 (94) 3424-34355. 
    
    :coords		: coordinates of protein 
                         (coords = parsePDB('PDBfile')) for more see ``getCoords`` method.
    :n_modes		: how many modes you want to include in your computations 
                         If ``None`` the number of 3*N-3 modes will be calculated.
    :selstr		: atom selection (defult ``calpha``)
    :cutoff		: cutoff distance (Ang) for ANM calculations.
    :gamma		: spring constant.
    :saveMatrix		: If ``True``, it will save data matrix to ``stiffnessMatrix.txt file`` 
                         and a mean value for each residue to ``stiffnessMatrix_mean.txt``.
    :saveMap		: If ``True`` Striffness Matrix map will be save as PNG file. 
                         To this option matplotlib library is required.
    
    ---> Usage example: calcStiffnessMatrix(coords, saveMap='True')
    
    """
    
    try:
        coords = coords.select(selstr) 
        coords_all = coords.select(selstr)
        numCalphas = coords.numAtoms() 
        coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                    coords.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be a Numpy array or an object '
                            'with `getCoords` method'
                            'Use: coords = parsePDB(PDBfile) or coords = parsePDB(fetchPDB("PDBcode", compressed=False))')

    if(n_modes=='None'):
        n_modes = int(3*numCalphas)
    else: 
        n_modes = int(n_modes) 
    
    SM  	  = SMBase()
    anm           = ANM('prot analysis')                        
    anm.buildHessian(coords, cutoff, gamma)
    anm.calcModes(n_modes, zeros=True)
    
    eigvals	  = anm.getEigvals().transpose()
    eigvecs	  = anm.getEigvecs().transpose()                            

    sm = SM.StiffnessMatrix(numCalphas, n_modes, coords, eigvals, eigvecs)
    meanV = SM.MeanValueOfStiffness(sm)

    if(saveMatrix == 'True'):
        np.savetxt('StiffnessMatrix.txt', sm, fmt='%.6f')
        
        MValfile = open('StiffnessMatrix_mean.txt','w')
        for nr_i, i in enumerate(meanV):
            MValfile.write(str(nr_i)+'\t'+str(i)+'\n')
        MValfile.close()
        
    if(saveMap == 'True'):
	SMrange = SM.getSMrange(sm)
	LOGGER.info('Maximum value of effective force constant is {0}.'
				.format(SMrange[1]))

        import matplotlib
        import matplotlib.pylab as plt
	
	matplotlib.rcParams['font.size'] = '16'
	fig = plt.figure(num=None, figsize=(10,8), dpi=100, facecolor='w')
	plt.imshow(sm, origin='lower')
	plt.axis([0, numCalphas, numCalphas, 0])
	plt.xlabel('residue', fontsize = '18')
	plt.ylabel('residue', fontsize = '18')
	cbar = plt.colorbar()
	cbar.set_label('N/m', rotation=90, fontsize = '18')
        plt.savefig('StiffnessMatrixMap.png', dpi=100)
        plt.show()    


def calcPairDeformationDist(coords, selstr1, selstr2, n_modes='None', selstr='calpha', cutoff=13., gamma=1., saveDataDefPair='True', saveDefPlot='False'):

    """This function calculates distribution of the deformations in distance (dij) 
    for a pair of amino acids. 
    It was described in: "Eran Eyal, Ivet Bahar, 
    Toward a Molecular Understanding of the Anisotropic Response of Proteins to External Forces: 
    Insights from Elastic Network Models, Biophysical Journal, 2008 (94) 3424-34355. 
    Eq. (10) and Fig. (2). 

    :coords		: coordinates of protein 
                          (coords = parsePDB('PDBfile')) for more see ``getCoords`` method.
    :selstr1		: number of first amino acid, number from PDB file (ResID).
    :selstr2		: number of second amino acid.
    :n_modes		: how many modes you want to include in your computations 
                          If ``None`` the number of 3*N-3 modes will be calculated.
    :selstr		: atom selection (defult ``calpha``).
    :cutoff		: cutoff distance (Ang) for ANM calculations.
    :gamma		: spring constant.   
    :saveDataDefPair	: If ``True``, it will save data to ``defpair_(selstr1)_(selstr2).dat file``.
    :saveDefPlot	: If ``True`` Plot dij(k) vs mode will be save as ``defpair_(selstr1)_(selstr2).png`` file.  
                          To this option matplotlib library is required.
    
    ---> Usage example: calcPairDeformationDist(coords, 3, 212, saveDataDefPair='False', saveDefPlot='True')
    
    """

    try:
        coords = coords.select(selstr)
        coords_all = coords.select(selstr)
        numCalphas = coords.numAtoms()
        coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                    coords.getCoords())

    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be a Numpy array or an object '
                            'with `getCoords` method'
                            'Use: coords = parsePDB(PDBfile) or coords = parsePDB(fetchPDB("PDBcode", compressed=False))')

    try:
        selstr1 = int(selstr1)
        selstr2 = int(selstr2)
    except TypeError:
        raise TypeError('selstr should be a residue number')

    if(n_modes == 'None'):
        n_modes = int(3*numCalphas)
    else: 
        n_modes = int(n_modes) 
    
    SM		  = SMBase()
    anm           = ANM('prot analysis')                        
    anm.buildHessian(coords, cutoff, gamma)
    anm.calcModes(n_modes, zeros=True)
    
    eigvals	  = anm.getEigvals()
    eigvecs	  = anm.getEigvecs()                            

    pairdef = SM.calcPairDeformation(coords, selstr1, selstr2, coords_all, numCalphas, eigvals, eigvecs, n_modes)

    if(saveDataDefPair == 'True'):
        DataFileName = "defpair_"+str(selstr1)+"_"+str(selstr2)+".dat"
        DataFile = open(DataFileName, 'w')
        for i in xrange(len(pairdef[0])):
            DataFile.write("{} {}\n".format(pairdef[0][i], pairdef[1][i]))
        DataFile.close()
        LOGGER.info('Data file has been saved.')
            
    if(saveDefPlot == 'True'):
        import matplotlib
        import matplotlib.pylab as plt
        
        matplotlib.rcParams['font.size'] = '16'
        fig = plt.figure(num=None, figsize=(12,8), dpi=100, facecolor='w')
        plt.plot(pairdef[0], pairdef[1])
        plt.xlabel('mode (k)', fontsize = '18')
        plt.ylabel('d$^k$' '($\AA$)', fontsize = '18')    
        plt.savefig('defpair_'+str(selstr1)+'_'+str(selstr2)+'.png', dpi=100)
        LOGGER.info('Plot has been saved.')
        plt.show()


def getVMDStiffness(coords, AAselstr, AArange, n_modes='None', selstr='calpha', cutoff=13., gamma=1., saveTCLfile='True', saveTextFile='True', loadToVMD='False'):
    
    """This function finding amino acids with specified AArange=[value_low, value_high] of effective force constant
    for particular amino acids AAselstr=[#aa1, #aa2] (from resid #aa1 to #aa2).
    Its creates TCL file for VMD program with founded amino acids and a text file.
    Effective force constant calculation based on: "Eran Eyal, Ivet Bahar, 
    Toward a Molecular Understanding of the Anisotropic Response of Proteins to External Forces: 
    Insights from Elastic Network Models, Biophysical Journal, 2008 (94) 3424-34355.  
     
    :coords		: coordinates of protein (parsePDB(pdf_file) or parse(fetchPDB(PDB_code)).
    :AAselstr		: range of amino acid structure using resid from PDB file ``[selstr1, selstr2]`` 
                            or ``[selstr1]`` for only one amino acid. 
    :AArange		: range of Kij - effective force constant [N/m], ``[min_Kvalue, max_Kvalue]``.
    :selstr		: atom selection (defult ``calpha``).
    :cutoff		: cutoff distance (Ang) for ANM calculations.
    :gamma		: spring constant.    
    :saveTCLfile	: If ``True`` TCL file will be save.
    :saveTextFile	: If ``True`` text file with pair of residues will be created. 
    :loadToVMD		: If ``True`` PDB file with coords and created TCL file will be loaded to VMD program.
                            VMD program is required.
    
    --> Usage example: getVMDStiffness(coords, [1,10], [0,8], loadToVMD='True'), 
        searching aa foe each amino acid from range 1-10 
        if they have Kij between 0-8 it will add it to TCL file.
        This function is prepared for single molecules. The numeration may change in a complex systems.
        It will be changed in the future.
    
    """    
    
    try:
        coords_all = coords
        coords = coords.select(selstr)
        numCalphas = coords.numAtoms()
        coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                    coords.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be a Numpy array or an object '
                            'with `getCoords` method'
                            'Use: coords = parsePDB(PDBfile) or coords = parsePDB(fetchPDB("PDBcode", compressed=False))')

    if(n_modes=='None'):
        n_modes = int(3*numCalphas)
    else: 
        n_modes = int(n_modes) 
    
    SM		  = SMBase()
    anm           = ANM('prot analysis')                        
    anm.buildHessian(coords, cutoff, gamma)
    anm.calcModes(n_modes, zeros=True)
    
    eigvals	  = anm.getEigvals().transpose()
    eigvecs	  = anm.getEigvecs().transpose()                                
    
    sm = SM.StiffnessMatrix(numCalphas, n_modes, coords, eigvals, eigvecs)
    vmdname = SM.VMDdeformationPair(AAselstr, AArange, coords, coords_all, sm, selstr)

    if vmdname != 'None':    
        if(saveTCLfile == 'False'):
            os.system('rm '+str(vmdname)+'.pdb '+str(vmdname)+'.tcl')
        else: pass

        if(saveTextFile == 'False'):
            os.system('rm '+str(vmdname)+'.txt')
            LOGGER.info('Text file has been saved.')
        
        if(loadToVMD == 'True'):
            from prody import pathVMD
            os.system(pathVMD()+" -e "+str(vmdname)+".tcl")
        else: pass
    else: pass

