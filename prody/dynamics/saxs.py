from prody import Atom
import sys
import numpy as np
from numba import jit
from .saxstools import calcSAXSNumeric
#import saxstools
from math import sqrt
from prody.utilities import saxsWater

__all__ = ['buildSolvShell', 'showSaxsProfiles', 
           'calcSaxsChi', 'parseSaxsData', 
           'calcSaxsPerModel','interpolateMode',
           'showChivsFrames', 'writeChivsFrames',
           'writeSaxsProfile']

WATER_BOX_SIZE = 119.7
NUM_WATER_ATOMS = 60656
water = saxsWater()
@jit
def min_dist_to_mol(X, Y, Z, x, y, z, atom_count, mol_type):
    mtype = 0
    dist = 0.0
    dist_x = 0.0
    dist_y = 0.0
    dist_z = 0.0
    min_dist = 100000000
    for i in range (0, atom_count):
        dist_x = (X[i]-x)
        dist_y = (Y[i]-y)
        dist_z = (Z[i]-z)
        #    dist=((dist_x*dist_x)+(dist_y*dist_y)+(dist_z*dist_z))
        dist = ((dist_x**2)+(dist_y**2)+(dist_z**2))
        if (dist < 9.0):
            min_dist = 0
            break

        if (dist < min_dist):
            min_dist = dist
            mtype = mol_type[i]
      
    return [mtype, min_dist]

@jit
def min_dist_to_mol_v4(X, Y, Z, x, y, z, atom_count, mol_type):
    mtype = 0
    dist = 0.0
    dist_x = 0.0
    dist_y = 0.0
    dist_z = 0.0
    min_dist = 100000000
    thickness = 3.0
    thicknessSqrd = 9.0
    for i in range (0, atom_count):
        dist_x = (X[i]-x)
        dist_y = (Y[i]-y)
        dist_z = (Z[i]-z)

        if (dist_x < thickness):
            if (dist_y < thickness):
                if (dist_z < thickness):
                    dist = ((dist_x*dist_x)+(dist_y*dist_y)+(dist_z*dist_z))    
                    if (dist < thicknessSqrd):
                        min_dist = 0
                        break

        dist = ((dist_x*dist_x)+(dist_y*dist_y)+(dist_z*dist_z))    
        if (dist < min_dist):
            min_dist = dist
            mtype = mol_type[i]
      
    return [mtype, min_dist]

@jit
def min_dist_to_mol_v5(X, Y, Z, x, y, z, atom_count, mol_type):
    mtype = 0
    dist = 0.0
    dist_x = 0.0
    dist_y = 0.0
    dist_z = 0.0
    min_dist = 100000000
    thickness = 3.0
    thicknessSqrd = 9.0
    b = np.array(x, y, z)  
    for i in range (0, atom_count):
        a = np.array([X[i], Y[i], Z[i]])

        dist = np.sum( (a - b)**2)
        if (dist < thicknessSqrd):
            min_dist = 0
            break

        if (dist < min_dist):
            min_dist = dist
            mtype = mol_type[i]
      
    return [mtype, min_dist]

def buildSolvShell(protein, W, X, Y, Z, closest_dist, thickness, wDNA, wRNA, wPROT, mol_type, MAX_ATOM):
    """
    This function build a 3 Angstrom water shell around a given biological 
    macromolecule. The solvation shell improves agreement of theoretical and 
    experimental Small Angle Scattering profiles. When all computations are 
    completed, a solvated structure file is written in PDB format.
    """
    total_atom = 0
    pdb_out_file = "cg_solvated.pdb"
    x = 0.0
    y = 0.0
    z = 0.0
    maxx = -10000
    maxy = -10000
    maxz = -10000
    minx = 10000
    miny = 10000
    minz = 10000

    mol_type_flag = 0
    flag = 0

    atom = [None]*MAX_ATOM
    cgatom = [None]*MAX_ATOM

    atom_count = 0
    RNA_count = 0
    DNA_count = 0
    PROT_count = 0
    closest_dist_plus_thickness = (closest_dist+thickness)
    closest_dist_plus_thickness_sqrd = closest_dist_plus_thickness*closest_dist_plus_thickness
    closest_dist_sqrd = closest_dist*closest_dist

    mtype = 0

    #mol_type_flag:  0-DNA, 1-RNA, 2-Protein
    for item in protein.iterAtoms():
        atomname = item.getName()
        resname = item.getResname()
        atomCoords = item.getCoords()

        if( (atomname == " O2'") or (atomname == " O2*") or (resname == "RNA") ): 
            mol_type_flag = 1
            if((RNA_count == 0) and (resname != "RNA") and (DNA_count == 1)):
                cgatom[atom_count-1] = "RNA"
                mol_type[atom_count-1] = 1
                W[atom_count-1] = wRNA
                RNA_count+=1
                DNA_count-=1

        if((atomname==" O5'") or (atomname==" O5*")): 
            if(mol_type_flag==1):
                cgatom[atom_count] = "RNA"
                RNA_count+=1

            else:
                cgatom[atom_count]="DNA";
                DNA_count+=1
                mol_type_flag = 0
            atom[atom_count]=" O5'";
            flag = 1

        if(((resname=="  A") or (resname==" DA") or
            (resname=="DA ") or (resname==" RA") or
            (resname=="RA ") or (resname=="ADE") )
           and
           ((atomname==" C5 ") or (atomname=="  C5"))):
            cgatom[atom_count] = "ADE"
            flag = 1
            atom[atom_count] = " C5 "

        if(((resname=="  C") or (resname==" DC") or
            (resname=="DC ") or (resname==" RC") or
            (resname=="RC ") or (resname=="CYT") ) 
           and
           ((atomname==" N3 ") or (atomname=="  N3"))):
            cgatom[atom_count] = "CYT"
            flag = 1
            atom[atom_count] = " N3 "

        if(((resname=="  G") or (resname==" DG") or
            (resname=="DG ") or (resname==" RG") or
            (resname=="RG ") or (resname=="GUA") ) 
           and
           ((atomname==" C4 ") or (atomname=="  C4"))):
            cgatom[atom_count] = "GUA"
            flag = 1
            atom[atom_count] = " C4 "

        if(((resname=="  U") or (resname==" DU") or
            (resname=="DU ") or (resname==" RU") or
            (resname=="RU ") or (resname=="URA") )
           and
           ((atomname==" N3 ") or (atomname=="  N3"))): 
            cgatom[atom_count] = "URA"
            flag = 1
            atom[atom_count] = " N3 "

        if(((resname=="  T") or (resname==" DT") or
            (resname=="DT ") or (resname==" RT") or
            (resname=="RT ") or (resname=="THY") ) 
           and
           ((atomname==" N3 ") or (atomname=="  N3"))):
            cgatom[atom_count] = "THY"
            flag = 1
            atom[atom_count] = " N3 "

        #Check if atom belongs to protein
        if((atomname==" CA ") or (atomname=="  CA") or
           (atomname=="CA  ") or (atomname=="CA")):

            cgatom[atom_count] = "   "
            if((resname=="GLY")): cgatom[atom_count] = "GLY"
            if((resname=="ALA")): cgatom[atom_count] = "ALA"
            if((resname=="VAL")): cgatom[atom_count] = "VAL"
            if((resname=="LEU")): cgatom[atom_count] = "LEU"
            if((resname=="ILE")): cgatom[atom_count] = "ILE"
            if((resname=="MET")): cgatom[atom_count] = "MET"
            if((resname=="PHE")): cgatom[atom_count] = "PHE"
            if((resname=="TRP")): cgatom[atom_count] = "TRP"
            if((resname=="PRO")): cgatom[atom_count] = "PRO"
            if((resname=="SER")): cgatom[atom_count] = "SER"
            if((resname=="THR")): cgatom[atom_count] = "THR"
            if((resname=="CYS")): cgatom[atom_count] = "CYS"
            if((resname=="TYR")): cgatom[atom_count] = "TYR"
            if((resname=="ASN")): cgatom[atom_count] = "ASN"
            if((resname=="GLN")): cgatom[atom_count] = "GLN"
            if((resname=="ASP")): cgatom[atom_count] = "ASP"
            if((resname=="GLU")): cgatom[atom_count] = "GLU"
            if((resname=="LYS")): cgatom[atom_count] = "LYS"
            if((resname=="ARG")): cgatom[atom_count] = "ARG"
            if((resname=="HIS")): cgatom[atom_count] = "HIS"
            if((resname=="HSD")): cgatom[atom_count] = "HIS"
            if((resname=="HYP")): cgatom[atom_count] = "PRO"
            PROT_count+=1
            mol_type_flag = 2
            flag = 1
            atom[atom_count] = " CA "

        if(flag==1):
            x = atomCoords[0]	  
            y = atomCoords[1]
            z = atomCoords[2]

            if(x > maxx): maxx = x
            if(y > maxy): maxy = y
            if(z > maxz): maxz = z

            if(x < minx): minx = x
            if(y < miny): miny = y
            if(z < minz): minz = z

            X[atom_count] = x
            Y[atom_count] = y
            Z[atom_count] = z

            mol_type[atom_count]=mol_type_flag
            if (mol_type_flag==0): W[atom_count] = wDNA
            if (mol_type_flag==1): W[atom_count] = wRNA
            if (mol_type_flag==2): W[atom_count] = wPROT

            atom_count+=1
            flag = 0

    pdb_flag = 0
    if(pdb_flag==1):
        fout_pdb = open(pdb_out_file, "w");
        for i in range (0, atom_count):
            print >>fout_pdb, "ATOM %6d %4s %3s %1d%4d    %8.3f%8.3f%8.3f" \
                %((i+1), atom[i], cgatom[i], mol_type[i],(i+1), X[i], Y[i], Z[i]) 


    water_count=1
    solvent_flag=1
    x_range = int((maxx-minx)/WATER_BOX_SIZE)+1
    y_range = int((maxy-miny)/WATER_BOX_SIZE)+1
    z_range = int((maxz-minz)/WATER_BOX_SIZE)+1
    if(solvent_flag == 1):
        for i in range(0, x_range):
            for j in range(0, y_range):
                for k in range(0, z_range):
                    for w in range(0, NUM_WATER_ATOMS):
                        x = water[w][0]+minx + (i*WATER_BOX_SIZE)
                        y = water[w][1]+miny + (j*WATER_BOX_SIZE)
                        z = water[w][2]+minz + (k*WATER_BOX_SIZE)
                        [mtype, min_dist] = min_dist_to_mol(X,Y,Z,x,y,z,atom_count,mol_type)

                        if (mtype==0): W[atom_count+water_count-1] = wDNA
                        if (mtype==1): W[atom_count+water_count-1] = wRNA
                        if (mtype==2): W[atom_count+water_count-1] = wPROT
                        cgatom[atom_count+water_count-1] = "SOL"
                        X[atom_count+water_count-1] = x
                        Y[atom_count+water_count-1] = y
                        Z[atom_count+water_count-1] = z

                        if(min_dist < closest_dist_sqrd): continue
                        if(min_dist > closest_dist_plus_thickness_sqrd): continue

                        if(pdb_flag ==1):
                            print >>fout_pdb,"ATOM %6d  OW  %3s %1d%4d    %8.3f%8.3f%8.3f" \
                            %(water_count, cgatom[atom_count+water_count-1], mtype, water_count, x, y, z)

                        water_count+=1

                        if(water_count==10000):
                            water_count = 1
    if(pdb_flag==1):
        fout_pdb.close()
    total_atom = atom_count + water_count - 1;

    return total_atom



#def get_f(resname, w, q):
#    f=0.0001
#
#    if(resname=="RNA"):
#        f = RS(q)
#    if(resname=="DNA"):
#        f = DS(q)
#    if(resname=="ADE"):
#        f = ADE(q)
#    if(resname=="GUA"):
#        f = GUA(q)
#    if(resname=="THY"):
#        f = THY(q)
#    if(resname=="URA"):
#        f = URA(q)
#    if(resname=="CYT"):
#        f = CYT(q)
#    
#    if(resname=="GLY"):
#        f = GLY(q)
#    if(resname=="ALA"):
#        f = ALA(q)
#    if(resname=="VAL"):
#        f = VAL(q)
#    if(resname=="ILE"):
#        f = ILE(q)
#    if(resname=="LEU"):
#        f = LEU(q)
#    if(resname=="MET"):
#        f = MET(q)
#    if(resname=="PHE"):
#        f = PHE(q)
#    if(resname=="TRP"):
#        f = TRP(q)
#    if(resname=="PRO"):
#        f = PRO(q)
#    if(resname=="SER"):
#        f = SER(q)
#    if(resname=="THR"):
#        f = THR(q)
#    if(resname=="CYS"):
#        f = CYS(q)
#    if(resname=="TYR"):
#        f = TYR(q)
#    if(resname=="ASN"):
#        f = ASN(q)
#    if(resname=="GLN"):
#        f = GLN(q)
#    if(resname=="ASP"):
#        f = ASP(q)
#    if(resname=="GLU"):
#        f = GLU(q)
#    if(resname=="LYS"):
#        f = LYS(q)
#    if(resname=="ARG"):
#        f = ARG(q)
#    if(resname=="HIS"):
#        f = HIS(q)
#    if(resname=="HSD"):
#        f = HIS(q)
#    
#    if(resname=="SOL"):
#        f = w*SOL(q)
#    if(resname=="TIP"):
#        f = w*SOL(q)
#    if(resname=="SPC"):
#        f = w*SOL(q)
#    #printf("%s %f %f %f\n", resname,w,q,f);
#    return f;
def numericResname(resname):
    f = 0

    if(resname=="RNA"):
        f = 1
    if(resname=="DNA"):
        f = 2
    if(resname=="ADE"):
        f = 3
    if(resname=="GUA"):
        f = 4
    if(resname=="THY"):
        f = 5
    if(resname=="URA"):
        f = 6
    if(resname=="CYT"):
        f = 7
    if(resname=="GLY"):
        f = 8
    if(resname=="ALA"):
        f = 9
    if(resname=="VAL"):
        f = 10
    if(resname=="ILE"):
        f = 11
    if(resname=="LEU"):
        f = 12
    if(resname=="MET"):
        f = 13
    if(resname=="PHE"):
        f = 14
    if(resname=="TRP"):
        f = 15
    if(resname=="PRO"):
        f = 16
    if(resname=="SER"):
        f = 17
    if(resname=="THR"):
        f = 18
    if(resname=="CYS"):
        f = 19
    if(resname=="TYR"):
        f = 20
    if(resname=="ASN"):
        f = 21
    if(resname=="GLN"):
        f = 22
    if(resname=="ASP"):
        f = 23
    if(resname=="GLU"):
        f = 24
    if(resname=="LYS"):
        f = 25
    if(resname=="ARG"):
        f = 26
    if(resname=="HIS"):
        f = 27
    if(resname=="HSD"):
        f = 28
    
    if(resname=="SOL"):
        f = 29
    if(resname=="TIP"):
        f = 30
    if(resname=="SPC"):
        f = 31
    #printf("%s %f %f %f\n", resname,w,q,f);
    return f;

def showSaxsProfiles(exp_data_file, model_data_file):
    """This function reads experimental and model data from txt files and
    plots them by using matplotlib. 
    """
    #Read experimental data
    #Read model data
    #Plot data
    from matplotlib import pyplot
    from pylab import genfromtxt
    mat0 = genfromtxt(exp_data_file)
    mat1 = genfromtxt(model_data_file)

    I_q_exp = mat0[:,1]
    sigma_q = mat0[:,2]
    I_q_model = mat1[:,1]

    part1 = np.sum((I_q_exp*I_q_model)/(sigma_q*sigma_q))
    part2 = np.sum((I_q_model*I_q_model)/(sigma_q*sigma_q))
    c = part1/part2
 
    offset = mat0[:,1][0] - mat1[:,1][0]
#    print c
    pyplot.plot(mat0[:,0], mat0[:,1], label = "Experimental")
#    pyplot.plot(mat1[:,0], c*mat1[:,1], label = "Model")
    pyplot.plot(mat1[:,0], mat1[:,1]+offset, label = "Model")
    pyplot.legend()
    pyplot.show()


def calcSaxsChi(q_exp, I_q_exp, sigma_q, q_model, I_q_model):
    """
    This function calculates Chi squared value between SAXS profile of an 
    experimental( or simulated) data and theoretical SAXS profile of a model.
    The model data can be obtained from an experimental macromolecular 
    structure, a molecular dynamics conformation or a normal mode.
    The formula is based on Equation 15 and Equation 16 of 
    Svergun, D., Barberato, C. & Koch, M. H. J. (1995). 
    J. Appl. Cryst. 28, 768-773.
    """
    #Read model data
#    q_model, I_q_model = np.loadtxt(model_I_q_file, unpack=True)

    #Check length of model data and compare it with experimental data set. 
    #Chech if model data was taken at experimental q values.
    #If not, interpolate

    #Ensure that q_exp[0] and q_model[0] points to the same value. 
    #Compare head of data sets to find c value.
    if(q_exp[0]==q_model[0]):
        offset = I_q_exp[0] - I_q_model[0]
#        part1=np.sum((I_q_exp*I_q_model)/(sigma_q*sigma_q))
#        part2=np.sum((I_q_model*I_q_model)/(sigma_q*sigma_q))
#        c=part1/part2

    else:
        print("Experimental and theoretical Q[0] values do not match!")
        sys.exit(-1)

    #Calculate Chi    
#    diff_array=((I_q_exp-c*I_q_model)/sigma_q)
    diff_array = ((I_q_exp - I_q_model - offset)/sigma_q)
    chi_sqrd = (np.sum(np.square(diff_array)))/len(q_exp)
    return sqrt(chi_sqrd)


#def parseSaxsData(I_q_file, simulated=False, isLogScale=True):
def parseSaxsData(I_q_file, **kwargs):
    """
    This function parses Small Angle X-ray Scattering data.
    ARGS:
    I_q_file:   This file is a plain text file with at least two
                columns, namely, q column and I_q column. If the 
                file is an experimental data file, there should be a third
                column of sigma_q values, namely, error bars. 
    simulated:  This parameter is a boolean parameter and its default 
                value is False. If the data
                is simulated SAXS data, you should make it True so 
                that it will assign 1.0 values to all sigma_q values.
              
    isLogScale: This parameter is also a boolean parameter and its 
                default value is False. prody_saxs module requires 
                the intensity (I_q) data to be in 
                logscale for all calculations. Therefore, if the 
                intensities are not in logscale, you should make 
                isLogScale=False, so that the function return the 
                input data in log10 scale. 
    """
    if kwargs is not None:
        simulated = kwargs.get('simulated', False)
        isLogScale = kwargs.get('isLogScale', True)

    #Check if data contains error bars or if it is simulated data!
    if(simulated==False):
        #Read experimental data
        q_exp, I_q_exp, sigma_q = np.loadtxt(I_q_file, unpack=True)
    elif(simulated==True):
        #Read simulated data
        q_exp, I_q_exp = np.loadtxt(I_q_file, unpack=True)

        #Since it is simulated data, assign 1.0 values to errorbars.
        sigma_q = np.ones(len(q_exp))
            
    #Check if data is in log scale or not.
    if(isLogScale==False):
        return (q_exp, np.log10(I_q_exp), np.log10(sigma_q))
    else:
        return (q_exp, I_q_exp, sigma_q)

#def calcSaxsPerModel(calphas, numCalphas, I, Q_exp):
def calcSaxsPerModel(calphas, I_model, Q_exp):
    ####Assign coordinates to a temporary array!#################################
    numCalphas=calphas.numAtoms()
    wDNA = 0.070
    wRNA = 0.126
    wPROT = 0.040
    thickness = 3.0
    closest_dist = 3.5
    delta = 0.0
#    total_atom=0
    pdb_flag = 1
    solvent_flag = 1
    total_atom = numCalphas
    MAX_ATOM = 100000
    W = np.zeros(MAX_ATOM)
    X = np.zeros(MAX_ATOM)
    Y = np.zeros(MAX_ATOM)
    Z = np.zeros(MAX_ATOM)
    mol_type = np.zeros(MAX_ATOM, dtype=int)
    cgatom_num = np.zeros(MAX_ATOM, dtype='int32')
    cgatom = calphas.getResnames()
    origCoords = calphas.getCoords()
    
    X[:len(origCoords[:,0])] = origCoords[:,0]
    Y[:len(origCoords[:,1])] = origCoords[:,1]
    Z[:len(origCoords[:,2])] = origCoords[:,2]

    #At first, solvate the protein:
    #print "@> Solvating the system"
    #start = timeit.timeit()
#    start = time.clock()
    total_atom = buildSolvShell(calphas, W, X, Y, Z, closest_dist, thickness, \
                              wDNA, wRNA, wPROT, mol_type, MAX_ATOM)
#    end= time.clock()
#    print "Solvation time: %f"%(end - start)

    #end = timeit.timeit()
    #print (end - start)
#    total_atom=buildSolvShell(calphas, W, X, Y, Z, closest_dist, thickness, \
#                              wDNA, wRNA, wPROT, mol_type, MAX_ATOM)
    #C version of this building solvation shell. I have to check which version works well. 
#    writePDB('model_no_water.pdb', calphas)
#    total_atom=saxstools.cgSolvateNumeric("model_no_water.pdb", X, Y, Z, W, \
    #                                      wDNA, wRNA, wPROT, \
    #                                      thickness, closest_dist, \
    #                                      pdb_flag, solvent_flag, MAX_ATOM)

    #print "@> Total number of atoms after solvation: %d" %total_atom
    
    for j in range(0, total_atom):
        if(j<numCalphas):
            cgatom_num[j] = numericResname(cgatom[j])
        else:
            cgatom_num[j] = numericResname("SOL")

    #Now, let reread solvated structure. And calculate SAXS profile of it with 
    #solvated_protein = parsePDB("cg_solvated.pdb")

    #You need a function to check experimental SAXS data. It should give Q_exp
    #and it should check whether experimental data is in log scale or not!
#    _saxs.calcSAXSNumeric(I, X, Y, Z, total_atom, cgatom_num, W, Q_exp, len(Q_exp))
    calcSAXSNumeric(I_model, X, Y, Z, total_atom, cgatom_num, W, Q_exp, len(Q_exp))
    ####Finish key part!#########################################################


def writeSaxsProfile(I_q, q_exp, sigma_q, filename):
    """ Write a theoretical or modified experimental 
    SAXS profile to a file.
    """
    
    combined = np.vstack([q_exp, I_q, sigma_q]).T
    np.savetxt(filename, combined)

    
def interpolateMode(calphas, mode, Q_exp, I_q_exp, sigma_q, max_chi,**kwargs):
  
    if kwargs is not None:
        args_numFrames = kwargs.get('numFrames', 20)
        args_scalCoeff = kwargs.get('scalCoeff', 3.0)
        args_out_pdb_file = kwargs.get('out_pdb_file', 'best_model.pdb')
#        args_out_saxs_file=kwargs.get('out_saxs_file', None)
    
    chi_list = []
    frames_list = []

    eigenvalue = mode.getEigval()
    eigenvector = np.transpose(mode.getEigvec())

    i = mode.getIndex()
    mod_num = None
    origCoords = calphas.getCoords()
    numCalphas = calphas.numAtoms()
    I_model = np.zeros(len(Q_exp))
    # setup toolbar
    sys.stdout.write("@> Calculating SAXS profiles for nonzero mode %d: " % (i+1))
    sys.stdout.write("[%s]" % (" " * (numFrames+1)))
    sys.stdout.flush()
    sys.stdout.write("\b" * (numFrames+2)) # return to start of line, after '['
    prody.LOGGER.timeit('_intplt_mode')
    invEigVal = (1.0/eigenvalue)
    for j in range((-numFrames/2), ((numFrames/2)+1)):
        coeff = j*rmsdScalingCoef*invEigVal*2.0/numFrames
        
        newCoords = calphas.getCoords().flatten()+(coeff*eigenvector)
        calphas.setCoords(newCoords.reshape((numCalphas, 3), order='C'))
        calcSaxsPerModel(calphas, numCalphas, I_model, Q_exp)
        chi = calcSaxsChi(Q_exp, I_q_exp, sigma_q, Q_exp, I_model)
        chi_list.append(chi)
        frames_list.append(j)
        if(chi < max_chi):
            max_chi = chi
            mod_num = i
            writePDB(out_pdb_file, calphas)
#           extendModel(calphas, 'calphas', protein)
#           writePDB('best_model.pdb', protein)

            #Reset coordinates to the original values
        calphas.setCoords(origCoords)
            
        sys.stdout.write('#')
        sys.stdout.flush()
    sys.stdout.write("\n")
    prody.LOGGER.report('SAXS profile calculations were performed in %2fs.', '_intplt_mode')
    
    return chi_list, frames_list


def showChivsFrames(chi_list, frames_list, numFrames=20):
    """Show Chi values vs Frames number to see if interpolating a mode reduces
    Chi value or not. 
    """

    from matplotlib import pyplot
    pyplot.xticks(fontsize = 16)
    pyplot.yticks(fontsize = 16)
    pyplot.grid()
    
    linestyles = ['-', '-.','--', ':']
    numModes=len(chi_list)/(numFrames+1)
    print("@> Number of modes in chi list is %d"%numModes)
    for i in range (0, numModes):
        pyplot.plot(frames_list[(i*(numFrames+1)):((i+1)*(numFrames+1))], \
                    chi_list[(i*(numFrames+1)):((i+1)*(numFrames+1))], \
                    linestyle = linestyles[i%4],\
                    label='Mode %d'%(i+1));
    
    pyplot.xlabel('Frame Number', fontsize = 16)
    pyplot.ylabel('$\chi$', fontsize = 18)
    pyplot.legend(loc='best')    
    
    pyplot.show();

def writeChivsFrames(chi_list, frames_list, outChiFile, numFrames=20):
    """Write Chi values vs Frames number to see if interpolating a mode reduces
    Chi value or not. 
    """

    from matplotlib import pyplot
    pyplot.xticks(fontsize = 16)
    pyplot.yticks(fontsize = 16)
    pyplot.grid()
    
    linestyles = ['-', '-.','--', ':']
    numModes=len(chi_list)/(numFrames+1)
    print("@> Number of modes in chi list is %d"%numModes)
    for i in range (0, numModes):
        pyplot.plot(frames_list[(i*(numFrames+1)):((i+1)*(numFrames+1))], \
                    chi_list[(i*(numFrames+1)):((i+1)*(numFrames+1))], \
                    linestyle = linestyles[i%4],\
                    label='Mode %d'%(i+1));
    

    pyplot.xlabel('Frame Number', fontsize = 16)
    pyplot.ylabel('$\chi$', fontsize = 18)
    pyplot.legend(loc='best')    
    
    pyplot.savefig(outChiFile);



def main():
    #A simple test of just solvation procedure
  
    protein = parsePDB(sys.argv[1])
    MAX_ATOM = 100000
    W = np.zeros(MAX_ATOM)
    X = np.zeros(MAX_ATOM)
    Y = np.zeros(MAX_ATOM)
    Z = np.zeros(MAX_ATOM)
    
    wDNA = 0.070
    wRNA = 0.126
    wPROT= 0.040
    mol_type = np.zeros(MAX_ATOM, dtype=int)
    thickness = 3.0
    closest_dist = 3.5
    
    buildSolvShell(protein, W, X, Y, Z, closest_dist, thickness, wDNA, wRNA, wPROT, mol_type, MAX_ATOM)

if __name__ == "__main__":
    main()
