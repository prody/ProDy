# -*- coding: utf-8 -*-
"""Perform Small Angle Scattering profile calculations of ANM modes 
   to enlighten solution structure of proteins. 
"""

from prody import *
import numpy as np
import os
import sys
import argparse
from math import *
import timeit
import time

from solvate import buildSolvShell
import saxstools

from matplotlib import pyplot

__all__ = ['prody_saxs', 'calcSAXSPerModel']

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
    f=0

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

def showSAXSProfiles(exp_data, model_data):
    #Read experimental data
    #Read model data
    #Plot data
    from matplotlib import pyplot;
    from pylab import genfromtxt;
    mat0 = genfromtxt(exp_data);
    mat1 = genfromtxt(model_data);
    pyplot.plot(mat0[:,0], mat0[:,1], label = "Experimental");
    pyplot.plot(mat1[:,0], mat1[:,1], label = "Model");
    pyplot.legend();
    pyplot.show();


def calcSaxsChi(q_exp, I_q_exp, sigma_q, q_model, I_q_model):
    """
    This function calculates Chi squared value between SAXS profile of an experimental/
    simulated data and SAXS profile of a model. The model data can be obtained
    from an experimental macromolecular structure, a molecular dynamics 
    conformation or a normal mode.
    """
    #Read model data
#    q_model, I_q_model = np.loadtxt(model_I_q_file, unpack=True)

    #Check length of model data and compare it with experimental data set. 
    #Chech if model data was taken at experimental q values.
    #If not, interpolate

    #Ensure that q_exp[0] and q_model[0] points to the same value. 
    #Compare head of data sets to find c value.
    if(q_exp[0]==q_model[0]):
        c=I_q_exp[0]/I_q_model[0]
    else:
        print "Experimental and theoretical Q[0] values do not match!"
        sys.exit(-1)

    #Calculate Chi    
    diff_array=((I_q_exp-c*I_q_model)/sigma_q)
    chi_sqrt=(np.sum(np.square(diff_array)))/len(q_exp)
    return sqrt(chi_sqrt)


def parseSaxsData(I_q_file, simulated=False, isLogScale=True):
    """
    This function parses Small Angle X-ray Scattering data.
    ARGS:
    I_q_file:   This file is a plain text file with at least two
                columns, namely, q column and I_q column. If the 
                file is an experimental data file, there should be a third
                column of sigma_q values, namely, error bars. 
    simulated:  This parameter is a boolean parameter. If the data
                is simulated SAXS data, you should make it True so 
                that it will assign 1.0 values to all sigma_q values.
              
    isLogScale: This parameter is also a boolean parameter. prody_saxs 
                module requires the intensity (I_q) data to be in 
                logscale for all calculations. Therefore, if the 
                intensities are not in logscale, you should make 
                isLogScale=False, so that the function return the 
                input data in log10 scale. 
    """
    #Check if data contains error bars or if it is simulated data!
    if(simulated==False):
        #Read experimental data
        q_exp, I_q_exp, sigma_q = np.loadtxt(I_q_file, unpack=True)
    elif(simulated==True):
        #Read simulated data
        q_exp, I_q_exp = np.loadtxt(I_q_file, unpack=True)

        #Since it is simulated data, assign 1.0 values to errorbars.
        sigma_q=np.ones(len(q_exp))
            
    #Check if data is in log scale or not.
    if(isLogScale==False):
        return (q_exp, np.log10(I_q_exp), np.log10(sigma_q))
    else:
        return (q_exp, I_q_exp, sigma_q)

def calcSAXSPerModel(calphas, numCalphas, I, Q_exp):
    ####Assign coordinates to a temporary array!#################################
    wDNA =0.070
    wRNA =0.126
    wPROT=0.040
    thickness=3.0
    closest_dist=3.5
    delta=0.0
#    total_atom=0
    pdb_flag=1
    solvent_flag=1
    total_atom=numCalphas
    MAX_ATOM=100000
    W=np.zeros(MAX_ATOM)
    X=np.zeros(MAX_ATOM)
    Y=np.zeros(MAX_ATOM)
    Z=np.zeros(MAX_ATOM)
    mol_type=np.zeros(MAX_ATOM, dtype=int)
    cgatom_num=np.zeros(MAX_ATOM, dtype='int32')
    cgatom=calphas.getResnames()
    origCoords=calphas.getCoords()
    
    X[:len(origCoords[:,0])]=origCoords[:,0]
    Y[:len(origCoords[:,1])]=origCoords[:,1]
    Z[:len(origCoords[:,2])]=origCoords[:,2]

    #At first, solvate the protein:
    #print "@> Solvating the system"
    #start = timeit.timeit()
#    start = time.clock()
    total_atom=buildSolvShell(calphas, W, X, Y, Z, closest_dist, thickness, \
                              wDNA, wRNA, wPROT, mol_type, MAX_ATOM)
#    end= time.clock()
#    print "Solvation time: %f"%(end - start)

    #end = timeit.timeit()
    #print (end - start)
#    total_atom=buildSolvShell(calphas, W, X, Y, Z, closest_dist, thickness, \
#                              wDNA, wRNA, wPROT, mol_type, MAX_ATOM)
    #C version of this building solvation shell. I have to check which version works well. 
    #start = timeit.timeit()
#    writePDB('model_no_water.pdb', calphas)
#    total_atom=saxstools.cgSolvateNumeric("model_no_water.pdb", X, Y, Z, W, \
    #                                      wDNA, wRNA, wPROT, \
    #                                      thickness, closest_dist, \
    #                                      pdb_flag, solvent_flag, MAX_ATOM)
#    end = timeit.timeit()
#    print (end - start)

    #print "@> Total number of atoms after solvation: %d" %total_atom
    
    for j in range(0, total_atom):
        if(j<numCalphas):
            cgatom_num[j]=numericResname(cgatom[j])
        else:
            cgatom_num[j]=numericResname("SOL")

    #Now, let reread solvated structure. And calculate SAXS profile of it with 
    #solvated_protein = parsePDB("cg_solvated.pdb")

    #You need a function to check experimental SAXS data. It should give Q_exp
    #and it should check whether experimental data is in log scale or not!
#    start = time.clock()
#    _saxs.calcSAXSNumeric(I, X, Y, Z, total_atom, cgatom_num, W, Q_exp, len(Q_exp))
    saxstools.calcSAXSNumeric(I, X, Y, Z, total_atom, cgatom_num, W, Q_exp, len(Q_exp))
    ####Finish key part!#########################################################
#    end= time.clock()
#    print (end - start)

def interpolateMode(calphas, mode,\
                    out_pdb_file,\
                    Q_exp, I_q_exp, sigma_q,\
                    max_chi,\
                    rmsdScalingCoef=3.0,\
                    numFrames=20):
    
    chi_list=[]
    frames_list=[]

    eigenvalue=mode.getEigval()
    eigenvector=np.transpose(mode.getEigvec())

    i=mode.getIndex()
    mod_num=None
    origCoords=calphas.getCoords()
    numCalphas=calphas.numAtoms()
    I_model=np.zeros(len(Q_exp))
    # setup toolbar
    sys.stdout.write("@> Calculating SAXS profiles for nonzero mode %d: " % (i+1))
    sys.stdout.write("[%s]" % (" " * (numFrames+1)))
    sys.stdout.flush()
    sys.stdout.write("\b" * (numFrames+2)) # return to start of line, after '['
    prody.LOGGER.timeit('_intplt_mode')
    invEigVal=(1.0/eigenvalue)
    for j in range((-numFrames/2), ((numFrames/2)+1)):
        coeff=j*rmsdScalingCoef*invEigVal*2.0/numFrames
        
        newCoords=calphas.getCoords().flatten()+(coeff*eigenvector)
        calphas.setCoords(newCoords.reshape((numCalphas, 3), order='C'))
        calcSAXSPerModel(calphas, numCalphas, I_model, Q_exp)
        chi=calcSaxsChi(Q_exp, I_q_exp, sigma_q, Q_exp, I_model)
        chi_list.append(chi)
        frames_list.append(j)
        if(chi<max_chi):
            max_chi=chi
            mod_num=i
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


def showChivsFrames(chi_list, frames_list, numFrames):
    numModes=len(chi_list)/(numFrames+1)
    print "@> Number of modes in chi list is %d"%numModes
    for i in range (0, numModes):
        pyplot.plot(frames_list[(i*(numFrames+1)):((i+1)*(numFrames+1))], \
                    chi_list[(i*(numFrames+1)):((i+1)*(numFrames+1))], label='Mode %d'%(i+1));
    
    pyplot.xlabel('Frame Number')
    pyplot.ylabel('Chi')
    pyplot.legend(loc='best')    
    
    pyplot.show();


#def main():
def prody_saxs():
    #Set default values for calculations.

    #0-Parse all arguments
    parser = argparse.ArgumentParser(description=\
    'Find the best mode fitting to a given SAXS profile.',\
    epilog=\
    'Example: python prody_saxs.py 4ake_chainA.pdb -s 1ake_chainA_saxs_w_yerrorbars.dat -o output.txt')
    parser.add_argument('pdb_file', nargs='?', type=str, \
                        help='Mandatory input pdb file.')
    
    parser.add_argument('-s', '--saxs', type=str, dest='saxs_file',\
                        help='Mandatory experimental/simulated SAXS profile.')

    parser.add_argument('-n', '--nmodes', type=int, dest='numModes',\
                        default=5, \
                        help='Number of nonzero modes to be used in calculations.')

    parser.add_argument('-c', '--coeff', type=float, dest='scalCoeff',\
                        default=3.0, \
                        help='''Scaling coefficient for interpolating a single mode. A positive value of larger than 1 is recommended.''')

    parser.add_argument('-p', '--out-pdb', \
                        type=str, dest='out_pdb_file', \
                        default='best_model.pdb', \
                        help='Output pdb file for the best model.')

    parser.add_argument('-o', '--out-saxs', type=str, dest='out_saxs_file', \
                        help='Output SAXS profile for the best model')

    parser.add_argument('-v', '--version', action='version', \
                        version='%(prog)s (version 0.1)')

    args = parser.parse_args()
    
#    print args.pdb_file
#    print args.saxs_file
#    print args.numModes
#    print args.scalCoeff
#    print args.out_pdb_file
#    print args.out_saxs_file

    
    #1-This module produces normal modes of a protein structure and 
    #by using anisotropic network model.    
    protein = parsePDB(args.pdb_file)
    calphas = protein.select('calpha')
    origCoords=calphas.getCoords()
    
    anm = ANM('ANM Analysis')
    anm.buildHessian(calphas, cutoff=15.0)

    numCalphas=calphas.numAtoms()

    modes=anm.calcModes(n_modes=args.numModes, zeros=False)

#    mode_ensemble=traverseMode(anm[0], calphas, n_steps=10, rmsd=3.0)
#    print mode_ensemble[0]
#    writePDB('traverseMode_0.pdb', mode_ensemble, csets=None, autoext=True)
#    sys.exit(-1)

    #Parse experimental/simulated SAXS data.
    Q_exp, I_q_exp, sigma_q=parseSaxsData(args.saxs_file, simulated=False, isLogScale=True)

    I_model=np.zeros(len(Q_exp))
    prody.LOGGER.info('Number of experimental data points=%.d'%len(Q_exp))

    #Calculate a SAXS profile for initial pdb file by using Fast-SAXS approach.
    prody.LOGGER.info('Solvating the system and calculating SAXS profile.')
    calcSAXSPerModel(calphas, numCalphas, I_model, Q_exp)

    #Get initial chi value between pdb_file and saxs_file
    max_chi=calcSaxsChi(Q_exp, I_q_exp, sigma_q, Q_exp, I_model)
    prody.LOGGER.info('Chi value between pdb file and experimental SAXS profile=%.3f'%max_chi)

    #A SAXS profile is produced for each model in a mode using Fast-SAXS approach.
    #Now, lets do it with python code rather than calling an external fast-saxs-pro program

    chi_overall=[]
    frames_overall=[]

    #All modes are interpolated in +/- directions. rmsdScalingCoef is scaling coefficient of interpolation.
    #  A positive value of larger than 1 is recommended.
    rmsdScalingCoef=3.0
    numFrames=20
    mod_num=None

    prody.LOGGER.timeit('_intplt_mode')    
    for i in range (0, args.numModes):
            eigenvalue=anm[i].getEigval()
            eigenvector=np.transpose(anm[i].getEigvec())

            # setup toolbar
            sys.stdout.write("@> Calculating SAXS profiles for nonzero mode %d: " % (i+1))
            sys.stdout.write("[%s]" % (" " * (numFrames+1)))
            sys.stdout.flush()
            sys.stdout.write("\b" * (numFrames+2)) # return to start of line, after '['

            invEigVal=(1.0/eigenvalue)
            for j in range((-numFrames/2), ((numFrames/2)+1)):
                #coeff=j*rmsdScalingCoef*invEigVal*2.0/numFrames
                coeff=j*(args.scalCoeff)*invEigVal*2.0/numFrames
                
                newCoords=calphas.getCoords().flatten()+(coeff*eigenvector)
                calphas.setCoords(newCoords.reshape((numCalphas, 3), order='C'))
                calcSAXSPerModel(calphas, numCalphas, I_model, Q_exp)
                chi=calcSaxsChi(Q_exp, I_q_exp, sigma_q, Q_exp, I_model)
                chi_overall.append(chi)
                frames_overall.append(j)
                if(chi<max_chi):
                    max_chi=chi
                    mod_num=i
                    writePDB(args.out_pdb_file, calphas)

                    #           extendModel(calphas, 'calphas', protein)
                    #           writePDB('best_model.pdb', protein)
                    
                #Reset coordinates to the original values
                calphas.setCoords(origCoords)
            
                sys.stdout.write('#')
                sys.stdout.flush()
            sys.stdout.write("\n")

    prody.LOGGER.report('SAXS profile calculations were performed in %2fs.', '_intplt_mode')
    
    showChivsFrames(chi_overall, frames_overall, numFrames)

    #The model with the lowest Chi value is written to a pdb file.
    prody.LOGGER.info('Chi value between the best model and the experimental SAXS data=%.3f'%np.amin(chi_overall))

    #    print np.argmin(chi_overall)
    best_model_all = parsePDB('best_model.pdb')
    best_model_calphas = best_model_all.select('calpha')
    calcSAXSPerModel(best_model_calphas, numCalphas, I_model, Q_exp)
            
if __name__ == "__main__":
    prody_saxs()

