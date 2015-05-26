import re

from numpy import linalg, sqrt, loadtxt, savetxt, argmax, floor

from prody import LOGGER, PY3K, dynamics, measure, ensemble, atomic
from prody.dynamics import ANM
from prody.proteins import parsePDB

from random import random

import os.path

from subprocess import call

from time import sleep

import multiprocessing as mp

import thread

try:
	range = xrange
except NameError:
	pass

global vmd
global namd


__all__ = ['perturb_structure', 'solvate_and_ionize', 'initial_minimization','align_two_pdbs', 'extract_final_CA', 'extract_tmd_target', 'run_TMD', 'minimization', 'RMSD_check', 'coMD']

def set_VMD_path():
	a = open('vmd.dat', 'w')
	call(["which","vmd"], stdout=a)
	a.close()
	with open ('vmd.dat', 'r') as a:
		vmd = a.read().replace('\n','')	
	call(["rm","vmd.dat"])

def set_NAMD_path():
	a = open('namd.dat', 'w')
	call(["which","namd2"], stdout=a)
	a.close()
	with open ('namd.dat', 'r') as a:
		namd = a.read().replace('\n','')	
	call(["rm","namd.dat"])

def perturb_structure(initial_pdb, final_pdbn, N):
	a=open(initial_pdb + '_' + final_pdbn + '.running', 'w')
	a.write('going\n')
	a.close()
	pdb = parsePDB(initial_pdb)
	final_pdb = parsePDB(final_pdbn)
	# Current Structure 
	pdb_ca = pdb.ca
	# ANM calculation based on current
	pdb_anm = ANM('pdb ca')
	# Build Hessian Matrix
	pdb_anm.buildHessian(pdb_ca)
	hessian = pdb_anm.getHessian()
	# SVD of Hessian matrix
	U, s, V = linalg.svd(hessian)
	# Cumulative sum vector preparation for metropolis sampling
	eigs = 1/sqrt(s)
	eigs = eigs / sum(eigs)
	eigscumsum = eigs.cumsum()
	# Target Structure
	final_pdb_ca = final_pdb.ca
	# Number of residues on protein structure 		
	size=pdb_ca.getResnums().shape[0]	
	# Difference between current and final structure 
	deviation = final_pdb_ca.getCoords() - pdb_ca.getCoords()
	# Cutoff to check the structure deviated a lot	
	stepcutoff = 0.5*(size**0.5)
	scale_devi = 0.1
	# Scale factor for estimation of energy
	scale_factor = sqrt(np.abs(scale_devi*np.max(s)))
	# counts for metropolis sampling
	count1 = 0 # Up-hill moves
	count2 = 0 # Accepted up-hill moves
	count3 = 0 # Down-hill moves
	# read MC parameter from file
	if os.path.isfile(initial_pdb + '_ratio.dat'):
		MCpara = loadtxt(initial_pdb + '_ratio.dat')
		accept_para = MCpara[4]
		if MCpara[1]>0.95:
			accept_para*=1.5
		elif MCpara[1]<0.85:
			accept_para/=1.5
		else:
			savetxt(initial_pdb + '_status.dat',[1])
	else:
		accept_para = 0.1
	# the best parameter is around 0.9 so that the parameters below than 0.85 and higher than 0.95 are not preferred and adjusted to limits.  
	
	# difference from the target structure is defined as the energy and the minimum is zero. 
	native_dist = buildDistMatrix(final_pdb_ca)
	Ep = 0
	dist = buildDistMatrix(pdb_ca)
	for i in range(size-1):
		for j in range(i+1, size):
			Ep += (native_dist[i][j]-dist[i][j])**2
	pdb_ca_ini = pdb_ca
	ensemble = Ensemble()
	ensemble_final = Ensemble()
	# MC Loop 
	for i in range(N):
		pdb_ca_temp = pdb_ca
		rand = random()	
		ID = argmax(rand<eigscumsum)		
		direction = 2*(random()>0.5)-1

		coords_temp = pdb_ca_temp.getCoords()
		
		coords_temp[0:,0] = coords_temp[0:,0] + direction * U[range(0,len(U),3),ID]
		coords_temp[0:,1] = coords_temp[0:,1] + direction * U[range(1,len(U),3),ID] 
		coords_temp[0:,2] = coords_temp[0:,2] + direction * U[range(2,len(U),3),ID] 
		pdb_ca_temp.setCoords(coords_temp)
		En = 0
		dist = buildDistMatrix(pdb_ca_temp)
		for i in range(size-1):
			for j in range(i+1, size):
				En += (native_dist[i][j]-dist[i][j])**2
		
		if Ep > En:
			count3 += 1
			pdb_ca = pdb_ca_temp
			Ep = En
		elif exp(-(En-Ep)*accept_para)>random():
			pdb_ca = pdb_ca_temp 
			count1 += 1
			count2 += 1
			Ep = En
		else:
			count1 += 1
		
		if calcRMSD(pdb_ca, pdb_ca_ini)	> stepcutoff: 
			break
		
		ensemble.addCoordset(pdb_ca.getCoords())
	
	ensemble_final.addCoordset(pdb_ca.getCoords())
	
	writeDCD(initial_pdb + '_cg.dcd', ensemble)
	writeDCD(initial_pdb + '_final_structure.dcd', ensemble_final)
	ratios = [count2*1.0/N, count2*1.0/count1, count2, N, accept_para ]
	savetxt(initial_pdb + '_ratio.dat', ratios)
	call(["rm", initial_pdb + '_' + final_pdbn + '.running'])

def solvate_and_ionize(pdb_id, chain_id):
	topology_file = 'top_all27_prot_lipid.rtf'
	a = open('ionize_script_' + pdb_id + '_' + chain_id + '.pgn','w')
	a.write("mol new " + pdb_id + ".pdb\nset pro [atomselect top \"protein and not altloc B and not hydrogen and chain " + chain_id + "\"]\n$pro writepdb prop.pdb\npackage require psfgen\ntopology " + topology_file + "\npdbalias residue HIS HSD\npdbalias atom ILE CD1 CD\nsegment R {pdb prop.pdb}\ncoordpdb prop.pdb R\nguesscoord\nwritepdb pro.pdb\nwritepsf pro.psf\npackage require solvate\nsolvate pro.psf pro.pdb -t 15 -o pro_wb\n")
	protein = parsePDB(pdb_id)
	pro1=protein.getHierView()
	seq=''
	for b in pro1.iterChains():
		seq= ''.join([seq,b.getSequence()])
	arg = seq.count('R')
	lys = seq.count('K')
	asp = seq.count('D')
	glu = seq.count('E')
	charge=arg+lys-asp-glu
	if charge>=0:
		na=5
		cl=na+charge
	else:
		charge *=-1
		cl=5
		na=cl+charge
	a.write("package require autoionize\nautoionize -psf pro_wb.psf -pdb pro_wb.pdb -ncl " + str(cl) + " -nna " + str(na) + " -o " + pdb_id + "_" + chain_id + "_ionized\nexit")
	a.close()
	a = open('min_max_center_' + pdb_id + '_' + chain_id + '.pgn','w')
	a.write("mol new " + pdb_id + ".pdb\nset everyone [atomselect top \"all and chain " + chain_id + "\"]\nset a [measure minmax $everyone]\nset b [measure center $everyone]\nset file [open " + pdb_id + "_" +chain_id + "_bound.dat w]\nputs $file [lindex $a 0]\nputs $file [lindex $a 1]\nputs $file $b\nexit")
	a.close()
	a = open('solv_ion_' + pdb_id + '_' + chain_id + '.log', 'w')
	call([vmd,"-dispdev","text","-e",'ionize_script_' + pdb_id + '_' + chain_id + '.pgn'], stdout=a)
	call([vmd,"-dispdev","text","-e",'min_max_center_' + pdb_id + '_' + chain_id + '.pgn'], stdout=a)
	call(['rm','ionize_script_' + pdb_id + '_' + chain_id + '.pgn'], stdout=a)
	call(['rm','min_max_center_' + pdb_id + '_' + chain_id + '.pgn'], stdout=a)
	a.close()
def initial_minimization(pdb_id, chain_id, minimization_length):
	structure = pdb_id + "_" + chain_id + '_ionized.psf'
	coords = pdb_id + "_" + chain_id + '_ionized.pdb'
	parameter = 'par_all27_prot_lipid.prm'
	bound_file = pdb_id + "_" + chain_id + '_bound.dat'
	bounds = loadtxt(bound_file)
	l = floor(abs(bounds[1,:]-bounds[0,:]+[16,16,16]))
	center = bounds[2,:]
	a = open('minimize' + pdb_id + '_' + chain_id + '.conf', 'w')
	a.write("structure " + structure + "\ncoordinates " + coords + "\nset temperature 310\nset outputname " + pdb_id + "_" + chain_id + "_minimized\nfirsttimestep 0\nparaTypeCharmm on\nparameters " + parameter + "\ntemperature $temperature\ncellBasisVector1 " + str(l[0]) + ",0,0\ncellBasisVector2 0," + str(l[1]) + ",0\ncellBasisVector3 0,0," + str(l[2]) + "\ncellOrigin " + str(center[0]) + "," + str(center[1]) + "," + str(center[2]) + "\nexclude scaled1-4\n1-4scaling 1.0\ncutoff 12.\nswitching on\nswitchdist 10.\npairlistdist 13.5\ntimestep 1.0\nrigidBonds none\nnonbondedFreq 1\nfullElectFrequency 1\nstepspercycle 5\nlangevin on\nlangevinDamping 5\nlangevinTemp $temperature\nlangevinHydrogen on\noutputname $outputname\noutputEnergies " + str(minimization_length) + "\noutputPressure " + str(minimization_length) + "\nrestartfreq " + str(minimization_length) + "\ndcdfreq " + str(minimization_length) + "\nxstFreq " + str(minimization_length) + "\nminimize " + str(minimization_length) + "\nreinitvels $temperature\n")
	a.close()
	a = open('minimize' + pdb_id + '_' + chain_id + '.log', 'w')
	num_cores = mp.cpu_count()/2
	call([namd,"+p"+str(num_cores),'minimize' + pdb_id + '_' + chain_id + '.conf'], stdout=a)
	call(["rm", 'minimize' + pdb_id + '_' + chain_id + '.conf'], stdout=a)
	a.close()

		
def align_two_pdbs(starting_rescoor, starting_psf, ending_rescoor, ending_psf, starting_name):
	a = open('alignment_' + starting_psf + '.pgn','w')
	a.write("mol delete all\nmol load psf " + starting_psf + "\nmol addfile " + starting_rescoor + "\nmol load psf " + ending_psf + "\nmol addfile " + ending_rescoor + "\nset s1 [atomselect 0 \"name CA\"]\nset s2 [atomselect 0 \"all\"]\nset s3 [atomselect 1 \"name CA\"]\nset trans_mat [measure fit $s1 $s3]\n$s2 move $trans_mat\n$s1 writepdb "+ starting_name + "\nexit" )
	a.close()
	a = open('align_' + starting_psf + '.log', 'w')
	call([vmd,"-dispdev","text","-e",'alignment_' + starting_psf + '.pgn'], stdout=a)
	call(["rm", 'alignment_' + starting_psf + '.pgn'], stdout=a)
	a.close()
	
def extract_final_CA(ending_rescoor, ending_psf, target_name):
	a = open("ca_extract_" + ending_psf + ".pgn", 'w')
	a.write("mol delete all\nmol load psf " + ending_psf + "\nmol addfile " + ending_rescoor + "\nset s1 [atomselect 0 \"name CA\"]\n$s1 writepdb "+ target_name + "\nexit")
	a.close()
	a = open('ca_extract_' + ending_psf + '.log', 'w')
	call([vmd,"-dispdev","text","-e","ca_extract_" + ending_psf + ".pgn"], stdout=a)
	call(["rm","ca_extract_" + ending_psf + ".pgn"], stdout=a)
	a.close()

def extract_tmd_target(starting_rescoor, starting_psf, final_pdb, final_dcd, adjust):
	a = open('tmd_target_' + starting_psf + '.pgn','w')
	a.write("mol delete all\nmol load psf " + starting_psf + "\nmol addfile " + starting_rescoor + "\nmol load pdb " + final_pdb + "\nmol addfile " + final_dcd + "\nset s1 [atomselect 0 \"name CA\"]\nset s2 [atomselect 0 \"all\"]\nset s3 [atomselect 1 \"name CA\"]\nset trans_mat [measure fit $s1 $s3]\n$s3 move $trans_mat\n$s1 set {x y z} [$s3 get {x y z}]\n $s2 set occupancy 0\n$s1 set occupancy 1\n$s2 writepdb "+ adjust + "\nexit" )
	a.close()
	a = open('etmd_' + starting_psf + '.log','w')
	call([vmd,"-dispdev","text","-e",'tmd_target_' + starting_psf + '.pgn'], stdout=a)
	call(["rm",'tmd_target_' + starting_psf + '.pgn'])
	a.close()

def run_TMD(pdb_id, chain_id,starting_rescoor, starting_resvel, starting_resxsc, adjust, MD_length):
	structure = pdb_id + '_' + chain_id + '_ionized.psf'
	coords = pdb_id + '_' + chain_id + '_ionized.pdb'
	parameter = 'par_all27_prot_lipid.prm'
	a = open('targeted_' + pdb_id + '_' + chain_id + '.conf', 'w')
	a.write("structure " + structure + "\ncoordinates " + coords + "\nset temperature 310\nset outputname " + pdb_id + "_" + chain_id + "_moved\nset restartname res\nbincoordinates " + starting_rescoor + "\nbinvelocities " + starting_resvel + "\nextendedSystem " + starting_resxsc + "\nfirsttimestep 0\nparaTypeCharmm on\nparameters " + parameter + "\nwrapWater on\nwrapAll on\nexclude scaled1-4\n1-4scaling 1.0\ncutoff 12.\nswitching on\nswitchdist 10.\npairlistdist 13.5\ntimestep 1.0\nrigidBonds none\nnonbondedFreq 1\nfullElectFrequency 1\nstepspercycle 5\nPME yes\nPMEGridSpacing 1.0\nlangevin on\nlangevinDamping 5\nlangevinTemp $temperature\nlangevinHydrogen on\nTMD on\nTMDk 200\nTMDOutputFreq " + str(MD_length) + "\nTMDFile " + adjust + "\nTMDFirstStep 0\nTMDLastStep " + str(MD_length) + "\noutputname $outputname\nrestartName $restartname\noutputEnergies " + str(MD_length) + "\noutputPressure " + str(MD_length) + "\nrestartfreq " + str(MD_length) + "\ndcdfreq " + str(MD_length) + "\nxstFreq " + str(MD_length) + "\nrun " + str(MD_length) + "\n")
	a.close()
	a = open('tmd_' + pdb_id + '_' + chain_id + '.log','w')
	num_cores = mp.cpu_count()/2
	call([namd,"+p"+str(num_cores),'targeted_' + pdb_id + '_' + chain_id + '.conf'],stdout=a)
	rm(["rm",'targeted_' + pdb_id + '_' + chain_id + '.conf'])
	a.close()	

def minimization(pdb_id, chain_id, starting_rescoor, starting_resvel, starting_resxsc, minimization_length):
	structure = pdb_id + '_' + chain_id + '_ionized.psf'
	coords = pdb_id + '_' + chain_id + '_ionized.pdb'
	parameter = 'par_all27_prot_lipid.prm'
	a = open('after_minimization_' + pdb_id + '_' + chain_id + '.conf', 'w')
	a.write("structure " + structure + "\ncoordinates " + coords + "\nset temperature 310\nset outputname " + pdb_id + "_" + chain_id  + "_minimized\nset restartname res\nbincoordinates " + starting_rescoor + "\nbinvelocities " + starting_resvel + "\nextendedSystem " + starting_resxsc + "\nfirsttimestep 0\nparaTypeCharmm on\nparameters " + parameter + "\nwrapWater on\nwrapAll on\nexclude scaled1-4\n1-4scaling 1.0\ncutoff 12.\nswitching on\nswitchdist 10.\npairlistdist 13.5\ntimestep 1.0\nrigidBonds none\nnonbondedFreq 1\nfullElectFrequency 1\nstepspercycle 5\nPME yes\nPMEGridSpacing 1.0\nlangevin on\nlangevinDamping 5\nlangevinTemp $temperature\nlangevinHydrogen on\noutputname $outputname\nrestartName $restartname\noutputEnergies " + str(minimization_length) + "\noutputPressure " + str(minimization_length) + "\nrestartfreq " + str(minimization_length) + "\ndcdfreq " + str(minimization_length) + "\nxstFreq " + str(minimization_length) + "\nminimize " + str(minimization_length) + "\nreinitvels $temperature\n")
	a.close()
	a = open('minimize_' + pdb_id + '_' + chain_id + '.log','w')
	num_cores = mp.cpu_count()/2
	call([namd,"+p"+str(num_cores),'after_minimization_' + pdb_id + '_' + chain_id + '.conf'], stdout=a)	
	a.close()

def RMSD_check(starting_pdb_id, starting_chain_id, starting_rescoor, ending_pdb_id, ending_chain_id, ending_rescoor):
	a = open('rmsd_check_' + starting_pdb_id + '_' + ending_pdb_id + '.pgn','w')
	a.write("mol delete all\nmol load psf " + starting_pdb_id + "_" + starting_chain_id + "_ionized.psf\nmol addfile " + starting_rescoor + "\nset sel1 [atomselect top \"name CA\"]\nset sel1a [atomselect top all]\nmol load psf " + ending_pdb_id + "_" + ending_chain_id + "_ionized.psf\nmol addfile " + ending_rescoor + "\nset sel2 [atomselect top \"name CA\"]\nset sel2a [atomselect top all]\nset trans_mat [measure fit $sel2 $sel1]\n$sel2a move $trans_mat\nset rmsd [measure rmsd $sel2 $sel1]\nset file [open rmsd.dat w]\nputs $file $rmsd\nexit\n")
	a.close()
	a = open('rmsd.log','w')
	call([vmd,"-dispdev","text","-e",'rmsd_check_' + starting_pdb_id + '_' + ending_pdb_id + '.pgn'], stdout=a)
	call(["rm", 'rmsd_check_' + starting_pdb_id + '_' + ending_pdb_id + '.pgn'])
	a.close()
	
def coMD(initial_pdb, initial_chain, final_pdb, final_chain, coMD_cycle=None, ANM_cycle=None, MD_length=None, minimization_length=None):
	import hwloc
	top = hwloc.Topology()
	topology.load()
	cpu_count = top.get_nbobjs_bey_type(hwloc.OBJ_CORE)
	if coMD_cycle == None:
		coMD_cycle = 1
	if ANM_cycle == None:
		ANM_cycle = 100000
	if MD_length == None:
		MD_length = 1000000
	if minimization_length == None:
		minimization_length = 2500
	set_VMD_path()
	set_NAMD_path()
	print "solvating\n"
	solvate_and_ionize(initial_pdb, initial_chain)
	solvate_and_ionize(final_pdb, final_chain)
	print "minimizing\n"
	thread.start_new_thread(initial_minimization, (initial_pdb, initial_chain, minimization_length, ))
	thread.start_new_thread(initial_minimization, (final_pdb, final_chain, minimization_length, ))
	#initial_minimization(initial_pdb, initial_chain, minimization_length)
	#initial_minimization(final_pdb, final_chain, minimization_length)
	while os.path.isfile('minimize' + initial_pdb + '_' + initial_chain + '.conf') or os.path.isfile('minimize' + final_pdb + '_' + final_chain + '.conf'):
		sleep(60)
	call(["cp",initial_pdb + "_" + initial_chain + "_minimized.dcd", initial_pdb + "_" + initial_chain + "_minimized0.dcd"])
	call(["cp",final_pdb + "_" + final_chain + "_minimized.dcd", final_pdb + "_" + final_chain + "_minimized0.dcd"])
	rmsd=[]
	RMSD_check(initial_pdb, initial_chain, initial_pdb + "_" + initial_chain + "_minimized.coor", final_pdb, final_chain, final_pdb + "_" + final_chain + "_minimized.coor")
	a = loadtxt('rmsd.dat')
	rmsd.append(a)
	for i in range(coMD_cycle):
		print "coMD cycle " + str(i+1) + "\n"
		print "aligning structures\n"
		align_two_pdbs(initial_pdb + "_" + initial_chain + "_minimized.coor", initial_pdb + "_" + initial_chain + "_ionized.psf", final_pdb + "_" + final_chain + "_minimized.coor", final_pdb + "_" + final_chain + "_ionized.psf" , initial_pdb + "_" + initial_chain + "_starting.pdb")
		align_two_pdbs(final_pdb + "_" + final_chain + "_minimized.coor", final_pdb + "_" + final_chain + "_ionized.psf", initial_pdb + "_" + initial_chain + "_minimized.coor", initial_pdb + "_" + initial_chain + "_ionized.psf" , final_pdb + "_" + final_chain + "_starting.pdb")  
		print "extract final CA's\n"
		extract_final_CA(final_pdb + "_" + final_chain + "_minimized.coor", final_pdb + "_" + final_chain + "_ionized.psf", initial_pdb + "_" + initial_chain + "_target.pdb")
		extract_final_CA(initial_pdb + "_" + initial_chain + "_minimized.coor", initial_pdb + "_" + initial_chain + "_ionized.psf", final_pdb + "_" + final_chain + "_target.pdb")  
		print "ANM perturbation\n"
		thread.start_new_thread(perturb_structure, (initial_pdb + "_" + initial_chain + "_starting.pdb", initial_pdb + "_" + initial_chain + "_target.pdb", ANM_cycle, ))
		thread.start_new_thread(perturb_structure, (final_pdb + "_" + final_chain + "_starting.pdb", final_pdb + "_" + final_chain + "_target.pdb", ANM_cycle, ))
		while os.path.isfile(initial_pdb + "_" + initial_chain + "_starting.pdb_", initial_pdb + "_" + initial_chain + "_target.pdb.running") or os.path.isfile(final_pdb + "_" + final_chain + "_starting.pdb_", final_pdb + "_" + final_chain + "_target.pdb.running"):
			sleep(60)
		print "TMD extraction\n"
		extract_tmd_target(initial_pdb + "_" + initial_chain + "_minimized.coor", initial_pdb + "_" + initial_chain + "_ionized.psf", final_pdb + "_" + final_chain + "_target.pdb", initial_pdb + "_" + initial_chain + "_starting.pdb_final_structure.dcd", initial_pdb + "_" + initial_chain + "_adjust.pdb")
		extract_tmd_target(final_pdb + "_" + final_chain + "_minimized.coor", final_pdb + "_" + final_chain + "_ionized.psf", initial_pdb + "_" + initial_chain + "_target.pdb", final_pdb + "_" + final_chain + "_starting.pdb_final_structure.dcd", final_pdb + "_" + final_chain + "_adjust.pdb")
		print "running TMD\n"
		run_TMD(initial_pdb, initial_chain, initial_pdb + "_" + initial_chain + "_minimized.coor", initial_pdb + "_" + initial_chain + "_minimized.vel", initial_pdb + "_" + initial_chain + "_minimized.xsc", initial_pdb + "_" + initial_chain + "_adjust.pdb", MD_length)
		run_TMD(final_pdb, final_chain, final_pdb + "_" + final_chain + "_minimized.coor", final_pdb + "_" + final_chain + "_minimized.vel", final_pdb + "_" + final_chain + "_minimized.xsc", final_pdb + "_" + final_chain + "_adjust.pdb", MD_length)
		call(["cp",initial_pdb + "_" + initial_chain + "_moved.dcd", initial_pdb + "_" + initial_chain + "_moved" + str(i+1) + ".dcd"])
		call(["cp",final_pdb + "_" + final_chain + "_moved.dcd", final_pdb + "_" + final_chain + "_moved" + str(i+1) + ".dcd"])
		print "minimization\n"
		minimization(initial_pdb, initial_chain, initial_pdb + "_" + initial_chain + "_moved.coor", initial_pdb + "_" + initial_chain + "_moved.vel", initial_pdb + "_" + initial_chain + "_moved.xsc", minimization_length)
		minimization(final_pdb, final_chain, final_pdb + "_" + final_chain + "_moved.coor", final_pdb + "_" + final_chain + "_moved.vel", final_pdb + "_" + final_chain + "_moved.xsc", minimization_length)
		call(["cp",initial_pdb + "_" + initial_chain + "_minimized.dcd", initial_pdb + "_" + initial_chain + "_minimized" + str(i+1) + ".dcd"])
		call(["cp",final_pdb + "_" + final_chain + "_minimized.dcd", final_pdb + "_" + final_chain + "_minimized" + str(i+1) + ".dcd"])
		RMSD_check(initial_pdb, initial_chain, initial_pdb + "_" + initial_chain + "_minimized.coor", final_pdb, final_chain, final_pdb + "_" + final_chain + "_minimized.coor")
		a = np.loadtxt('rmsd.dat')
		print "RMSD at the end of cycle " + str(i+1) + " is " + str(a) + ".\n" 
		rmsd.append(a)
		if rmsd[-1]<1.5 or rmsd[-2]-rmsd[-1]<0.15:
			break
			

	
		
	
	
			
