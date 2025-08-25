__author__ = 'Carlos Ventura'
__email__ = ['carlosventura@stonybrook.edu']

#!/usr/bin/python
import os
import subprocess
import shutil
import glob
import sys
import ast
from pathlib import Path
import importlib

current_dir = Path(__file__).resolve().parent

package_dir = current_dir / "DruGUI-script"

if package_dir.is_dir() and str(package_dir) not in sys.path:
    sys.path.insert(0, str(package_dir)) 

try:
    druggability = importlib.import_module("druggability")
except ImportError as e:
    raise ImportError(
        f"Could not import 'druggability' from {package_dir}. "
        f"({e})"
    )

from numpy import array, ceil, histogramdd, arange
from prody.proteins.compare import matchAlign
from prody.proteins.pdbfile import parsePDB, writePDB
from prody.measure.measure import calcCenter
from prody.measure.transform import moveAtoms, wrapAtoms
from prody.trajectory.dcdfile import DCDFile, parseDCD
from prody.trajectory.psffile import parsePSF
from prody.trajectory.trajectory import Trajectory

from druggability.grid import OpenDX

__all__ = ['drugui_prepare', 'drugui_analysis']

def drugui_data(vmd_executable):
    """Prepares the nessecary files needed for both preparing and analyzing druggability simulations"""


    os.chdir('/Users/carlosventura/Desktop/prody_drugui/ProDy/prody/drugui/DruGUI-script')
    Druggability_path = os.getcwd()
    global PROBEDATA
    global PROBETYPES
    global PROBETOPPAR

    probedata = f"""
    set PROBEDATA [dict create]
    set PROBETYPES [dict create core "Core probes"]
    dict set PROBETYPES polar "Polar probes"
    dict set PROBETYPES hydrophobe "Hydrophobes"
    dict set PROBETYPES negative "Negatively charged"
    dict set PROBETYPES positive "Positively charged"
    dict set PROBETYPES ring5 "5-membered rings"
    dict set PROBETYPES ring6 "6-membered rings"
    set PROBETOPPAR [dict create PBDA "probe2.top probe.prm"]
    dict set PROBETOPPAR CGenFF "top_all36_cgenff.rtf par_all36_cgenff.prm"
    set PACKAGEPATH {Druggability_path}
    """
    probedata += """
    foreach {key} [dict keys $PROBETYPES] {
    dict set PROBEDATA $key ""
    }

    set inf [open [file join $PACKAGEPATH "probeV2.dat"] r]
    foreach line [split [read -nonewline $inf] \\n] {
        if {[llength $line] < 3 || [string equal -length 1 $line "#"]} {continue}
        set resi [lindex $line 0]
        set key [lindex $line 1]
        if {![dict exists $PROBEDATA $resi]} {
            dict set PROBEDATA $resi default 0
            dict set PROBEDATA $resi alias ""
            dict set PROBEDATA $resi atomnames ""
            dict set PROBEDATA $resi charge 0
            dict set PROBEDATA $resi source ""
        }
        if {$key == "default" && [lindex $line 2] > 0} {
            dict lappend PROBEDATA defaults $resi
        }
        if {$key == "type" && [lindex $line 2] > 0} {
            dict lappend PROBEDATA [lindex $line 2] $resi
        }
        dict set PROBEDATA $resi $key [lrange $line 2 end]
        if {$key == "alias" && [expr [llength [dict get $PROBEDATA $resi $key]] % 2] != 0} {
            error "Problem in aliases of $resi!"
        }
     }
    close $inf

    set ipronames "C2 H21 C1 H11 H12 H13 C3 C4 C5 C6 H31 H32 H33 OH2 HO2"
    dict for {src toppar} $PROBETOPPAR {
        set inf [open [file join $PACKAGEPATH [lindex $toppar 0]] r]
        puts "$inf"
        set resi "____"
        set prev "____"
        foreach line [split [read -nonewline $inf] \\n] {
            if {[dict exists $PROBEDATA $resi] && [string range $line 0 3] == "ATOM"} {
                dict set PROBEDATA $resi atomnames "[dict get $PROBEDATA $resi atomnames] [lindex $line 1]"
            } elseif {[string range $line 0 3] == "RESI"} {
              set prev $resi
              set resi [string trim [lindex $line 1]]
              if {[dict exists $PROBEDATA $resi]} {
                dict set PROBEDATA $resi charge [lindex $line 2]
                dict set PROBEDATA $resi source $src
              }
              if {[dict exists $PROBEDATA $prev]} {
                set thisnames [dict get $PROBEDATA $prev atomnames]
                set remove 0
                set aliases [dict get $PROBEDATA $prev alias]
                foreach {name alias} $aliases {
                    if {[lsearch $ipronames $name] == -1} {
                        error "Invalid $prev alias, $name is not found in IPRO!"
                        #set remove 1
                    } elseif {[lsearch $thisnames $alias] == -1} {
                      error "Invalid $prev alias, $alias is not found in $prev!"
                    #set remove 1
                    }
                    #if {$remove} {
                    #  set PROBEDATA [dict remove $PROBEDATA $prev]
                    #  break
                    #}
                }
            }
        }
     }
    close $inf
    }

    set probe_data [open "probedata.txt" w]
    puts $probe_data $PROBEDATA
    close $probe_data

    set probetypes [open "probetypes.txt" w]
    puts $probetypes $PROBETYPES
    close $probetypes

    set probetoppar [open "probetoppar.txt" w]
    puts $probetoppar $PROBETOPPAR
    close $probetoppar
    exit
    """
    with open('probedata.tcl', 'w') as tcl_file:
        tcl_file.write(probedata)
    vmd = vmd_executable
    subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'probedata.tcl'])
    with open('probedata.txt', 'r') as output_file:
            PROBEDATA = output_file.read().strip()

def drugui_prepare(pdb, psf, **kwargs):
    """Prepares protein systems for druggability simulations. Both, a pdb and psf file are required to prepare druggability simulations.
    DruGUI needs your computer's VMD executable to solvate, neutralize, and add probes to the system.
    
    prefix: The name the prepared files will be given
    Probes: Whether the system will have probes. Systems without probes will have a different equilibration and minimization process.
    Systems without probes will be solvated and neutralized similarily to druggability simulations with probes.
    probes_and_percent: The probe and probe composition the user can customized for druggability simulations.
    solvent_padding: The minimum solvent padding in all directions
    boundary_padding: The minimum distance between water/probes and solute
    neutralize: Add counter ions to make the system neutral
    lipids: The system is solvated by considering the lipid bilayer's xy plane
    write_conf: Writes configuration files to run druggability simulations using NAMD
    n_sims: Number of independent simulations
    sim_length: The length of individual simulations
    outdir_location: Files will be written in specified directory
    constrain: Which atoms to constrain in equilibration steps
    vmd: VMD executable needed to prepare the protein system for druggability simulations
    """
    
    global PROBEDATA
    global PROBETYPES
    protein_pdb = pdb
    protein_psf = psf
    prefix = kwargs.pop('prefix', "md")
    Probes = kwargs.pop('Probes', True)
    probes_and_percent = kwargs.pop("probes_and_percent", 
                                    {"IPRO": 16, "IMID": 14, "ACTT": 14, "ACAM": 14, "IPAM": 14, "IBTN": 14, "BENZ": 14})
    solvent_padding = kwargs.pop('solvent_padding', 15)
    boundary_padding = kwargs.pop('boundary_padding', 2.4)
    neutralize = kwargs.pop('neutralize', True)
    lipids = kwargs.pop('lipids', False)
    write_conf = kwargs.pop('write_conf', True)
    n_sims = kwargs.pop('n_sims', 4)
    sim_length = kwargs.pop('sim_length', 40)
    outdir_location = kwargs.pop('outdir_location', "")
    constrain = kwargs.pop("constrain", "heavy")
    vmd = kwargs.pop("vmd","")
    additional_parameters = kwargs.pop("additional_parameters", [])
    drugui_data(vmd)

    if package_dir.is_dir():
        os.chdir(package_dir)
    else:
        raise FileNotFoundError(f"DruGUI-script not found at {package_dir}")

    Druggability_path = os.getcwd()

    probepsf = os.path.join(Druggability_path, "probe.psf")
    probepdb = os.path.join(Druggability_path, "probe.pdb")
    probetop = os.path.join(Druggability_path, "probe2.top")
    probeprm = os.path.join(Druggability_path, "probe.prm")
    cgenfftop = os.path.join(Druggability_path, "top_all36_cgenff.rtf")
    cgenffprm = os.path.join(Druggability_path, "par_all36_cgenff.prm")
    probebox = 62.3572
    probekey = "name OH2"
    intermediate = os.path.join(Druggability_path, "intermediate")
    drugui_title = "Druggability GUI Version 2.0"

    if not (os.path.exists(probepsf) and os.path.exists(probepdb) and os.path.exists(probetop) and os.path.exists(probeprm) and os.path.exists(cgenfftop) and os.path.exists(cgenffprm)):
        raise ValueError("ERROR", f"One of the probe PSF, PDB, TOP, or PRM files is not found in {Druggability_path}")
                         
    if protein_pdb == "" or protein_psf == "":
        raise ValueError("ERROR", "Both PSF and PDB files must be specified.")
    
    if prefix == "":
        raise ValueError("ERROR", "Prefix must not be left blank.")

    if vmd == "":
        raise ValueError("ERROR", "The location of your VMD executable is needed to setup druggability simulations.")

    par_files_list = [
    "probe.prm", "par_all36_cgenff.prm", "par_all27_prot_lipid_na.inp",
    "par_all36_lipid.prm", "par_all36_prot.prm", "par_all36_carb.prm",
    "par_all36_na.prm", "toppar_water_ions.str"
    ] 

    if additional_parameters:
        par_files_list.extend(additional_parameters)

    parameterfiles = par_files_list

    log = open(f"{prefix}.log",'w')
    with open(f"{prefix}.log",'a') as log:
        log.write(f"{drugui_title}\n")
        log.write(f"Input PSF: {protein_pdb}.\n")
        log.write(f"Input PDB: {protein_psf}.\n")
        log.write(f"Output directory: {outdir_location}.\n")
        log.write(f"Output prefix: {prefix}.\n")
        log.write(f"Intermediate = {intermediate}\n")

        if Probes == True:
            opts = list(probes_and_percent.keys())
            opts_percentage = list(probes_and_percent.values())
            log.write(f"{opts}\n")
            log.write(f"Probe compostion: \n")
            opts_dict = dict(zip(opts, opts_percentage))

            for probe in opts:
                if probe not in PROBEDATA:
                    raise ValueError("ERROR", f"Unknown probe '{probe}' not found in PROBEDATA.")
                
            for probe_percentage in opts_percentage:
                if probe_percentage < 0:
                    raise ValueError("ERROR", "Probe percentages must be a positive number")

            percent_total = sum(opts_percentage)

            if percent_total != 100:
                raise ValueError("ERROR", "Probe percentages must sum up to 100.\nCurrent total is "+ str(percent_total))

            tcl_opts = "set opts [dict create " + " ".join(f"{key} \"{value}\"" for key, value in opts_dict.items()) + "]"
            probe_analysis = f"""
            set PROBEDATA {{{PROBEDATA}}}
            set PROBETOPPAR {{{PROBETOPPAR}}}
            set PROBETYPES {{{PROBETYPES}}}
            {tcl_opts}
            set logfile [open probe_analysis.log" a]

            """
            probe_analysis += """
            set percent_charge [open "percent_charge.txt" w]

            dict for {key info} $PROBEDATA {
                    if {[dict exists $opts "$key"]} {
                        dict with info {
                            set percentage [dict get $opts "$key"]
                            if {![string is digit $percentage]} {
                            error "\"$key $percentage\" is not valid, probe percentages must be positive integers."
                            }
                            dict unset opts "$key"
                            if {$source == "CGenFF"} {
                                set general 1
                            }           
                            dict set probePercent $key $percentage
                            incr probeTotal $percentage
                            set holder "$key $percentage% ($name; $charge e)"
                            puts $percent_charge $holder
                        }
                    }
                }
            close $percent_charge
            exit
            """
            with open('probe_analysis.tcl', 'w') as tcl_file:
                tcl_file.write(probe_analysis)
            subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'probe_analysis.tcl'])
            with open('percent_charge.txt', 'r') as percent_charge_file:
                probe_information = percent_charge_file.read().strip()
            log.write(f"{probe_information} \n")
        else :
            log.write("Probe molecules are not added to the system.")

        if Probes == True:
            i = 0
            n = len(probes_and_percent)
            probePercent = {}
            while i<n:
                probe = opts[i]
                percentage = opts_percentage[i]
                probePercent[probe] = []
                probePercent[probe].append(int(percentage))
                i+=1

        if neutralize == True:
            log.write("System is neutralized by adding counter ions.\n")
        else :
            log.write("System is not neutralized by adding counter ions.\n")
        
        log.write(f"Minimum solvent padding is set to {solvent_padding} Å.\n")
        log.write(f"Minimum solvent boundary is to {boundary_padding} Å.\n")

        if lipids == False:
            log.write(f"System does not contain lipid bilayer.\n")
        else :
            lipids == True
            log.write(f"System does contain lipid bilayer.\n")

        log.write(f"NAMD configuration files for {n_sims} independent {sim_length} ns simuation(s) are written.\n")
    
        #set padding for the system and probes
        padding_x = solvent_padding
        padding_y = solvent_padding
        padding_z = solvent_padding
        padx = 0
        pady = 0
        padz = 0

        if Probes == True:
            padx = 5
            pady = 5
            padz = 5
        else :
            solvent_pad = solvent_padding

        if lipids == True:
            padding_x = 0
            padding_y = 0
            padding_z = solvent_padding
            padx = -3
            pady = -3
            padz = 5
        
        solvate_options = f"{protein_psf} {protein_pdb} -o {intermediate}"
        
        if Probes == True:
            solvate_options += f" -b {boundary_padding} -x {padding_x + padx} -y {padding_y + pady} -z {padding_z + padz}"
            solvate_options += f" +x {padding_x + padx} +y {padding_y + pady} +z {padding_z + padz}"
            solvate_options += f' -spdb "{probepdb}" -spsf "{probepsf}" -stop "{probetop}" -ks "{probekey}" -ws {probebox}\n'

            log.write(f"Command solvate: {solvate_options}")
                
            #Solvating the system with VMD
            tcl_script = f"""
            package require solvate
            solvate {solvate_options}
            set all [atomselect top "protein"]
            set proteincharge [vecsum [$all get charge]]
            set pc_output [open "pc_output.txt" w]
            puts $pc_output $proteincharge
            close $pc_output
            exit
            """
            with open('solvate_script.tcl', 'w') as tcl_file:
                tcl_file.write(tcl_script)
            subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'solvate_script.tcl'])
            with open('pc_output.txt', 'r') as output_file:
                proteincharge = output_file.read().strip()
        else :
            if lipids == True:
                
                solvate_options += f" -b {boundary_padding} -x {padding_x + padx} -y {padding_y + pady} -z {padding_z + padz}"
                solvate_options += f" +x {padding_x + padx} +y {padding_y + pady} +z {padding_z + padz}\n"

                log.write(f"Command solvate: {solvate_options}")
                #Solvating the system with VMD
                tcl_script = f"""
                package require solvate
                solvate {solvate_options}
                set all [atomselect top "protein"]
                set proteincharge [vecsum [$all get charge]]
                set pc_output [open "pc_output.txt" w]
                puts $pc_output $proteincharge
                close $pc_output
                exit
                 """
                with open('solvate_script.tcl', 'w') as tcl_file:
                    tcl_file.write(tcl_script)
                subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'solvate_script.tcl'])
                with open('pc_output.txt', 'r') as output_file:
                    proteincharge = float(output_file.read().strip())
            else:
                solvate_options += f" -b {boundary_padding} -t {solvent_pad}\n"
                log.write(f"Command solvate: {solvate_options}")
                tcl_script = f"""
                package require solvate
                solvate {solvate_options}
                set all [atomselect top "protein"]
                set proteincharge [vecsum [$all get charge]]
                set pc_output [open "pc_output.txt" w]
                puts $pc_output $proteincharge
                close $pc_output
                exit
                """
                with open('solvate_script.tcl', 'w') as tcl_file:
                    tcl_file.write(tcl_script)
                subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'solvate_script.tcl'])
                with open('pc_output.txt', 'r') as output_file:
                    proteincharge = float(output_file.read().strip())

        #Neutraling the system with VMD
        if neutralize == True:
            if Probes == False:
                log.write(f"Ionization: System has a total charge of {proteincharge} electrons.\n")
                if proteincharge > 0:
                    nna = 0
                    ncl = round(proteincharge)
                    log.write(f"Ionization: {ncl} chloride ions will be added.\n")
                else :
                    ncl = 0
                    nna = round(proteincharge)
                    log.write(f"Ionization: {nna} sodium ions will be added.\n")
                ionization_script = f"""
                package require autoionize
                autoionize -psf {intermediate}.psf -pdb {intermediate}.pdb -o {prefix} -from 5 -between 5 -ncl {ncl} -nna {nna} -seg ION
                exit
                 """
                with open('ionization_script.tcl', 'w') as tcl_file:
                    tcl_file.write(ionization_script)
                subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'ionization_script.tcl'])

        if Probes == True:
        #===============================================================================
        # START - ADJUST RELATIVE NUMBER of WATER and PROBE molecules
        # minimum and maximum coordinates for PROTEIN is calculated.
        # if a molecule other than a PROTEIN is of interest, change selection in the next line.
            protein = 'atomselect top "not water and not segid ION \\"WT.*\\""'
            minmaxP = f"measure minmax $protein"
            minP = "lindex $minmaxP 0"
            maxP = "lindex $minmaxP 1"
            minPx = "lindex $minP 0"
            minPy = "lindex $minP 1"
            minPz = "lindex $minP 2"
            maxPx = "lindex $maxP 0"
            maxPy = "lindex $maxP 1"
            maxPz = "lindex $maxP 2"

            padx = int(padding_x)
            pady = int(padding_y)
            padz = int(padding_z)

            if lipids == False:
                selWstring = f"water and name OH2 and x > [expr $minPx-{padx}] and y > [expr $minPy-{pady}] and z > [expr $minPz-{padz}] and x < [expr $maxPx+{padx}] and y < [expr $maxPy+{pady}] and z < [expr $maxPz+{padz}]"
                
            else :
                sel = 'set sel [atomselect top "lipid"]'
                minmaxL = "measure minmax $sel"
                minLz = "expr [lindex [lindex $minmaxL 0] 2] + 10"
                maxLz = "expr [lindex [lindex $minmaxL 0] 2] - 10"
                selWstring = f"water and name OH2 and x > [expr $minPx - {padx}] and y > [expr $minPy - {pady}] and x < [expr $maxPx+{padx}] and y < [expr $maxPy+{pady}] and ((z < [expr $maxPz+{padz}] and z > {maxLz}) or (z < {minLz} and z > [expr $minPz-{padz}]))"
                
            # select waters in the box of requested size

            if lipids == False:
                selWater = f'atomselect top "{selWstring}"'
                nWater = "$selWater num"
                water_script = f"""
                mol load psf {intermediate}.psf
                mol load pdb {intermediate}.pdb
                set protein [{protein}]
                set minmaxP [{minmaxP}]
                set minP [{minP}]
                set maxP [{maxP}]
                set minPx [{minPx}]
                set minPy [{minPy}]
                set minPz [{minPz}]
                set maxPx [{maxPx}]
                set maxPy [{maxPy}]
                set maxPz [{maxPz}]
                set selWater [{selWater}]
                set nWater [{nWater}]
                set indicesWater [$selWater get index]
                set output_file [open "nwater_output.txt" w]
                puts $output_file $nWater
                close $output_file
                exit
                """
            else :
                selWater = f'atomselect top "{selWstring}"'
                nWater = "$selWater num"
                water_script = f"""
                mol load psf {intermediate}.psf
                mol load pdb {intermediate}.pdb
                set protein [{protein}]
                set minmaxP [{minmaxP}]
                set minP [{minP}]
                set maxP [{maxP}]
                set minPx [{minPx}]
                set minPy [{minPy}]
                set minPz [{minPz}]
                set maxPx [{maxPx}]
                set maxPy [{maxPy}]
                set maxPz [{maxPz}]
                {sel}
                set minmaxL [{minmaxL}]
                $sel delete
                set minLz [{minLz}]
                set maxLz [{maxLz}]
                set selWater [{selWater}]
                set nWater [{nWater}]
                set indicesWater [$selWater get index]
                set output_file [open "nwater_output.txt" w]
                puts $output_file $nWater
                close $output_file
                exit
                """
            with open('water_script.tcl', 'w') as tcl_file:
                tcl_file.write(water_script)
            subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'water_script.tcl'])
            with open('nwater_output.txt', 'r') as output_file:
                nWater = int(output_file.read().strip())
            log.write(f"Number of waters: {nWater}\n")

            # THE RATIO OF WATER TO PROBE MOLECULE (INTEGERS)
            # - 20 is an ideal ratio. It has worked fine in the test cases.
            # - Less than 20 leads to underestimates.
            if lipids == False:
                selPstring = f"resname IPRO and name OH2 and x > [expr $minPx-{padx}] and y > [expr $minPy-{pady}] and x < [expr $maxPx+{padx}] and y < [expr $maxPy+{padz}] and z < [expr $maxPz+{padz}] and z > [expr $minPz-{padz}]"
            else :
                selPstring = f"resname IPRO and name OH2 and x > [expr $minPx-{padx}] and y > [expr $minPy-{pady}] and x < [expr $maxPx+{padx}] and y < [expr $maxPy+{padz}] and ((z < [expr $maxPz+{padz}] and z > $maxLz) or (z < $minLz and z > [expr $minPz-{padz}]))"
            selPstring
            water2probeRatio = 20

            def calc_gcd(p,q):
            #calculate the greater common divisor
                while q != 0:
                    p, q = q, p % q
                return abs(p)
            gcd = 0
            for key, value in probePercent.items():
                for v in value:
                    gcd = calc_gcd(gcd, v)
            modWater = (water2probeRatio * 100) / gcd
            log.write(f"Number of waters must be multiples of {modWater}\n")
            if nWater % modWater == 0:
                howManyMoreWater = 0
            else :
                howManyMoreWater = int(modWater - (nWater % modWater))

            log.write(f"Change in number of waters: {howManyMoreWater}\n")

            if howManyMoreWater:
                pad = 0.1 
                if lipids == False:
                    addWater = f'[atomselect top "water and name OH2 and exwithin {pad} of index $indicesWater"]'
                else :
                    addWater = f'[atomselect top "water and name OH2 and exwithin {pad} of index $indicesWater and (z > $maxLz or z < $minLz)"]'

                if lipids == False:
                    add_water_script = f"""
                    package require psfgen
                    mol load psf {intermediate}.psf
                    mol load pdb {intermediate}.pdb
                    set protein [{protein}]
                    set minmaxP [{minmaxP}]
                    set minP [{minP}]
                    set maxP [{maxP}]
                    set minPx [{minPx}]
                    set minPy [{minPy}]
                    set minPz [{minPz}]
                    set maxPx [{maxPx}]
                    set maxPy [{maxPy}]
                    set maxPz [{maxPz}]
                    set selWater [{selWater}]
                    set nWater {nWater}
                    set indicesWater [$selWater get index]
                    set addWater {addWater}
                    set pad 0.1
                    set padding {boundary_padding}
                    set howManyMoreWater {howManyMoreWater}
                    """
                    add_water_script += """
                    while {[$addWater num] < $howManyMoreWater && $pad < [expr {$padding + 5}] } {
                        $addWater delete
                        set pad [expr {$pad + 0.1}]
                        set addWater [atomselect top "water and name OH2 and exwithin $pad of index $indicesWater"]
                    } 
                    set indicesWater "$indicesWater [lrange [$addWater get index] 0 [expr $howManyMoreWater - 1]]"
                    $addWater delete
                    set numWater [llength $indicesWater]
                    set output_file [open "nwater_output.txt" w]
                    puts $output_file $numWater
                    close $output_file
                    set indiceswater_file [open "indiceswater_ouput.txt" w]
                    puts $indiceswater_file $indicesWater
                    close $indiceswater_file
                    exit
                    """
                else :
                    add_water_script = f"""
                    package require psfgen
                    mol load psf {intermediate}.psf
                    mol load pdb {intermediate}.pdb
                    set protein [{protein}]
                    set minmaxP [{minmaxP}]
                    set minP [{minP}]
                    set maxP [{maxP}]
                    set minPx [{minPx}]
                    set minPy [{minPy}]
                    set minPz [{minPz}]
                    set maxPx [{maxPx}]
                    set maxPy [{maxPy}]
                    set maxPz [{maxPz}]
                    {sel}
                    set minmaxL [{minmaxL}]
                    $sel delete
                    set minLz [{minLz}]
                    set maxLz [{maxLz}]
                    set selWater [{selWater}]
                    set nWater {nWater}
                    set indicesWater [$selWater get index]
                    set addWater {addWater}
                    set pad 0.1
                    set padding {boundary_padding}
                    set howManyMoreWater {howManyMoreWater}
                    """
                    add_water_script += """
                    while {[$addWater num] < $howManyMoreWater && $pad < [expr {$padding + 5}] } {
                        $addWater delete
                        set pad [expr {$pad + 0.1}]
                        set addWater [atomselect top "water and name OH2 and exwithin $pad of index $indicesWater and (z > $maxLz or z < $minLz)"]
                    } 
                    set indicesWater "$indicesWater [lrange [$addWater get index] 0 [expr $howManyMoreWater - 1]]"
                    $addWater delete
                    set numWater [llength $indicesWater]
                    set output_file [open "nwater_output.txt" w]
                    puts $output_file $numWater
                    close $output_file
                    set indiceswater_file [open "indiceswater_ouput.txt" w]
                    puts $indiceswater_file $indicesWater
                    close $indiceswater_file
                    exit
                    """
                with open('add_water_script.tcl', 'w') as tcl_file:
                    tcl_file.write(add_water_script)
                subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'add_water_script.tcl'])
                with open('nwater_output.txt', 'r') as output_file:
                    numWater = int(output_file.read().strip())
                with open('indiceswater_ouput.txt', 'r') as waterindice_file:
                    indicesWater = waterindice_file.read().strip()
                log.write(f"Final number of waters: {numWater}\n")
            
            #Select Probes
            if lipids == False:
                selProbe = f'[atomselect top "{selPstring}"]'
                nProbe = "$selProbe num"

                probe_script = f"""
                mol load psf {intermediate}.psf
                mol load pdb {intermediate}.pdb
                set protein [{protein}]
                set minmaxP [{minmaxP}]
                set minP [{minP}]
                set maxP [{maxP}]
                set minPx [{minPx}]
                set minPy [{minPy}]
                set minPz [{minPz}]
                set maxPx [{maxPx}]
                set maxPy [{maxPy}]
                set maxPz [{maxPz}]
                set selProbe {selProbe}
                set nProbe [{nProbe}]
                set indicesProbe [$selProbe get index]
                set output_file [open "nprobe_output.txt" w]
                puts $output_file $nProbe
                close $output_file
                set indices_file [open "indices_output.txt" w]
                puts $indices_file $indicesProbe
                close $indices_file
                set Probelength [llength $indicesProbe]
                exit
                """
            else :
                selProbe = f'[atomselect top "{selPstring}"]'
                nProbe = "$selProbe num"

                probe_script = f"""
                mol load psf {intermediate}.psf
                mol load pdb {intermediate}.pdb
                set protein [{protein}]
                set minmaxP [{minmaxP}]
                set minP [{minP}]
                set maxP [{maxP}]
                set minPx [{minPx}]
                set minPy [{minPy}]
                set minPz [{minPz}]
                set maxPx [{maxPx}]
                set maxPy [{maxPy}]
                set maxPz [{maxPz}]
                {sel}
                set minmaxL [{minmaxL}]
                $sel delete
                set minLz [{minLz}]
                set maxLz [{maxLz}]
                set selProbe {selProbe}
                set nProbe [{nProbe}]
                set indicesProbe [$selProbe get index]
                set output_file [open "nprobe_output.txt" w]
                puts $output_file $nProbe
                close $output_file
                set indices_file [open "indices_output.txt" w]
                puts $indices_file $indicesProbe
                close $indices_file
                set Probelength [llength $indicesProbe]
                exit
                """

            with open('probe_sel_script.tcl', 'w') as tcl_file:
                tcl_file.write(probe_script)
            subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'probe_sel_script.tcl'])
            with open('indices_output.txt', 'r') as output_file1:
                indicesProbe = output_file1.read().strip()
            with open('nprobe_output.txt', 'r') as output_file:
                nProbe = int(output_file.read().strip())
            log.write(f"Number of probes: {nProbe}\n")

            #Selection of final probe number
            howManyMoreProbe = int((numWater / water2probeRatio) - nProbe)
            log.write(f"Change in number of probes: {howManyMoreProbe}\n")

            if howManyMoreProbe > 0:
                pad = 5.0
                if lipids == False:
                    addProbe = f'[atomselect top "resname IPRO and name OH2 and exwithin {pad} of index {nProbe}"]'
                else :
                    addProbe = f'[atomselect top "resname IPRO and name OH2 and exwithin {pad} of index {nProbe} and (z > $maxLz or z < $minLz)"]'
                
                if lipids == False:
                    add_probe_script = f"""
                    mol load psf {intermediate}.psf
                    mol load pdb {intermediate}.pdb
                    set protein [{protein}]
                    set minmaxP [{minmaxP}]
                    set minP [{minP}]
                    set maxP [{maxP}]
                    set minPx [{minPx}]
                    set minPy [{minPy}]
                    set minPz [{minPz}]
                    set maxPx [{maxPx}]
                    set maxPy [{maxPy}]
                    set maxPz [{maxPz}]
                    set selProbe {selProbe}
                    set nProbe [$selProbe num]
                    set indicesProbe [$selProbe get index]
                    set padding {solvent_padding}
                    set howManyMoreProbe {howManyMoreProbe}
                    set addProbe {addProbe}
                    """
                    add_probe_script += """
                    if {$howManyMoreProbe > 0} {
                        set pad 5.0
                        while {[$addProbe num] < $howManyMoreProbe && $pad < [expr {$padding + 5}] } {
                            $addProbe delete
                            set pad [expr {$pad + 0.25}]
                            set addProbe [atomselect top "resname IPRO and name OH2 and exwithin $pad of index $indicesProbe"]
                        }
                        set indicesProbe "$indicesProbe [lrange [$addProbe get index] 0 [expr $howManyMoreProbe - 1]]"
                        $addProbe delete
                    } elseif {$howManyMoreProbe < 0} {
                        set indicesProbe [lrange $indicesProbe 0 end+$howManyMoreProbe]
                    }
                    set numProbe [llength $indicesProbe]
                    set indicesProbe [lsort $indicesProbe]
                    set output_file [open "nprobe_output.txt" w]
                    puts $output_file $numProbe
                    close $output_file
                    set indices_file [open "indices_files.txt" w]
                    puts $indices_file $indicesProbe
                    close $indices_file
                    exit
                    """
                else :
                    add_probe_script = f"""
                    mol load psf {intermediate}.psf
                    mol load pdb {intermediate}.pdb
                    set protein [{protein}]
                    set minmaxP [{minmaxP}]
                    set minP [{minP}]
                    set maxP [{maxP}]
                    set minPx [{minPx}]
                    set minPy [{minPy}]
                    set minPz [{minPz}]
                    set maxPx [{maxPx}]
                    set maxPy [{maxPy}]
                    set maxPz [{maxPz}]
                    {sel}
                    set minmaxL [{minmaxL}]
                    $sel delete
                    set minLz [{minLz}]
                    set maxLz [{maxLz}]
                    set selProbe {selProbe}
                    set nProbe [$selProbe num]
                    set indicesProbe [$selProbe get index]
                    set padding {boundary_padding}
                    set howManyMoreProbe {howManyMoreProbe}
                    set addProbe {addProbe}
                    """
                    add_probe_script += """
                    if {$howManyMoreProbe > 0} {
                        set pad 5.0
                        while {[$addProbe num] < $howManyMoreProbe && $pad < [expr {$padding + 5}] } {
                            $addProbe delete
                            set pad [expr {$pad + 0.25}]
                            set addProbe [atomselect top "resname IPRO and name OH2 and exwithin $pad of index $indicesProbe"]
                        }
                        set indicesProbe "$indicesProbe [lrange [$addProbe get index] 0 [expr $howManyMoreProbe - 1]]"
                        $addProbe delete
                    } elseif {$howManyMoreProbe < 0} {
                        set indicesProbe [lrange $indicesProbe 0 end+$howManyMoreProbe]
                    }
                    set numProbe [llength $indicesProbe]
                    set indicesProbe [lsort $indicesProbe]
                    set output_file [open "nprobe_output.txt" w]
                    puts $output_file $numProbe
                    close $output_file
                    set indices_file [open "indices_files.txt" w]
                    puts $indices_file $indicesProbe
                    close $indices_file
                    exit
                    """
                with open('add_probe_script.tcl', 'w') as tcl_file:
                    tcl_file.write(add_probe_script)
                subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'add_probe_script.tcl'])
                with open('nprobe_output.txt', 'r') as output_file:
                    numProbe = int(output_file.read().strip())
                with open('indices_files.txt', 'r') as output_file1:
                    indicesProbe = output_file1.read().strip()
                log.write(f"Final number of probes: {numProbe}\n")
                log.write(f"System contains {numWater} water and {numProbe} probe molecules\n")

            #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            # END - ADJUST RELATIVE NUMBER of WATER and PROBE molecules

            # WRITE PDB files for SOLVENT and IONS
            # PSFGEN

            
            #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            # START - PROBE RENUMBER and MUTATE
            # Renumber probe molecule resids starting from 1
            # This is useful when mutating probe molecules
            log.write(f"System contains the following probes: \n")
            probecharge = 0
            tcl_probepercent = "set probePercent [dict create " + " ".join(f"{key} \"{value}\"" for key, value in opts_dict.items()) + "]"

            probe_charge = f"""
            set probeCount [dict create]
            set howmanylist [list]
            set probeidlist [list]
            set aliaslist [list]
            set nProbe {numProbe}
            set PROBEDATA {{{PROBEDATA}}}
            set probecharge {probecharge}
            {tcl_probepercent}
            """
            probe_charge +="""
            set probe_numcharge [open "probe_numcharge.txt" w]
            dict for {key value} $probePercent {
                    set howmanyPROB [::tcl::mathfunc::int [expr $nProbe * $value / 100.0]]
                    dict set probeCount $key $howmanyPROB
                    set holder "$howmanyPROB [dict get $PROBEDATA $key name] molecules ($key)"
                    puts $probe_numcharge $holder
                    lappend probeidlist $key
                    lappend howmanylist $howmanyPROB
                    lappend aliaslist [dict get $PROBEDATA $key alias]
                    set charge [::tcl::mathfunc::int [dict get $PROBEDATA $key charge]]
                    if {$charge} {
                        incr probecharge [expr $howmanyPROB * $charge]
                    }
                }
            close $probe_numcharge
            set output [open "probe_charge.txt" w]
            puts $output $probecharge
            close $output

            set howmanyoutput [open "howmanylist.txt" w]
            puts $howmanyoutput $howmanylist
            close $howmanyoutput

            set probeidoutput [open "probeidlist.txt" w]
            puts $probeidoutput $probeidlist
            close $probeidoutput

            set aliasoutput [open "aliaslist.txt" w]
            puts $aliasoutput $aliaslist
            close $aliasoutput

            exit
            """
            with open('probe_charge.tcl', 'w') as tcl_file:
                tcl_file.write(probe_charge)
            subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'probe_charge.tcl'])
            with open('probe_numcharge.txt', 'r') as probe_numcharge_file:
                probe_amount = probe_numcharge_file.read().strip()
            log.write(f"{probe_amount} \n")
            with open('probe_charge.txt', 'r') as probe_charge_file:
                probecharge = float(probe_charge_file.read().strip())
            with open('howmanylist.txt', 'r') as howmanylist_file:
                howmanylist = howmanylist_file.read().strip()
            with open('probeidlist.txt', 'r') as probeidlist_file:
                probeidlist = probeidlist_file.read().strip()
            with open('aliaslist.txt', 'r') as aliaslist_file:
                aliaslist = aliaslist_file.read().strip()

            totalcharge = (probecharge + float(proteincharge))
            log.write(f"The system has a total charge of {totalcharge} electrons.\n")
            
            n_ions = 0

            if totalcharge > 0:
                n_ions = round(totalcharge)
                log.write(f"Ionization: {n_ions} chloride ions will be added.\n")
            else :
                n_ions = (round(totalcharge) * -1)
                log.write(f"Ionization: {n_ions} sodium ions will be added.\n")

            probe_script = f"""
            package require psfgen
            set intermediate {intermediate}
            set prefix {prefix}
            mol load psf intermediate.psf
            mol load pdb intermediate.pdb
            psfcontext [psfcontext create]
            topology {cgenfftop}
            topology {probetop}
            readpsf {protein_psf}
            set sel [atomselect top "not segid ION \\"WT.*\\""]
            $sel writepdb $intermediate.pdb 
            $sel delete
            coordpdb $intermediate.pdb
            set residProbe 1
            set indicesProbe {{{indicesProbe}}}
            set PROBEDATA [dict create {PROBEDATA}]
            {tcl_probepercent}
            set nProbe {numProbe}
            set probeidlist {{{probeidlist}}}
            set howmanylist {{{howmanylist}}}
            set aliaslist {{{aliaslist}}}
            set totalcharge {totalcharge}
            set n_ions {n_ions}
            set indicesWater {{{indicesWater}}}
            set neutral {{{neutralize}}}
            """
            probe_script += """
            foreach indexProbe $indicesProbe {
                set sel [atomselect top "same residue as index $indexProbe"]
                $sel set resid $residProbe
                $sel set chain " X"
                incr residProbe
                $sel delete
            }
            if {![dict exist $probePercent 'IPRO'] || [dict get $probePercent 'IPRO'] < 100} {
                set residProbe 1
                while {$residProbe <= $nProbe} {
                    set whichProbe [::tcl::mathfunc::int [expr rand() * [llength $probeidlist]]]
                    if {[lindex $howmanylist $whichProbe] > 0} {
                        set sel [atomselect top "chain  X and resid $residProbe"]
                        $sel set resname [lindex $probeidlist $whichProbe]
                        $sel delete
                        foreach {old new} [lindex $aliaslist $whichProbe] {
                        set sel [atomselect top "chain  X and resid $residProbe and name $old"]
                        $sel set name $new
                        $sel delete
                        }
                        incr residProbe
                        lset howmanylist $whichProbe [expr [lindex $howmanylist $whichProbe] -1]
                    }
                }
                set selstr [list]
                dict for {key value} $probePercent {
                    set info [dict get $PROBEDATA $key]
                    dict with info {
                    lappend selstr "(resname $key and name $atomnames)"
                    }
                }
                set selstr [join $selstr " or "]
                set sel [atomselect top "(same residue as index $indicesProbe) and ($selstr)"]
                } else {
                    set sel [atomselect top "same residue as index $indicesProbe"]
                }
            $sel writepdb $intermediate.pdb
            $sel delete
            set residProbe 1
            segment PROB { pdb $intermediate.pdb }
            coordpdb $intermediate.pdb PROB

            if {$totalcharge > 0} {
                set n_ions [expr round($totalcharge)]
                set ion_name "CLA"
                set ion_resname "CLA"
                set ncl $n_ions
                set nna 0
            } else  {
                set n_ions [expr -1 * round($totalcharge)]
                set ion_name "SOD"
                set ion_resname "SOD"
                set ncl 0
                set nna $n_ions
            }

            set sel [atomselect top "segid \\"WT.*\\""]
            set segidWTs [lsort -unique [$sel get segid]]
            $sel delete
            foreach segidWT $segidWTs {
                set sel [atomselect top "segid $segidWT and index $indicesWater"]
                # While at it, renumber water molecule resids starting from 1 for each segment
                set residWater 1
                foreach indexWater [$sel get index] {
                    set sel [atomselect top "same residue as index $indexWater"]
                    $sel set resid $residWater
                    incr residWater
                    $sel delete
                }
                set sel [atomselect top "segid $segidWT and (same residue as index $indicesWater)"]
                $sel writepdb $intermediate.pdb
                segment $segidWT {pdb $intermediate.pdb}
                coordpdb $intermediate.pdb $segidWT
                $sel delete
            }

            guesscoord

            writepsf $intermediate.psf
            writepdb $intermediate.pdb

            if {$neutral}{
                package require autoionize
                autoionize -psf $intermediate.psf -pdb $intermediate.pdb -o $prefix -from 5 -between 5 -ncl $ncl -nna $nna -seg ION
            }
                
            exit
            """
            with open('probe_script.tcl', 'w') as tcl_file:
                tcl_file.write(probe_script)
            subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'probe_script.tcl'])

        log.write(f"Output: Structural data is written into {prefix}.psf file. \n")
        log.write(f"Output: Structural data is written into {prefix}.pdb file. \n")

        if constrain:
            constrain_script = f"""
            mol load psf {prefix}.psf
            mol load pdb {prefix}.pdb
            set prefix {prefix}

            set all [atomselect top "all"]
            $all set beta 0
            $all set occupancy 0
            # protein heavy atoms BETA 1
            $all delete
            set constrain {constrain}
            """
            constrain_script += """
            if {$constrain == "heavy"} {
            set protein [atomselect top "noh and not water and not segid PROB ION \\"WT.*\\""]
            set protein_num [$protein num]
            }
            set output_file [open "protein_num_constrain.txt" w]
            puts $output_file $protein_num
            close $output_file
            $protein set beta 1
            $protein delete
            # alpha carbons OCCUPANCY 1
            set protein [atomselect top "protein and name CA and not segid PROB ION \\"WT.*\\""]
            $protein set occupancy 1
            set geomcent [measure center $protein]
            $protein delete
            set all [atomselect top "all"]
            $all writepdb $prefix.pdb
            $all delete
            set geomcent_output [open "geomcent_output.txt" w]
            puts $geomcent_output $geomcent
            close $geomcent_output
            exit
            """
            with open('constrain_script.tcl', 'w') as tcl_file:
                tcl_file.write(constrain_script)
            subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'constrain_script.tcl'])
            with open('protein_num_constrain.txt', 'r') as output_file:
                protein_num = int(output_file.read().strip())
            with open('geomcent_output.txt', 'r') as geo_output_file:
                geomcent = geo_output_file.read().strip()
            log.write(f"Constraints: {protein_num} heavy atoms are constrained in equilibration\n")

        #XSC File generation
        xsc_script = f"""
        mol load psf {prefix}.psf
        mol load pdb {prefix}.pdb

        set selWater [atomselect top "noh water"]
        set minmaxW [measure minmax $selWater]
        set minW [lindex $minmaxW 0]
        set maxW [lindex $minmaxW 1]
        set minWx [lindex $minW 0]
        set minWy [lindex $minW 1]
        set minWz [lindex $minW 2]
        set maxWx [lindex $maxW 0]
        set maxWy [lindex $maxW 1]
        set maxWz [lindex $maxW 2]
        set geomcent {{{geomcent}}}

        set desired_density 0.62

        set total_mass [vecsum [[atomselect top "all"] get mass]]
        set dimScale [::tcl::mathfunc::pow [expr $total_mass / $desired_density / ($maxWx - $minWx) / ($maxWy - $minWy) / ($maxWz - $minWz)] [expr 1.0 / 3.0]]
        set xLength [expr ($maxWx - $minWx)*$dimScale]
        set yLength [expr ($maxWy - $minWy)*$dimScale]
        set zLength [expr ($maxWz - $minWz)*$dimScale]
            
        set xsc_file [open "{prefix}.xsc" w]
        puts $xsc_file "# NAMD extended system configuration output file"
        puts $xsc_file "#\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z"
        puts $xsc_file "0 $xLength  0  0 0 $yLength  0 0  0 $zLength [join $geomcent]" 
        close $xsc_file

        set atomtotal [[atomselect top all] num]
        set num_atom [open "numatom_ouput.txt" w]
        puts $num_atom $atomtotal
        close $num_atom 

        set mass_total [open "masstotal_output.txt" w]
        puts $mass_total $total_mass
        close $mass_total

        set density [expr $total_mass / $xLength / $yLength / $zLength]
        set dens [open "density_output.txt" w]
        puts $dens $density
        close $dens
        exit
        """

        log.write(f"Output: Extended system coordinates are written into {prefix}.xsc file. \n")

        with open('xsc_script.tcl', 'w') as tcl_file:
            tcl_file.write(xsc_script)
        subprocess.run([f'{vmd}', '-dispdev', 'text', '-e', 'xsc_script.tcl'])
        with open('numatom_ouput.txt', 'r') as output_file:
            numAtom = int(output_file.read().strip())
        log.write(f"Statistics: System contains {numAtom} atoms\n")
        with open('masstotal_output.txt', 'r') as output_file:
            massTotal = float(output_file.read().strip())
        log.write(f"Statistics: Mass of the system is {massTotal} amu\n")
        with open('density_output.txt', 'r') as output_file:
            density = float(output_file.read().strip())
        log.write(f"Statistics: Density of the system is {density} amu/A^3\n")

        parfolder = outdir_location + f'/{prefix}/parameters'
        par_lines = []
        par_lines.append("paraTypeCharmm  on")

        par_files = [] 

        for par_file in parameterfiles:
            par_files = os.path.basename(par_file)
            par_files_location = os.path.join(Druggability_path, par_files)
            par_files_path = os.path.join(parfolder,par_files)

            if not os.path.exists(parfolder):
                os.makedirs(parfolder)

            if not os.path.exists(par_files_path):
                shutil.copy(par_files_location, parfolder)
            par_lines.append(f"parameters      ../parameters/{par_files}")

        par_lines = "\n".join(par_lines)
        log.write("Simulation: Parameter files are copied into /parameter folder.\n")
        log.close()

        if write_conf:
            nonbonded = []
            nonbonded.append("cutoff          10.0")
            nonbonded.append("switching       on")
            nonbonded.append("switchdist      8.0")  
            nonbonded.append("pairlistdist    12.0")
            nonbonded.append("margin          1.0")
            nonbonded.append("exclude         scaled1-4")
            nonbonded = "\n".join(nonbonded)

            minfix = "min"
            min = open(f"min.conf",'w')
            min.write(f"coordinates     ./{prefix}.pdb\n")
            min.write(f"structure       ./{prefix}.psf\n")
            min.write(f"{par_lines}\n")
            min.write(f"outputname      {minfix}\n")
            min.write(f"binaryoutput    no\n")
            min.write(f"restartname     {minfix}\n")
            min.write(f"restartfreq     10000\n")
            min.write(f"timestep        1.0\n")
            min.write(f"{nonbonded}\n")
            min.write(f"temperature     0\n")
            min.write(f"constraints     on\n")
            min.write(f"consref         ./{prefix}.pdb\n")
            min.write(f"conskfile       ./{prefix}.pdb\n")
            min.write(f"conskcol        B\n")
            min.write(f"constraintScaling  1.0\n")
            min.write(f"PME             yes\n")
            min.write(f"PMEGridSpacing  1.0\n")
            min.write(f"extendedSystem  ./{prefix}.xsc\n")
            min.write(f"wrapWater       on\n")
            min.write(f"wrapAll         on\n")
            if lipids == True:
                min.write(f"minimize        20000\n")
            else :
                min.write(f"minimize        10000\n")
            min.close()

            eq1 = open(f"eq1.conf",'w')
            eq1.write(f"coordinates     ./min.coor\n")
            eq1.write(f"structure       ./{prefix}.psf\n")
            eq1.write(f"{par_lines}\n")
            eq1.write(f"outputname      eq1\n")
            eq1.write(f"binaryoutput    no\n")
            eq1.write(f"restartname     eq1\n")
            eq1.write(f"restartfreq     2000\n")
            eq1.write(f"binaryrestart   no\n")
            eq1.write(f"DCDfreq         2000\n")
            eq1.write(f"DCDfile         eq1.dcd\n")
            eq1.write(f"outputEnergies  2000\n")
            if lipids == True:
                eq1.write(f"timestep        1.0\n")
            else :
                eq1.write(f"timestep        2.0\n")
            eq1.write(f"fullElectFrequency 2\n")
            eq1.write(f"nonbondedFreq   1\n")
            eq1.write(f"rigidBonds      all\n")
            eq1.write(f"{nonbonded}\n")
            eq1.write(f"temperature     0\n")
            eq1.write(f"constraints     on\n")
            eq1.write(f"consref         ./min.coor\n")
            eq1.write(f"conskfile       ./{prefix}.pdb\n")
            eq1.write(f"conskcol        B\n")
            eq1.write(f"constraintScaling  0.5\n")
            eq1.write(f"PME             yes\n")
            eq1.write(f"PMEGridSpacing  1.0\n")
            eq1.write(f"langevin        on\n")
            eq1.write(f"langevinTemp    100\n")
            eq1.write(f"langevinDamping 5\n")
            eq1.write(f"langevinHydrogen off\n")
            eq1.write(f"useGroupPressure yes\n")
            eq1.write(f"useFlexibleCell  no\n")
            eq1.write(f"useConstantArea  no\n")
            eq1.write(f"useConstantRatio no\n")
            eq1.write(f"langevinPiston   on\n")
            eq1.write(f"langevinPistonTarget  1.01325\n")
            eq1.write(f"langevinPistonPeriod  100.0\n")
            eq1.write(f"langevinPistonDecay   50.0\n")
            eq1.write(f"langevinPistonTemp    100\n")
            eq1.write(f"extendedSystem    min.xsc\n")
            eq1.write(f"wrapWater       on\n")
            eq1.write(f"wrapAll         on\n")
            eq1.write(f"reinitvels      100\n")
            eq1.write("for {set T 100} {$T < 300} {incr T 10} {\n")
            eq1.write("    langevinTemp      $T;\n")
            eq1.write("    run              1000;\n")
            eq1.write("}\n")
            eq1.write(f"langevinTemp    300\n")
            eq1.write(f"run           40000;")
            eq1.close()

            if Probes == True:
                eq2 = open(f"eq2.conf",'w')
                eq2.write(f"coordinates     eq1.coor\n")
                eq2.write(f"structure       ./{prefix}.psf\n")
                eq2.write(f"{par_lines}\n")
                eq2.write(f"outputname      eq2\n")
                eq2.write(f"binaryoutput    no\n")
                eq2.write(f"restartname     eq2\n")
                eq2.write(f"restartfreq     2000\n")
                eq2.write(f"binaryrestart   no\n")
                eq2.write(f"DCDfreq         2000\n")
                eq2.write(f"DCDfile         eq2.dcd\n")
                eq2.write(f"outputEnergies  2000\n")
                if lipids == True:
                    eq2.write(f"timestep        1.0\n")
                else :
                    eq2.write(f"timestep        2.0\n")
                eq2.write(f"fullElectFrequency 2\n")
                eq2.write(f"nonbondedFreq   1\n")
                eq2.write(f"rigidBonds      all\n")
                eq2.write(f"{nonbonded}\n")
                eq2.write(f"velocities      eq1.vel\n")
                eq2.write(f"constraints     on\n")
                eq2.write(f"consref         ./min.coor\n")
                eq2.write(f"conskfile       ./{prefix}.pdb\n")
                eq2.write(f"conskcol        B\n")
                eq2.write(f"constraintScaling  1.0\n")
                eq2.write(f"PME             yes\n")
                eq2.write(f"PMEGridSpacing  1.0\n")
                eq2.write(f"langevin        on\n")
                eq2.write(f"langevinDamping 5\n")
                eq2.write(f"langevinHydrogen off\n")
                eq2.write(f"extendedSystem  eq1.xsc\n")
                eq2.write(f"wrapWater       on\n")
                eq2.write(f"wrapAll         on\n")
                eq2.write("for {set T 300} {$T < 600} {incr T  10} {\n")
                eq2.write("    langevinTemp      $T;\n")
                eq2.write("    run              1000;\n")
                eq2.write("}\n")
                eq2.write(f"langevinTemp    600\n")
                eq2.write(f"run           300000;")
                eq2.write("for {set T 570} {$T >= 300} {incr T -30} {\n")
                eq2.write(f"    langevinTemp     $T;\n")
                eq2.write(f"	   run             1000;\n")
                eq2.write("}\n")
                eq2.close()

                eq3 = open(f"eq3.conf",'w')
                eq3.write(f"coordinates     eq2.coor\n")
                eq3.write(f"structure       ./{prefix}.psf\n")
                eq3.write(f"{par_lines}\n")
                eq3.write(f"outputname      eq3\n")
                eq3.write(f"binaryoutput    no\n")
                eq3.write(f"restartname     eq3\n")
                eq3.write(f"restartfreq     2000\n")
                eq3.write(f"binaryrestart   no\n")
                eq3.write(f"DCDfreq         2000\n")
                eq3.write(f"DCDfile         eq3.dcd\n")
                eq3.write(f"outputEnergies  2000\n")
                if lipids == True:
                    eq3.write(f"timestep        1.0\n")
                else :
                    eq3.write(f"timestep        2.0\n")
                eq3.write(f"fullElectFrequency 2\n")
                eq3.write(f"nonbondedFreq   1\n")
                eq3.write(f"rigidBonds      all\n")
                eq3.write(f"{nonbonded}\n")
                eq3.write(f"velocities      eq2.vel\n")
                eq3.write(f"PME             yes\n")
                eq3.write(f"PMEGridSpacing  1.0\n")
                eq3.write(f"langevin        on\n")
                eq3.write(f"langevinTemp    300\n")
                eq3.write(f"langevinDamping 5\n")
                eq3.write(f"langevinHydrogen off\n")
                eq3.write(f"useGroupPressure        yes\n")
                eq3.write(f"useFlexibleCell         no\n")
                eq3.write(f"useConstantArea         no\n")
                eq3.write(f"useConstantRatio        no\n")
                eq3.write(f"langevinPiston          on\n")
                eq3.write(f"langevinPistonTarget    1.01325\n")  
                eq3.write(f"langevinPistonPeriod    100.0\n")     
                eq3.write(f"langevinPistonDecay     50.0\n")  
                eq3.write(f"langevinPistonTemp      300.0\n")                                               
                eq3.write(f"extendedSystem  eq2.xsc\n")
                eq3.write(f"wrapWater       on\n")
                eq3.write(f"wrapAll         on\n")
                eq3.write(f"run                  300000")
                eq3.close()
                
                sim_step = int(sim_length) * 500000

                sim = open(f"sim.conf",'w')
                sim.write(f"coordinates     eq3.coor\n")
                sim.write(f"structure       ./{prefix}.psf\n")
                sim.write(f"{par_lines}\n")
                sim.write(f"outputname      sim\n")
                sim.write(f"binaryoutput    no\n")
                sim.write(f"restartname     sim\n")
                sim.write(f"restartfreq     2000\n")
                sim.write(f"binaryrestart   no\n")
                sim.write(f"DCDfreq         2000\n")
                sim.write(f"DCDfile         sim.dcd\n")
                sim.write(f"outputEnergies  2000\n")
                if lipids == True:
                    sim.write(f"timestep        2.0\n")
                else :
                    sim.write(f"timestep        2.0\n")
                sim.write(f"fullElectFrequency 2\n")
                sim.write(f"nonbondedFreq   1\n")
                sim.write(f"rigidBonds      all\n")
                sim.write(f"{nonbonded}\n")
                sim.write(f"velocities      eq3.vel\n")
                sim.write(f"PME             yes\n")
                sim.write(f"PMEGridSpacing  1.0\n")
                sim.write(f"langevin        on\n")
                sim.write(f"langevinTemp    300\n")
                sim.write(f"langevinDamping 5\n")
                sim.write(f"langevinHydrogen off\n")
                sim.write(f"useGroupPressure        yes\n")
                sim.write(f"useFlexibleCell         no\n")
                sim.write(f"useConstantArea         no\n")
                sim.write(f"useConstantRatio        no\n")
                sim.write(f"langevinPiston          on\n")
                sim.write(f"langevinPistonTarget    1.01325\n")  
                sim.write(f"langevinPistonPeriod    100.0\n")     
                sim.write(f"langevinPistonDecay     50.0\n")  
                sim.write(f"langevinPistonTemp      300.0\n")                                               
                sim.write(f"extendedSystem  eq3.xsc\n")
                sim.write(f"wrapWater       on\n")
                sim.write(f"wrapAll         on\n")
                sim.write(f"run                  {sim_step}")
                sim.close()

            sh_file = open(f"{prefix}.sh",'w')
            sh_file.write("namd2 +p8 min.conf > min.log\n")
            sh_file.write("namd2 +p8 eq1.conf > eq1.log\n")
            if Probes == True:
                sh_file.write("namd2 +p8 eq2.conf > eq2.log\n")
                sh_file.write("namd2 +p8 eq3.conf > eq3.log\n")
            sh_file.write("namd2 +p8 sim.conf > sim.log\n")
            sh_file.close()

        if Probes == True:
            conf_files = ["min.conf", "eq1.conf", "eq2.conf", "eq3.conf", "sim.conf",f"{prefix}.sh", f"{prefix}.psf", f"{prefix}.pdb",f"{prefix}.xsc", f"{prefix}.log"]
        else :
            conf_files = ["min.conf", "eq1.conf", "sim.conf",f"{prefix}.sh", f"{prefix}.psf", f"{prefix}.pdb",f"{prefix}.xsc", f"{prefix}.log"]
            
        ztrj_files = ["0chlcmpall.tcl", "1a.sh"]

        output_location = os.path.join(outdir_location, prefix)
        con_folder = os.path.join(output_location, 'simulation_run')
        total_num_sims = n_sims + 1
        sim_num = 1

        while sim_num < total_num_sims:
            final_folder = f"{con_folder}{sim_num}"

            os.makedirs(final_folder, exist_ok = True)

            for files in conf_files:
                shutil.copy(files, final_folder)
            
            ztrj_folder = f'{final_folder}/ztrj'
            os.makedirs(ztrj_folder, exist_ok = True)
            
            for files in ztrj_files:
                shutil.copy(files, ztrj_folder)  

            sim_num += 1

        remove_files = ['*.pdb', '*.psf', '*.log', '*.xsc', '*.tcl', '*.txt', '*.conf', '*.sh']
        important_files = ['probe.pdb', 'probe.psf', '0chlcmpall.tcl', '1a.sh']
        for items in remove_files:
            files = glob.glob(items)
            for file in files:
                if file not in important_files:
                    shutil.os.remove(file)

def drugui_analysis(pdb, psf, dcds, **kwargs):
    """Analyzes druggability simulations to determine potential druggable sites of a protein.
    
    prefix: The name the prepared files will be given
    outdir_location: Files will be written in specified directory
    selection: Which part of the protein will be align during analysis
    grid_spacing: The size of a grid element, along X, Y, and Z dimensions
    contact_distance: Distance between probe atom and protein that will be used for grid calculation
    align: How molecules will be wrapped in druggability calculations
    temperature: Temperature of the system in the productive simulation
    merge_radius: Probe merge radius in angstroms
    n_frames: Number of frames used in determining the grid data
    n_probes: Number of probe binding hotspots to merge to make a drug-size solution
    delta_g: Probe binding free energy to determine binding hotspots
    min_n_probes: Minimum number of hotspots in an acceptable drug-size solution
    low_affinity: Lowest affinity to report a solution in micromolar units
    max_charge: Maximum absolute charge to accept solutions
    n_solutions: Number of solutions to report in each distinct potential binding site
    n_charged: Maximum number of charged hotspots in a solution
    probes: Selected probes used for druggability simulation
    """
    protein_pdb = pdb
    protein_psf = psf
    dcds_raw = dcds
    prefix = kwargs.pop('prefix', 'dg')
    outdir_location = kwargs.pop('outdir_location', "")
    selection = kwargs.pop('selection', "noh and protein")
    grid_spacing = kwargs.pop('grid_spacing', 0.5)
    contact_distance = kwargs.pop('contact_distance', 4.0)
    align = kwargs.pop('align', 'calpha')
    temperature = kwargs.pop('temperature', 300.)
    merge_radius = kwargs.pop('merge_radius', 5.6)
    n_frames = kwargs.pop('n_frames', 1)
    n_probes = kwargs.pop('n_probes', 7)
    delta_g = kwargs.pop('delta_g', -1.0)
    min_n_probes = kwargs.pop('min_n_probes', 6)
    low_affinity = kwargs.pop('low_affinity', 10)
    max_charge = kwargs.pop('max_charge', 2)
    n_solutions = kwargs.pop('n_solutions', 3)
    n_charged = kwargs.pop('n_charged', 3)
    probes = kwargs.pop('probes', ['IPAM', 'IPRO', 'ACTT', 'IBUT', 'ACAM', 'IMID', 'BENZ'])

    verbose = 'info'

    def buildGrids(prefix, pdb_psf_dcds, probes, align = align, protein = selection, contacti = contact_distance, resolution = grid_spacing, savedcd=False):

        if len(pdb_psf_dcds) > 1:
            reference = parsePDB(pdb_psf_dcds[0][0]).select(align).copy()
        else:
            reference = None
        
        UNITCELL = []
        DCDOUT = {}
        nframe = 0

        for p in probes:
            DCDOUT[p] = DCDFile(prefix + '_' + p + '.dcd', 'w')

        from prody import startLogfile
        startLogfile(prefix + '_grid.log')

        pdb, psf = pdb_psf_dcds[0][:2]
        dcds = pdb_psf_dcds[0][2:]

        pdb = parsePDB(pdb, AtomGroup=parsePSF(psf), long_resname = True)
        palign = pdb.select(align)
        if reference is not None:
            matchAlign(palign, reference)
        pcenter = calcCenter(palign)

        # make sure all probe names select some residues
        from prody import plog
        probe_selstr = 'noh and resname'
        for p in probes:
            sel = pdb.select('noh and resname ' + p)
            if sel is None:
                continue
                #raise ValueError('probe ' + p + ' is not found in the system')
            hv = sel.getHierView()
            n = hv.numResidues()
            res = next(hv.iterResidues())
            writePDB(prefix + '_' + p + '.pdb', res, long_resname = True)
            plog(str(n) + ' copies of ' + p + ' is found.')
            probe_selstr += ' ' + p

        # start trajectories for reading
        plog('Opening DCD files for reading.')
        dcd = Trajectory(dcds[0]) 
        for fn in dcds[1:]:
            dcd.addFile(fn)
        plog(str(len(dcd)) + ' frames from ' + str(len(dcds)) +  
            ' file(s) will be evaluated.')
        
        dcd.link(pdb)
        # make alignment selection, calling `frame.superpose` will align frames 
        # based on this selection
        dcd.setAtoms(palign)

        # make a probe selection
        PRBSEL = pdb.select('noh and ' + probe_selstr)
        PRBIDX = PRBSEL.getIndices()
        # make a copy of probe selection, it will be faster to make contact search
        # in this copy after copying the coordinates
        PROBES = PRBSEL.copy().toAtomGroup()

        pcontact = pdb.select(protein)
        writePDB(prefix + '_protein_heavyatoms.pdb', pcontact)

        from prody import LOGGER
        LOGGER.progress('Evaluating frames:', len(dcd))
        if savedcd:
            dcdout = DCDFile(pdb.getTitle() + '_aligned_wrapped.dcd', 'w')

        for i, frame in enumerate(dcd):
            #if i % 20 != 0:
            #    continue
            # first center the frame
            moveAtoms(palign, to=pcenter, ag=True)

            # wrap probes that are out of the box
            unitcell = frame.getUnitcell()[:3]
            UNITCELL.append(unitcell)
            coords = pdb._coords[0]
            coords[PRBIDX] = wrapAtoms(coords[PRBIDX], unitcell, pcenter)  
            
            # superpose frame coordinates onto selected atoms
            frame.superpose()
            
            if savedcd:
                dcdout.write(coords)
            
            # identify probes in contact and write them in DCD file
            PROBES.setCoords(PRBSEL.getCoords())
            cont = PROBES.select(f'same residue as within {contacti} of pcontact', 
                                 pcontact=pcontact)
            
            if cont:
                for res in cont.getHierView().iterResidues():
                    DCDOUT[res.getResname()].write(res._getCoords()) 
            nframe += 1
            LOGGER.update(i)
        
        dcd.close()

        for p in probes: DCDOUT.pop(p).close()
    
        UNITCELL = array(UNITCELL).max(0)

        offset = pcenter - UNITCELL / 2.
        length = UNITCELL
        n_bins = ceil(length / resolution).astype(int)
        bins = [arange(n) * resolution + offset[i] for i, n in enumerate(n_bins)]
        probe_grids = []
        for p in probes:    
            fn = prefix + '_' + p
            if not os.path.getsize(fn + '.dcd'):
                plog('No ' + p + ' found in contact with the protein.')
                continue
            e = parseDCD(fn + '.dcd')
            c = calcCenter(e._getCoordsets())
            garray = histogramdd(c, bins)
            grid = OpenDX()
            grid.filename = fn + '.dx'
            grid.name = fn
            grid.spacing = array([resolution, resolution, resolution])
            grid._origin = offset
            grid.offset = offset
            grid.array = garray[0] / nframe
            grid.shape = grid.array.shape
            grid._comments = ['# comment', 'object']
            grid.write(fn + '.dx')
            probe_grids.append((p, fn + '.dx'))

        from prody import closeLogfile
        closeLogfile(prefix + '_grid.log')
        return probe_grids

    def calcDruggability(prefix, probe_grids, **kwargs):



        dia = druggability.DIA(prefix, workdir=outdir_location, verbose=verbose)
        # check parameters for their correctness
        dia.set_parameters(temperature=kwargs.get('temperature', temperature)) # K (productive simulation temperature)
        dia.set_parameters(delta_g=kwargs.get('delta_g', delta_g)) # kcal/mol (probe binding hotspots with lower values will be evaluated)
        dia.set_parameters(n_probes=kwargs.get('n_probes', n_probes)) # (number of probes to be merged to determine achievable affinity of a potential site)
        dia.set_parameters(min_n_probes=kwargs.get('min_n_probes', min_n_probes)) # (minimum number of probes to be merged for an acceptable soltuion)
        dia.set_parameters(merge_radius=kwargs.get('merge_radius', merge_radius)) # A (distance within which two probes will be merged)
        dia.set_parameters(low_affinity=kwargs.get('low_affinity', low_affinity)) # microMolar (potential sites with affinity better than this value will be reported)
        dia.set_parameters(n_solutions=kwargs.get('n_solutions', n_solutions)) # (number of drug-size solutions to report for each potential binding site)
        dia.set_parameters(max_charge=kwargs.get('max_charge', max_charge)) # (maximum absolute total charge, where total charge is occupancy weighted sum of charges on probes)
        dia.set_parameters(n_charged=kwargs.get('n_charged', n_charged)) # (maximum number of charged hotspots in a solution)
        dia.set_parameters(n_frames=kwargs.get('n_frames', n_frames)) # number of frames (if volmap was used (.dx), 1 is correct)
    
        for probe_type, grid_file in probe_grids: dia.add_probe(probe_type, grid_file)

            # Do the calculations
        dia.perform_analysis()
        dia.pickle()
            # Evaluate a ligand. Be sure that the ligand bound structure is superimposed
            # onto PROTEIN_heavyatoms.pdb
        ligand = kwargs.get('ligand', None)
        if ligand:
            dia.evaluate_ligand(ligand)

    def evalLigandSite(prefix, ligand, radius=1.5, delta_g=-0.5):

        dia = pickler(os.path.join(prefix, prefix + '.dso.gz'))
        dia.evaluate_ligand(ligand, radius=radius, delta_g=delta_g)

    
    # output file and folder names will start with the following
    prefix_n = prefix
    os.chdir(outdir_location)
    
    try:
        dcds_list = ast.literal_eval(dcds_raw)
        if not isinstance(dcds_list, (list, tuple)):
            raise ValueError("Parsed object is not a list or tuple")
    except Exception:
        dcds_list = dcds_raw.split()

    probes = list(probes)
    pdb_psf_dcds = [[
    protein_pdb,
    protein_psf,
    *dcds_list
    ]]

    # build grids
    buildGrids(prefix_n, pdb_psf_dcds, probes, savedcd= True)

    # DRUGGABILITY
    parameters = {}
    parameters['temperature'] = temperature    # K (productive simulation temperature)
    parameters['delta_g'] = delta_g        # kcal/mol (probe binding hotspots with lower values will be evaluated)
    parameters['n_probes'] = n_probes         # (number of probes to be merged to determine achievable affinity of a potential site)
    parameters['min_n_probes'] = min_n_probes    # (minimum number of probes to be merged for an acceptable soltuion)
    parameters['merge_radius'] = merge_radius    # A (distance within which two probes will be merged)
    parameters['low_affinity'] = low_affinity     # microMolar (potential sites with affinity better than this value will be reported)
    parameters['n_solutions'] = n_solutions       # (number of drug-size solutions to report for each potential binding site)
    parameters['max_charge'] = max_charge        # (maximum absolute total charge, where total charge is occupancy weighted sum of charges on probes)
    parameters['n_charged'] = n_charged         # (maximum number of charged hotspots in a solution)
    parameters['n_frames'] = n_frames          # number of frames (if volmap was used (.dx), 1 is correct)
    
    # Probe grid files are automatically determined based on prefix
    from glob import glob
    probe_grids = [(os.path.splitext(fn)[0].split('_')[-1], fn) 
                    for fn in glob(prefix_n + '_*.dx')]  
                     
    calcDruggability(prefix_n, probe_grids, **parameters)
    
    # LIGAND SITE
    # Evaluate a ligand. Be sure that the ligand bound structure is superimposed
    # onto PROTEIN_heavyatoms.pdb
    #evalLigandSite(prefix, 'ligand.pdb', radius=1.5, delta_g=-0.5)

PROBETOPPAR = {
        "PBDA": "probe2.top probe.prm",
        "CGenff": "top_all36_cgenff.rtf par_all36_cgenff.prm"
    }

PROBETYPES = {
        "core": "Core probes",
        "polar": "Polar probes",
        "hydrophobe": "Hydrophobes",
        "negative": "Negatively charged",
        "positive": "Positively charged",
        "ring5": "5-membered rings",
        "ring6": "6-membered rings"
    }
