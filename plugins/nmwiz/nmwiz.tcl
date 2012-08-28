# NMWiz: Normal Mode Visualization, Animation, and Plotting
# 
# University of Illinois Open Source License
# Copyright 2010-2011 Ahmet Bakan
# All rights reserved.
# 
# Developed by:		
#       Ahmet Bakan
# 			http://www.pitt.edu/~ahb12/
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the Software), to deal with 
# the Software without restriction, including without limitation the rights to 
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software is furnished to 
# do so, subject to the following conditions:
# 
# Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimers.
# 
# Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimers in the documentation 
# and/or other materials provided with the distribution.
# 
# Neither the names of Theoretical and Computational Biophysics Group, 
# University of Illinois at Urbana-Champaign, nor the names of its contributors 
# may be used to endorse or promote products derived from this Software without 
# specific prior written permission.
# 
# THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
# THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
# OTHER DEALINGS WITH THE SOFTWARE.

package require exectool
package require multiplot
package provide nmwiz 1.2


proc sign x {expr {($x>0) - ($x<0)}}

# List of NMWiz functions
#   ::NMWiz::initGUI                 make main window
#   ::NMWiz::loadNMD                 load NMD file

#   ::nmguiX::deleteMolecules        delete molecules storing coordinates/graphics/animations/selection
#   ::nmguiX::prepareSelmol          prepare the molecule where selections are displayed
#   ::nmguiX::getPDBLines            return PDB files for the molecule
#   ::nmguiX::clearSelection         turn selection representation of and clear labels


namespace eval ::NMWiz:: {
  namespace export nmwizgui
  namespace export initialize

  variable guicount -1
  variable tmpdir
  variable titles [list]
  variable preserview 1
  variable plothandles [list] 
  #variable namespaces [list]
  #variable nmwizguis [list]
  variable platform $tcl_platform(platform) 
  switch $platform {
    unix {
      set tmpdir "/tmp" ;  # or even $::env(TMPDIR), at times.
    } macintosh {
      set tmpdir $::env(TRASH_FOLDER)  ;# a better place?
    } default {
      set tmpdir [pwd]
      catch {set tmpdir $::env(TMP)}
      catch {set tmpdir $::env(TEMP)}
    }
  }
  
  proc showHelp {context} {
    set windowname nmwizhelp
    if {[winfo exists .$windowname] == 0} {
      set log [toplevel ".$windowname"]
      wm title $log "NMWiz Help"
      wm resizable $log 1 1
      incr logcount

      text $log.text -bg White -bd 2 -font Courier \
        -yscrollcommand ".$windowname.vscr set"
      scrollbar $log.vscr -command ".$windowname.text yview"
      pack $log.text -side left -fill both -expand 1  
      pack $log.vscr -side right -fill y
    } else {
      set log .$windowname
    }
    $log.text configure -state normal -wrap word
    $log.text delete 1.0 end
      
    if {$context == "wizard"} {
      $log.text insert end "ProDy Interface\n"
      $log.text insert end "===============\n\n"
      $log.text insert end "\nActive Mode\n"
      $log.text insert end "--------------\n\n"
      $log.text insert end "Select the active mode for which you want to draw arrows or make an animation. "
      $log.text insert end "Direction of arrows depicting the normal mode can be changed using +/- button. "
      $log.text insert end "Arrows can be drawn along both directions by changing the options Mode Graphics Options panel. "
      $log.text insert end "The selected color effects both arrow graphics and square fluctuation plots."
      $log.text insert end "\n\n**RMSD**\n\n"
      $log.text insert end "The RMSD corresponding to the displacement described by the arrows is displayed. User can change the RMSD value to rescale the arrows. "
      $log.text insert end "The scaling factor that produces the specified RMSD is printed to the VMD console (along with the magnitude of the mode provided in NMD file). "
      $log.text insert end "\n\n**Selection**\n\n"
      $log.text insert end "Selection entry allows the user to display arrows for a subset of atoms.\n\n"
      $log.text insert end "*TIP*: If the arrow graphics are too crowded or the display is slow, draw arrows for an evenly spaced subset of residues, e.g try 'name CA and residue % 4 == 0', which will draw an arrow for every fourth residue."
      $log.text insert end "\n\n\n"
      $log.text insert end "Mode Graphics\n"
      $log.text insert end "--------------\n\n"
      $log.text insert end "Id of the molecule that contains the arrow graphics of the active mode is shown in parentheses.\n\n"
      $log.text insert end "Buttons:\n\n"
      $log.text insert end " * Draw: draw/redraw arrow graphics for the active mode\n"
      $log.text insert end " * Clean: remove most recently drawn arrow graphics\n"
      $log.text insert end " * Hide/Show: hide/show most recently drawn arrow graphics\n"
      $log.text insert end " * Options: show/hide arrow graphics option panel\n"
      $log.text insert end "\nOptions:\n\n"
      $log.text insert end "User can change arrow graphics properties and how NMWiz behave upon such changes in this panel.\n"
      $log.text insert end "\nBy default:\n\n"
      $log.text insert end " * arrow graphics are set to automatically change when graphics properties are changed by the user\n"
      $log.text insert end " * current graphics are hidden the active mode is changed\n"
      $log.text insert end "\nOptionally:\n\n"
      $log.text insert end " * arrows can be drawn in both directions to look like a double headed arrow\n"
      $log.text insert end " * arrows shorter than a length (A) can be hidden\n"
      $log.text insert end "\nAdditionally, user can change:\n\n"
      $log.text insert end " * width of the arrow cylinder\n"
      $log.text insert end " * width/height of the arrow head code\n"
      $log.text insert end " * graphics material and resolution"
      $log.text insert end "\n\n\n"
      $log.text insert end "Animations\n"
      $log.text insert end "----------\n\n"
      $log.text insert end "Id of the molecule that contains the most recently generated animation is shown in parentheses.\n\n"
      $log.text insert end "Buttons:\n\n"
      $log.text insert end " * Draw: animate fluctuations along the active mode\n"
      $log.text insert end " * Play : play/pause the animation\n"
      $log.text insert end " * Hide : hide/show the animation\n"
      $log.text insert end " * Options: show/hide animation option panel\n"
      $log.text insert end "\nOptions:\n\n"
      $log.text insert end "User can elect automatic generation and continuous play of animations when the active mode changes. User can also select the number of frames in the animation."
      $log.text insert end "\n\n\n"
      $log.text insert end "Plotting\n"
      $log.text insert end "--------\n\n"
      $log.text insert end "Id of the molecule for displaying selected residues is shown in parentheses.\n\n"
      $log.text insert end "Buttons:\n\n"
      $log.text insert end " * Plot: plot squared-fluctuations along the active mode\n"
      $log.text insert end " * Clear: clear all selections and selected atom labels\n"
      $log.text insert end " * Hide/Show: hide/show the selected residues\n"
      $log.text insert end " * Options: change plotting options\n"
      $log.text insert end "\n\n\n"
      $log.text insert end "Molecule Representations\n"
      $log.text insert end "------------------------\n\n"
      $log.text insert end "Id of the molecule that contains the structure is shown in parentheses.\n\n"
      $log.text insert end "Buttons:\n\n"
      $log.text insert end " * Update: update molecule representation\n"
      $log.text insert end " * Focus: reset view to focus on the structure\n"
      $log.text insert end " * Hide/Show: hide/show strudture\n"
      $log.text insert end " * Options: change molecular system representation\n"
      $log.text insert end "\nOptions:\n\n"
      $log.text insert end "User can select the representation and coloring scheme. User can change the molecule representation settings manually, by setting 'Show structure as' to 'Custom'.\n\n"      
      $log.text insert end "Structure can be colored based on the `Mobility` of the residues in the active mode, based on 'Bfactors' that came in NMD file, or based on residue/atom 'Index'.\n\n"
      $log.text insert end "In addition to the standard representations (e.g. Tube/Trace/Licorice), structure can be represented as an elastic network. Color scale method and midpoint can be used to adjust mobility and Bfactors based coloring."
      $log.text insert end "User can set the cutoff distance, width of dynamic bonds, and node spheres. Note that changing the cutoff distance distance only affects representation, not the precalculated normal mode data.\n\n"
      $log.text insert end "*TIP*: When visualizing a large system, display molecule at lower resolutions and/or try displaying fewer atoms if all atoms are displayed."
    } elseif {$context == "prody"} {
      $log.text insert end "ProDy Interface\n"
      $log.text insert end "===============\n\n"
      $log.text insert end "ProDy interface allows users to perform the following calculations for molecules loaded in VMD:\n\n"
      $log.text insert end "* Anisotropic Network Model (ANM)\n"
      $log.text insert end "* Gaussian Network Model (GNM)\n"
      $log.text insert end "* Principal Component Analysis (PCA) a.k.a. Essential Dynamics Analysis (EDA)\n\n\n"
      $log.text insert end "Atom Selection\n"
      $log.text insert end "--------------\n\n"
      $log.text insert end "First thing you need to do is selecting the molecule and specifying the atoms that you want to include in the calculations. "
      $log.text insert end "If you do not see all molecules in the menu, click 'Update'."
      $log.text insert end "\n\n\n"
      $log.text insert end "ProDy Job Settings\n"
      $log.text insert end "------------------\n\n"
      $log.text insert end "Specify the calculation type and output options in this panel. "
      $log.text insert end "Coordinate data for selected atoms and the NMD data after calculations will be written into the 'Output directory'. "
      $log.text insert end "All output files will named after 'Output filename'.\n\n"
      $log.text insert end "**ProDy Scripts**\n\n"
      $log.text insert end "NMWiz will try to find the path to Python executable ('python' or 'python.exe') and ProDy script ('prody'). For this to work, both of these files must be in you PATH environment variable. If they are not found, you will be prompted to specify the path."
      $log.text insert end "\n\n\n"
      $log.text insert end "ANM/GNM Settings\n"
      $log.text insert end "----------------\n\n"
      $log.text insert end "Specify the following:\n\n"
      $log.text insert end " * number of modes to be calculated\n"
      $log.text insert end " * index of the frame (coordinate set) to be used in calculations\n" 
      $log.text insert end " * cutoff distance\n" 
      $log.text insert end " * force constant"
      $log.text insert end "\n\n\n"
      $log.text insert end "PCA/EDA Settings\n"
      $log.text insert end "----------------\n\n"
      $log.text insert end "Note that for PCA/EDA calculations molecule must have multiple frames. Specify the range of frames to be used in calculations. For large systems, prefer to write coordinates in DCD format to gain IO speed and save save disk space."
    } elseif {$context == "compare"} {
      $log.text insert end "Structure Comparison\n"
      $log.text insert end "====================\n\n"
      $log.text insert end "This interface allows the user to compare two molecules (or two frames of the same molecule) loaded in VMD. It can be used to align structures and draw deformation vector."
      $log.text insert end "Follow these steps to compare structures:\n\n"
      $log.text insert end "1) Load molecules into VMD\n\n"
      $log.text insert end "Select the molecules that you want to compare. If you don't see the molecule, click 'Update' button. "
      $log.text insert end "\n\n"
      $log.text insert end "2) Select atoms and specify frames\n\n"
      $log.text insert end "Make atom selections for each molecule (or frame). Note that the number of selected atoms must be the same. "
      $log.text insert end "\n\n"
      $log.text insert end "3) Align molecules (or frames)\n\n"
      $log.text insert end "Before deformation vector is calculated, you need to align the molecules (or frames) to have a meaningful depiction of structural change. "
      $log.text insert end "\n\n"
      $log.text insert end "Finally, click 'Calculate' button to generate depiction of structural change and NMWiz GUI."
    } elseif {$context == "frommolecule"} {
      $log.text insert end "From Molecule\n"
      $log.text insert end "=============\n\n"
      $log.text insert end "This interface allows the user to analyze normal mode data present in file formats that are recognized by VMD. "
      $log.text insert end "Follow these steps to analyze your data:\n\n"
      $log.text insert end "1) Load data into VMD\n\n"
      $log.text insert end "Normal mode data can be retrieved from a molecule with multiple frames in VMD. "
      $log.text insert end "The molecule must contain both coordinate and normal mode data. "
      $log.text insert end "First, you need to load the coordinate data as a new molecule and then load the normal mode data into the same molecule."
      $log.text insert end "\n\n"
      $log.text insert end "2) Select the molecule and atoms\n\n"
      $log.text insert end "Select the molecule with normal mode data. If you don't see the molecule, click 'Update' button. "
      $log.text insert end "You can also select a subset of atoms for which you want to display normal mode data graphics."
      $log.text insert end "\n\n"
      $log.text insert end "3) Specify data frames\n\n"
      $log.text insert end "Frames that contain coordinate and normal mode data must be specified. "
      $log.text insert end "Note that '0' is the index of the very first frame, and 'end' can be used to specify the last frame of the molecule. "
      $log.text insert end "\n\n"
      $log.text insert end "Finally, click load button to instantiate NMWiz window for selected data."
      $log.text insert end "\n\n\n"
      $log.text insert end "*TIP*: If normal mode data is calculated for all atom data for a large molecular system, select backbone or carbon alpha atoms "
      $log.text insert end "for more responsive visual analysis experience, e.g enter 'name CA' as the selection string."
    } elseif {$context == "main"} {
      $log.text insert end "NMWiz Main\n"
      $log.text insert end "==========\n\n\n"
      $log.text insert end "Load Normal Mode Data\n"
      $log.text insert end "---------------------\n\n"
      $log.text insert end "Main interface allows user to load data into NMWiz in two ways: "
      $log.text insert end "\n\n"
      $log.text insert end "**Load NMD File**\n\n"
      $log.text insert end "If you have an NMD file, click 'Load NMD file' button to select the file. "
      $log.text insert end "The contents will be loaded and a Wizard window associated with the data will appear."
      $log.text insert end "\n\n"
      $log.text insert end "**From Molecule**\n\n"
      $log.text insert end "Alternatively, when normal mode data is present in a file format that is recognized by VMD, "
      $log.text insert end "load the files into VMD as a molecule and click 'From Molecule' button. A window will appear to "
      $log.text insert end "facilitate selection of normal mode data from a molecule with multiple frames.\n\n"
      $log.text insert end "Note that, data obtained from a molecule can be saved in NMD format from the Main window and NMD files can be parsed with ProDy for further analysis."
      $log.text insert end "\n\n\n"
      $log.text insert end "Perform NMA Calculations\n"
      $log.text insert end "------------------------\n\n"
      $log.text insert end "You can use NMWiz to perform NMA calculations via ProDy for molecules loaded in VMD. "
      $log.text insert end "Click 'ProDy Interface' and follow the instructions therein for ANM, GNM, and PCA (EDA) calculations."
      $log.text insert end "\n\n\n"
      $log.text insert end "Settings and Options\n"
      $log.text insert end "--------------------\n\n"
      $log.text insert end "NMWiz saves some user settings in your home folder. These settings can be changed using 'Settings' window."
      $log.text insert end "\n\n"
      $log.text insert end "**Preserve View**\n\n"
      $log.text insert end "When NMWiz loads data, VMD will shift focus to the new molecule. Check this to preserve the current view when loading a new dataset."
    }
      
    $log.text yview moveto 0
    $log.text configure -state disabled
  }
  
  # Called by nmwiz_tk function
  # Makes the Main Window
  proc initGUI {} {
    variable w
    variable platform  
    if [winfo exists .nmwizgui] {
      wm deiconify .nmwizgui
      raise .nmwizgui
      return 
    }
    set w [toplevel .nmwizgui]
    wm title $w "NMWiz 1.0 - Main"
    wm resizable $w 0 0

    set wmf [frame $w.mainframe -bd 2]
    
    grid [button $wmf.loadnmd -width 20 -text "Load NMD File" -command {
      set tempfile [tk_getOpenFile \
        -filetypes {{"NMD files" { .nmd .NMD }} {"Text files" { .txt .TXT }} {"All files" *}}]
        if {![string equal $tempfile ""]} {::NMWiz::loadNMD $tempfile}}] \
      -row 3 -column 0 -columnspan 2 -sticky we

    grid [button $wmf.fromol -width 20 -text "From Molecule" -command ::NMWiz::initFromMolecule] \
      -row 5 -column 0 -columnspan 2 -sticky we

    grid [button $wmf.prody -width 20 -text "ProDy Interface" -command ::NMWiz::initProdyGUI] \
      -row 6 -column 0 -columnspan 2 -sticky we

    grid [button $wmf.compare -width 20 -text "Structure Comparison" -command ::NMWiz::initStrComp] \
      -row 7 -column 0 -columnspan 2 -sticky we

   
    grid [button $wmf.showhelp -text "Help" \
        -command {::NMWiz::showHelp main}] \
      -row 8 -column 0 -sticky we
    grid [button $wmf.website -text "Website" \
        -command "vmd_open_url http://www.csb.pitt.edu/NMWiz/"] \
      -row 8 -column 1 -sticky we
    
    if {[molinfo num] > 0} {
      set ::NMWiz::preserview 1
    }
    
    grid [checkbutton $wmf.preserview -text " preserve current view" \
        -variable ::NMWiz::preserview] \
      -row 10 -column 0 -columnspan 3 -sticky w

    #pack $wmf.options -side top -fill x -expand 1
    pack $wmf -side top -fill x -expand 1

    #if {$::NMWiz::guicount > -1} {    
    #  for {set i 0} {$i <= $::NMWiz::guicount} {incr i} {
    #    set ns "::nmgui$i"
    #    if {[namespace exists $ns]} {
        
        foreach ns [namespace children :: "nmdset*"] {
          set wgf [labelframe $w.{[string range $ns 2 end]}frame -text "[subst $${ns}::title]" -bd 2]
          grid [button $wgf.show -text "GUI" \
              -command "${ns}::nmwizgui" ] \
            -row 0 -column 0 -sticky we
          grid [button $wgf.remove -text "Remove" \
              -command "lset ::NMWiz::titles $::NMWiz::guicount NONE; pack forget $wgf; ${ns}::deleteMolecules; namespace delete $ns; destroy .[string range $ns 2 end]"] \
            -row 0 -column 1 -sticky we
          grid [button $wgf.save -text "Save" \
              -command "::NMWiz::writeNMD $ns"] \
            -row 0 -column 2 -sticky we
          pack $wgf -side top -fill x -expand 1
        }
    #    }
    #  }
    #}
  }
  
  
  variable nmwizColors "blue red gray orange yellow tan green white pink \
cyan purple black yellow2 yellow3 green2 green3 \
cyan2 cyan3 blue2 blue3 violet magenta magenta2 red2 red3 orange2 \
orange3"


  variable outputdir [pwd]
  variable defaultColor "yellow3"
  
  variable prodyMolecule
  variable prodyMolid -1
  variable prodyNFrames 0
  variable prodySelstr "protein and name CA or nucleic and name P C4' C2"
  variable prodySelAtoms 0
  variable prodyScript "ANM"
  variable prodyPrefix ""
  variable prodyRmCoords 0
  variable prodyPCAfile "DCD"
  variable prodyPCAAligned 0
  variable prodyTask ""
  variable prodyFrame 0
  variable prodyCutoff 15
  variable prodyExtend none
  variable prodyGNMCutoff 10
  variable prodyGamma 1
  variable prodyNModes 10
  variable prodyFirstFrame 0
  variable prodySkipFrame 1 
  variable prodyLastFrame end

  variable fromolMolecule
  variable fromolMolid -1
  variable fromolNFrames 0
  variable fromolSelstr "all"
  variable fromolSelAtoms 0
  variable fromolFrame 0
  variable fromolFirstFrame 1
  variable fromolLastFrame end
  
  variable strcompRefSelstr "name CA and protein"
  variable strcompTarSelstr "name CA and protein"
  variable strcompRefMol ""
  variable strcompTarMol ""
  variable strcompRefFrame 0
  variable strcompTarFrame 0
  variable strcompRefid -1
  variable strcompTarid -1
  variable strcompRefN 0
  variable strcompTarN 0
  
  proc initStrComp {} {
    variable strcompGUI
    # If already initialized, just turn on
    if [winfo exists .nmwizstrcomp] {
      wm deiconify .nmwizstrcomp
      raise .nmwizstrcomp
      return 
    }    
    set strcompGUI [toplevel .nmwizstrcomp]
    wm title $strcompGUI "NMWiz - Structure Comparison"
    wm resizable $strcompGUI 0 0

    # Main frame (molecule and selection)
    set wmf [labelframe $strcompGUI.mainFrame -text "Molecule Selection" -bd 2]
    grid [label $wmf.molRLabel -text "Reference"] \
      -row 1 -column 2 -sticky w
    grid [label $wmf.molTLabel -text "Target"] \
      -row 1 -column 3 -sticky w

    grid [label $wmf.molLabel -text "Molecule:"] \
      -row 2 -column 1 -sticky w
    grid [frame $wmf.refFrame] \
      -row 2 -column 2 -sticky ew
    tk_optionMenu $wmf.refFrame.list ::NMWiz::strcompRefMol "" 
    grid [frame $wmf.tarFrame] \
      -row 2 -column 3 -sticky ew
    tk_optionMenu $wmf.tarFrame.list ::NMWiz::strcompTarMol "" 
    grid [button $wmf.molUpdate -text "Update" \
        -command ::NMWiz::strcompUpdateMolList] \
      -row 2 -column 4 -sticky ew
      
    grid [label $wmf.selstrLabel -text "Selection:"] \
      -row 5 -column 1 -sticky w
    grid [entry $wmf.selstrRefEntry -width 12 -textvariable ::NMWiz::strcompRefSelstr] \
      -row 5 -column 2 -sticky ew
    grid [entry $wmf.selstrTarEntry -width 12 -textvariable ::NMWiz::strcompTarSelstr] \
      -row 5 -column 3 -sticky ew
    grid [button $wmf.selUpdate -text "Select" \
        -command ::NMWiz::strcompUpdateMolinfo] \
      -row 5 -column 4 -sticky ew
      
    grid [label $wmf.frameLabel -text "Frame:"] \
      -row 7 -column 1 -sticky w
    grid [entry $wmf.frameRefEntry -width 4 -textvariable ::NMWiz::strcompRefFrame] \
      -row 7 -column 2 -sticky w
    grid [entry $wmf.frameTarEntry -width 4 -textvariable ::NMWiz::strcompTarFrame] \
      -row 7 -column 3 -sticky w

    grid [label $wmf.selinfoLbl -text "Information:"] \
      -row 8 -column 1 -sticky w
    grid [label $wmf.selinfoRefLabel -text ""] \
      -row 8 -column 2 -sticky w
    grid [label $wmf.selinfoTarLabel -text ""] \
      -row 8 -column 3 -sticky w
    grid [label $wmf.selinfoRMSDLabel -text ""] \
      -row 8 -column 4 -sticky w
      
    grid [button $wmf.showHelp -text "Help" \
        -command {::NMWiz::showHelp compare}] \
      -row 12 -column 1 -sticky we
    grid [button $wmf.rmsdUpdate -text "RMSD" \
        -command ::NMWiz::strcompRMSD] \
      -row 12 -column 2 -sticky ew
    grid [button $wmf.align -text "Align" \
        -command ::NMWiz::strcompAlign] \
      -row 12 -column 3 -sticky ew
    grid [button $wmf.prodySubmit -text "Calculate" \
        -command ::NMWiz::calcDeform] \
      -row 12 -column 4 -sticky we
      
    pack $wmf -side top -fill x -expand 1
    ::NMWiz::strcompUpdateMolList
    ::NMWiz::strcompUpdateMolinfo


  }

  proc initFromMolecule {} {
    variable fromolGUI
    # If already initialized, just turn on
    if [winfo exists .nmwizfromol] {
      wm deiconify .nmwizfromol
      raise .nmwizfromol
      return 
    }    
    set fromolGUI [toplevel .nmwizfromol]
    wm title $fromolGUI "NMWiz - From Molecule"
    wm resizable $fromolGUI 0 0

    # Main frame (molecule and selection)
    set wmf [labelframe $fromolGUI.mainFrame -text "Molecule Selection" -bd 2]
    grid [label $wmf.molLabel -text "Molecule:"] \
      -row 2 -column 1 -sticky w
    grid [frame $wmf.molFrame] \
      -row 2 -column 2 -sticky ew
    tk_optionMenu $wmf.molFrame.list ::NMWiz::fromolMolecule "" 
    grid [button $wmf.molUpdate -text "Update" \
        -command ::NMWiz::fromolUpdateMolList] \
      -row 2 -column 3 -sticky ew
      
    grid [label $wmf.molinfoLbl -text "Information:"] \
      -row 3 -column 1 -sticky w
    grid [label $wmf.molinfoLabel -text ""] \
      -row 3 -column 2 -columnspan 2 -sticky w
    
    grid [label $wmf.selstrLabel -text "Selection:"] \
      -row 5 -column 1 -sticky w
    grid [entry $wmf.selstrEntry -width 20 -textvariable ::NMWiz::fromolSelstr] \
      -row 5 -column 2 -sticky ew
    grid [button $wmf.selUpdate -text "Select" \
        -command ::NMWiz::fromolUpdateSelection] \
      -row 5 -column 3 -sticky ew
      
    grid [label $wmf.selinfoLbl -text "Information:"] \
      -row 6 -column 1 -sticky w
    grid [label $wmf.selinfoLabel -text ""] \
      -row 6 -column 2 -columnspan 2 -sticky w

    grid [label $wmf.frameLabel -text "Coordinate frame:"] \
      -row 7 -column 1 -sticky w
    grid [entry $wmf.frameEntry -width 4 -textvariable ::NMWiz::fromolFrame] \
      -row 7 -column 2 -sticky w
    
    grid [label $wmf.firstLabel -text "First mode frame:"] \
      -row 8 -column 1 -sticky w
    grid [entry $wmf.firstEntry -width 4 -textvariable ::NMWiz::fromolFirstFrame] \
      -row 8 -column 2 -sticky w
    
    grid [label $wmf.lastLabel -text "Last mode frame:"] \
      -row 10 -column 1 -sticky w
    grid [entry $wmf.lastEntry -width 4 -textvariable ::NMWiz::fromolLastFrame] \
      -row 10 -column 2 -sticky w
      
    grid [button $wmf.showHelp -text "Help" \
        -command {::NMWiz::showHelp frommolecule}] \
      -row 12 -column 1 -sticky we
    grid [button $wmf.prodySubmit -text "Load data from molecule" \
        -command ::NMWiz::fromMolecule] \
      -row 12 -column 2 -columnspan 2 -sticky we
      
    pack $wmf -side top -fill x -expand 1
    ::NMWiz::fromolUpdateMolList
    ::NMWiz::fromolUpdateMolinfo

  }

  proc initProdyGUI {} {
    variable prodyGUI
    # If already initialized, just turn on
    if [winfo exists .nmwizprody] {
      wm deiconify .nmwizprody
      raise .nmwizprody
      return 
    }    
    set prodyGUI [toplevel .nmwizprody]
    wm title $prodyGUI "NMWiz - ProDy Interface"
    wm resizable $prodyGUI 0 0
    
    # Main frame (molecule and selection)
    set wmf [labelframe $prodyGUI.mainFrame -text "Atom Selection" -bd 2]
    grid [label $wmf.molLabel -text "Molecule:"] \
      -row 2 -column 1 -sticky w
    grid [frame $wmf.molFrame] \
      -row 2 -column 2 -sticky ew
    tk_optionMenu $wmf.molFrame.list ::NMWiz::prodyMolecule "" 
    grid [button $wmf.molUpdate -text "Update" \
        -command ::NMWiz::prodyUpdateMolList] \
      -row 2 -column 3 -sticky ew
    
    grid [label $wmf.molinfoLbl -text "Information:"] \
      -row 3 -column 1 -sticky w
    grid [label $wmf.molinfoLabel -text ""] \
      -row 3 -column 2 -columnspan 2 -sticky w
    
    grid [label $wmf.selstrLabel -text "Selection:"] \
      -row 5 -column 1 -sticky w
    grid [entry $wmf.selstrEntry -width 20 -textvariable ::NMWiz::prodySelstr] \
      -row 5 -column 2 -sticky ew
    grid [button $wmf.selUpdate -text "Select" \
        -command ::NMWiz::prodyUpdateSelection] \
      -row 5 -column 3 -sticky ew
      
    grid [label $wmf.selinfoLbl -text "Information:"] \
      -row 6 -column 1 -sticky w
    grid [label $wmf.selinfoLabel -text ""] \
      -row 6 -column 2 -columnspan 2 -sticky w
    
    pack $wmf -side top -fill x -expand 1
    ::NMWiz::prodyUpdateMolList
    ::NMWiz::prodyUpdateMolinfo
     
    # ProDy job frame
    set wf [labelframe $prodyGUI.jobFrame -text "ProDy Job Settings" -bd 2]
      
    grid [label $wf.scriptLabel -text "ProDy job:"] \
      -row 7 -column 1 -sticky w
    grid [frame $wf.scriptFrame] \
      -row 7 -column 2 -columnspan 2 -sticky ew
    tk_optionMenu $wf.scriptFrame.list ::NMWiz::prodyTask "ANM calculation" 
    $wf.scriptFrame.list.menu delete 0 last
    foreach script "ANM GNM PCA" {
      $wf.scriptFrame.list.menu add radiobutton -label "$script calculation" \
          -variable ::NMWiz::prodyTask \
          -command "set ::NMWiz::prodyScript $script; ::NMWiz::prodyChangeTask; ::NMWiz::prodyUpdatePrefix"
      incr counter  
    }
    pack $wf.scriptFrame.list -side left -anchor w -fill x
    variable prodyTask "ANM calculation"

    grid [label $wf.outdLabel -text "Output directory:"] \
      -row 8 -column 1 -sticky w
    grid [entry $wf.outdEntry -width 16 -textvariable ::NMWiz::outputdir] \
      -row 8 -column 2 -sticky ew
    grid [button $wf.outdBrowse -text "Browse" \
        -command {
      set tempdir [tk_chooseDirectory -initialdir $::NMWiz::outputdir ]
        if {![string equal $tempdir ""]} {set ::NMWiz::outputdir $tempdir}
        }] \
      -row 8 -column 3 -sticky ew

    grid [label $wf.filepLabel -text "Output filename:"] \
      -row 9 -column 1 -sticky w
    grid [entry $wf.filepEntry -width 20 -textvariable ::NMWiz::prodyPrefix] \
      -row 9 -column 2 -columnspan 2 -sticky we

    grid [checkbutton $wf.rmfileEntry -text " remove coordinate file upon job completion" \
        -variable ::NMWiz::prodyRmCoords] \
      -row 10 -column 1 -columnspan 3 -sticky w

    pack $wf -side top -fill x -expand 1
    

    # ANM frame
    set wf [labelframe $prodyGUI.anmFrame -text "ANM Settings" -bd 2]
    grid [label $wf.modesLabel -text "Number of modes:"] \
      -row 6 -column 1 -sticky w
    grid [entry $wf.modesEntry -width 4 -textvariable ::NMWiz::prodyNModes] \
      -row 6 -column 2 -sticky w   
    pack $wf -side top -fill x -expand 1
    
    grid [label $wf.frameLabel -text "Frame number:"] \
      -row 6 -column 3 -sticky w
    grid [entry $wf.frameEntry -width 4 -textvariable ::NMWiz::prodyFrame] \
      -row 6 -column 4 -sticky w

    grid [label $wf.cutoffLabel -text "Cutoff distance (A):"] \
      -row 8 -column 1 -sticky w
    grid [entry $wf.cutoffEntry -width 4 -textvariable ::NMWiz::prodyCutoff] \
      -row 8 -column 2 -sticky w

    grid [label $wf.gammaLabel -text "Force constant:"] \
      -row 8 -column 3 -sticky w
    grid [entry $wf.gammaEntry -width 4 -textvariable ::NMWiz::prodyGamma] \
      -row 8 -column 4 -sticky w    

    grid [label $wf.extendLabel -text "Extend model to:"] \
      -row 9 -column 1 -sticky w
    grid [frame $wf.extendto] \
      -row 9 -column 2 -columnspan 3 -sticky w    
    radiobutton $wf.extendto.all -text "all atoms"  -value "all" -variable ::NMWiz::prodyExtend
    radiobutton $wf.extendto.bb -text "backbone"  -value "bb" -variable ::NMWiz::prodyExtend
    radiobutton $wf.extendto.none -text "none"  -value "none" -variable ::NMWiz::prodyExtend
    pack $wf.extendto.all $wf.extendto.bb $wf.extendto.none -side left


    # GNM frame
    set wf [labelframe $prodyGUI.gnmFrame -text "GNM Settings" -bd 2]
    grid [label $wf.modesLabel -text "Number of modes:"] \
      -row 6 -column 1 -sticky w
    grid [entry $wf.modesEntry -width 4 -textvariable ::NMWiz::prodyNModes] \
      -row 6 -column 2 -sticky w   
    
    grid [label $wf.frameLabel -text "Frame number:"] \
      -row 6 -column 3 -sticky w
    grid [entry $wf.frameEntry -width 4 -textvariable ::NMWiz::prodyFrame] \
      -row 6 -column 4 -sticky w

    grid [label $wf.cutoffLabel -text "Cutoff distance (A):"] \
      -row 8 -column 1 -sticky w
    grid [entry $wf.cutoffEntry -width 4 -textvariable ::NMWiz::prodyGNMCutoff] \
      -row 8 -column 2 -sticky w

    grid [label $wf.gammaLabel -text "Force constant:"] \
      -row 8 -column 3 -sticky w
    grid [entry $wf.gammaEntry -width 4 -textvariable ::NMWiz::prodyGamma] \
      -row 8 -column 4 -sticky w    
      
    grid [label $wf.extendLabel -text "Extend model to:"] \
      -row 9 -column 1 -sticky w
    grid [frame $wf.extendto] \
      -row 9 -column 2 -columnspan 3 -sticky w    
    radiobutton $wf.extendto.all -text "all atoms"  -value "all" -variable ::NMWiz::prodyExtend
    radiobutton $wf.extendto.bb -text "backbone"  -value "bb" -variable ::NMWiz::prodyExtend
    radiobutton $wf.extendto.none -text "none"  -value "none" -variable ::NMWiz::prodyExtend
    pack $wf.extendto.all $wf.extendto.bb $wf.extendto.none -side left


    # PCA frame
    set wf [labelframe $prodyGUI.pcaFrame -text "PCA (EDA) Settings" -bd 2]
    
    grid [label $wf.modesLabel -text "Number of modes:"] \
      -row 6 -column 1 -sticky w
    grid [entry $wf.modesEntry -width 4 -textvariable ::NMWiz::prodyNModes \
      ] -row 6 -column 2 -sticky w
      
    grid [label $wf.skipLabel -text "Frame stride:"] \
      -row 6 -column 3 -sticky w
    grid [entry $wf.skipEntry -width 4 -textvariable ::NMWiz::prodySkipFrame \
          -validate all -validatecommand ::NMWiz::calcOutputSize \
      ] -row 6 -column 4 -sticky w

    grid [label $wf.firstLabel -text "First frame:"] \
      -row 8 -column 1 -sticky w
    grid [entry $wf.firstEntry -width 4 -textvariable ::NMWiz::prodyFirstFrame \
          -validate all -validatecommand ::NMWiz::calcOutputSize \
      ] -row 8 -column 2 -sticky w    

    grid [label $wf.lastLabel -text "Last frame:"] \
      -row 8 -column 3 -sticky w
    grid [entry $wf.lastEntry -width 4 -textvariable ::NMWiz::prodyLastFrame \
          -validate all -validatecommand ::NMWiz::calcOutputSize \
      ] -row 8 -column 4 -sticky w

    grid [label $wf.filetypeLabel -text "Trajectory type:"] \
      -row 10 -column 1 -sticky w
    grid [frame $wf.filetypeFrame] \
      -row 10 -column 2 -columnspan 3 -sticky ew
    tk_optionMenu $wf.filetypeFrame.list ::NMWiz::prodyPCAfile "DCD" 
    $wf.filetypeFrame.list.menu delete 0 last
    foreach script "DCD PDB" {
      $wf.filetypeFrame.list.menu add radiobutton -label "$script" \
          -variable ::NMWiz::prodyPCAfile \
          -command "::NMWiz::calcOutputSize"
      incr counter  
    }
    checkbutton $wf.filetypeFrame.alignedEntry -text " aligned" \
        -variable ::NMWiz::prodyPCAAligned
    
    pack $wf.filetypeFrame.list $wf.filetypeFrame.alignedEntry -side left \
      -anchor w -fill x
    variable prodyPCAfiletype "DCD"

    grid [label $wf.infolabel -text "Trajectory info:"] \
      -row 11 -column 1 -sticky w
    grid [label $wf.infoentry -text "" ] \
      -row 11 -column 2 -columnspan 3 -sticky w

    grid [label $wf.extendLabel -text "Extend model to:"] \
      -row 12 -column 1 -sticky w
    grid [frame $wf.extendto] \
      -row 12 -column 2 -columnspan 3 -sticky w    
    radiobutton $wf.extendto.all -text "all atoms"  -value "all" -variable ::NMWiz::prodyExtend
    radiobutton $wf.extendto.bb -text "backbone"  -value "bb" -variable ::NMWiz::prodyExtend
    radiobutton $wf.extendto.none -text "none"  -value "none" -variable ::NMWiz::prodyExtend
    pack $wf.extendto.all $wf.extendto.bb $wf.extendto.none -side left


    # Submit button
    set wf [frame $prodyGUI.submitFrame -bd 2]
      
    grid [button $wf.showHelp -text "Help" \
        -command {::NMWiz::showHelp prody}] \
      -row 0 -column 0 -sticky we
    grid [button $wf.prodySubmit -text "Submit Job" \
        -command ::NMWiz::prodySubmitJob] \
      -row 0 -column 1 -sticky we
    grid [button $wf.prodyWebsite -text "ProDy Website" \
        -command "vmd_open_url http://www.csb.pitt.edu/ProDy"] \
      -row 0 -column 2 -sticky we
    pack $wf -side top -fill x -expand 1
         
  }
  
  proc prodyUpdateMolList {} {
    variable prodyGUI
    set wf $prodyGUI.mainFrame
    $wf.molFrame.list.menu delete 0 last
    set counter 0
    variable prodyMolid -1
    foreach id [molinfo list] {
      if {[molinfo $id get numatoms] > 0 && [molinfo $id get numframes] > 0} {
        if {$counter == 0} {
          variable prodyMolid $id
        }
        $wf.molFrame.list.menu add radiobutton -label "[::NMWiz::cleanMolName $id] ($id)" \
            -variable ::NMWiz::prodyMolecule \
            -command "set ::NMWiz::prodyMolid $id; ::NMWiz::prodyUpdateMolinfo"
        incr counter  
      }
    }
    pack $wf.molFrame.list -side left -anchor w -fill x
    variable prodyMolid
    if {$prodyMolid > -1} {
      variable prodyMolecule "[::NMWiz::cleanMolName $prodyMolid] ($prodyMolid)"
    }
    ::NMWiz::prodyUpdateMolinfo
    ::NMWiz::calcOutputSize
  }
  
  proc calcOutputSize {} {
  
    variable prodyScript
    if {$prodyScript != "PCA"} { return 1}
     
    set numframes [molinfo $::NMWiz::prodyMolid get numframes]
    
    if {!([string is digit $::NMWiz::prodyFirstFrame] && $::NMWiz::prodyFirstFrame >= 0 && 
        $::NMWiz::prodyFirstFrame < $numframes)} { return 1}
    set first $::NMWiz::prodyFirstFrame
    
    if {!([string is digit $::NMWiz::prodySkipFrame] && $::NMWiz::prodySkipFrame > 0 && 
        $::NMWiz::prodySkipFrame < $numframes)} {return 1}
     set skip $::NMWiz::prodySkipFrame
    
    if {!($::NMWiz::prodyLastFrame == "end" || ([string is digit $::NMWiz::prodyLastFrame]
       && $::NMWiz::prodyLastFrame > 0 && $::NMWiz::prodyLastFrame < $numframes))} {return 1}
      
    if {$::NMWiz::prodyLastFrame == "end"} {
      set last [expr $numframes - 1]
    } else {
      set last $::NMWiz::prodyLastFrame
    }
    
    set count 0 
    for {set i $first} {$i <= $last} {incr i $skip} {
      incr count
    }
    if {$::NMWiz::prodyPCAfile == "DCD"} {
      set size [expr ($count * (56 + ($::NMWiz::prodySelAtoms + 2) * 12. ) + 276.)/ 1048576 ]
    } else {
      set size [expr ($count * ($::NMWiz::prodySelAtoms * 79. + 4) + 71.) / 1048576]
    }
    .nmwizprody.pcaFrame.infoentry configure -text [format "$count frames (~%.2f MB)" $size]
    return 1
  }

  proc strcompUpdateMolList {} {
    variable strcompGUI
    set wf $strcompGUI.mainFrame
    $wf.refFrame.list.menu delete 0 last
    $wf.tarFrame.list.menu delete 0 last
    set counter 0
    variable strcompRefid -1
    variable strcompTarid -1
    foreach id [molinfo list] {
      if {[molinfo $id get numatoms] > 0 && [molinfo $id get numframes] > 1} {
        if {$counter == 0} {
          variable strcompRefid $id
          variable strcompTarid $id
        }
        $wf.refFrame.list.menu add radiobutton -label "[::NMWiz::cleanMolName $id] ($id)" \
            -variable ::NMWiz::strcompRefMolecule \
            -command "set ::NMWiz::strcompRefid $id; ::NMWiz::strcompUpdateMolinfo"
        $wf.tarFrame.list.menu add radiobutton -label "[::NMWiz::cleanMolName $id] ($id)" \
            -variable ::NMWiz::strcompTarMolecule \
            -command "set ::NMWiz::strcompTarid $id; ::NMWiz::strcompUpdateMolinfo"
        incr counter  
      }
    }
    pack $wf.refFrame.list -side left -anchor w -fill x
    pack $wf.tarFrame.list -side left -anchor w -fill x
    if {$strcompTarid > -1} {
      variable strcompTarMol "[::NMWiz::cleanMolName $strcompTarid] ($strcompTarid)"
      variable strcompRefMol "[::NMWiz::cleanMolName $strcompRefid] ($strcompRefid)"
    }
    ::NMWiz::strcompUpdateMolinfo
  }
  
  proc strcompUpdateMolinfo {} {
    variable strcompRefid
    variable strcompTarid
    if {$strcompRefid > -1} {
      set ::NMWiz::strcompRefMol "[::NMWiz::cleanMolName $strcompRefid] ($strcompRefid)"
      set ref [atomselect $strcompRefid "$::NMWiz::strcompRefSelstr" frame $::NMWiz::strcompRefFrame]
      set ::NMWiz::strcompRefN [$ref num]
      .nmwizstrcomp.mainFrame.selinfoRefLabel configure \
        -text "$::NMWiz::strcompRefN selected"
      $ref delete
    } else {
      set ::NMWiz::strcompRefMol ""
      .nmwizstrcomp.mainFrame.selinfoRefLabel configure \
        -text "0 selected"
    }
    if {$strcompTarid > -1} {
      set ::NMWiz::strcompTarMol "[::NMWiz::cleanMolName $strcompTarid] ($strcompTarid)"
      set tar [atomselect $strcompTarid "$::NMWiz::strcompTarSelstr" frame $::NMWiz::strcompTarFrame]
      set ::NMWiz::strcompTarN [$tar num]
      .nmwizstrcomp.mainFrame.selinfoTarLabel configure \
        -text "$::NMWiz::strcompTarN selected"
      $tar delete
    } else {
      set ::NMWiz::strcompTarMol ""
      .nmwizstrcomp.mainFrame.selinfoTarLabel configure \
        -text "0 selected"
    }
    
  }
  proc strcompAlign {} {
    set ref [atomselect $::NMWiz::strcompRefid "$::NMWiz::strcompRefSelstr" frame $::NMWiz::strcompRefFrame]
    set tar [atomselect $::NMWiz::strcompTarid "$::NMWiz::strcompTarSelstr" frame $::NMWiz::strcompTarFrame]
    if {[$ref num] == [$tar num]} {
      set all [atomselect $::NMWiz::strcompTarid "all" frame $::NMWiz::strcompTarFrame]
      $all move [measure fit $tar $ref]
      .nmwizstrcomp.mainFrame.selinfoRMSDLabel configure \
        -text "RMSD = [format %.2f [measure rmsd $ref $tar]] A"
      $all delete
    } else {      
      .nmwizstrcomp.mainFrame.selinfoRMSDLabel configure \
        -text "Length mismatch"
    }
    $tar delete
    $ref delete
  }  
  proc strcompRMSD {} {
    set ref [atomselect $::NMWiz::strcompRefid "$::NMWiz::strcompRefSelstr" frame $::NMWiz::strcompRefFrame]
    set tar [atomselect $::NMWiz::strcompTarid "$::NMWiz::strcompTarSelstr" frame $::NMWiz::strcompTarFrame]
    if {[$ref num] == [$tar num]} {
      .nmwizstrcomp.mainFrame.selinfoRMSDLabel configure \
        -text "RMSD = [format %.2f [measure rmsd $ref $tar]] A"
    } else {
      .nmwizstrcomp.mainFrame.selinfoRMSDLabel configure \
        -text "Length mismatch"
    }
    $tar delete
    $ref delete
  }


  proc fromolUpdateMolList {} {
    variable fromolGUI
    set wf $fromolGUI.mainFrame
    $wf.molFrame.list.menu delete 0 last
    set counter 0
    variable fromolMolid -1
    foreach id [molinfo list] {
      if {[molinfo $id get numatoms] > 0 && [molinfo $id get numframes] > 1} {
        if {$counter == 0} {
          variable fromolMolid $id
        }
        $wf.molFrame.list.menu add radiobutton -label "[::NMWiz::cleanMolName $id] ($id)" \
            -variable ::NMWiz::fromolMolecule \
            -command "set ::NMWiz::fromolMolid $id; ::NMWiz::fromolUpdateMolinfo"
        incr counter  
      }
    }
    pack $wf.molFrame.list -side left -anchor w -fill x
    variable fromolMolid
    if {$fromolMolid > -1} {
      variable fromolMolecule "[::NMWiz::cleanMolName $fromolMolid] ($fromolMolid)"
    } 
    ::NMWiz::fromolUpdateMolinfo
  }

  
  proc prodyCheckMolecule {} {
    if {[lsearch [molinfo list] $::NMWiz::prodyMolid] > -1 && [molinfo $::NMWiz::prodyMolid get numframes] > 0} {
      return 1
    } 
    ::NMWiz::prodyUpdateMolList
    return 0
  }
  
  proc fromolCheckMolecule {} {
    if {[lsearch [molinfo list] $::NMWiz::fromolMolid] > -1 && [molinfo $::NMWiz::fromolMolid get numframes] > 1} {
      return 1
    }
    ::NMWiz::fromolUpdateMolList
    return 0
  }
  
  proc prodyUpdateMolinfo {} {
    variable prodyMolid
    if {$prodyMolid > -1} {
      variable prodyMolecule "[::NMWiz::cleanMolName $prodyMolid] ($prodyMolid)"
      set ::NMWiz::prodyNFrames [molinfo $::NMWiz::prodyMolid get numframes]
      .nmwizprody.mainFrame.molinfoLabel configure \
        -text "[molinfo $::NMWiz::prodyMolid get numatoms] atoms, $::NMWiz::prodyNFrames frames"
      ::NMWiz::prodyUpdateSelection
      ::NMWiz::prodyUpdatePrefix
    } else {
      set ::NMWiz::prodyNFrames 0
      variable fromolMolecule ""
      .nmwizprody.mainFrame.molinfoLabel configure \
        -text "Load a molecule and click Update."
      .nmwizprody.mainFrame.selinfoLabel configure \
        -text "Load a molecule and click Update."
    }
  }
  
  proc fromolUpdateMolinfo {} {
    variable fromolMolid
    if {$fromolMolid > -1} {
      variable fromolMolecule "[::NMWiz::cleanMolName $fromolMolid] ($fromolMolid)"
      set ::NMWiz::fromolNFrames [molinfo $::NMWiz::fromolMolid get numframes]
      .nmwizfromol.mainFrame.molinfoLabel configure \
        -text "[molinfo $::NMWiz::fromolMolid get numatoms] atoms, $::NMWiz::fromolNFrames frames"
      ::NMWiz::fromolUpdateSelection
    } else {
      set ::NMWiz::fromolMolecule ""
      set ::NMWiz::fromolNFrames 0
      .nmwizfromol.mainFrame.molinfoLabel configure \
        -text "Load a molecule and click Update."
      .nmwizfromol.mainFrame.selinfoLabel configure \
        -text "Load a molecule and click Update."
        
    }
  }
  
  proc prodyUpdatePrefix {} {
    if {[::NMWiz::prodyCheckMolecule]} {
      set prefix [::NMWiz::cleanMolName $::NMWiz::prodyMolid]
      
      if {[string range $prefix [expr [string length $prefix] - 4] end] == ".pdb"} {
        set prefix [string range $prefix 0 [expr [string length $prefix] - 5]]
      }
      if {$::NMWiz::prodyScript == "ANM"} {
        set ::NMWiz::prodyPrefix "$prefix\_anm"
      } elseif {$::NMWiz::prodyScript == "GNM"} {
        set ::NMWiz::prodyPrefix "$prefix\_gnm"
      } else {
        set ::NMWiz::prodyPrefix "$prefix\_pca"
      }
    }
  }
  
  proc cleanMolName {molid} {
    set name "[molinfo $molid get name]"
    foreach char {"\}" "\{"} { 
      set first [string first $char $name] 
      while {$first > -1} {
        set name [string replace $name $first $first ""] 
        set first [string first $char $name]
      }
    }
    foreach char {" " "."} { 
      set first [string first $char $name] 
      while {$first > -1} {
        set name [string replace $name $first $first "_"] 
        set first [string first $char $name]
      }
    }
    return $name
  }
  
  proc prodyUpdateSelection {} {
    ::NMWiz::prodyCheckMolecule
    variable prodyMolid
    variable prodySelstr
    set sel [atomselect $prodyMolid $prodySelstr]
    variable prodySelAtoms [$sel num]
    $sel delete
    variable prodyGUI
    $prodyGUI.mainFrame.selinfoLabel configure \
      -text "$prodySelAtoms atoms are selected"
    ::NMWiz::calcOutputSize  
  }
  
  proc fromolUpdateSelection {} {
    ::NMWiz::fromolCheckMolecule
    variable fromolMolid
    variable fromolGUI
    if {$fromolMolid > -1} {
      
      variable fromolSelstr
      set sel [atomselect $fromolMolid $fromolSelstr]
      variable fromolSelAtoms [$sel num]
      $sel delete
      $fromolGUI.mainFrame.selinfoLabel configure \
        -text "$fromolSelAtoms atoms are selected"
    } else {
      $fromolGUI.mainFrame.selinfoLabel configure \
        -text "Load a molecule and click Update."
    }
  }
  
  proc prodyChangeTask {} {
    variable prodyGUI
    variable prodyScript
    if {$prodyScript == "ANM"} {
      pack forget $prodyGUI.gnmFrame
      pack forget $prodyGUI.pcaFrame
      pack forget $prodyGUI.submitFrame
      pack $prodyGUI.anmFrame -side top -fill x -expand 1
      pack $prodyGUI.submitFrame -side top -fill x -expand 1
    } elseif {$prodyScript == "GNM"} {
      pack forget $prodyGUI.pcaFrame
      pack forget $prodyGUI.anmFrame
      pack forget $prodyGUI.submitFrame
      pack $prodyGUI.gnmFrame -side top -fill x -expand 1
      pack $prodyGUI.submitFrame -side top -fill x -expand 1
    } else {
      pack forget $prodyGUI.anmFrame
      pack forget $prodyGUI.gnmFrame
      pack forget $prodyGUI.submitFrame
      pack $prodyGUI.pcaFrame -side top -fill x -expand 1
      pack $prodyGUI.submitFrame -side top -fill x -expand 1
      ::NMWiz::calcOutputSize
    }
  }
  
  proc prodySubmitJob {} {
    
    set ::NMWiz::pybin [::ExecTool::find -interactive -description "Python executable" python]
    set ::NMWiz::prody [::ExecTool::find -interactive -description "ProDy script" prody]
  
    if {$::NMWiz::prodySelAtoms == 0} {
      tk_messageBox -type ok -title "ERROR" \
        -message "You need to make an atom selection before you can submit a job."
      return 
    }
    if {!([string is digit $::NMWiz::prodyNModes] && $::NMWiz::prodyNModes > 0)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Number of modes must be a number larger than 0."
      return 
    }
    if {![file isdirectory $::NMWiz::outputdir]} {
      tk_messageBox -type ok -title "ERROR" \
        -message "$::NMWiz::outputdir is not a valid directory."
      return 
    }
    
    if {$::NMWiz::prodyScript == "ANM"} {
      ::NMWiz::prodySubmitANMjob
    } elseif {$::NMWiz::prodyScript == "GNM"} {
      ::NMWiz::prodySubmitGNMjob
    } else {
      ::NMWiz::prodySubmitPCAjob
    }
  }
  proc prodySubmitANMjob {} {
    set n_frames [molinfo $::NMWiz::prodyMolid get numframes]
    if {!([string is digit $::NMWiz::prodyFrame] && $::NMWiz::prodyFrame >= 0 && 
        $::NMWiz::prodyFrame < $n_frames)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Frame number must be an integer from 0 to [expr $n_frames - 1]."
      return 
    }
    if {!([string is double $::NMWiz::prodyCutoff] && $::NMWiz::prodyCutoff > 4.5)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Cutoff distance (A) must be a number greater than 4.5."
      return 
    }
    if {!([string is double $::NMWiz::prodyGamma] && $::NMWiz::prodyGamma > 0)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Force constant must be a positive number."
      return 
    }
    set pdbfn [file join $::NMWiz::outputdir $::NMWiz::prodyPrefix.pdb]
    if {$::NMWiz::prodyExtend == "none"} {
      set sel [atomselect $::NMWiz::prodyMolid $::NMWiz::prodySelstr]
    } elseif {$::NMWiz::prodyExtend == "all"} {
      set sel [atomselect $::NMWiz::prodyMolid "same residue as ($::NMWiz::prodySelstr)"]
    } elseif {$::NMWiz::prodyExtend == "bb"} {
      set sel [atomselect $::NMWiz::prodyMolid "$::NMWiz::prodySelstr or (backbone and same residue as ($::NMWiz::prodySelstr))"]
    }    
    
    $sel frame $::NMWiz::prodyFrame
    $sel writepdb $pdbfn
    $sel delete
    
    set prefix [file join $::NMWiz::outputdir $::NMWiz::prodyPrefix]    
    if {$::NMWiz::prodyExtend == "none"} {
      vmdcon -info "Executing: $::NMWiz::pybin $::NMWiz::prody anm --quiet -s all -o \"$::NMWiz::outputdir\" -p \"$prefix\" -n $::NMWiz::prodyNModes -c $::NMWiz::prodyCutoff -g $::NMWiz::prodyGamma \"$pdbfn\""
      set status [exec $::NMWiz::pybin $::NMWiz::prody anm --quiet -s all -o "$::NMWiz::outputdir" -p "$prefix" -n $::NMWiz::prodyNModes -c $::NMWiz::prodyCutoff -g $::NMWiz::prodyGamma "$pdbfn"]
      set nmdfile "$prefix.nmd"
    } else {
      vmdcon -info "Executing: $::NMWiz::pybin $::NMWiz::prody anm --quiet -s \"$::NMWiz::prodySelstr\" -o \"$::NMWiz::outputdir\" -p \"$prefix\" -n $::NMWiz::prodyNModes -c $::NMWiz::prodyCutoff -g $::NMWiz::prodyGamma -t $::NMWiz::prodyExtend \"$pdbfn\""
      set status [exec $::NMWiz::pybin $::NMWiz::prody anm --quiet -s "$::NMWiz::prodySelstr" -o "$::NMWiz::outputdir" -p "$prefix" -n $::NMWiz::prodyNModes -c $::NMWiz::prodyCutoff -g $::NMWiz::prodyGamma -t $::NMWiz::prodyExtend "$pdbfn"]
      set nmdfile "$prefix\_extended_$::NMWiz::prodyExtend.nmd"
    }

    if {$status != -1} {
      tk_messageBox -type ok -title "INFO" \
        -message "ProDy ANM calculation is finished and results are being loaded."
      ::NMWiz::loadNMD $nmdfile
      if {$::NMWiz::prodyRmCoords} {
        file delete -force $pdbfn
      }
    }  else {
      tk_messageBox -type ok -title "ERROR" \
        -message "An error occured."
      file delete -force $pdbfn
    }
  }  
  proc prodySubmitGNMjob {} {
    set n_frames [molinfo $::NMWiz::prodyMolid get numframes]
    if {!([string is digit $::NMWiz::prodyFrame] && $::NMWiz::prodyFrame >= 0 && 
        $::NMWiz::prodyFrame < $n_frames)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Frame number must be an integer from 0 to [expr $n_frames - 1]."
      return 
    }
    if {!([string is double $::NMWiz::prodyGNMCutoff] && $::NMWiz::prodyGNMCutoff > 4.5)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Cutoff distance (A) must be a number greater than 4.5."
      return 
    }
    if {!([string is double $::NMWiz::prodyGamma] && $::NMWiz::prodyGamma > 0)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Force constant must be a positive number."
      return 
    }
    set pdbfn [file join $::NMWiz::outputdir $::NMWiz::prodyPrefix.pdb]
    if {$::NMWiz::prodyExtend == "none"} {
      set sel [atomselect $::NMWiz::prodyMolid $::NMWiz::prodySelstr]
    } elseif {$::NMWiz::prodyExtend == "all"} {
      set sel [atomselect $::NMWiz::prodyMolid "same residue as ($::NMWiz::prodySelstr)"]
    } elseif {$::NMWiz::prodyExtend == "bb"} {
      set sel [atomselect $::NMWiz::prodyMolid "$::NMWiz::prodySelstr or (backbone and same residue as ($::NMWiz::prodySelstr))"]
    }
    $sel frame $::NMWiz::prodyFrame
    $sel writepdb $pdbfn
    $sel delete
    
    set prefix [file join $::NMWiz::outputdir $::NMWiz::prodyPrefix]
    if {$::NMWiz::prodyExtend == "none"} {
      vmdcon -info "Executing: $::NMWiz::pybin $::NMWiz::prody gnm --quiet -s all -o \"$::NMWiz::outputdir\" -p \"$prefix\" -n $::NMWiz::prodyNModes -c $::NMWiz::prodyGNMCutoff -g $::NMWiz::prodyGamma \"$pdbfn\""
      set status [exec $::NMWiz::pybin $::NMWiz::prody gnm --quiet -s all -o "$::NMWiz::outputdir" -p "$prefix" -n $::NMWiz::prodyNModes -c $::NMWiz::prodyGNMCutoff -g $::NMWiz::prodyGamma "$pdbfn"]
      set nmdfile "$prefix.nmd"
    } else {
      vmdcon -info "Executing: $::NMWiz::pybin $::NMWiz::prody gnm --quiet -s \"$::NMWiz::prodySelstr\" -o \"$::NMWiz::outputdir\" -p \"$prefix\" -n $::NMWiz::prodyNModes -c $::NMWiz::prodyGNMCutoff -g $::NMWiz::prodyGamma -t $::NMWiz::prodyExtend \"$pdbfn\""
      set status [exec $::NMWiz::pybin $::NMWiz::prody gnm --quiet -s "$::NMWiz::prodySelstr" -o "$::NMWiz::outputdir" -p "$prefix" -n $::NMWiz::prodyNModes -c $::NMWiz::prodyGNMCutoff -g $::NMWiz::prodyGamma -t $::NMWiz::prodyExtend "$pdbfn"]
      set nmdfile "$prefix\_extended_$::NMWiz::prodyExtend.nmd"
    }

    if {$status != -1} {
      tk_messageBox -type ok -title "INFO" \
        -message "ProDy GNM calculation is finished and results are being loaded."
      ::NMWiz::loadNMD $nmdfile
      if {$::NMWiz::prodyRmCoords} {
        file delete -force $pdbfn
      } 
    }  else {
      tk_messageBox -type ok -title "ERROR" \
        -message "An error occured."
      file delete -force $pdbfn
    }
  }  
  proc prodySubmitPCAjob {} {
    if {$::NMWiz::prodyNFrames < 2} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Selected molecule must have more than 1 frames for PCA calculations."
      return       
    }
    if {!([string is digit $::NMWiz::prodyFirstFrame] && $::NMWiz::prodyFirstFrame >= 0 && 
        $::NMWiz::prodyFirstFrame < [molinfo $::NMWiz::prodyMolid get numframes])} {
      tk_messageBox -type ok -title "ERROR" \
        -message "First frame must be a number and must be in the valid range."
      return 
    }
    if {!([string is digit $::NMWiz::prodySkipFrame] && $::NMWiz::prodySkipFrame > 0 && 
        $::NMWiz::prodySkipFrame < [molinfo $::NMWiz::prodyMolid get numframes])} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Frame stride must be a positive number and must be in the valid range."
      return 
    }
    if {!($::NMWiz::prodyLastFrame == "end" || ([string is digit $::NMWiz::prodyLastFrame]
       && $::NMWiz::prodyLastFrame > 0 && $::NMWiz::prodyLastFrame < [molinfo $::NMWiz::prodyMolid get numframes]))} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Last frame may be \"end\" or a number in the valid range."
      return 
    }
    set sel [atomselect $::NMWiz::prodyMolid $::NMWiz::prodySelstr]
    set pdbfn [file join $::NMWiz::outputdir $::NMWiz::prodyPrefix.[string tolower $::NMWiz::prodyPCAfile]]
    set end $::NMWiz::prodyLastFrame 
    if {$end == "end"} {
      set end [expr $::NMWiz::prodyNFrames -1]
    }
    set prefix [file join $::NMWiz::outputdir $::NMWiz::prodyPrefix]
    #vmdcon -info "animate write pdb $pdbfn beg $::NMWiz::prodyFirstFrame end $end skip $::NMWiz::prodySkipFrame waitfor all sel $sel $::NMWiz::prodyMolid"
    if {$::NMWiz::prodyPCAfile == "DCD" | $::NMWiz::prodyExtend != "none"} {
      set nwritten [animate write dcd $pdbfn beg $::NMWiz::prodyFirstFrame end $end skip $::NMWiz::prodySkipFrame waitfor all sel $sel $::NMWiz::prodyMolid]
      if {$::NMWiz::prodyExtend != "none"} {
        $sel delete
        if {$::NMWiz::prodyExtend == "all"} {
          set sel [atomselect $::NMWiz::prodyMolid "same residue as ($::NMWiz::prodySelstr)"]
        } elseif {$::NMWiz::prodyExtend == "bb"} {
          set sel [atomselect $::NMWiz::prodyMolid "$::NMWiz::prodySelstr or (backbone and same residue as ($::NMWiz::prodySelstr))"]
        }
        $sel writepdb $prefix.pdb
      }
    } else {
      set nwritten [animate write pdb $pdbfn beg $::NMWiz::prodyFirstFrame end $end skip $::NMWiz::prodySkipFrame waitfor all sel $sel $::NMWiz::prodyMolid]
    }
    vmdcon -info "$nwritten frames are written as $pdbfn"
    
    $sel delete

    if {$::NMWiz::prodyPCAAligned} {
      if {$::NMWiz::prodyExtend == "none"} {
        vmdcon -info "Executing: $::NMWiz::pybin $::NMWiz::prody pca --quiet --aligned -s all -o \"$::NMWiz::outputdir\" -p \"$prefix\" -n $::NMWiz::prodyNModes \"$pdbfn\""
        set status [exec $::NMWiz::pybin $::NMWiz::prody pca --quiet --aligned -s all -o "$::NMWiz::outputdir" -p "$prefix" -n $::NMWiz::prodyNModes "$pdbfn"]
        set nmdfile "$prefix.nmd"
      } else {
        vmdcon -info "Executing: $::NMWiz::pybin $::NMWiz::prody pca --quiet --aligned -s \"$::NMWiz::prodySelstr\" -o \"$::NMWiz::outputdir\" -p \"$prefix\" -n $::NMWiz::prodyNModes -t $::NMWiz::prodyExtend --pdb \"$prefix.pdb\" \"$pdbfn\""
        set status [exec $::NMWiz::pybin $::NMWiz::prody pca --quiet --aligned -s "$::NMWiz::prodySelstr" -o "$::NMWiz::outputdir" -p "$prefix" -n $::NMWiz::prodyNModes -t $::NMWiz::prodyExtend --pdb "$prefix.pdb" "$pdbfn"]
        set nmdfile "$prefix\_extended_$::NMWiz::prodyExtend.nmd"
      }
    } else {
      if {$::NMWiz::prodyExtend == "none"} {
        vmdcon -info "Executing: $::NMWiz::pybin $::NMWiz::prody pca --quiet -s all -o \"$::NMWiz::outputdir\" -p \"$prefix\" -n $::NMWiz::prodyNModes \"$pdbfn\""
        set status [exec $::NMWiz::pybin $::NMWiz::prody pca --quiet -s all -o "$::NMWiz::outputdir" -p "$prefix" -n $::NMWiz::prodyNModes "$pdbfn"]
        set nmdfile "$prefix.nmd"
      } else {
        vmdcon -info "Executing: $::NMWiz::pybin $::NMWiz::prody pca --quiet -s \"$::NMWiz::prodySelstr\" -o \"$::NMWiz::outputdir\" -p \"$prefix\" -n $::NMWiz::prodyNModes -t $::NMWiz::prodyExtend --pdb \"$prefix.pdb\" \"$pdbfn\""
        set status [exec $::NMWiz::pybin $::NMWiz::prody pca --quiet -s "$::NMWiz::prodySelstr" -o "$::NMWiz::outputdir" -p "$prefix" -n $::NMWiz::prodyNModes -t $::NMWiz::prodyExtend --pdb "$prefix.pdb" "$pdbfn"]
        set nmdfile "$prefix\_extended_$::NMWiz::prodyExtend.nmd"
      }    }
    
    if {$status != -1} {
      tk_messageBox -type ok -title "INFO" \
        -message "ProDy PCA calculation is finished and results are being loaded."
      ::NMWiz::loadNMD $nmdfile 
      if {$::NMWiz::prodyRmCoords} {
        file delete -force $pdbfn
        if {$::NMWiz::prodyExtend == "none"} {
          file delete -force $prefix.pdb
        }
      }  
    }  else {
      tk_messageBox -type ok -title "ERROR" \
        -message "An error occured."
      file delete -force $pdbfn
      if {$::NMWiz::prodyExtend == "none"} {
        file delete -force $prefix.pdb
      }
    }
  }

  proc writeNMD {ns} {
    
    set tempfile [tk_getSaveFile -filetypes {{"NMD files" { .nmd .NMD }} {"All files" *}}]
    if {$tempfile == ""} {return}
    
    set fl [open $tempfile w]
    puts $fl "nmwiz_load $tempfile"
    puts $fl "name [subst $${ns}::title]"
    puts $fl "atomnames [subst $${ns}::atomnames]"
    puts $fl "resnames [subst $${ns}::resnames]"
    puts $fl "resids [subst $${ns}::resids]"
    puts $fl "chainids [subst $${ns}::chainids]"
    puts $fl "bfactors [subst $${ns}::bfactors]"
    puts $fl "coordinates [subst $${ns}::coordinates]"
    set indices [subst $${ns}::indices] 
    set lengths [subst $${ns}::lengths] 
    set modes [subst $${ns}::modes]
    for {set i 0} {$i < [llength $indices]} {incr i} {
      puts $fl "mode [lindex $indices $i] [lindex $lengths $i] [lindex $modes $i]"
    }
    close $fl    
  }
  
  proc loadNMD {fn} {
    if {![file isfile $fn]} {
      vmdcon -err "$fn is not a valid NMD file path."
      return
    }
    variable filename $fn
    vmdcon -info "NMWiz: Parsing file $filename"
    # Parse the file, and make sure coordinates are stored in the file
    #variable namespaces
    #variable nmwizguis
    
    set nmdfile [open $filename]
    set nmdlist [list]
    set coordinates 0 
    set n_dims 0
    while {[gets $nmdfile nmdline] != -1} { 
      if {[lindex $nmdline 0] == "coordinates"} {
        if {[expr [llength $nmdline] % 3] != 1} {
          tk_messageBox -type ok -title "ERROR" \
            -message "Length of the coordinate array in $filename must be a\
                      multiple of 3. An array of length\
                      [llength $coordinates] is provided."
          return
        }
        set coordinates [lrange $nmdline 1 end]
        set n_atoms [expr [llength $coordinates] / 3]
      }
      lappend nmdlist $nmdline
    }
    if {$coordinates == 0} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Coordinate data was not found in the input file.\
                  NMD files must contain system coordinate data."
      return
    }
    
    # Evaluate each line
    set title ""
    set atomnames ""
    set resnames ""
    set resids ""
    set bfactors ""
    set chainids ""
    set modes [list]
    set modelength 0
    set modecounter 0
    set dof 0
    foreach nmdline $nmdlist {
      switch -exact [lindex $nmdline 0] {
        name -
        title {
          set title [lrange $nmdline 1 end]
        }
        coordinates {
        }
        atomnames -
        names {
          if {[llength $nmdline] != [expr $n_atoms + 1]} {
            vmdcon -info "NMWiz WARNING: Length of atomnames array must be $n_atoms, not [expr [llength $nmdline] -1]."
          } else {
            set atomnames [lrange $nmdline 1 end]
          }
        }
        resnames {
          if {[llength $nmdline] != [expr $n_atoms + 1]} {
            vmdcon -info "NMWiz WARNING: Length of resnames array must be $n_atoms, not [expr [llength $nmdline] -1]."
          } else {
            set resnames [lrange $nmdline 1 end]
          }
        }
        chainids -
        chids -
        chains {
          if {[llength $nmdline] != [expr $n_atoms + 1]} {
            vmdcon -info "NMWiz WARNING: Length of chainids array must be $n_atoms, not [expr [llength $nmdline] -1]."
          } else {
            set chainids [lrange $nmdline 1 end]
          }
        }
        resids -
        resnums {
          if {[llength $nmdline] != [expr $n_atoms + 1]} {
            vmdcon -info "NMWiz WARNING: Length of resids array must be $n_atoms, not [expr [llength $nmdline] -1]."
          } else {
            set resids [lrange $nmdline 1 end]
          }
        }
        bfactors -
        betas {
          if {[llength $nmdline] != [expr $n_atoms + 1]} {
            vmdcon -info "NMWiz WARNING: Length of bfactors array must be $n_atoms, not [expr [llength $nmdline] -1]."
          } else {
            set bfactors [lrange $nmdline 1 end]
          }      
        }
        mode {
          incr modecounter
          set l [expr [llength $nmdline] - 1]
          if {$n_dims == 0} {
            set modelength $l
            if {$l >= [llength $coordinates]} {
              set diff [expr $l - [llength $coordinates]] 
              if {$diff > 2} {
                tk_messageBox -type ok -title "ERROR" \
                  -message "First mode line is not formatted correctly."
                return
              }
              set dof [llength $coordinates]
              set n_dims 3
              vmdcon -info "NMWiz INFO: File contains a 3D model."
            } else {
              set diff [expr $l - $n_atoms]
              if {$diff > 2} {
                tk_messageBox -type ok -title "ERROR" \
                  -message "First mode line is not formatted correctly."
                return
              }
              set dof $n_atoms
              set n_dims 1
              vmdcon -info "NMWiz INFO: File contains a 1D model."
            }
          }
          if {$modelength > 0 && $l != $modelength} {
            vmdcon -info "NMWiz WARNING: Mode line $modecounter does not have the same length as the first mode line."
          } else {
            switch -exact [expr $l - $dof] {
              0 {
                lappend modes [lrange $nmdline 1 end]
                lappend indices [llength $modes]
                lappend lengths 1
              }
              1 {
                lappend modes [lrange $nmdline 2 end]
                if {[string is integer [lindex $nmdline 1]]} {
                  lappend indices [lindex $nmdline 1]
                  lappend lengths 1
                } else {
                  lappend lengths [lindex $nmdline 1]
                  lappend indices [llength $modes]
                }
              }
              2 {
                lappend modes [lrange $nmdline 3 end]
                if {[string is integer [lindex $nmdline 1]]} {
                  lappend indices [lindex $nmdline 1]
                  lappend lengths [lindex $nmdline 2]
                } else {
                  lappend indices [lindex $nmdline 2]
                  lappend lengths [lindex $nmdline 1]
                }
              } 
              default {
                vmdcon -info "NMWiz WARNING: Mode data was not understood. Line starts with [lrange $nmdline 0 4]."
              }
            } 
          }
        }
        default {
          vmdcon -info "NMWiz WARNING: Unrecognized line starting with \"[lindex $nmdline 0]\""
        }
      }
      
    } 
    
    if {[llength $modes] == 0} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Mode data was not found in the input file."
      return
    }

    if {$title == ""} {
      set title "Untitled ($::NMWiz::guicount)"
      vmdcon -info "NMWiz INFO: Dataset is named as \"$title\"."
    }

    set ns [::NMWiz::makeNMWizGUI]
    ${ns}::initialize $coordinates $modes $n_dims $title $lengths $indices $atomnames $resnames $resids $chainids $bfactors
    ${ns}::nmwizgui
    ::NMWiz::appendGUIcontrols $ns
    return ${ns}::handle
  }

  proc calcDeform {} {
    variable strcompRefN
    variable strcompTarN
    variable strcompRefid
    variable strcompTarid
    variable strcompRefFrame 
    variable strcompTarFrame
    variable strcompRefMol 
    variable strcompTarMol
    variable strcompRefSelstr 
    variable strcompTarSelstr
    if {$strcompRefid == $strcompTarid && $strcompRefFrame == $strcompTarFrame} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Reference and target cannot be the same frame of the same molecule."
      return
    }
    if {$strcompRefN == 0 || $strcompRefN != $strcompTarN} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Reference and target selections must have same number of atoms."
      return
    }
    set title "Deformation $strcompRefMol -> $strcompTarMol"
    set ref [atomselect $strcompRefid $strcompRefSelstr frame $strcompRefFrame]
    set coordinates [concat {*}[$ref get {x y z}]]
    set modes [list]
    set lengths [list]
    set indices [list]
    set tar [atomselect $strcompTarid $strcompTarSelstr frame $strcompTarFrame]
    set target [concat {*}[$tar get {x y z}]]
    $tar delete
    set m [vecsub $target $coordinates]
    set len [veclength $m]
    lappend indices 1
    lappend modes [vecscale $m [expr 1 / $len]]
    lappend lengths $len
    set atomnames [$ref get name]
    set resnames [$ref get resname]
    set resids [$ref get resid]
    set chainids [$ref get chain]
    set bfactors [$ref get beta]
    $ref delete
    
    set ns [::NMWiz::makeNMWizGUI]
    ${ns}::initialize $coordinates $modes 3 $title $lengths $indices $atomnames $resnames $resids $chainids $bfactors
    ${ns}::nmwizgui
    ::NMWiz::appendGUIcontrols $ns
    
  }

  proc fromMolecule {} {  
    variable fromolMolecule
    variable fromolMolid
    variable fromolSelstr
    variable fromolFrame
    variable fromolFirstFrame
    variable fromolLastFrame
    set lastframe $fromolLastFrame
    if {$lastframe == "end"} {
      set lastframe [expr [molinfo $fromolMolid get numframes] - 1]
    }
    if {$fromolFrame >= $fromolFirstFrame && $fromolFrame <= $lastframe} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Coordinate frame cannot be one of the mode frames."
      return
    }
    set title $fromolMolecule
    set n_dims 3
    set sel [atomselect $fromolMolid $fromolSelstr frame $fromolFrame]
    set coordinates [concat {*}[$sel get {x y z}]]
    set modes [list]
    set lengths [list]
    set indices [list]
    for {set i $fromolFirstFrame} {$i <= $lastframe} {incr i} {
      $sel frame $i
      set m [concat {*}[$sel get {x y z}]]
      set length [veclength $m]
      lappend lengths $length  
      lappend modes [vecscale $m [expr 1 / $length]]
      lappend indices $i
    }
    $sel frame $fromolFrame
    set atomnames [$sel get name]
    set resnames [$sel get resname]
    set resids [$sel get resid]
    set chainids [$sel get chain]
    set bfactors [$sel get beta]
    
    $sel delete
    
    set ns [::NMWiz::makeNMWizGUI]
    ${ns}::initialize $coordinates $modes $n_dims $title $lengths $indices $atomnames $resnames $resids $chainids $bfactors
    ${ns}::nmwizgui
    ::NMWiz::appendGUIcontrols $ns
    
  }
  
  proc appendGUIcontrols {ns} {
    
    if [winfo exists .nmwizgui] {
      set w .nmwizgui
      set wgf [labelframe $w.{[string range $ns 2 end]}frame -text "[subst $${ns}::title]" -bd 2]
      
      grid [button $wgf.show -text "GUI" \
          -command "${ns}::nmwizgui" ] \
        -row 0 -column 0 -sticky we
      grid [button $wgf.remove -text "Remove" \
          -command "lset ::NMWiz::titles $::NMWiz::guicount NONE; pack forget $wgf; ${ns}::deleteMolecules; namespace delete $ns; destroy .[string range $ns 2 end]"] \
        -row 0 -column 1 -sticky we
      grid [button $wgf.save -text "Save" \
          -command "::NMWiz::writeNMD $ns"] \
        -row 0 -column 2 -sticky we

      pack $wgf -side top -fill x -expand 1
    
    }
  }
  
  proc array2xyz {arr} {
    set coords [list]
    foreach {x y z} $arr {
      lappend coords "$x $y $z"
    }
    return $coords
  }
  
  proc makeNMWizGUI {} {

    variable guicount
    incr guicount
    set ns "::nmdset$guicount"
    
    namespace eval $ns {
      variable tempfn ".[string range [namespace current] 2 end].pdb"
      variable w
      variable name
      variable molid -1
      variable arrid -1
      variable arridlist [list]
      variable animid -1
      variable animidlist [list]
      variable stopped 1 
      variable scalearrows 1
      variable sense +1
      variable color $::NMWiz::defaultColor
      variable colorlist [list]
      variable materials on
      variable material "HardPlastic"
      variable material_protein "HardPlastic"
      variable resolution 10
      variable resolution_protein 10
      variable ndim 3
      variable bothdirections 0
      variable porcupine 0
      
      #GNM option
      variable msformode "Mobility"
      
      variable selstr "all"
      variable selrep 0
      
      variable autoupdate 1
      variable autoanimate 0
      
      variable scalearrows_list [list]
      variable hide_shorter_list [list] 
      variable cylinder_radius_list [list]
      variable cone_radius_list [list]
      variable cone_height_list [list]
      variable resolution_list [list]
      variable material_list [list]

      variable hide_shorter 0.0 
      variable cylinder_radius 0.3
      variable cone_radius 0.5
      variable cone_height 1.0
      
      variable showproteinas "Tube"
      variable tuberadius 0.4
      variable bondradius 0.3   
      variable spherescale 0.6
      variable cutoffdistance 8.0
      variable proteincolor "Mobility" 
      variable color_scale_method [colorinfo scale method]
      variable color_scale_midpoint [colorinfo scale midpoint]
      
      variable nframes 50
      
      variable betalist {}
      variable betamin
      variable betamax
      
      variable hideprev 1
      variable autoplay 1

      variable arropt 0
      variable anmopt 0
      variable prtopt 0
      variable pltopt 0

      variable selid -1
      variable overplot 0
      variable selectscale 0.8
      variable plotwidth 800
      variable plotheight 600
      variable linewidth 1
      variable mradius 2
      variable dash "-"
      variable lornol "lines"
      variable marker "circle"
      variable plothandles [list]
      
      
      proc handle { args } {
        set cmd [lindex $args 0]
        
        if {![llength $cmd]} {
          vmdcon -info "nmdhandle commands: getcoords, setcoords, getmode, \
setmode, getlen, setlen, addmode"
          return
        }
        
        if {$cmd=="getcoords"} {
          variable coordinates
          return $coordinates
          
        } elseif {$cmd=="setcoords"} {
          set data [lindex $args 1]
          if {![llength $data]} {
            vmdcon -err "coordinate array is not provided"
            return
          }
          variable coordinates
          if {[llength $data] != [llength $coordinates]} {
            vmdcon -err "length of coordinate array is incorrect"
            return
          }
          set coordinates $data
          list xcoords
          list ycoords
          list zcoords
          foreach {x y z} $data {
            lappend xcoords $x 
            lappend ycoords $y 
            lappend zcoords $z
          }
          variable molid
          set ns [namespace current]
          set sel [atomselect $molid all]
          $sel set x $xcoords
          $sel set y $ycoords
          $sel set z $zcoords
          $sel delete
          ${ns}::updateProtRep $molid
          ${ns}::drawArrows
          
        } elseif {$cmd=="getmode"} {
          set index [lindex $args 1]
          if {![llength $index]} {
            vmdcon -err "mode index is not provided"
            return
          }
          if {![string is digit $index]} {
            vmdcon -err "mode index is not valid"
            return
          }
          variable indices
          set index [lsearch $indices $index] 
          if {$index < 0} {
            vmdcon -err "mode index is not valid"
            return            
          }
          variable modes
          return [lindex $modes $index] 
          
        } elseif {$cmd=="setmode"} {
          set index [lindex $args 1]
          if {![llength $index]} {
            vmdcon -err "mode index is not provided"
            return
          }
          if {![string is digit $index]} {
            vmdcon -err "mode index is not valid"
            return
          }
          set idx $index
          variable indices
          set index [lsearch $indices $index] 
          if {$index < 0} {
            vmdcon -err "mode index is not valid"
            return            
          }
          set data [lindex $args 2]
          if {![llength $data]} {
            vmdcon -err "mode array is not provided"
            return
          }
          variable modes
          if {[llength $data] != [llength [lindex $modes $index]]} {
            vmdcon -err "length of mode array is incorrect"
            return
          }
          lset modes $index $data
          variable activemode
          set activemode $idx 
          set ns [namespace current]
          ${ns}::changeMode
          ${ns}::drawArrows
        } elseif {$cmd=="nummodes"} {
          variable numofmodes
          return $numofmode
        
        } elseif {$cmd=="numatoms"} {
          variable n_atoms
          return $n_atoms

        } elseif {$cmd=="addmode"} {
          set data [lindex $args 1]
          if {![llength $data]} {
            vmdcon -err "mode array is not provided"
            return
          }
          variable modes
          if {[llength $data] != [llength [lindex $modes 0]]} {
            vmdcon -err "length of mode array is incorrect"
            return
          }
          lappend modes $data
          variable indices
          set index [expr [lindex $indices end] + 1]
          lappend indices $index 
          variable lengths
          lappend lengths 1

          variable w
          set ns [namespace current]
          $w.active_mode.active.list.menu add radiobutton -label $index \
              -variable ${ns}::activemode \
              -command "${ns}::changeMode;"
          $w.graphics_options.copyfrom.list.menu add radiobutton -label $index \
              -variable ${ns}::copyfrom

          variable arridlist          
          lappend arridlist -1
          variable animidlist 
          lappend animidlist -1
          variable colorlist
          variable color
          lappend colorlist $color
          variable hide_shorter_list
          variable hide_shorter
          lappend hide_shorter_list $hide_shorter
          variable cylinder_radius_list 
          variable cylinder_radius
          lappend cylinder_radius_list $cylinder_radius
          variable cone_radius_list
          variable cone_radius
          lappend cone_radius_list $cone_radius
          variable cone_height_list
          variable cone_height
          lappend cone_height_list $cone_height
          variable resolution_list
          variable resolution
          lappend resolution_list $resolution
          variable material_list
          variable material
          lappend material_list $material
          
          variable numofmodes
          incr numofmodes
          
          variable rmsd_list
          variable n_atoms
          set rmsd 2.0
          set one_over_root_of_n_atoms [expr 1 / $n_atoms ** 0.5]
          variable scalearrows_list
          lappend scalearrows_list [expr $rmsd / $one_over_root_of_n_atoms / [veclength $data]]
          variable rmsd_list
          lappend rmsd_list $rmsd
          
          set ns [namespace current]
          ${ns}::nmwizgui
          variable activemode
          set activemode $index 
          ${ns}::changeMode
          #${ns}::drawArrows
          
        } elseif {$cmd=="getlen"} {
          set index [lindex $args 1]
          if {![llength $index]} {
            vmdcon -err "mode index is not provided"
            return
          }
          variable indices
          set index [lsearch $indices $index] 
          if {$index < 0} {
            vmdcon -err "mode index is not valid"
            return            
          }
          variable lengths
          return [lindex $lengths $index] 

        } elseif {$cmd=="setlen"} {
          set index [lindex $args 1]
          if {![llength $index]} {
            vmdcon -err "mode index is not provided"
            return
          }
          if {![string is digit $index]} {
            vmdcon -err "mode index is not valid"
            return
          }
          variable indices
          set index [lsearch $indices $index] 
          if {$index < 0} {
            vmdcon -err "mode index is not valid"
            return            
          }

          set data [lindex $args 2]
          if {![llength $data]} {
            vmdcon -err "mode scalar is not provided"
            return
          }
          if {![string is double $data]} {
            vmdcon -err "mode scalar is not valid"
            return
          }
          variable lengths
          lset lengths $index $data 
        } else {
          vmdcon -err "$cmd is not a valid nmwiz command"
          return     
        }
        
        return
      }
      
      proc initialize {xyz m d t l i an rn ri ci bf} {
        variable arridlist
        variable animidlist
        variable colorlist
        variable coordinates $xyz
        variable modes $m
        variable ndim $d
        variable title $t
        variable lengths $l
        variable indices $i
        variable copyfrom [lindex $i 0]
        variable atomnames $an
        variable resnames $rn
        variable chainids $ci
        variable resids $ri
        variable bfactors $bf
        variable plotrids [list]
        variable bfactormin
        variable bfactormax
        variable n_dims 0
        variable n_atoms [expr [llength $coordinates] / 3]
        
        if {[lsearch $::NMWiz::titles $title] > -1} {
          set title "$title ($::NMWiz::guicount)"
        }
        lappend ::NMWiz::titles $title

        
        if {$atomnames == ""} {
          set atomnames [string repeat "CA " $n_atoms]
          vmdcon -info "NMWiz INFO: All atom names are set as \"CA\"."
        } 
        #else {
        #  variable showproteinas
        #  foreach an $atomnames {
        #    if {$an != "CA"} {
        #      set showproteinas "Licorice"
        #      break
        #    }
        #  }
        #}
        if {$resnames == ""} {
          set resnames [string repeat "GLY " $n_atoms]
          vmdcon -info "NMWiz INFO: All residue names are named set as \"GLY\"."
        }
        if {$chainids == ""} {
          set chainids [string repeat "A " $n_atoms]
          vmdcon -info "NMWiz INFO: All chain identifiers are set as \"A\"."
        }
        if {[llength $resids] == 0} {
          for {set i 1} {$i <= $n_atoms} {incr i} {lappend resids $i}
          vmdcon -info "NMWiz INFO: Residues are numbered starting from 1."
        }

        foreach i [lrange $resids 0 end-1] j [lrange $resids 1 end] {
          if {$i >= $j} {
            for {set i 1} {$i <= $n_atoms} {incr i} {lappend plotrids $i}
            break
          }
        }
        if {[llength $plotrids] == 0} {
          set plotrids $ri
          vmdcon -info "NMWiz INFO: Residue numbers will be used for plotting."
        }

        if {[llength $bfactors] == 0} {
          vmdcon -info "NMWiz INFO: Experimental bfactors were not found in the data file."
          for {set i 1} {$i <= $n_atoms} {incr i} {lappend bfactors 0.0}
          set bfactormin 0.0
          set bfactormax 0.0
        } else {
          set bfactormin [lindex $bfactors 0]
          set bfactormax [lindex $bfactors 0]
          foreach i $bfactors {
            if {$i < $bfactormin} {set bfactormin $i}
            if {$i > $bfactormax} {set bfactormax $i}
          }
        }
        
        variable prefix $title
        while {[string first " " $prefix] > -1} {
          set i [string first " " $prefix]
          set prefix [string replace $prefix $i $i "_"]
        }
        while {[string first "." $prefix] > -1} {
          set i [string first "." $prefix]
          set prefix [string replace $prefix $i $i "_"]
        }

        variable numofmodes [llength $modes] 
        variable activeindex 0
        variable drawlength [lindex $lengths $activeindex]
        variable drawlengthstr [format "%.1f" [lindex $lengths $activeindex]]
        variable activemode [lindex $indices $activeindex]
        variable betalist
        foreach atnm $atomnames {
          lappend betalist 0
        } 
        variable length [lindex $lengths 0]
        variable scalearrows_list
        variable rmsd 2
        variable rmsdprev 0
        variable activeindexprev 0
        variable rmsd_list
        set one_over_root_of_n_atoms [expr 1 / $n_atoms ** 0.5]
        foreach len $lengths mode $modes {
            set l [::tcl::mathfunc::abs $len] 
            lappend scalearrows_list [expr $rmsd / $one_over_root_of_n_atoms / [veclength $mode] / $l]
            lappend rmsd_list $rmsd  
        }
        variable scalearrows
        set scalearrows [lindex $scalearrows_list 0]
        variable scalarprev 0
                
        variable arridlist
        variable animidlist
        variable colorlist
        
        variable hide_shorter_list 
        variable cylinder_radius_list
        variable cone_radius_list
        variable cone_height_list
        variable resolution_list
        variable material_list
        
        variable hide_shorter 
        variable cylinder_radius
        variable cone_radius
        variable cone_height
        variable resolution
        variable material
        set shift 0
        for {set i 0} {$i < [llength $modes]} {incr i} {
          lappend arridlist -1
          lappend animidlist -1
          set curcolor [lindex $::NMWiz::nmwizColors [expr ($i + $shift + $::NMWiz::guicount) % [llength $::NMWiz::nmwizColors]]]
          if {$curcolor == $::NMWiz::defaultColor} {
            incr shift 
            set curcolor [lindex $::NMWiz::nmwizColors [expr ($i + $shift + $::NMWiz::guicount) % [llength $::NMWiz::nmwizColors]]]
          }
          lappend colorlist $curcolor
          lappend hide_shorter_list $hide_shorter 
          lappend cylinder_radius_list $cylinder_radius
          lappend cone_radius_list $cone_radius
          lappend cone_height_list $cone_height
          lappend resolution_list $resolution
          lappend material_list $material
        }
        if {$::NMWiz::guicount == 0} {
          lset colorlist 0 $::NMWiz::defaultColor  
        }
        variable color [lindex $colorlist 0]
      }
      
      #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      proc deleteMolecules {} {
        variable molid
        variable animidlist
        variable selid
        variable arridlist
        if {[lsearch [molinfo list] $molid] > -1} {
          mol delete $molid
        }
        foreach arrid $arridlist {
          if {[lsearch [molinfo list] $arrid] > -1} {
            mol delete $arrid
          }
        }
        foreach animid $animidlist {
          if {[lsearch [molinfo list] $animid] > -1} {
            mol delete $animid
          }
        }
        if {[lsearch [molinfo list] $selid] > -1} {
          mol delete $selid
        }
      }
      
      proc Plot {ns} {
        
        #puts [subst $${ns}::overplot]
        set plothandle 0
        if {[subst $${ns}::overplot]} {
          for {set i [expr [llength [subst $${ns}::plothandles]] - 1]} {$i >= 0} {incr i -1} {
            set temphandle [lindex [subst $${ns}::plothandles] $i]
            if {[namespace exists [string range $temphandle 0 end-12]]} {
              set plothandle $temphandle
              break
            } else {
              set ${ns}::plothandles [lrange [subst $${ns}::plothandles] 0 $i]
            }
          }
          
        }
        if {$plothandle != 0} {
          $plothandle add \
            [subst $${ns}::plotrids] [subst $${ns}::betalist] \
            -title "[subst $${ns}::title] square fluctuations" \
            -linewidth [subst $${ns}::linewidth] \
            -legend "Mode [subst $${ns}::activemode]" \
            -[subst $${ns}::lornol] -linecolor [subst $${ns}::color] \
            -xsize [subst $${ns}::plotwidth] -ysize [subst $${ns}::plotheight] \
            -radius [subst $${ns}::mradius] \
            -fillcolor [subst $${ns}::color] -marker [subst $${ns}::marker] \
            -xlabel "Atom/Residue #" \
            -callback $ns\::highlight \
            -plot
             
        } else {
          lappend ${ns}::plothandles [multiplot \
            -x [subst $${ns}::plotrids] -y [subst $${ns}::betalist] \
            -title "[subst $${ns}::title] square fluctuations" \
            -linewidth [subst $${ns}::linewidth] \
            -legend "Mode [subst $${ns}::activemode]" \
            -[subst $${ns}::lornol] -linecolor [subst $${ns}::color] \
            -xsize [subst $${ns}::plotwidth] -ysize [subst $${ns}::plotheight] \
            -radius [subst $${ns}::mradius] \
            -fillcolor [subst $${ns}::color] -marker [subst $${ns}::marker] \
            -xlabel "Atom/Residue #" \
            -callback $ns\::highlight \
            -plot]
        }
        vmdcon -info "Plot handle: [lindex [subst $${ns}::plothandles] end]"
        #-dash [subst $${ns}::dash] \
      }
      
      proc prepareSelmol {} {
        variable selid
        variable molid
        # make sure selmol exists
        if {$selid == -1 || [lsearch [molinfo list] $selid] == -1} {
          set selid [mol new]
          variable title
          variable w
          mol rename $selid "$title selections"
          $w.draw_arrows.plot_label configure -text "Plotting ($selid):"
        }
        # make sure coordinates are loaded
        if {[molinfo $selid get numframes] == 0} {
          variable molid
          set currentview [molinfo $molid get {rotate_matrix center_matrix scale_matrix global_matrix}]

          variable coordinates
          variable tempfn
          set outfile [open [file join $::NMWiz::tmpdir $tempfn] w]
          foreach line [[namespace current]::getPDBLines $coordinates] {
            puts $outfile $line
          } 
          close $outfile
          mol addfile [file join $::NMWiz::tmpdir $tempfn] molid $selid
          
          foreach id [molinfo list] {
            molinfo $id set {rotate_matrix center_matrix scale_matrix global_matrix} $currentview
          }
        }
        
        variable resids
        if {[molinfo $selid get numreps] != [llength $resids]} {
          for {set i [molinfo $selid get numreps]} {$i >= 0} {incr i -1} {
            mol delrep $i $selid
          }
          for {set i 0} {$i < [llength $resids]} {incr i} {
            mol addrep $selid
            mol modstyle $i $selid VDW
            mol showrep $selid $i off
          }
        }
        mol top $molid
      }

      proc clearSelection {} {
        variable selid
        if {$selid > -1 && [lsearch [molinfo list] $selid] > -1} {
          for {set i [expr [molinfo $selid get numreps] - 1]} {$i >= 0} {incr i -1} {
            mol showrep $selid $i off
          }
        }
        set labels [label list Atoms]
        for {set i [expr [llength $labels] - 1]} {$i >= 0} {incr i -1} {
          if {[lindex [lindex [lindex $labels $i] 0] 0] == $selid} {
            label delete Atoms $i
          }
        }
      }

      proc highlight {args} {

        [namespace current]::loadCoordinates
        set resid [lindex $args 0] 
        set y [lindex $args 1]
        set color [lindex $args 2]
        variable plotrids
        set which [lsearch $plotrids $resid]
        if {$which == -1} {return 0}
        [namespace current]::prepareSelmol
        
        variable selid
        variable selection
        variable selectscale
        variable resolution
        variable material
        variable chainids
        variable resnames
        variable resids

        label add Atoms $selid/$which

        #set i [molinfo $selid get numreps]
        #mol addrep $selid
        if {[mol showrep $selid $which]} {
          mol showrep $selid $which off
          vmdcon -info "Deselected [lindex $chainids $which]:[lindex $resnames $which][lindex $resids $which]"
        } else {
          vmdcon -info "Selected [lindex $chainids $which]:[lindex $resnames $which][lindex $resids $which]"
          mol showrep $selid $which on
          mol modstyle $which $selid VDW $selectscale $resolution
          mol modmaterial $which $selid $material
          mol modselect $which $selid "index $which"
          mol modcolor $which $selid ColorID [lsearch "blue red gray orange yellow tan silver green white pink cyan purple lime mauve ochre iceblue black yellow2 yellow3 green2 green3 cyan2 cyan3 blue2 blue3 violet violet2 magenta magenta2 red2 red3 orange2 orange3" $color]
        }
      }

      proc updateProtRep {targetid} {
        variable molid

        if {[lsearch [molinfo list] $targetid] == -1} {
          if {$targetid == $molid} {
            [namespace current]::loadCoordinates
          } else {
            return 0
          }
        }
        variable showproteinas
        if {$showproteinas == "Custom" && $targetid == $molid} {
          return
        }
        variable tuberadius
        variable bondradius
        variable cutoffdistance
        variable resolution_protein
        variable material_protein
        variable betamin
        variable betamax
        variable spherescale
        variable proteincolor
        variable bfactors
        variable bfactormin
        variable bfactormax
        variable betalist
        variable selstr
        variable selrep
        variable color_scale_method
        variable color_scale_midpoint
        
        for {set i [molinfo $targetid get numreps]} {$i >= 0} {incr i -1} {
          mol delrep $i $targetid
        }
        
        set all [atomselect $targetid "all"]
        set midpoint 0.5
        if {$proteincolor == "Bfactors" | $proteincolor == "Mobility"} {
          if {$proteincolor == "Bfactors"} {
            $all set beta $bfactors
          } elseif {$proteincolor == "Mobility"} {
            $all set beta $betalist
          }
        }
        $all delete
        
        variable msformode
        if {$showproteinas == "Network"} {
          mol addrep $targetid
          mol modstyle 0 $targetid DynamicBonds $cutoffdistance $bondradius $resolution_protein
          mol modmaterial 0 $targetid $material_protein
          mol addrep $targetid
          mol modstyle 1 $targetid VDW $spherescale $resolution_protein
          mol modmaterial 1 $targetid $material_protein
          if {$proteincolor == "Mobility"} {
            if {$msformode == "Eigenvector"} {
                mol modcolor 0 $targetid Beta

                mol addrep $targetid
                mol modstyle 2 $targetid DynamicBonds $cutoffdistance $bondradius $resolution_protein
                mol modmaterial 2 $targetid $material_protein
    
                mol modcolor 1 $targetid ColorID 0
                mol modselect 1 $targetid "beta <= 0"
                mol modcolor 2 $targetid ColorID 0
                mol modselect 2 $targetid "beta <= 0"
                
                mol addrep $targetid
                mol modstyle 3 $targetid VDW $spherescale $resolution_protein
                mol modmaterial 3 $targetid $material_protein
                mol addrep $targetid
                mol modstyle 4 $targetid DynamicBonds $cutoffdistance $bondradius $resolution_protein
                mol modmaterial 4 $targetid $material_protein
                
                mol modcolor 3 $targetid ColorID 1
                mol modselect 3 $targetid "beta > 0"
                mol modcolor 4 $targetid ColorID 1
                mol modselect 4 $targetid "beta > 0"


            } else { 
              mol modcolor 0 $targetid Beta
              mol scaleminmax $targetid 0 $betamin $betamax 
              mol modcolor 1 $targetid Beta
              mol scaleminmax $targetid 1 $betamin $betamax
            } 
          } elseif {$proteincolor == "Bfactors"} {
            mol modcolor 0 $targetid Beta
            mol scaleminmax $targetid 0 $bfactormin $bfactormax 
            mol modcolor 1 $targetid Beta
            mol scaleminmax $targetid 1 $bfactormin $bfactormax
          } else {
            mol modcolor 0 $targetid $proteincolor
            mol modcolor 1 $targetid $proteincolor
          }
          if {$selrep} {
            mol modselect 0 $targetid $selstr
            mol modselect 1 $targetid $selstr
          }
        } else {
          if {$msformode == "Eigenvector"} { 
            set n_reps 2 
          } else {
            set n_reps 1
          }
          for {set i 0} {$i < $n_reps} {incr i} {
              mol addrep $targetid
              switch $showproteinas {
                "Ribbons" {
                  mol modstyle $i $targetid $showproteinas 0.3 $resolution_protein 2
                }
                "NewRibbons" {
                  mol modstyle $i $targetid $showproteinas 0.3 $resolution_protein 3
                }
                "Cartoon" {
                  mol modstyle $i $targetid $showproteinas 2.1 $resolution_protein 5
                }
                "NewCartoon" {
                  mol modstyle $i $targetid $showproteinas 0.3 $resolution_protein 4.1
                }
                "CPK" {
                  mol modstyle $i $targetid $showproteinas 1.0 0.3 $resolution_protein $resolution_protein
                }
                "VDW" {
                  mol modstyle $i $targetid $showproteinas $spherescale $resolution_protein
                }
                "Lines" {
                  mol modstyle $i $targetid $showproteinas
                }
                "Licorice" {
                  mol modstyle $i $targetid $showproteinas $tuberadius $resolution_protein $resolution_protein
                }
                default {
                  mol modstyle $i $targetid $showproteinas $tuberadius $resolution_protein
                }            
              }
          }
         
          mol modmaterial 0 $targetid $material_protein
          if {$selrep} {
            mol modselect 0 $targetid $selstr
          }
          if {$proteincolor == "Mobility"} {
            variable msformode
            if {$msformode == "Eigenvector"} {
                mol modcolor 0 $targetid ColorID 0
                mol modselect 0 $targetid "beta <= 0"
                mol modcolor 1 $targetid ColorID 1
                mol modselect 1 $targetid "beta > 0"
            } else {          
                mol modcolor 0 $targetid Beta
                mol scaleminmax $targetid 0 $betamin $betamax
            }
          } elseif {$proteincolor == "Bfactors"} {
            mol modcolor 0 $targetid Beta
            mol scaleminmax $targetid 0 $bfactormin $bfactormax
          } else {
            mol modcolor 0 $targetid $proteincolor
          }
        }
      }

      proc calcMSF {} {
        variable molid
        variable activemode
        variable indices
        variable lengths
        variable scalearrows
        variable modes
        variable animid
        variable material
        variable selstr
        variable ndim
        variable color_scale_method
        variable color_scale_midpoint

        variable length
        set mode [lindex $modes [lsearch $indices $activemode]]
        variable msformode
        if {$msformode == "Mobility"} {
          set mode [vecscale [expr $length * $length] [vecmul $mode $mode]]  
        } else {
          set mode [vecscale $length $mode]
        }

        set index 0
        variable betalist {}
        variable betamin 10000
        variable betamax -10000
        if {$ndim == 3} {
          foreach {mx my mz} $mode {
            set beta [expr $mx + $my + $mz]
            lappend betalist $beta
            if {$beta < $betamin} {set betamin $beta}
            if {$beta > $betamax} {set betamax $beta}
            incr index
          }
        } else {
          foreach mx $mode {
            set beta [expr $mx]
            lappend betalist $beta
            if {$beta < $betamin} {set betamin $beta}
            if {$beta > $betamax} {set betamax $beta}
            incr index
          }
        }
        set all [atomselect $molid "all"] 
        $all set beta $betalist
        $all delete 
        [namespace current]::updateProtRep $molid
      }
      
      proc autoUpdate {} {
        variable autoupdate
        if {$autoupdate} {
          [namespace current]::drawArrows
        }
      }
      
      proc copySettings {} {
        
        variable indices
        variable copyfrom 
        set from [lsearch $indices $copyfrom]
        variable hide_shorter_list 
        variable cylinder_radius_list
        variable cone_radius_list
        variable cone_height_list
        variable resolution_list
        variable material_list
      
        variable hide_shorter [lindex $hide_shorter_list $from]
        variable cylinder_radius [lindex $cylinder_radius_list $from]
        variable cone_radius [lindex $cone_radius_list $from]
        variable cone_height [lindex $cone_height_list $from]
        variable resolution [lindex $resolution_list $from]
        variable material [lindex $material_list $from]
        [namespace current]::autoUpdate
      }
      
      proc evalRMSD {} {
        variable rmsd
        if {![string is double $rmsd]} {
          set rmsd 2
        } elseif {$rmsd <= 0} {
          set rmsd 0.1
        }
        variable scalearrows
        variable scalarprev
        variable lengths
        variable modes
        variable activemode
        variable indices
        variable n_atoms
        set whichmode [lsearch $indices $activemode]
        variable length [lindex $lengths $whichmode]
        set scalearrows [expr [sign $scalearrows] * $rmsd / $length / [veclength [lindex $modes $whichmode]] * $n_atoms ** 0.5 ]
        set scalarprev $scalearrows
        #vmdcon -info "Mode $activemode is scaled by [format %.2f $scalearrows](x[format %.2f $length]) for [format %.2f $rmsd] A RMSD."
        variable rmsdprev $rmsd
      }

      proc evalScale {} {
        variable scalearrows
        variable rmsd
        if {![string is double $scalearrows]} {
          set rmsd 2
          [namespace current]::evalRMSD
          return
        }
        variable lengths
        variable modes
        variable activemode
        variable indices
        variable n_atoms
        set whichmode [lsearch $indices $activemode]
        variable length [lindex $lengths $whichmode]
        variable scaleprev $scalearrows
        variable rmsdprev
        set rmsd [expr [sign $scalearrows] * $scalearrows * $length * [veclength [lindex $modes $whichmode]] / $n_atoms ** 0.5]
        set rmsdprev $rmsd
        #vmdcon -info "Mode $activemode is scaled by [format %.2f $scalearrows](x[format %.2f $length]) for [format %.2f $rmsd] A RMSD."
      }

      proc drawArrows {} {
        variable rmsd
        variable rmsdprev
        variable activeindex
        variable activeindexprev
        variable scalearrows 
        variable scalarprev
        if {$rmsdprev != $rmsd || $activeindex != $activeindexprev} {
          [namespace current]::evalRMSD
        } elseif {$scalarprev != $scalearrows || $activeindex != $activeindexprev} {
          [namespace current]::evalScale
        }
        variable color
        variable material
        variable resolution
        variable coordinates
        variable sense
        variable arrid
        variable activemode
        variable indices
        variable modes
        variable title
        variable prefix
        variable lengths
        variable molid
        variable w
        variable selstr
        
        variable hide_shorter 
        variable cylinder_radius
        variable cone_radius
        variable cone_height
        
        set whichmode [lsearch $indices $activemode] 

        if {[lsearch [molinfo list] $molid] == -1} {
          [namespace current]::loadCoordinates
        }

        if {[lsearch [molinfo list] $arrid] == -1} {
          set arrid [mol new]
        } else {
          graphics $arrid delete all
        }
        graphics $arrid color $color
        graphics $arrid materials on
        graphics $arrid material $material
        #set length [lindex $lengths $whichmode]
        variable length
        set mode [vecscale [expr $length * $scalearrows] [lindex $modes $whichmode]]
        
        variable bothdirections
        variable porcupine
        set sel [atomselect $molid $selstr]
        foreach index [$sel get index] {
          set from [expr $index * 3]
          set to  [expr $from + 2]
          set xyz [lrange $coordinates $from $to ] 
          set v [lrange $mode $from $to ] 
          if {$hide_shorter < [veclength $v]} {
            set temp [vecadd $xyz $v]
            if {$porcupine} {
              graphics $arrid cone $xyz $temp radius $cone_radius resolution $resolution
            } else { 
              graphics $arrid cylinder $xyz $temp radius $cylinder_radius resolution $resolution
              set temp2 [vecadd $temp [vecscale $v [expr $cone_height / [veclength $v]]]]
              graphics $arrid cone $temp $temp2 radius $cone_radius resolution $resolution
            }
            if {$bothdirections} {
              set v [vecscale -1 $v]
              set temp [vecadd $xyz $v]
              if {$porcupine} {
                graphics $arrid cone $xyz $temp radius $cone_radius resolution $resolution
              } else { 
                graphics $arrid cylinder $xyz $temp radius $cylinder_radius resolution $resolution
                set temp2 [vecadd $temp [vecscale $v [expr $cone_height / [veclength $v]]]]
                graphics $arrid cone $temp $temp2 radius $cone_radius resolution $resolution
              }
            }
          }
        }
        $sel delete

        mol rename $arrid "$title mode $activemode arrows"
        set currentview [molinfo $molid get {rotate_matrix center_matrix scale_matrix global_matrix}]
        display resetview
        foreach id [molinfo list] {
          molinfo $id set {rotate_matrix center_matrix scale_matrix global_matrix} $currentview
        }
        [namespace current]::calcMSF
        $w.draw_arrows.arrowbuttons_label configure -text "Mode ($arrid):"
                
        variable arridlist
        lset arridlist $whichmode $arrid
      }

      proc getPDBLines {coords} {
        variable atomnames
        variable resnames
        variable chainids
        variable resids
        variable betalist
        set pdblines ""
        set i 0
        foreach an $atomnames rn $resnames ci $chainids ri $resids {x y z} $coords b $betalist {
          incr i
          if {[string length $an] < 4} {
            #lappend pdblines [format "ATOM  %5d  %-3s %-4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" \
            #                     $i  $an $rn  $ci $ri $x $y $z 1.0 $b]
            lappend pdblines [format "ATOM  %5d  %-3s %-4s%1s%4d    %8.3f%8.3f%8.3f" \
                                 $i  $an [string range $rn 0 4] $ci $ri $x $y $z]
          } else {
            set an [string range $an 0 3]
            lappend pdblines [format "ATOM  %5d %4s %-4s%1s%4d    %8.3f%8.3f%8.3f" \
                                 $i  $an [string range $rn 0 4] $ci $ri $x $y $z]
          }
        }
        return $pdblines
      }

      proc locateCoordinates {} {
        variable pdbfile
        variable molid
        if {$molid > -1} {
          tk_messageBox -type ok -title "WARNING" \
            -message "Coordinate data from $pdbfile is already loaded."
          return
        }
        
        set tempfile [tk_getOpenFile \
          -filetypes {{"PDB files" { .pdb .PDB }} {"All files" *}}]
        if {![string equal $tempfile ""]} { 
          set pdbfile $tempfile 
        }
        [namespace current]::loadCoordinates
      }


      proc Animate {} {
        variable nframes
        variable activemode
        variable coordinates
        variable prefix
        variable lengths
        variable indices
        variable scalearrows
        variable modes
        variable animid
        variable title
        variable molid
        variable autoplay
        variable w
        variable tempfn
        variable selstr
        variable selrep
        variable betalist
        variable betamin
        variable betamax
        #puts [namespace current]
        set whichmode [lsearch $indices $activemode] 
        animate pause
        set animfn [file join $::NMWiz::tmpdir $tempfn]
        [namespace current]::loadCoordinates
        variable molid
        set sel [atomselect $molid "all"]
        $sel writepdb $animfn
        $sel delete
        
        set currentview [molinfo $molid get {rotate_matrix center_matrix scale_matrix global_matrix}]
        if {[lsearch [molinfo list] $animid] == -1} {
          set animid [mol new $animfn]
        } else {
          animate delete beg 0 end -1 skip 0 $animid
          mol addfile $animfn waitfor all $animid
        }
        mol delrep 0 $animid
        mol off $animid

        mol top $animid
        $w.draw_arrows.animbuttons_label configure -text "Animation ($animid):"
        mol rename $animid "$title mode $activemode animation"
        foreach id [molinfo list] {
          molinfo $id set {rotate_matrix center_matrix scale_matrix global_matrix} $currentview
        }
        for {set i 1} {$i <= $nframes} {incr i} {
          
        }

        variable length
        #set length [lindex $lengths [lsearch $indices $activemode]]
        set mode [vecscale [expr $length * [::tcl::mathfunc::abs $scalearrows]] [lindex $modes [lsearch $indices $activemode]]]
        set coords [vecadd $coordinates $mode]
        set mode [::NMWiz::array2xyz [vecscale $mode [expr  -2.0 / $nframes]]]
        set sel [atomselect $animid "all"]
        $sel set {x y z} [::NMWiz::array2xyz $coords]
        for {set i 1} {$i <= $nframes} {incr i} {
          animate dup frame [expr $i - 1] $animid
          $sel frame $i
          $sel lmoveby $mode

        }
        $sel delete
        [namespace current]::updateProtRep $animid 
        mol off $animid
        if {$autoplay} {
          animate speed 0.96
          animate style rock
          animate forward
        }
        mol on $animid
        if {$selrep} {
          mol modselect 0 $animid $selstr
        }
        eval "\$[namespace current]::w.draw_arrows.animbuttons_showhide configure -text Hide"
        eval "\$[namespace current]::w.draw_arrows.animbuttons_stop configure -text Pause"
        set [namespace current]::stopped 0
        
        variable animidlist
        lset animidlist $whichmode $animid

      }


      proc loadCoordinates {} {
        variable molid
        if {[lsearch [molinfo list] $molid] != -1} {
          return 0
        }
        variable coordinates
        variable title
        variable w
        variable tempfn        
        set outfile [open [file join $::NMWiz::tmpdir $tempfn] w]
        foreach line [[namespace current]::getPDBLines $coordinates] {
          puts $outfile $line
        } 
        close $outfile
        set preserve 0
        if {[molinfo num] > 0 && $::NMWiz::preserview} {
          set currentview [molinfo [lindex [molinfo list] 0] get {rotate_matrix center_matrix scale_matrix global_matrix}]
          set preserve 1
        }

        set molid [mol new [file join $::NMWiz::tmpdir $tempfn]]

        $w.draw_arrows.protbuttons_label configure -text "Molecule ($molid):"
        mol rename $molid "$title coordinates"
        [namespace current]::calcMSF

        if {$preserve} {
          foreach id [molinfo list] {
            molinfo $id set {rotate_matrix center_matrix scale_matrix global_matrix} $currentview
          }
        }
      }

      proc checkCoordinates {} {
        variable molid
        variable modes
        variable coordinates
        set sel [atomselect $molid "all"]
        if {[$sel num] != [llength [lindex $modes 0]]} {
          mol delete $molid
          set molid -1
          if {"ok" == [tk_messageBox -type okcancel -title "ERROR" \
              -message "[[atomselect $molid all] num] atoms are loaded. Coordinate data file must contain [llength [lindex $modes 0]] atoms. Please locate the correct file."]} {
            [namespace current]::locateCoordinates
          } 
        } else {
          set $all [atomselect $molid all]
          set coordinates [$all get {x y z}]
          $all delete
        }
        $sel delete
      }
      
      proc changeColor {} {
        variable color
        variable colorlist
        variable indices
        variable activemode
        lset colorlist [lsearch $indices $activemode] $color
      }
      
      proc prevMode {} {
        variable activeindex
        if {$activeindex > 0} {
          variable activemode
          variable indices
          set activemode [lindex $indices [expr $activeindex - 1]];
          [namespace current]::changeMode;
        }
      }
      
      proc nextMode {} {
        variable activeindex
        variable numofmodes
        if {$activeindex < [expr $numofmodes - 1]} {
          variable activemode
          variable indices
          set activemode [lindex $indices [expr $activeindex + 1]];
          [namespace current]::changeMode;
        }
      }
    
      proc changeMode {} {
        [namespace current]::loadCoordinates
        variable w
        variable activemode
        variable indices
        variable arrid
        variable arridlist
        variable animid
        variable animidlist
        variable color
        variable colorlist
        variable drawlengthstr
        variable lengths
        variable activeindex 
        variable ndim
        set inactiveindex $activeindex 
        set activeindex [lsearch $indices $activemode]
        variable length [lindex $lengths $activeindex]
        if {$ndim == 3} {
          variable scalearrows_list
          variable rmsd_list
          variable hide_shorter_list 
          variable cylinder_radius_list
          variable cone_radius_list
          variable cone_height_list
          variable resolution_list
          variable material_list
          
          variable scalearrows
          variable hide_shorter 
          variable cylinder_radius
          variable cone_radius
          variable cone_height
          variable resolution
          variable material
          variable rmsd
          
          lset scalearrows_list $inactiveindex $scalearrows
          lset hide_shorter_list $inactiveindex $hide_shorter  
          lset cylinder_radius_list $inactiveindex $cylinder_radius
          lset cone_radius_list $inactiveindex $cone_radius
          lset cone_height_list $inactiveindex $cone_height
          lset resolution_list $inactiveindex $resolution
          lset material_list $inactiveindex $material
          lset rmsd_list $inactiveindex $rmsd

          set drawlengthstr [format "%.1f" [lindex $lengths $activeindex]];
          
          set scalearrows [lindex $scalearrows_list $activeindex]
          set hide_shorter [lindex $hide_shorter_list $activeindex]
          set cylinder_radius [lindex $cylinder_radius_list $activeindex]
          set cone_radius [lindex $cone_radius_list $activeindex]
          set cone_height [lindex $cone_height_list $activeindex]
          set resolution [lindex $resolution_list $activeindex]
          set material [lindex $material_list $activeindex]
          set rmsd [lindex $rmsd_list $activeindex]

      
          variable hideprev
          if {$hideprev} {
            if {$arrid > -1 && [lsearch [molinfo list] $arrid] > -1} {
              mol off $arrid
            }
          }
          if {$animid > -1 && [lsearch [molinfo list] $animid] > -1} {
            mol off $animid
          }
          
          set which [lsearch $indices $activemode]
          set arrid [lindex $arridlist $which]      
          set animid [lindex $animidlist $which]

          
          set color [lindex $colorlist $which]

          if {$arrid > -1 && [lsearch [molinfo list] $arrid] > -1} {
            mol on $arrid
          } else {
            [namespace current]::drawArrows
          }
          [namespace current]::calcMSF
          
          $w.draw_arrows.arrowbuttons_showhide configure -text Hide

          if {$animid > -1 && [lsearch [molinfo list] $animid] > -1} {
            mol on $animid
            mol top $animid
          } else {
            variable autoanimate  
            if {$autoanimate} {
              [namespace current]::Animate
            }
          }
          $w.draw_arrows.animbuttons_showhide configure -text Hide
          $w.draw_arrows.animbuttons_stop configure -text Play
        } else {
          [namespace current]::calcMSF          
        }
        variable activeindexprev $activeindex
        
      }
      
      #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
      
      proc nmwizgui {} {
        variable w
        set ns [namespace current]
        if [winfo exists .[string range $ns 2 end]] {
          wm deiconify .[string range $ns 2 end]
          raise .[string range $ns 2 end]
          return
        }
        set w [toplevel .[string range $ns 2 end]]
        variable title

        wm title $w "NMWiz - $title"
        wm resizable $w 0 0
        set wam [labelframe $w.active_mode -text "$title" -bd 2] 

        variable ndim
        grid [label $wam.active_label -text "Active mode:"] \
            -row 0 -column 1 -sticky w

        grid [frame $wam.active] \
            -row 0 -column 2 -columnspan 3 -sticky w

        tk_optionMenu $wam.active.list ${ns}::activemode 0
        $wam.active.list.menu delete 0
        variable indices
        variable lengths
        foreach index $indices length $lengths {
          $wam.active.list.menu add radiobutton -label $index \
              -variable ${ns}::activemode \
              -command "${ns}::changeMode;"
        }
        
        button $wam.active.prev -text "<=" -command "${ns}::prevMode"
        button $wam.active.negate -text "+/-" -command "set ${ns}::scalearrows \[expr - \$${ns}::scalearrows]; ${ns}::autoUpdate"  
        button $wam.active.next -text "=>" -command "${ns}::nextMode"
        
        variable ndim
        if {$ndim == 3} {
          tk_optionMenu $wam.active.color ${ns}::color "blue"
          $wam.active.color.menu delete 0
          foreach acolor "blue red gray orange yellow tan green white pink \
        cyan purple black yellow2 yellow3 green2 green3 \
        cyan2 cyan3 blue2 blue3 violet magenta magenta2 red2 red3 orange2 \
        orange3" {
            $wam.active.color.menu add radiobutton -label $acolor \
                -variable ${ns}::color \
                -command "${ns}::changeColor; ${ns}::autoUpdate; "
          }
          pack $wam.active.list $wam.active.prev $wam.active.negate $wam.active.next $wam.active.color -side left -anchor w -fill x
        } else {
          tk_optionMenu $wam.active.color ${ns}::msformode "Mobility"
          $wam.active.color.menu delete 0
          foreach acolor "Mobility Eigenvector" {
            $wam.active.color.menu add radiobutton -label $acolor \
                -variable ${ns}::msformode \
                -command "${ns}::calcMSF;"
          }
          pack $wam.active.list $wam.active.prev $wam.active.next $wam.active.color -side left -anchor w -fill x
        }
        
        #blue red gray orange yellow tan silver green white pink cyan purple lime mauve ochre iceblue black yellow2 yellow3 green2 green3 cyan2 cyan3 blue2 blue3 violet violet2 magenta magenta2 red2 red3 orange2 orange3

        if {$ndim == 3} {
          
          grid [label $wam.scale_label -text "Scale by:"] \
            -row 1 -column 1 -sticky w
          grid [frame $wam.scale_frame] \
            -row 1 -column 2 -columnspan 3 -sticky w
          entry $wam.scale_frame.entry -width 4 -state readonly -textvariable ${ns}::length
          label $wam.scale_frame.times -text "x"
          entry $wam.scale_frame.scalar -width 6 -textvariable ${ns}::scalearrows
          button $wam.scale_frame.incr1 -text "+1" -command \
            "set ${ns}::scalearrows \[expr \$${ns}::scalearrows + 1]; ${ns}::autoUpdate"
          button $wam.scale_frame.incr5 -text "+5" -command \
            "set ${ns}::scalearrows \[expr \$${ns}::scalearrows + 5]; ${ns}::autoUpdate"
          button $wam.scale_frame.decr5 -text "-5" -command \
            "set ${ns}::scalearrows \[expr \$${ns}::scalearrows - 5]; ${ns}::autoUpdate"
          button $wam.scale_frame.decr1 -text "-1" -command \
            "set ${ns}::scalearrows \[expr \$${ns}::scalearrows - 1]; ${ns}::autoUpdate"
          pack $wam.scale_frame.entry $wam.scale_frame.times $wam.scale_frame.scalar \
            $wam.scale_frame.incr1 $wam.scale_frame.incr5 \
            $wam.scale_frame.decr5 $wam.scale_frame.decr1 \
            -side left -anchor w -fill x
            
          grid [label $wam.rmsd_label -text "RMSD (A):"] \
            -row 2 -column 1 -sticky w
          grid [frame $wam.adjust_frame] \
            -row 2 -column 2 -columnspan 3 -sticky w
          entry $wam.adjust_frame.entry -width 4 -textvariable ${ns}::rmsd
          button $wam.adjust_frame.incr1 -text "+0.1" -command \
            "set ${ns}::rmsd \[expr \$${ns}::rmsd + 0.1]; ${ns}::autoUpdate"
          button $wam.adjust_frame.incr5 -text "+0.5" -command \
            "set ${ns}::rmsd \[expr \$${ns}::rmsd + 0.5]; ${ns}::autoUpdate"
          button $wam.adjust_frame.decr5 -text "-0.5" -command \
            "set ${ns}::rmsd \[expr \$${ns}::rmsd - 0.5]; ${ns}::autoUpdate"
          button $wam.adjust_frame.decr1 -text "-0.1" -command \
            "set ${ns}::rmsd \[expr \$${ns}::rmsd - 0.1]; ${ns}::autoUpdate"
          pack $wam.adjust_frame.entry \
            $wam.adjust_frame.incr1 $wam.adjust_frame.incr5 \
            $wam.adjust_frame.decr5 $wam.adjust_frame.decr1 \
            -side left -anchor w -fill x

          grid [label $wam.selstr_label -text "Selection:"] \
            -row 5 -column 1 -sticky w
          grid [entry $wam.selstr_entry \
            -textvariable ${ns}::selstr] \
            -row 5 -column 2 -columnspan 2 -sticky we
          grid [button $wam.selstr_draw -width 4 -text "Redraw" \
              -command ${ns}::drawArrows] \
            -row 5 -column 4 -sticky we
        }
        grid [button $wam.showmain -text "Main" \
            -command nmwiz_tk] \
          -row 8 -column 1 -sticky we
        grid [button $wam.save -text "Save" \
            -command "::NMWiz::writeNMD $ns"] \
          -row 8 -column 2 -sticky we
        grid [button $wam.remove -text "Remove" \
            -command "lset ::NMWiz::titles $::NMWiz::guicount NONE; pack forget .nmwizgui.{[string range $ns 2 end]}frame; ${ns}::deleteMolecules; namespace delete $ns; destroy .[string range $ns 2 end]"] \
          -row 8 -column 3 -sticky we
        grid [button $wam.showhelp -text "Help" \
            -command {::NMWiz::showHelp wizard}] \
          -row 8 -column 4 -sticky we

        pack $wam -side top -ipadx 10 -ipady 5 -fill x -expand 1

        set wda [labelframe $w.draw_arrows -text "Actions" -bd 2]
        
        if {$ndim == 3} {
          grid [label $wda.arrowbuttons_label -text "Mode:"] \
            -row 5 -column 0 -sticky w
          grid [button $wda.arrowbuttons_draw -text "Draw" \
              -command ${ns}::drawArrows] \
            -row 5 -column 2 -sticky ew
          grid [button $wda.arrowbuttons_clean -text "Clean" \
              -command "if {\[lsearch \[molinfo list] \$${ns}::arrid] != -1} {graphics \$${ns}::arrid delete all}"] \
            -row 5 -column 3
          grid [button $wda.arrowbuttons_showhide -text "Hide" \
              -command "if {\[lsearch \[molinfo list] \$${ns}::arrid] != 0} { if {\[molinfo \$${ns}::arrid get displayed]} {mol off \$${ns}::arrid;\
                        \$${ns}::w.draw_arrows.arrowbuttons_showhide configure -text Show} else {mol on \$${ns}::arrid;\
                        \$${ns}::w.draw_arrows.arrowbuttons_showhide configure -text Hide} }"] \
            -row 5 -column 4 -sticky ew
          grid [button $wda.arrowbuttons_options -text "Options" \
              -command "if {\$${ns}::arropt} {pack forget \$${ns}::w.graphics_options;\
                        set ${ns}::arropt 0; \$${ns}::w.draw_arrows.arrowbuttons_options configure -relief raised} else {pack \$${ns}::w.graphics_options -side top -ipadx 10 -ipady 5 -fill x -expand 1;\
                        set ${ns}::arropt 1; \$${ns}::w.draw_arrows.arrowbuttons_options configure -relief sunken}"] \
            -row 5 -column 5 -sticky ew

          grid [label $wda.animbuttons_label -text "Animation:"] \
            -row 6 -column 0 -sticky w
          grid [button $wda.animbuttons_animate -text "Make" \
              -command ${ns}::Animate] \
            -row 6 -column 2 -sticky ew
          grid [button $wda.animbuttons_stop -text "Play" \
              -command "if {\$${ns}::animid == -1 || \[lsearch \[molinfo list] \$${ns}::animid] == -1} {${ns}::Animate} else {if {\$${ns}::stopped} {mol top \$${ns}::animid; animate forward; \$${ns}::w.draw_arrows.animbuttons_stop configure -text Pause; set ${ns}::stopped 0} else {animate pause; \$${ns}::w.draw_arrows.animbuttons_stop configure -text Play; set ${ns}::stopped 1}}"] \
            -row 6 -column 3 -sticky ew
          grid [button $wda.animbuttons_showhide -text "Hide" \
              -command "if {\$${ns}::animid > -1 && \[lsearch \[molinfo list] \$${ns}::animid] > -1} {if {\[molinfo \$${ns}::animid get displayed]} {animate pause; mol off \$${ns}::animid; \$${ns}::w.draw_arrows.animbuttons_showhide configure -text Show} else {mol on \$${ns}::animid; \$${ns}::w.draw_arrows.animbuttons_showhide configure -text Hide; animate forward}}"] \
            -row 6 -column 4 -sticky ew
          grid [button $wda.animbuttons_options -text "Options" \
              -command "if {\$${ns}::anmopt} {pack forget \$${ns}::w.animation_options; set ${ns}::anmopt 0; \$${ns}::w.draw_arrows.animbuttons_options configure -relief raised} else {pack \$${ns}::w.animation_options -side top -ipadx 10 -ipady 5 -fill x -expand 1; set ${ns}::anmopt 1; \$${ns}::w.draw_arrows.animbuttons_options configure -relief sunken}"] \
            -row 6 -column 5 -sticky ew
        }

        grid [label $wda.plot_label -text "Plotting:"] \
          -row 8 -column 0 -sticky w
        grid [button $wda.plot_plot -text "Plot" \
            -command "${ns}::Plot ${ns}"] \
          -row 8 -column 2 -sticky ew
        grid [button $wda.plot_clear -text "Clear" \
            -command "${ns}::clearSelection"] \
          -row 8 -column 3 -sticky ew
        grid [button $wda.plot_showhide -text "Hide" \
            -command "if {\$${ns}::selid > -1 && \[lsearch \[molinfo list] \$${ns}::selid] > -1} {if {\[molinfo \$${ns}::selid get displayed]} {mol off \$${ns}::selid; \$${ns}::w.draw_arrows.plot_showhide configure -text Show} else {mol on \$${ns}::selid; \$${ns}::w.draw_arrows.plot_showhide configure -text Hide}}"] \
          -row 8 -column 4 -sticky ew
        grid [button $wda.plot_options -text "Options" \
            -command "if {\$${ns}::pltopt} {pack forget \$${ns}::w.plotting_options; set ${ns}::pltopt 0; \$${ns}::w.draw_arrows.plot_options configure -relief raised} else {pack \$${ns}::w.plotting_options -side top -ipadx 10 -ipady 5 -fill x -expand 1; set ${ns}::pltopt 1; \$${ns}::w.draw_arrows.plot_options configure -relief sunken}"] \
          -row 8 -column 5 -sticky ew
         
        ##-command "if {\$${ns}::pltopt} {pack forget \$${ns}::w.animation_options; set ${ns}::pltopt 0; \$${ns}::w.draw_arrows.plot_options configure -relief raised} else {pack \$${ns}::w.plot_options -side top -ipadx 10 -ipady 5 -fill x -expand 1; set ${ns}::pltopt 1; \$${ns}::w.draw_arrows.plot_options configure -relief sunken}"] \

        grid [label $wda.protbuttons_label -text "Molecule:"] \
          -row 9 -column 0 -sticky w
        grid [button $wda.prt_update -text "Update" \
            -command "${ns}::loadCoordinates; ${ns}::updateProtRep \$${ns}::molid"] \
          -row 9 -column 2 -sticky ew
        grid [button $wda.protbuttons_focus -text "Focus" \
            -command "${ns}::loadCoordinates; mol top \$${ns}::molid; display resetview"] \
          -row 9 -column 3  -sticky ew
        grid [button $wda.protbuttons_showhide -text "Hide" \
            -command "${ns}::loadCoordinates; if {\[molinfo \$${ns}::molid get displayed]} {mol off \$${ns}::molid; \$${ns}::w.draw_arrows.protbuttons_showhide configure -text Show;} else {mol on \$${ns}::molid; \$${ns}::w.draw_arrows.protbuttons_showhide configure -text Hide;}"] \
          -row 9 -column 4 -sticky ew
        grid [button $wda.protbuttons_repoptions -text "Options" \
            -command "if {\$${ns}::prtopt} {pack forget \$${ns}::w.prograph_options; set ${ns}::prtopt 0; \$${ns}::w.draw_arrows.protbuttons_repoptions configure -relief raised} else {pack \$${ns}::w.prograph_options -side top -ipadx 10 -ipady 5 -fill x -expand 1; set ${ns}::prtopt 1; \$${ns}::w.draw_arrows.protbuttons_repoptions configure -relief sunken}"] \
          -row 9 -column 5 -sticky ew

        pack $wda -side top -fill x -expand 1

        set wgo [labelframe $w.graphics_options -text "Mode Graphics Options" -bd 2]
        
        grid [checkbutton $wgo.auto_check -text " auto update graphics" \
            -variable ${ns}::autoupdate] \
          -row 0 -column 1 -sticky w

        grid [checkbutton $wgo.hideprev_check -text " auto hide inactive mode" \
            -variable ${ns}::hideprev] \
          -row 0 -column 2 -sticky w

        grid [checkbutton $wgo.both_check -text " draw in both directions" \
            -variable ${ns}::bothdirections] \
          -row 7 -column 1 -sticky w

        grid [checkbutton $wgo.porcupine_check -text " use porcupine style" \
            -variable ${ns}::porcupine] \
          -row 7 -column 2 -sticky w

        grid [label $wgo.hide_label -text "Draw if longer than:"] \
          -row 9 -column 1 -sticky w
        grid [frame $wgo.hide_frame] \
          -row 9 -column 2 -sticky w
        entry $wgo.hide_frame.entry -width 4 -textvariable ${ns}::hide_shorter
        button $wgo.hide_frame.decr -text "-0.5" \
          -command "set ${ns}::hide_shorter \[::tcl::mathfunc::abs \[expr \$${ns}::hide_shorter - 0.5]]; ${ns}::autoUpdate"
        button $wgo.hide_frame.incr -text "+0.5" \
          -command "set ${ns}::hide_shorter \[::tcl::mathfunc::abs \[expr \$${ns}::hide_shorter + 0.5]]; ${ns}::autoUpdate"
        label $wgo.hide_frame.angstrom -text "A"
        pack $wgo.hide_frame.entry $wgo.hide_frame.decr $wgo.hide_frame.incr \
          $wgo.hide_frame.angstrom -side left -anchor w -fill x

        grid [label $wgo.cylinder_label -text "Arrow cylinder radius:"] \
          -row 10 -column 1 -sticky w
        grid [frame $wgo.cylinder_frame] \
          -row 10 -column 2 -sticky w
        entry $wgo.cylinder_frame.entry -width 4 -textvariable ${ns}::cylinder_radius
        button $wgo.cylinder_frame.decr -text "-0.1" \
          -command "set ${ns}::cylinder_radius \[::tcl::mathfunc::abs \[expr \$${ns}::cylinder_radius - 0.1]]; ${ns}::autoUpdate"
        button $wgo.cylinder_frame.incr -text "+0.1" \
          -command "set ${ns}::cylinder_radius \[::tcl::mathfunc::abs \[expr \$${ns}::cylinder_radius + 0.1]]; ${ns}::autoUpdate"
        label $wgo.cylinder_frame.angstrom -text "A"
        pack $wgo.cylinder_frame.entry $wgo.cylinder_frame.decr \
          $wgo.cylinder_frame.incr $wgo.cylinder_frame.angstrom \
          -side left -anchor w -fill x

        grid [label $wgo.coner_label -text "Arrow/porcupine cone radius:"] \
          -row 11 -column 1 -sticky w
        grid [frame $wgo.coner_frame] \
          -row 11 -column 2 -sticky w
        entry $wgo.coner_frame.entry -width 4 -textvariable ${ns}::cone_radius
        button $wgo.coner_frame.decr -text "-0.1" \
          -command "set ${ns}::cone_radius \[::tcl::mathfunc::abs \[expr \$${ns}::cone_radius - 0.1]]; ${ns}::autoUpdate"
        button $wgo.coner_frame.incr -text "+0.1" \
          -command "set ${ns}::cone_radius \[::tcl::mathfunc::abs \[expr \$${ns}::cone_radius + 0.1]]; ${ns}::autoUpdate"
        label $wgo.coner_frame.angstrom -text "A"
        pack $wgo.coner_frame.entry $wgo.coner_frame.decr $wgo.coner_frame.incr \
          $wgo.coner_frame.angstrom -side left -anchor w -fill x

        grid [label $wgo.coneh_label -text "Arrow cone height:"] \
          -row 12 -column 1 -sticky w
        grid [frame $wgo.coneh_frame] \
          -row 12 -column 2 -sticky w
        entry $wgo.coneh_frame.entry -width 4 -textvariable ${ns}::cone_height
        button $wgo.coneh_frame.decr -text "-0.2" \
          -command "set ${ns}::cone_height \[::tcl::mathfunc::abs \[expr \$${ns}::cone_height - 0.2]]; ${ns}::autoUpdate"
        button $wgo.coneh_frame.incr -text "+0.2" \
          -command "set ${ns}::cone_height \[::tcl::mathfunc::abs \[expr \$${ns}::cone_height + 0.2]]; ${ns}::autoUpdate"
        label $wgo.coneh_frame.angstrom -text "A"
        pack $wgo.coneh_frame.entry $wgo.coneh_frame.decr $wgo.coneh_frame.incr \
          $wgo.coneh_frame.angstrom -side left -anchor w -fill x

        grid [label $wgo.material_label -text "Graphics material:"] \
          -row 20 -column 1 -sticky w
        grid [frame $wgo.material_frame] \
          -row 20 -column 2 -sticky w
        tk_optionMenu $wgo.material_frame.list ${ns}::material "Opaque"
        $wgo.material_frame.list.menu delete 0
        foreach mtrl "Opaque Transparent BrushedMetal Diffuse Ghost Glass1 Glass2 Glass3 Glossy HardPlastic MetallicPastel Steel Translucent Edgy EdgyShiny EdgyGlass Goodsell AOShiny AOChalky AOEdgy" {
          $wgo.material_frame.list.menu add radiobutton -label $mtrl \
              -variable ${ns}::material \
              -command "${ns}::autoUpdate"    
        }
        pack $wgo.material_frame.list -side left -anchor w -fill x

        grid [label $wgo.resolution_label -text "Graphics resolution:"] \
          -row 21 -column 1 -sticky w
        grid [frame $wgo.resolution_frame] \
          -row 21 -column 2 -sticky w
        tk_optionMenu $wgo.resolution_frame.list ${ns}::resolution 6 
        $wgo.resolution_frame.list.menu delete 0
        foreach resol "6 10 15 20 25 30 35 40 45 50" {
          $wgo.resolution_frame.list.menu add radiobutton -label $resol \
              -variable ${ns}::resolution \
              -command "${ns}::autoUpdate"  
        } 
        pack $wgo.resolution_frame.list -side left -anchor w -fill x
        
        grid [label $wgo.copyfrom_label -text "Copy settings from mode:"] \
          -row 22 -column 1 -sticky w
        grid [frame $wgo.copyfrom] \
          -row 22 -column 2 -sticky w
        tk_optionMenu $wgo.copyfrom.list ${ns}::copyfrom 0
        $wgo.copyfrom.list.menu delete 0
        variable indices
        variable lengths
        foreach index $indices length $lengths {
          $wgo.copyfrom.list.menu add radiobutton -label $index \
              -variable ${ns}::copyfrom
        }
        button $wgo.copyfrom.copy -text "Copy" -command "${ns}::copySettings"
        pack $wgo.copyfrom.list $wgo.copyfrom.copy -side left -anchor w -fill x
        

        set wpgo [labelframe $w.prograph_options -text "Molecule Representations" -bd 2]
        
        grid [checkbutton $wpgo.selstr_check -text " show only selected atoms" \
            -variable ${ns}::selrep -command "${ns}::autoUpdate"] \
          -row 0 -column 1 -columnspan 2 -sticky w
        
        grid [label $wpgo.protas_label -text "Representation:"] \
          -row 13 -column 1 -sticky w
        grid [frame $wpgo.protas_frame] \
          -row 13 -column 2 -sticky w
        tk_optionMenu $wpgo.protas_frame.list ${ns}::showproteinas "Tube"
        $wpgo.protas_frame.list.menu delete 0
        $wpgo.protas_frame.list.menu add radiobutton -label "Custom" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "Network" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "Lines" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "VDW" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "CPK" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "Licorice" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "Trace" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "Tube" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "Ribbons" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "NewRibbons" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "Cartoon" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "NewCartoon" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        pack $wpgo.protas_frame.list -side left -anchor w -fill x

        grid [label $wpgo.procolor_label -text "Color scheme:"] \
          -row 14 -column 1 -sticky w
        grid [frame $wpgo.procolor_frame] \
          -row 14 -column 2 -sticky w
        tk_optionMenu $wpgo.procolor_frame.list ${ns}::proteincolor "Mobility"
        $wpgo.procolor_frame.list.menu delete 0
        $wpgo.procolor_frame.list.menu add radiobutton -label "Mobility" -variable ${ns}::proteincolor -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.procolor_frame.list.menu add radiobutton -label "Bfactors" -variable ${ns}::proteincolor -command "${ns}::updateProtRep \$${ns}::molid"
        foreach acolor "Index Chain ResName ResType" {
          $wpgo.procolor_frame.list.menu add radiobutton -label $acolor \
              -variable ${ns}::proteincolor \
              -command "${ns}::updateProtRep \$${ns}::molid"
        }
        pack $wpgo.procolor_frame.list -side left -anchor w -fill x


        grid [label $wpgo.cutoffdistance_label -text "Network cutoff distance:"] \
          -row 16 -column 1 -sticky w
        grid [frame $wpgo.cutoffdistance_frame] \
          -row 16 -column 2 -sticky w
        entry $wpgo.cutoffdistance_frame.entry -width 4 -textvariable ${ns}::cutoffdistance
        button $wpgo.cutoffdistance_frame.decr -text "-1.0" \
          -command "set ${ns}::cutoffdistance \[::tcl::mathfunc::abs \[expr \$${ns}::cutoffdistance - 1.0]]; ${ns}::updateProtRep \$${ns}::molid"
        button $wpgo.cutoffdistance_frame.incr -text "+1.0" \
          -command "set ${ns}::cutoffdistance \[::tcl::mathfunc::abs \[expr \$${ns}::cutoffdistance + 1.0]]; ${ns}::updateProtRep \$${ns}::molid"
        label $wpgo.cutoffdistance_frame.angstrom -text "A"
        pack $wpgo.cutoffdistance_frame.entry $wpgo.cutoffdistance_frame.decr $wpgo.cutoffdistance_frame.incr \
          $wpgo.cutoffdistance_frame.angstrom -side left -anchor w -fill x

        grid [label $wpgo.nodescale_label -text "Scale node spheres:"] \
          -row 17 -column 1 -sticky w
        grid [frame $wpgo.nodescale_frame] \
          -row 17 -column 2 -sticky w
        entry $wpgo.nodescale_frame.entry -width 4 -textvariable ${ns}::spherescale
        button $wpgo.nodescale_frame.decr -text "-0.1" \
          -command "set ${ns}::spherescale \[::tcl::mathfunc::abs \[expr \$${ns}::spherescale - 0.1]]; ${ns}::updateProtRep \$${ns}::molid"
        button $wpgo.nodescale_frame.incr -text "+0.1" \
          -command "set ${ns}::spherescale \[::tcl::mathfunc::abs \[expr \$${ns}::spherescale + 0.1]]; ${ns}::updateProtRep \$${ns}::molid"
        pack $wpgo.nodescale_frame.entry $wpgo.nodescale_frame.decr $wpgo.nodescale_frame.incr \
          -side left -anchor w -fill x
          
        grid [label $wpgo.bondradius_label -text "Network bond radius:"] \
          -row 18 -column 1 -sticky w
        grid [frame $wpgo.bondradius_frame] \
          -row 18 -column 2 -sticky w
        entry $wpgo.bondradius_frame.entry -width 4 -textvariable ${ns}::bondradius
        button $wpgo.bondradius_frame.decr -text "-0.1" \
          -command "set ${ns}::bondradius \[::tcl::mathfunc::abs \[expr \$${ns}::bondradius - 0.1]]; ${ns}::updateProtRep \$${ns}::molid"
        button $wpgo.bondradius_frame.incr -text "+0.1" \
          -command "set ${ns}::bondradius \[::tcl::mathfunc::abs \[expr \$${ns}::bondradius + 0.1]]; ${ns}::updateProtRep \$${ns}::molid"
        label $wpgo.bondradius_frame.angstrom -text "A"
        pack $wpgo.bondradius_frame.entry $wpgo.bondradius_frame.decr $wpgo.bondradius_frame.incr \
          $wpgo.bondradius_frame.angstrom -side left -anchor w -fill x


        grid [label $wpgo.tuberadius_label -text "Tube/licorice radius:"] \
          -row 19 -column 1 -sticky w
        grid [frame $wpgo.tuberadius_frame] \
          -row 19 -column 2 -sticky w
        entry $wpgo.tuberadius_frame.entry -width 4 -textvariable ${ns}::tuberadius
        button $wpgo.tuberadius_frame.decr -text "-0.1" \
          -command "set ${ns}::tuberadius \[::tcl::mathfunc::abs \[expr \$${ns}::tuberadius - 0.1]]; ${ns}::updateProtRep \$${ns}::molid"
        button $wpgo.tuberadius_frame.incr -text "+0.1" \
          -command "set ${ns}::tuberadius \[::tcl::mathfunc::abs \[expr \$${ns}::tuberadius + 0.1]]; ${ns}::updateProtRep \$${ns}::molid"
        label $wpgo.tuberadius_frame.angstrom -text "A"
        pack $wpgo.tuberadius_frame.entry $wpgo.tuberadius_frame.decr $wpgo.tuberadius_frame.incr \
          $wpgo.tuberadius_frame.angstrom -side left -anchor w -fill x


        grid [label $wpgo.material_label -text "Graphics material:"] \
          -row 20 -column 1 -sticky w
        grid [frame $wpgo.material_frame] \
          -row 20 -column 2 -sticky w
        tk_optionMenu $wpgo.material_frame.list ${ns}::material_protein "Opaque"
        $wpgo.material_frame.list.menu delete 0
        foreach mtrl "Opaque Transparent BrushedMetal Diffuse Ghost Glass1 Glass2 Glass3 Glossy HardPlastic MetallicPastel Steel Translucent Edgy EdgyShiny EdgyGlass Goodsell AOShiny AOChalky AOEdgy" {
          $wpgo.material_frame.list.menu add radiobutton -label $mtrl \
              -variable ${ns}::material_protein \
              -command "${ns}::updateProtRep \$${ns}::molid"    
        }
        pack $wpgo.material_frame.list -side left -anchor w -fill x

        grid [label $wpgo.resolution_label -text "Graphics resolution:"] \
          -row 21 -column 1 -sticky w
        grid [frame $wpgo.resolution_frame] \
          -row 21 -column 2 -sticky w
        tk_optionMenu $wpgo.resolution_frame.list ${ns}::resolution_protein 6 
        $wpgo.resolution_frame.list.menu delete 0
        foreach resol "6 10 15 20 25 30 35 40 45 50" {
          $wpgo.resolution_frame.list.menu add radiobutton -label $resol \
              -variable ${ns}::resolution_protein \
              -command "${ns}::updateProtRep \$${ns}::molid"  
        } 
        pack $wpgo.resolution_frame.list -side left -anchor w -fill x

        grid [label $wpgo.csm_label -text "Color scale method:"] \
          -row 22 -column 1 -sticky w
        grid [frame $wpgo.csm_frame] \
          -row 22 -column 2 -sticky w
        tk_optionMenu $wpgo.csm_frame.list ${ns}::color_scale_method [colorinfo scale method]
        $wpgo.csm_frame.list.menu delete 0
        foreach mtrl "RWB BWR RGryB BGryR RGB BGR RWG GWR GWB BWG BlkW WBlk" {
          $wpgo.csm_frame.list.menu add radiobutton -label $mtrl \
              -variable ${ns}::color_scale_method \
              -command "color scale method $mtrl;"
        }
        pack $wpgo.csm_frame.list -side left -anchor w -fill x
        
        grid [label $wpgo.csmp_label -text "Color scale midpoint:"] \
          -row 23 -column 1 -sticky w
        grid [frame $wpgo.csmp_frame] \
          -row 23 -column 2 -sticky w
        tk_optionMenu $wpgo.csmp_frame.list ${ns}::color_scale_midpoint 0.5 
        $wpgo.csmp_frame.list.menu delete 0
        foreach mtrl "0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0" {
          $wpgo.csmp_frame.list.menu add radiobutton -label $mtrl \
              -variable ${ns}::color_scale_midpoint \
              -command "color scale midpoint $mtrl;"
        }
        pack $wpgo.csmp_frame.list -side left -anchor w -fill x

        set wao [labelframe $w.animation_options -text "Animation Options" -bd 2]
        
        grid [checkbutton $wao.auto_check -text " auto animate" \
            -variable ${ns}::autoanimate] \
          -row 0 -column 1 -columnspan 2 -sticky w
        
        grid [checkbutton $wao.autoplay_check -text " continuous autoplay" -variable ${ns}::autoplay] \
          -row 0 -column 3 -columnspan 2 -sticky w
          
        grid [label $wao.nframes_label -text "Number of frames:"] \
          -row 9 -column 1 -sticky w
        grid [frame $wao.nframes_frame] \
          -row 9 -column 2 -columnspan 3 -sticky w
        entry $wao.nframes_frame.entry -width 4 -textvariable ${ns}::nframes
        button $wao.nframes_frame.decr5 -text "-5" \
          -command "set ${ns}::nframes \[expr \$${ns}::nframes - 5]"
        button $wao.nframes_frame.incr5 -text "+5" \
          -command "set ${ns}::nframes \[expr \$${ns}::nframes + 5]"
        pack $wao.nframes_frame.entry $wao.nframes_frame.decr5 \
          $wao.nframes_frame.incr5 -side left -anchor w -fill x

        set wpo [labelframe $w.plotting_options -text "Plotting Options" -bd 2]

        grid [checkbutton $wpo.overplot_check -text " add new plot to most recent MultiPlot" \
            -variable ${ns}::overplot] \
          -row 0 -column 1 -columnspan 2 -sticky w
  
        if {$ndim == 1} {        
          grid [label $wpo.plotcolor_label -text "Line color:"] \
            -row 0 -column 5 -sticky w
            
          grid [frame $wpo.pcf] \
            -row 0 -column 6 -sticky w
          tk_optionMenu $wpo.pcf.color ${ns}::color "blue"
            $wpo.pcf.color.menu delete 0
            foreach acolor "blue red gray orange yellow tan green white pink \
          cyan purple black yellow2 yellow3 green2 green3 \
          cyan2 cyan3 blue2 blue3 violet magenta magenta2 red2 red3 orange2 \
          orange3" {
              $wpo.pcf.color.menu add radiobutton -label $acolor \
                  -variable ${ns}::color
            }
          pack $wpo.pcf.color -side left -anchor w -fill x
        }
        grid [label $wpo.plotwidth_label -text "Plot width:"] \
          -row 1 -column 1 -sticky w
        grid [entry $wpo.plotwidth_entry -width 4 -textvariable ${ns}::plotwidth] \
          -row 1 -column 2 -sticky w
        grid [label $wpo.spacing_label -text "  "] \
          -row 1 -column 3 -sticky w

        grid [label $wpo.plotheight_label -text "Plot height:"] \
          -row 1 -column 5 -sticky w
        grid [entry $wpo.plotheight_entry -width 4 -textvariable ${ns}::plotheight] \
          -row 1 -column 6 -sticky w

        grid [label $wpo.line_label -text "Lines:"] \
          -row 3 -column 1 -sticky w
        grid [frame $wpo.line_frame] \
          -row 3 -column 2 -sticky w
        tk_optionMenu $wpo.line_frame.list ${ns}::lornol "lines"
        $wpo.line_frame.list.menu delete 0
        foreach lnl "lines nolines" {
          $wpo.line_frame.list.menu add radiobutton -label $lnl \
              -variable ${ns}::lornol 
        }
        pack $wpo.line_frame.list -side left -anchor w -fill x  

        grid [label $wpo.linewidth_label -text "Line width:"] \
          -row 3 -column 5 -sticky w
        grid [entry $wpo.linewidth_entry -width 4 -textvariable ${ns}::linewidth] \
          -row 3 -column 6 -sticky w
          
        grid [label $wpo.marker_label -text "Marker:"] \
          -row 5 -column 1 -sticky w
        grid [frame $wpo.marker_frame] \
          -row 5 -column 2 -sticky w
        tk_optionMenu $wpo.marker_frame.list ${ns}::marker "circle"
        $wpo.marker_frame.list.menu delete 0
        foreach mrkr "none point circle square" {
          $wpo.marker_frame.list.menu add radiobutton -label $mrkr \
              -variable ${ns}::marker     
        }
        pack $wpo.marker_frame.list -side left -anchor w -fill x
        
        grid [label $wpo.radius_label -text "Marker size:"] \
          -row 5 -column 5 -sticky w
        grid [entry $wpo.radius_entry -width 4 -textvariable ${ns}::mradius] \
          -row 5 -column 6 -sticky w

        #grid [button $wpo.dash_help -text "?" \
        #    -command {tk_messageBox -type ok -title "HELP" \
        #      -message "Draw dashed lines."}] \
        #  -row 7 -column 0 -sticky w
        #grid [label $wpo.dash_label -text "Dashed lines:"] \
        #  -row 7 -column 1 -sticky w
        #grid [frame $wpo.dash_frame] \
        #  -row 7 -column 2 -sticky w
        #tk_optionMenu $wpo.dash_frame.list ${ns}::dash "-"
        #$wpo.dash_frame.list.menu delete 0
        #foreach dsh "no - , . _" {
        #  $wpo.dash_frame.list.menu add radiobutton -label $dsh \
        #      -variable ${ns}::dash
        #}
        #pack $wpo.dash_frame.list -side left -anchor w -fill x  
        
        ${ns}::loadCoordinates
        if {$ndim == 3} {
          ${ns}::drawArrows
        }

        return $w
      }

    }

    #lappend namespaces $ns
    #lappend nmwizguis [string range $ns 2 end]
    return $ns
  }
}

proc nmwiz_tk {} {
  ::NMWiz::initGUI
  return $::NMWiz::w
}

proc nmwiz_load {filename} {
  return [nmwiz load $filename]
} 

proc nmwiz { args } {
  set cmd [lindex $args 0]
  if {![llength $cmd]} {
    vmdcon -info "nmwiz commands: load, list, main"
    return
  }
  if {$cmd=="list"} {
    set handles {}
    foreach ns [namespace children :: "nmdset*"] { 
      lappend handles [subst $ns]::handle
    }
    return $handles
  } elseif {$cmd=="main"} {
    nmwiz_tk
  } elseif {$cmd=="load"} {
    set fn [lindex $args 1]
    if {![llength $fn]} { 
      vmdcon -err "a .nmd filename needs to be specified"
      return     
    }
    if {![file isfile $fn]} {
      vmdcon -err "$fn is not a valid filename"
      return     
    }
    return [::NMWiz::loadNMD $fn]
  } else {
    vmdcon -err "$cmd is not a valid nmwiz command"
    return     
  }
}
