# 2020-05-12 Ji Young Lee

set struc     ./t0.pdb                       
set psffi     XXX
set trajec    RRR                        
set freq      SSTEP
set start     FSTEP 

mol load psf $psffi pdb $struc
animate read dcd $trajec beg $start end -1 skip $freq waitfor all

####wrap
package require pbctools
pbc wrap -compound fragment -center com -centersel "protein" -sel " SSS " -all

####alignment
set ref_molid [molinfo top get id]
set traj_molid [molinfo top get id]

set calpha "alpha carbon"
set CAreference [atomselect $ref_molid $calpha frame 0]
set CAcompare   [atomselect $traj_molid $calpha]
set System      [atomselect $traj_molid " SSS "]

set num_frames [molinfo $traj_molid get numframes]

for {set f 1} {$f < $num_frames} {incr f} {
        molinfo $traj_molid set frame $f
        $CAcompare frame $f
        $System    frame $f
        set trans_mat [measure fit $CAcompare $CAreference]
        $System    move $trans_mat
}
####alignment

animate write pdb pp0.pdb beg 0 end 0 sel [atomselect top " SSS "]
animate write psf pp0.psf beg 0 end 0 sel [atomselect top " SSS "]
animate write pdb pp1.pdb beg 1 end 1 sel [atomselect top " SSS "]
animate write psf pp1.psf beg 1 end 1 sel [atomselect top " SSS "]
animate write dcd pp.dcd  beg 1 end -1 waitfor all sel [atomselect top " SSS "]

exit
