# 2020-05-12 Ji Young Lee

###Check
Struc=../md.pdb
PSFfi=' ..\/md.psf '
Traje=' ..\/sim.dcd '
Selec='all     '
Step=1  
Start=0
Head='CRYST1   85.328   96.807   93.491  90     90     90'
###Check

sed '1,1d' $Struc  >  __test   ## Check!!!
echo "$Head "      >  t0.pdb
cat __test         >> t0.pdb

cp 0chlcmpall.tcl ___ligb.tcl
sed -e "s/XXX/$PSFfi/g" -e "s/RRR/$Traje/g" -e "s/SSS/$Selec/g" -e "s/FSTEP/$Start/g" -e "s/SSTEP/$Step/g" ___ligb.tcl  > __ligb.tcl

vmd -dispdev text -e __ligb.tcl    
rm __*
