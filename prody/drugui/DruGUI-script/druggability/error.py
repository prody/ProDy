"""Druggability Index Analysis"""

from prody import *

pdb = "/Users/carlosventura/Desktop/druggability/ampa.case.study/newcgenff/analysis/ztrj/pp0.pdb"
prefix = 'dg'
psf="/Users/carlosventura/Desktop/druggability/ampa.case.study/newcgenff/analysis/ztrj/pp0.psf"
pdb = parsePDB(pdb, ag=parsePSF(psf))
align = 'calpha'
palign = pdb.select(align)
pcenter = calcCenter(palign)
probes = ['IPA', 'PRO', 'ACE', 'IBU', 'ACE', 'IMI', 'BEN']

probe_selstr = 'noh and resname'

dcd = "/Users/carlosventura/Desktop/druggability/ampa.case.study/newcgenff/analysis/ztrj/pp.dcd" 

for p in probes:
            sel = pdb.select('resname ' + p)
            print(sel)
            if sel is None:
                continue
                #raise ValueError('probe ' + p + ' is not found in the system')
            hv = sel.getHierView()
            n = hv.numResidues()
            res = next(hv.iterResidues())
            writePDB(prefix + '_' + p + '.pdb', res)
            plog(str(n) + ' copies of ' + p + ' is found.')
            probe_selstr += ' ' + p

print(probe_selstr)
PRBSEL = pdb.select(probe_selstr)
PRBIDX = PRBSEL.getIndices()
PROBES = PRBSEL.copy()


pcontact = pdb.select('protein')
writePDB(prefix + '_protein_heavyatoms.pdb', pcontact)

UNITCELL = []

from prody import LOGGER
LOGGER.progress('Evaluating frames:', len(dcd))

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
    PROBES._setCoords(PRBSEL._getCoords())
    cont = PROBES.select('same residue as within 4 of pcontact', 
                                 pcontact=pcontact)
            
    if cont:
        for res in cont.getHierView().iterResidues():
                DCDOUT[res.getResname()].write(res._getCoords()) 
    nframe += 1
    LOGGER.update(i)
            
dcd.close()