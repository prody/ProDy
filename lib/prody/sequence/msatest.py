from os.path import isfile
from time import time

from numpy import *
from prody import *
from msatools import *

filename = '/home/abakan/Downloads/PF00069_full.sth'
if not isfile('/home/abakan/Downloads/PF00069_full.sth'):
    filename = 'msa_Cys_knot.slx' 

msa = MSAFile(filename)

t=time();
arr = zeros((80000, msa.numResidues()), '|S1')

print 'empty array', time()-t 

t=time()
titles = parseSelex(filename, arr)
print 'parse array', len(titles), time()-t 
#for title in titles:
#    print title
