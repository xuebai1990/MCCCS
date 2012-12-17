#!/usr/bin/env python

import sys
import convert2pdb

#>>> molSpec specifies the atom name of each bead
#   The example here is a system of TIP4P water, methanol, and ethanol
molSpec=[['O','H','H','M'],['O','H','C'],['O','H','C','C']]
#>>> printMask indicates what beads to include in the PDB file
#   Here atom M in the TIP4P model is not output
printMask=[[1,1,1,0],[1,1,1],[1,1,1,1]]
convert2pdb.cfg2pdb(sys.argv,molSpec,printMask,nbox=3)
