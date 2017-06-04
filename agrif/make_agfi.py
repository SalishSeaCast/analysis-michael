"""
Build AGRIF_FixedGrids.in by reading the values in the nesting tools' namelist
"""
import sys
from nemo_cmd.namelist import namelist2dict

namelistfile = sys.argv[1]
nst = namelist2dict(namelistfile)['nesting'][0]
imin,imax = nst['imin'],nst['imax']
jmin,jmax = nst['jmin'],nst['jmax']
rho,rhot = nst['rho'],nst['rhot']

print("Writing AGRIF_FixedGrids.in")
with open('AGRIF_FixedGrids.in','w') as f:
    f.write("1\n")
    f.write("{} {} {} {} {} {} {} {}\n".format(imin,imax,jmin,jmax,rho,rho,rhot,rhot))
    f.write("0\n\n")
    f.write("# number of children per parent\n")
    f.write("# imin imax jmin jmax spacerefx spacerefy timerefx timerefy\n")
    f.write("# [all coordinates are relative to each parent grid!]\n")
