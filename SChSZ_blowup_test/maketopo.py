
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw import topotools
from clawpack.geoclaw.util import get_remote_file
from numpy import *

def get_topo():
    file_name = "etopo4min120E65W65S65N.asc"
    url = "ptha@homer.u.washington.edu:CC/topo/" + file_name
    get_remote_file(url, verbose=True)

def makeqinit():
    """
    Create qinit data file
    """
    nxpoints = 101
    nypoints = 101
    xlower = 205.
    xupper = 217.
    yupper = -14.
    ylower = -22.
    outfile= "hump.xyz"     
    topotools.topo1writer(outfile,qinit,xlower,xupper,ylower,yupper,nxpoints,nypoints)


def qinit(x,y):
    """
    Gaussian hump:
    """
    from numpy import where
    ze = -((x-216.)**2 + (y+15)**2)/5.
    z = where(ze>-10., 0.5*exp(ze), 0.)
    return z

if __name__=='__main__':
    makeqinit()
    get_topo()