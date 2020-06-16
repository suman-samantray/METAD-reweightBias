#!/opt/local/bin/python
from __future__ import print_function
import copy, numpy, scipy, scipy.interpolate, sys

#sys.path.append('/Users/davidcheung/Dropbox/Programs/Analysis/CollectiveVariables')
sys.path.append('/Users/sumansamantray/Research/ab16-22/code/static-bias')
#sys.path.append('/Users/davidcheung/Dropbox/Programs/Analysis/FES')
sys.path.append('/Users/sumansamantray/Research/ab16-22/code/ptMetad')

import CollectiveVariables
import interpolate_FES

R=8.3143/4184.0
beta=1.0/R/300.0

## get SS file from command line
inf=open(sys.argv[1],'r')

## read in free energy grid and create interpolation
#nah_grid,dh_grid,fes_grid=interpolate_FES.read_fes_grid(sys.argv[2])
cgc_grid,dh_grid,fes_grid=interpolate_FES.read_fes_grid(sys.argv[2])
#fes_spline=scipy.interpolate.RectBivariateSpline(nah_grid,dh_grid,fes_grid)
fes_spline=scipy.interpolate.RectBivariateSpline(cgc_grid,dh_grid,fes_grid)

## get numbers of residues from command line
nRes=int(sys.argv[3])

## initialise stuff
SStypes=['H','G','I','E','T','B','C']
SScount={'C':0.0,'T':0.0,'E':0.0,'B':0.0,'H':0.0,'G':0.0,'I':0.0}
SScountRes=[copy.deepcopy(SScount) for i in range(nRes)]
SScountResUnwgt=[copy.deepcopy(SScount) for i in range(nRes)]
norm=0.0

## loop over saved datasets
for line in inf.readlines():

    ## read in data for this timeset
    data=line.split()

    ## extract values of CVs used in biasing
    #nah=float(data[-2])
    cgc=float(data[-2])
    dh=float(data[-1])

    ## get weight for this timeset
    fe=fes_spline.ev([cgc],[dh])[0]
    prob=numpy.exp(-beta*fe)
    norm+=prob

    ## get secondary structure information
    ssdata=data[1]

    ## accumulate
    for iRes in range(nRes):
        SScountRes[iRes][ssdata[iRes]]+=prob
        SScountResUnwgt[iRes][ssdata[iRes]]+=1.0

## write out final results
for iRes in range(nRes):
    #print ("%4d " % (iRes+1), end='')
    print ("%4d " % (iRes+1))
    for sst in SStypes:
        SScountRes[iRes][sst]/=norm
        #print ("%s %8.3f %8.3f" % (sst,SScountRes[iRes][sst],SScountResUnwgt[iRes][sst]), end='')
        print ("%s %8.3f %8.3f" % (sst,SScountRes[iRes][sst],SScountResUnwgt[iRes][sst]))
    print ()
