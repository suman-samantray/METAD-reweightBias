#!/opt/local/bin/python
from __future__ import print_function
import sys, numpy, scipy, scipy.interpolate,matplotlib,matplotlib.pyplot

#sys.path.append('/Users/davidcheung/Dropbox/Programs/Analysis/CollectiveVariables')
sys.path.append('/Users/sumansamantray/Research/ab16-22/code/static-bias')
#sys.path.append('/Users/davidcheung/Dropbox/Programs/Analysis/FES')
sys.path.append('/Users/sumansamantray/Research/ab16-22/code/ptMetad')

import CollectiveVariables
import interpolate_FES

def anint(x):
    xlo=numpy.floor(x)
    xhi=xlo+1.0
    if (xhi-x)<(x-xlo):
        return xhi
    else:
        return xlo

R=8.3143/4184.0
beta=1.0/R/300.0

inf=open(sys.argv[1],'r')
ouf=open(sys.argv[3],'w')
nRes=int(sys.argv[4])

rama=numpy.zeros((361,361))
rama_r=numpy.zeros((nRes-2,361,361))
cnt=0.0
cnt_r=numpy.zeros((nRes-2))
## read in free energy grid and create interpolation
#nah_grid,dh_grid,fes_grid=interpolate_FES.read_fes_grid(sys.argv[2])
cgc_grid,dh_grid,fes_grid=interpolate_FES.read_fes_grid(sys.argv[2])
#fes_spline=scipy.interpolate.RectBivariateSpline(nah_grid,dh_grid,fes_grid)
fes_spline=scipy.interpolate.RectBivariateSpline(cgc_grid,dh_grid,fes_grid)

i=0
for line in inf.readlines():
    data=line.split()
    #nah=float(data[5])
    cgc=float(data[8])
    dh=float(data[7])


    fe=fes_spline.ev([dh],[cgc])[0]
    prob=numpy.exp(-beta*fe)

    for iRes in range(nRes-2): 
        phi=float(data[10+2*iRes])
        psi=float(data[11+2*iRes])

        phi_bin=int(anint(phi+180.0))
        psi_bin=int(anint(psi+180.0))

        rama[psi_bin,phi_bin]+=prob
        rama_r[iRes,psi_bin,phi_bin]+=prob
        cnt+=prob
        cnt_r[iRes]+=prob

rama/=cnt
for iRes in range(nRes-2):
    rama_r[iRes,:,:]/=cnt_r[iRes]

sum_ll=0.0
sum_ul=0.0
sum_lu=0.0
sum_uu=0.0
for iphi in range(361):
    for ipsi in range(361):
        phi=iphi-180.0
        psi=ipsi-180.0
        if phi<0.0 and psi<0.0:
            sum_ll+=rama[ipsi,iphi]
        elif phi<0.0 and psi>=0.0:
            sum_ul+=rama[ipsi,iphi]
        elif phi>=0.0 and psi<0.0:
            sum_lu+=rama[ipsi,iphi]
        else:
            sum_uu+=rama[ipsi,iphi]
        print (phi,psi,rama[iphi,ipsi],file=ouf)
nn=sum_ll+sum_ul+sum_lu+sum_uu
sum_ll/=nn
sum_ul/=nn
sum_lu/=nn
sum_uu/=nn
print ("All   %8.3f  %8.3f  %8.3f  %8.3f" % (sum_ll,sum_ul,sum_lu,sum_uu))


for iRes in range(nRes-2):
    sum_ll=0.0
    sum_ul=0.0
    sum_lu=0.0
    sum_uu=0.0
    for iphi in range(361):
        for ipsi in range(361):
            phi=iphi-180.0
            psi=ipsi-180.0
            print (phi,psi,rama_r[iRes,iphi,ipsi],file=ouf)
            if phi<0.0 and psi<0.0:
                sum_ll+=rama_r[iRes,ipsi,iphi]
            elif phi<0.0 and psi>=0.0:
                sum_ul+=rama_r[iRes,ipsi,iphi]
            elif phi>=0.0 and psi<0.0:
                sum_lu+=rama_r[iRes,ipsi,iphi]
            else:
                sum_uu+=rama_r[iRes,ipsi,iphi]
    nn=sum_ll+sum_ul+sum_lu+sum_uu
    sum_ll/=nn
    sum_ul/=nn
    sum_lu/=nn
    sum_uu/=nn
    print ("%3d   %8.3f  %8.3f  %8.3f  %8.3f" % (iRes+2,sum_ll,sum_ul,sum_lu,sum_uu))
