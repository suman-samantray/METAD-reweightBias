#!/opt/local/bin/python
from __future__ import print_function
import numpy, pandas, sys, scipy, scipy.interpolate

sys.path.append('/Users/sumansamantray/Research/ab16-22/code/ptMetad')

import interpolate_FES

R=8.3143/4184.0
beta=1.0/R/300.0
RT=R*300.


cvf=open(sys.argv[2],'r')

cv_data=cvf.readlines()
ncv=len(cv_data)

t_cv=numpy.zeros((ncv))
cgc_t=numpy.zeros((ncv))
dh_t=numpy.zeros((ncv))
for i,line in enumerate(cv_data):
    data=line.split()
    t_cv[i]=float(data[0])
    cgc_t[i]=float(data[5])
    dh_t[i]=float(data[7])

## read in free energy grid and create interpolation
#nah_grid,dh_grid,fes_grid=interpolate_FES.read_fes_grid(sys.argv[3])
cgc_grid,dh_grid,fes_grid=interpolate_FES.read_fes_grid(sys.argv[3])

fmin=fes_grid.min()
fes_grid-=fmin
fes_spline=scipy.interpolate.RectBivariateSpline(cgc_grid,dh_grid,fes_grid)

## read in h-bond data
hbond_data=pandas.read_csv(sys.argv[1])

nres=int(sys.argv[5])
dt=float(sys.argv[6])


g=hbond_data.groupby('time')
timesteps=list(g.groups.keys())

timesteps.sort()
nn=int(timesteps[-1]/dt)
nstep=nn+1

nhbonds_time=numpy.zeros((nstep))
nhbonds_acc_res_time=numpy.zeros((nstep,nres))
nhbonds_don_res_time=numpy.zeros((nstep,nres))

nhbond_ave=0.0
nhbond_sqr=0.0
nhbond_acc_res_ave=numpy.zeros((nres))
nhbond_acc_res_sqr=numpy.zeros((nres))
nhbond_don_res_ave=numpy.zeros((nres))
nhbond_don_res_sqr=numpy.zeros((nres))

nhbond_pair=numpy.zeros((nres,nres))

prob=1.0
norm=0.0
for gg in g:


    istep=int(gg[0]/dt)

    ## get weight for this timeset
    fe=fes_spline.ev([cgc_t[istep]],[dh_t[istep]])[0]
    prob=numpy.exp(-beta*fe)
    norm+=prob

    nhbond=0.0
    nhbond_acc=numpy.zeros((nres))
    nhbond_don=numpy.zeros((nres))

    for donor_id,acceptor_id in zip(gg[1]['donor_resid'],gg[1]['acceptor_resid']):
        nhbond+=1.0
        nhbond_acc[acceptor_id-1]+=1.0
        nhbond_don[donor_id-1]+=1.0
        nhbond_pair[acceptor_id-1,donor_id-1]+=prob
        nhbond_pair[donor_id-1,acceptor_id-1]+=prob

    nhbond_ave+=prob*nhbond
    nhbond_sqr+=prob*nhbond**2
    nhbond_acc_res_ave+=prob*nhbond_acc
    nhbond_acc_res_sqr+=prob*nhbond_acc**2
    nhbond_don_res_ave+=prob*nhbond_don
    nhbond_don_res_sqr+=prob*nhbond_don**2
    print (istep,nhbond)

nhbond_ave/=norm
nhbond_sqr/=norm
nhbond_std=numpy.sqrt(nhbond_sqr-nhbond_ave**2)
nhbond_acc_res_ave/=norm
nhbond_acc_res_sqr/=norm
nhbond_acc_res_std=numpy.sqrt(nhbond_acc_res_sqr-nhbond_acc_res_ave**2)
nhbond_don_res_ave/=norm
nhbond_don_res_sqr/=norm
nhbond_don_res_std=numpy.sqrt(nhbond_don_res_sqr-nhbond_don_res_ave**2)

ouf=open(sys.argv[4],'w')
print ('### All  %8.3f  %8.3f' % (nhbond_ave,nhbond_std),file=ouf)
for ires in range(nres):
    print ('%4d %8.3f  %8.3f  %8.3f  %8.3f' % (ires+1,nhbond_acc_res_ave[ires],nhbond_acc_res_std[ires],nhbond_don_res_ave[ires],nhbond_don_res_std[ires]),file=ouf)

nhbond_pair/=norm
for ires in range(nres):
    for jres in range(ires+1,nres):
        if nhbond_pair[ires,jres]>0.5:
            print ("%4d %4d %8.3f" % (ires+1,jres+1,nhbond_pair[ires,jres]))
