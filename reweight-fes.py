#!/opt/local/bin/python

import sys, numpy, scipy, scipy.interpolate

## change as needed
#sys.path.append('/Users/davidcheung/Dropbox/Programs/Analysis/FES')
sys.path.append('/Users/sumansamantray/Research/ab16-22/code/ptMetad')

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


## read in free enregy grid and create interpolation
z_grid,cv_grid,fes_grid=interpolate_FES.read_fes_grid(sys.argv[2])
fes_spline=scipy.interpolate.RectBivariateSpline(z_grid,cv_grid,fes_grid)

## initialise values
## put in code to initialise averages of quantities that will be calculated
## for each need to average the quantity and the square of the quantity (which will be used for
## calculating standard deviations)
rg_ave=0.0
rg_sqr=0.0
major_ave=0.0
major_sqr=0.0
median_ave=0.0
median_sqr=0.0
minor_ave=0.0
minor_sqr=0.0
alh_ave=0.0
alh_sqr=0.0
three_ten_alh_ave=0.0
three_ten_alh_sqr=0.0
dih_ave=0.0
dih_sqr=0.0
cga_ave=0.0
cga_sqr=0.0
saltb_ave=0.0
saltb_sqr=0.0
cnt1=0.0
cnt2=0.0
cnt3=0.0
cnt4=0.0
cnt5=0.0
cnt6=0.0
cnt7=0.0
cnt8=0.0
cnt9=0.0
## loop over saved data sets
for line in inf.readlines():

    data=line.split()

    ## get values of collective variables for reweighting
    ## assume it is the first two - change if needed
    cv1=float(data[7])
    cv2=float(data[8])

    ## get values of other quantities of interest
    rg=float(data[1])
    major=float(data[2])
    median=float(data[3])
    minor=float(data[4])
    alh=float(data[5])
    three_ten_alh=float(data[6])
    dih=float(data[7])
    cga=float(data[8])
    saltb=float(data[9])

    ## get value of free energy surface and reweighting factor
    fe=fes_spline.ev(z_grid,cv_grid)[0]
    prob=numpy.exp(-beta*fe)

    ## update averages
    rg_ave+=prob*rg
    rg_sqr+=prob*rg**2
    cnt1+=prob
    major_ave+=prob*major
    major_sqr+=prob*major**2
    cnt2+=prob
    median_ave+=prob*median
    median_sqr+=prob*median**2
    cnt3+=prob
    minor_ave+=prob*minor
    minor_sqr+=prob*minor**2
    cnt4+=prob
    alh_ave+=prob*alh
    alh_sqr+=prob*alh**2
    cnt5+=prob
    three_ten_alh_ave+=prob*three_ten_alh
    three_ten_alh_sqr+=prob*three_ten_alh**2
    cnt6+=prob
    dih_ave+=prob*dih
    dih_sqr+=prob*dih**2
    cnt7+=prob
    cga_ave+=prob*cga
    cga_sqr+=prob*cga**2
    cnt8+=prob
    saltb_ave+=prob*saltb
    saltb_sqr+=prob*saltb**2
    cnt9+=prob

## get final averages and standard deviations
rg_ave/=cnt1
rg_sqr/=cnt1
rg_std=numpy.sqrt(rg_sqr-rg_ave**2)
major_ave/=cnt2
major_sqr/=cnt2
major_std=numpy.sqrt(major_sqr-major_ave**2)
median_ave/=cnt3
median_sqr/=cnt3
median_std=numpy.sqrt(median_sqr-median_ave**2)
minor_ave/=cnt4
minor_sqr/=cnt4
minor_std=numpy.sqrt(minor_sqr-minor_ave**2)
alh_ave/=cnt5
alh_sqr/=cnt5
alh_std=numpy.sqrt(alh_sqr-alh_ave**2)
three_ten_alh_ave/=cnt6
three_ten_alh_sqr/=cnt6
three_ten_alh_std=numpy.sqrt(three_ten_alh_sqr-three_ten_alh_ave**2)
dih_ave/=cnt7
dih_sqr/=cnt7
dih_std=numpy.sqrt(dih_sqr-dih_ave**2)
cga_ave/=cnt8
cga_sqr/=cnt8
cga_std=numpy.sqrt(cga_sqr-cga_ave**2)
saltb_ave/=cnt9
saltb_sqr/=cnt9
saltb_std=numpy.sqrt(saltb_sqr-saltb_ave**2)
## print out results
print "--CVs--weighted_average--standard_deviation--"
print "Rg %8.3f %8.3f" % (rg_ave,rg_std)
print "val_max %8.3f %8.3f" % (major_ave,major_std)
print "val_mid %8.3f %8.3f" % (median_ave,median_std)
print "val_min %8.3f %8.3f" % (minor_ave,minor_std)
print "nh_alpha %8.3f %8.3f" % (alh_ave,alh_std)
print "nh_three_ten %8.3f %8.3f" % (three_ten_alh_ave,three_ten_alh_std)
print "dih_offset %8.3f %8.3f" % (dih_ave,dih_std)
print "cg_contacts %8.3f %8.3f" % (cga_ave,cga_std)
print "sb_contacts %8.3f %8.3f" % (saltb_ave,saltb_std)