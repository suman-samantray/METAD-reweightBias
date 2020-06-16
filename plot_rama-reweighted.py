#!/opt/local/bin/python

from matplotlib import rc, rcParams
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Plsztino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Plsztino']})
rc('text', usetex=True)
#rc('text', dvipnghack=True)
rcParams['text.latex.preamble'] = [r'\usepackage{siunitx}',r'\sisetup{detect-all}',r'\usepackage{helvet}',r'\usepackage{sansmath}',r'\sansmath']

presentation=False
paper=True
if presentation:
    ticksize=36
    labelsize=36
    linethickness=2
    xmin=0.2
    ymin=0.23
    dx=0.6
    dy=0.67
if paper:
    ticksize=24
    labelsize=24
    linethickness=2
    xmin=0.15
    ymin=0.15
    dx=0.75
    dy=0.75


import numpy, matplotlib, matplotlib.colors, matplotlib.pyplot, gzip
import matplotlib.gridspec
import matplotlib.patches as patches

import sys

def make_file_name(root,letter,serial):
    if serial<10:
        return root+'_'+letter+"%1s" % (serial)+'.png'
    if serial<100:
        return root+'_'+letter+"%2s" % (serial)+'.png'

inf=open(sys.argv[1],'r')
seq=sys.argv[3]
nres=len(seq)

## read in data for overall Ramachandran plot
phi=numpy.zeros((361,361))
psi=numpy.zeros((361,361))
rama=numpy.zeros((361,361))
norm=0.0
for iphi in range(361):
    for ipsi in range(361):
        data=inf.readline().split()
        phi[ipsi,iphi]=float(data[0])
        psi[ipsi,iphi]=float(data[1])
        rama[ipsi,iphi]=float(data[2])
        norm+=float(data[2])

rama/=norm
rama=numpy.transpose(rama)

## Make overall Ramachandran plot and save to file
fig=matplotlib.pyplot.figure(figsize=(10,8))
cs=matplotlib.pyplot.contourf([float(iPsi) for iPsi in range(-180,181)],[float(iPhi) for iPhi in range(-180,181)],rama,[0.0,0.00025,0.00050,0.00075,0.001],extend='both',cmap='plasma_r')
#cs=matplotlib.pyplot.contourf(rama)
#,[0.0,0.00025,0.00050,0.00075,0.001],extend='both',cmap='plasma_r')
# ca=matplotlib.pyplot.gca()
# ca.add_patch(patches.Rectangle((-122,-84),70,70,fill=False))
# ca.add_patch(patches.Rectangle((-150,95),80,80,fill=False))
# ca.add_patch(patches.Rectangle((14,52),70,70,fill=False))
# ca.add_patch(patches.Rectangle((-180,50),130,130,fill=False))
# ca.add_patch(patches.Rectangle((150,50),30,130,fill=False))
# ca.add_patch(patches.Rectangle((-180,-180),130,30,fill=False))
cb=matplotlib.pyplot.colorbar()

matplotlib.pyplot.xlabel('$\\phi$ / degrees',name='helvetica-oblique',size=24)
matplotlib.pyplot.xticks([-180,-90,0,90,180],['-180','-90','0','90','180'],size=18)
matplotlib.pyplot.ylabel('$\\psi$ / degrees',name='helvetica-oblique',size=24)
matplotlib.pyplot.yticks([-180,-90,0,90,180],['-180','-90','0','90','180'],size=18)
matplotlib.pyplot.savefig(sys.argv[2]+'.png')
matplotlib.pyplot.close()

## Loop over residues
for ires in range(1,nres-1):

    ## Read in Ramachandran plot for this residue
    phi=numpy.zeros((361,361))
    psi=numpy.zeros((361,361))
    rama=numpy.zeros((361,361))
    norm=0.0
    for iphi in range(361):
        for ipsi in range(361):
            data=inf.readline().split()
            phi[ipsi,iphi]=float(data[0])
            psi[ipsi,iphi]=float(data[1])
            rama[ipsi,iphi]=float(data[2])
            norm+=float(data[2])

    rama/=norm
    rama=numpy.transpose(rama)

    fig=matplotlib.pyplot.figure(figsize=(10,8))
    cs=matplotlib.pyplot.contourf([float(iPsi) for iPsi in range(-180,181)],[float(iPhi) for iPhi in range(-180,181)],rama,[0.0,0.00025,0.00050,0.00075,0.001],extend='both',cmap='plasma_r')
    #cs=matplotlib.pyplot.contourf(rama)
    #,[0.0,0.00025,0.00050,0.00075,0.001],extend='both',cmap='plasma_r')
    # ca=matplotlib.pyplot.gca()
    # ca.add_patch(patches.Rectangle((-122,-84),70,70,fill=False))
    # ca.add_patch(patches.Rectangle((-150,95),80,80,fill=False))
    # ca.add_patch(patches.Rectangle((14,52),70,70,fill=False))
    # ca.add_patch(patches.Rectangle((-180,50),130,130,fill=False))
    # ca.add_patch(patches.Rectangle((150,50),30,130,fill=False))
    # ca.add_patch(patches.Rectangle((-180,-180),130,30,fill=False))
    cb=matplotlib.pyplot.colorbar()

    matplotlib.pyplot.xlabel('$\\phi$ / degrees',name='helvetica-oblique',size=24)
    matplotlib.pyplot.xticks([-180,-90,0,90,180],['-180','-90','0','90','180'],size=18)
    matplotlib.pyplot.ylabel('$\\psi$ / degrees',name='helvetica-oblique',size=24)
    matplotlib.pyplot.yticks([-180,-90,0,90,180],['-180','-90','0','90','180'],size=18)
    matplotlib.pyplot.savefig(make_file_name(sys.argv[2],seq[ires],ires+1))
    matplotlib.pyplot.close()
