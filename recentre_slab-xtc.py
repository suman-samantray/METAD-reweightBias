#!/opt/local/bin/python
#from __future__ import print_function
import numpy, os, sys
import MDAnalysis

## recenter_slab
##
## python program for recentring molecule to centre of box for AWI
##

groFile=sys.argv[1]
xtcFile=sys.argv[2:-1]
outFile=sys.argv[-1]

## check if output file already exists and exit if is does
if os.path.isfile(outFile):
    print "output file %s already exists" % (outFile)
    print "cowardly refusing to overwrite"
    sys.exit()


u=MDAnalysis.Universe(groFile,xtcFile)

## get protein and water selections
protein=u.select_atoms('protein')
water=u.select_atoms('name OW HW1 HW2')

nAtoms=len(u.atoms)
nprotein_atom=len(protein)
nwater_atom=len(water)
nwater_mol=int(nwater_atom/3)

print (nprotein_atom,nwater_atom,nwater_mol,3*nwater_mol)


with MDAnalysis.Writer(outFile,nAtoms) as w:
    for ts in u.trajectory:

        lx=ts.dimensions[0]
        ly=ts.dimensions[1]
        lz=ts.dimensions[2]


        ## get centre of mass of water
        com_water=water.centroid()
        print (com_water)

        ## get shift in z-direction (from water)
        dz=com_water[2]-lz/2.0

        ## shift so water slab in centre of z-dimension
        for i in range(len(u.atoms)):
            p=u.atoms[i].position
            p[2]-=dz
            u.atoms[i].position=p

        ## get centre of mass of water
        com_water=water.centroid()

        ## get centre of mass of protein
        com_protein=protein.centroid()
        dx=com_protein[0]-lx/2.0
        dy=com_protein[1]-ly/2.0
        print (com_protein,dx,dy)

        ## shift protein
        for i in range(nprotein_atom):
            p=u.atoms[i].position
            p[0]-=dx
            p[1]-=dy
            u.atoms[i].position=p

        com_protein=protein.centroid()

        ## shift water molecules
        for imol in range(nwater_mol):
            iO=nprotein_atom+3*imol
            iH1=iO+1
            iH2=iO+2
            dH1=u.atoms[iH1].position-u.atoms[iO].position
            dH2=u.atoms[iH2].position-u.atoms[iO].position

            rO=u.atoms[iO].position
            rO[0]-=dx
            if rO[0]<0.0:
                rO[0]+=lx
            elif rO[0]>lx:
                rO[0]-=lx
            rO[1]-=dy
            if rO[1]<0.0:
                rO[1]+=ly
            elif rO[1]>ly:
                rO[1]-=ly

            u.atoms[iO].position=rO
            rH1=rO+dH1
            u.atoms[iH1].position=rH1
            rH2=rO+dH2
            u.atoms[iH2].position=rH2


    w.write(u)
    #sys.exit()
