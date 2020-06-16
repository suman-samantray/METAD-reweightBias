#!/opt/local/bin/python

import sys, numpy, MDAnalysis, os, subprocess
#sys.path.append('/Users/davidcheung/Dropbox/Programs/Analysis/CollectiveVariables')
sys.path.append('/Users/sumansamantray/Research/ab16-22/code/static-bias')
import CollectiveVariables

stride_command="/Users/sumansamantray/Stride/stride"

dt=10.0
words=[]
ssdata=""

groFile=sys.argv[1]
xtcFile=sys.argv[2]
ff=sys.argv[3]

u=MDAnalysis.Universe(groFile,xtcFile)

protein=u.select_atoms('protein')
nprotein_atoms=len(protein)

nres=len(protein.residues)

## find alpha-helix hydrogen bonds
#alpha_hbonds=CollectiveVariables.get_alpha_hbonds(protein,ff=ff)
cgamma_atoms=CollectiveVariables.get_cgamma_atoms(protein)
## find phi and psi angles
phi_angles,psi_angles=CollectiveVariables.get_phi_psi(protein)

for i,ts in enumerate(u.trajectory):

    timestep=i*dt
    pdbFile='tmp.pdb'

    ## get number of alpha-helix H-bonds
    #nh_alpha=CollectiveVariables.sum_alpha_hbonds(protein,alpha_hbonds,pbc=True,boxx=ts.dimensions)
    ## get number of Cgamma contacts
    cgamma_contacts=CollectiveVariables.sum_cgamma_contacts(protein,cgamma_atoms,pbc=True,boxx=ts.dimensions)

    ## get dihedral angle correlation function
    dih_offset=CollectiveVariables.beta_dih_offset(phi_angles,psi_angles)

    W=MDAnalysis.Writer(pdbFile,nprotein_atoms)
    W.write(protein)
    W.close()

    try:
        output=subprocess.check_output([stride_command,pdbFile])
    except:
        continue
    words=output.decode('utf-8').splitlines()
    
    for word in words:
        if word[:3]=='STR':
            ssdata=''
            for letter in word[10:10+nres]:
                if letter==' ':
                    ssdata=ssdata+'C'
                else:
                    ssdata=ssdata+letter
    #print ("%12.3f %s %12.6f %12.6f " % (timestep,ssdata,nh_alpha,dih_offset))
    print ("%12.3f %s %12.6f %12.6f " % (timestep,ssdata,cgamma_contacts,dih_offset))
    subprocess.call(["rm",pdbFile])

    #sys.exit()