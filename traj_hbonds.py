#!/opt/local/bin/python

import sys,numpy, pandas
import MDAnalysis,MDAnalysis.analysis.hbonds

## read in simulation trajectory
u=MDAnalysis.Universe(sys.argv[1],sys.argv[2])


## make selections
pro_select="protein"
back_select="backbone"

## make selections
protein=u.select_atoms(pro_select)
backbone=u.select_atoms(back_select)


h=MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u,'protein','protein',detect_hydrogens='heuristic')
h.run()
h.generate_table()
data_frame=pandas.DataFrame.from_records(h.table)
data_frame.to_csv(sys.argv[3]+'_hbonds-protein.csv')


hb=MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u,'backbone','backbone',detect_hydrogens='heuristic')
hb.run()
hb.generate_table()
data_frame=pandas.DataFrame.from_records(hb.table)
data_frame.to_csv(sys.argv[3]+'_hbonds-backbone.csv')
