#!/usr/bin/env python
"""
HYDROID (HYDroxyl-Radical fOotprinting Interpretation for DNA)
Example 1, centromeric nucleosome of yeast reconstituted on a well-positioning sequence, Shaytan et al. 2017

HYDROIDpred, Step 2:
Calculating H-SASA profiles (deoxyribose hydrogen atoms solvent accessible surface area)
from a strucutre of a protein DNA complex.
FreeSASA software is used to perform the calculation of
solven accessible surface area (SASA) of deoxyribose hydrogen atoms.
SASA areas of these hydrogen atoms are summed up for every nucleotide.


See get_DNA_H_SASA inline help for explanation of parameters.

"""

import os
from Bio.Seq import Seq

from hydroid.HYDROIDpred import get_DNA_H_SASA

out_path="results"
str_path="data/structures"
try:
	os.mkdir(out_path)
except:
	pass

#Set DNA sequence
TS_seq=Seq('GGAGATACCCGGTGCTAAGGCCGCTTAATTGGTCGTAGCAAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCTACCGCGTTTTAACCGCCAATAGGATTACTTACTAGTCTCTAGGCACGTGTAAGATATATACATCCT')
BS_seq=TS_seq.reverse_complement()

prof_data=[
	{'prof_name':'scCSE4_601TA_TS','pdb_file':'scCSE4_601TA_nucl_H_Ndist.pdb','chain':'I','resids':range(-71,72),'seq':TS_seq,'vdw_set':'charmm36-rmin'},
	{'prof_name':'scCSE4_601TA_BS','pdb_file':'scCSE4_601TA_nucl_H_Ndist.pdb','chain':'J','resids':range(-71,72),'seq':BS_seq,'vdw_set':'charmm36-rmin'}
	]


#Common code:
#will calculate H-SASA profiles from PDB structure with hydrogen atoms.
#adjust n_threads to the number of CPU cores in your computer for optimal performance
###################################
for p in prof_data:
	get_DNA_H_SASA(os.path.join(str_path,p['pdb_file']),os.path.join(out_path,p['prof_name']+'_H-SASA.csv'),\
		chain=p['chain'],resids=p['resids'],seq=p['seq'],probe_radius=1.4,slicen=200,vdw_set=p['vdw_set'],\
		Hcontrib=[1.0]*7,n_threads=6,verbose=False)

