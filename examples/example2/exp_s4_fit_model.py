#!/usr/bin/env python
"""
HYDROID (HYDroxyl-Radical fOotprinting Interpretation for DNA)
Example 2, chicken nucleosomes reconstituted on a well-positioning sequence, Morozov et al. NAR 2009

HYDROIDexp, Stage 4:
Fitting a model based on a sum of Gaussian/Lorentzian functions representing each band peak.

In this stage a model will be fit to the original HRF profile curve,
this allows to extract unbaised DNA cleavage intensities for every DNA position.

"""

from Bio import SeqIO
from os import mkdir
import sys

from hydroid.HYDROIDexp import fit_peaks

lane_profile_file="data/lane_profiles.xls"
lane_config_file="data/lane_config.csv"
out_path="results"
try:
	mkdir(out_path)
except:
	pass

#Read DNA sequence from file via biopython
TS_seq=SeqIO.parse(open("data/DNA_seq.fasta"),'fasta').next().seq
BS_seq=TS_seq.reverse_complement()

lane_data=[
{'lane_name':'gg_601_TS_a','seq':TS_seq,'label':'five_prime','constraint':'dSIGMA>=0'},
{'lane_name':'gg_601_TS_b','seq':TS_seq,'label':'five_prime','constraint':'dSIGMA>=0'},
{'lane_name':'gg_601_BS_a','seq':BS_seq,'label':'five_prime','constraint':'SAFA'}, # dSIGMA>=0 solutions do not behave well on the right end (can be seen by visual inspection)
{'lane_name':'gg_601_BS_b','seq':BS_seq,'label':'five_prime','constraint':'dSIGMA>=0'}
]


#Common code:
# will iterate through the lanes,
# perform fitting and export a png plot and a csv file.
# For interactive exploration of the plot set graphshow=True
# To use precalculated results set csvfilein to the same as csvfileout
###################################
for s in lane_data:
	fit_peaks(lane_profile_file,lane_config_file,DNAseq=s['seq'],labeled_end=s['label'],lane_name=s['lane_name'],out_path=out_path,\
		peaktype='Gaussian',fitting_constraint=s['constraint'],maxfev=50000,graphshow=False,plotcontrib=True,\
		csvfileout=s['lane_name']+'_fitted_intensities.csv',pngfileout=s['lane_name']+'_fitted_intensities.png',\
		csvfilein=None,Nauxvirtpeaks=2)#s['lane_name']+'_fitted_intensities.csv')
