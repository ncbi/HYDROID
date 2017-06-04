#!/usr/bin/env python
"""
HYDROID (HYDroxyl-Radical fOotprinting Interpretation for DNA)
Example 1, centromeric nucleosome of yeast reconstituted on a well-positioning sequence, Shaytan et al. 2017

HYDROIDexp, Step 4:
Fitting a model based on a sum of Gaussian/Lorentzian functions representing each band peak.

In this step a model will be fit to the original HRF profile curve,
this allows to extract unbaised DNA cleavage intensities for every DNA position.

"""

from Bio import SeqIO
from os import mkdir
import sys
sys.path.insert(0,'../..')

from HYDROID.HYDROIDexp import fit_peaks

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
{'lane_name':'scCSE4_601TA_BS','seq':BS_seq,'label':'three_prime'},
{'lane_name':'scCSE4_601TA_TS','seq':TS_seq,'label':'three_prime'}
]


#Common code:
# will iterate through the lanes,
# perform fitting and export a png plot and a csv file.
# For interactive exploration of the plot set graphshow=True
# To use precalculated results set csvfilein to the same as csvfileout
###################################
for s in lane_data:
	fit_peaks(lane_profile_file,lane_config_file,DNAseq=s['seq'],labeled_end=s['label'],lane_name=s['lane_name'],out_path=out_path,\
		peaktype='Gaussian',fitting_constraint='dSIGMA>=0',maxfev=50000,graphshow=False,plotcontrib=True,\
		csvfileout=s['lane_name']+'_fitted_intensities.csv',pngfileout=s['lane_name']+'_fitted_intensities.png',\
		csvfilein=None)#s['lane_name']+'_fitted_intensities.csv')
