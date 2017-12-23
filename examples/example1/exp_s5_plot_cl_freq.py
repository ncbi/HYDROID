#!/usr/bin/env python
"""
HYDROID (HYDroxyl-Radical fOotprinting Interpretation for DNA)
Example 1, centromeric nucleosome of yeast reconstituted on a well-positioning sequence, Shaytan et al. 2017

HYDROIDexp, Step 5:
Plotting the extracted DNA cleavage frequencies along the seqeunce in a graphically compelling way.

"""


from Bio import SeqIO
import os
import sys

from hydroid.HYDROIDexp import plot_prof_on_seq

out_path="results"
try:
	os.mkdir(out_path)
except:
	pass

#Read DNA sequence from file via biopython
TS_seq=SeqIO.parse(open("data/DNA_seq.fasta"),'fasta').next().seq
BS_seq=TS_seq.reverse_complement()

lane_data=[
	{'file':'scCSE4_601TA_TS_fitted_intensities.csv','seq':TS_seq,'pngfileout':'scCSE4_601TA_TS_cl_freq_profile.png','title':'scCSE4_601TA, Top Strand, Cleavage frequency profile along DNA'},
	{'file':'scCSE4_601TA_BS_fitted_intensities.csv','seq':BS_seq,'pngfileout':'scCSE4_601TA_BS_cl_freq_profile.png','title':'scCSE4_601TA, Bottom Strand, Cleavage frequency profile along DNA'},
	]

#Common code:
# will iterate through the lanes and plot the profile values along DNA sequence
# For interactive exploration of the plot set graphshow=True
###################################
for s in lane_data:
	plot_prof_on_seq(os.path.join(out_path,s['file']),DNAseq=s['seq'],\
		graphshow=True,pngfileout=os.path.join(out_path,s['pngfileout']),title=s['title'],\
		prof_columns='Intensity',seq_column="Site",\
		colorb={'A':'#0b0','T':'#00b','G':'#fff','C':'#fff'},colorf={'A':'#fafafa','T':'#fafafa','G':'#000','C':'#000'})


