#!/usr/bin/env python
"""
HYDROID (HYDroxyl-Radical fOotprinting Interpretation for DNA)
Example 1, centromeric nucleosome of yeast reconstituted on a well-positioning sequence, Shaytan et al. 2017

HYDROIDpred, Stage 3:
Plotting the extracted theoretical DNA cleavage freqiencies (H-SASA profiles)
along the DNA sequence.


"""


from Bio.Seq import Seq
import os,sys,tempfile
import pandas as pd


from hydroid.HYDROIDexp import plot_prof_on_seq,simulate_gel

out_path="results"
try:
	os.mkdir(out_path)
except:
	pass

#Set DNA sequence
TS_seq=Seq('GGAGATACCCGGTGCTAAGGCCGCTTAATTGGTCGTAGCAAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCTACCGCGTTTTAACCGCCAATAGGATTACTTACTAGTCTCTAGGCACGTGTAAGATATATACATCCT')
BS_seq=TS_seq.reverse_complement()


lane_data=[
	{'file':'scCSE4_601TA_TS_H-SASA.csv','seq':TS_seq,'pngfileout':'scCSE4_601TA_TS_H-SASA.png',
	'title':'scCSE4_601TA, Top Strand, H-SASA profile'},
	{'file':'scCSE4_601TA_BS_H-SASA.csv','seq':BS_seq,'pngfileout':'scCSE4_601TA_BS_H-SASA.png',
	'title':'scCSE4_601TA, Bottom Strand, H-SASA profile'},
	]

#Common code:
# will iterate through the lanes and plot the profile values along DNA sequence
# For interactive exploration of the plot set graphshow=True
###################################
for s in lane_data:
    plot_prof_on_seq(os.path.join(out_path,s['file']),DNAseq=s['seq'],\
        graphshow=True,pngfileout=os.path.join(out_path,s['pngfileout']),title=s['title'],\
        prof_columns='H-SASA',seq_column="Site",\
        colorb={'A':'#0b0','T':'#00b','G':'#fff','C':'#fff'},colorf={'A':'#fafafa','T':'#fafafa','G':'#000','C':'#000'})
    #plotting gels
    dataframe=pd.read_csv(os.path.join(out_path,s['file']),comment='#')
    data=dataframe['H-SASA'].as_matrix()
    simulate_gel(data,10,model='ogston',
                    pngfileout=os.path.join(out_path,s['pngfileout'][:-4]+'_simulated.png'),graphshow=True,title=s['title'])




