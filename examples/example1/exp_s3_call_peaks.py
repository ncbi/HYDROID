#!/usr/bin/env python
"""
HYDROID (HYDroxyl-Radical fOotprinting Interpretation for DNA)
Example 1, centromeric nucleosome of yeast reconstituted on a well-positioning sequence, Shaytan et al. 2017

HYDROIDexp, Step 3:
Calling peaks (assigning DNA sequence positions to individual peaks)

Peak calling is usually accomplished by comparing footprinting lanes to sequencing lanes or other data.

In this step the user is allowed to compare Maxam-Gilbert and HRF profiles on same plots
and to interactively specify the location of any band peak on the DNA sequence.
Sequence calling for all other peaks will be done automatically.

The data is saved to data/lane_config.csv by pressing the Save button.
It is advised to check data/lane_config.csv afterwards,
any saved records will be written at the end of the file and override any previous records.

"""


from Bio import SeqIO

import sys

from hydroid.HYDROIDexp import call_peaks_interactive

lane_profile_file="data/lane_profiles.xls"
lane_config_file="data/lane_config.csv"
#Read DNA seqeunce from file via biopython
TS_seq=SeqIO.parse(open("data/DNA_seq.fasta"),'fasta').next().seq
BS_seq=TS_seq.reverse_complement()

lane_sets=[
{'footprinting_profile':'scCSE4_601TA_BS','helper_profiles':['GA_601TA_BS','CT_601TA_BS'],'seq':BS_seq,'label':'three_prime'},
{'footprinting_profile':'scCSE4_601TA_TS','helper_profiles':['GA_601TA_TS','CT_601TA_TS'],'seq':TS_seq,'label':'three_prime'}
]


#Common code:
# will iterate through the lanes sets opening an interactive window which will show
# footprinting profile overlayed with helper profiles,
# and will allow to interactively assign any identified peak to a position on the DNA seqeunce
###################################
for s in lane_sets:
	call_peaks_interactive(lane_profile_file,lane_config_file,DNAseq=s['seq'],labeled_end=s['label'],lane_name=s['footprinting_profile'],helper_prof_names=s['helper_profiles'])



