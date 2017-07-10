#!/usr/bin/env python
"""
HYDROID (HYDroxyl-Radical fOotprinting Interpretation for DNA)
Example 2, chicken nucleosomes reconstituted on a well-positioning sequence, Morozov et al. NAR 2009

HYDROIDexp, Step 3:
Calling peaks (assigning DNA sequence positions to individual peaks)

Peak calling is usually accomplished by comparing footprinting lanes to sequencing lanes or other data.

In this step the user is allowed to compare PCR profiles with dideoxynucleotide triphosphates
and HRF profiles on same plots
and to interactively specify the location of any band peak on the DNA sequence.
Sequence calling for all other peaks will be done automatically.

The data is saved to data/lane_config.csv by pressing the Save button.
It is advised to check data/lane_config.csv afterwards,
any saved records will be written at the end of the file and override any previous records.

"""

from Bio import SeqIO

import sys
sys.path.insert(0,'../..')

from HYDROID.HYDROIDexp import call_peaks_interactive

lane_profile_file="data/lane_profiles.xls"
lane_config_file="data/lane_config.csv"
#Read DNA seqeunce from file via biopython
TS_seq=SeqIO.parse(open("data/DNA_seq.fasta"),'fasta').next().seq
BS_seq=TS_seq.reverse_complement()

lane_sets=[
{'footprinting_profile':'gg_601_TS_a','helper_profiles':['ddG_601_TS','ddT_601_TS'],'seq':TS_seq,'label':'five_prime'},
{'footprinting_profile':'gg_601_TS_b','helper_profiles':['ddG_601_TS','ddT_601_TS'],'seq':TS_seq,'label':'five_prime'},
{'footprinting_profile':'gg_601_BS_a','helper_profiles':['ddA_601_BS','ddG_601_BS'],'seq':BS_seq,'label':'five_prime'},
{'footprinting_profile':'gg_601_BS_b','helper_profiles':['ddA_601_BS','ddG_601_BS'],'seq':BS_seq,'label':'five_prime'}
]


#Common code:
# will iterate through the lanes sets opening an interactive window which will show
# footprinting profile overlayed with helper profiles,
# and will allow to interactively assign any identified peak to a position on the DNA seqeunce
###################################
for s in lane_sets:
	call_peaks_interactive(lane_profile_file,lane_config_file,DNAseq=s['seq'],labeled_end=s['label'],lane_name=s['footprinting_profile'],helper_prof_names=s['helper_profiles'])



