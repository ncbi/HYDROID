#!/usr/bin/env python
"""
HYDROID (HYDroxyl-Radical fOotprinting Interpretation for DNA)
Example 2, chicken nucleosomes reconstituted on a well-positioning sequence, Morozov et al. NAR 2009

HYDROIDexp, Step 2:
assigning initial peaks positions to gel band peaks on 1-D lane profiles

The automatic peak identification algorithm is guided by parameters set in data/lane_config.csv
(see inside for a list of parameters).
HYDROIDexp functions triggered by this stript allow for interactive manual adjustment of the majority of
these paramters and saving them back to file (press Save button).
Positions of band peaks identified by the semi-automatic algorithm are maked by asterisks and change interactively
upon adjustments of the algorithm paramters.

As a result of this step it is important to make sure that every peak in the region of interest of
HRF lane profiles (specified by leftlim an rightlim paramters) is maked by exactly on asterisk.

For the Maxam-Guilbert lanes exact identification of peaks is not important,
however, this step helps to set the alignpos parameter - a characteristic position on a profile that
will be used in the next step to align HRF and Maxam-Guilbert profiles for sequence calling.

It is advised to check data/lane_config.csv afterwards,
any saved records will be written at the end of the file and override any previous records.


"""
import sys
sys.path.insert(0,'../..')

from HYDROID.HYDROIDexp import assign_peaks_interactive

lane_profile_file="data/lane_profiles.xls"
lane_config_file="data/lane_config.csv"
lane_names=['gg_601_TS_a','gg_601_TS_b','ddG_601_TS','ddT_601_TS',\
			'ddA_601_BS','ddG_601_BS','gg_601_BS_a','gg_601_BS_b']


#Common code:
# will iterate through the lanes opening an interactive window which will show
# results of automatic peak identification.
# Parameters of the automatic peak identification should be interactively adjusted
# until locations of all peaks are identified correctly
# This mainly applies to the OH-footprinting profiles that will be latter quantified,
# sequencing profiles might be left as is.
###################################
for LN in lane_names:
	assign_peaks_interactive(lane_profile_file,lane_config_file,LN)



