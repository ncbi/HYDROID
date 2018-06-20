#!/usr/bin/env python
"""
HYDROID (HYDroxyl-Radical fOotprinting Interpretation for DNA)
Example 2, chicken nucleosomes reconstituted on a well-positioning sequence, Morozov et al. NAR 2009

Comparison stage:
Comparing experimental profiles from different gel lanes
as well as similarity between top and bottom strand profiles due
to pseudosymmetry of nucleosome

"""

from Bio.Seq import Seq
import os,sys,tempfile
import pandas as pd


from hydroid.HYDROIDexp import plot_prof_on_seq


out_path="results"
try:
	os.mkdir(out_path)
except:
	pass


#601-DNA sequence as in experiment 147bp long

TS_seq=Seq('CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCACATATATACATCCTAT')
BS_seq=TS_seq.reverse_complement()

prof_datasets=[
	{
	'name':'Experimental cleavage frequencies from different gel lanes, gg-601 nucleosomes, Top Strand',
	'pngfileout':'exp_compar_TS.png',
	'zero_at':74, #move zero on the plot to this position
	'seq':TS_seq, #common sequence for all profiles
	#to align data from different files we set
	#start_on_seq - parameter, 1-based numbering position
	# on 'seq' determines where porfile from corresponding file starts (including any NaNs)
	'rescale':'every',
	'profiles':[
	{'file':'gg_601_TS_a_fitted_intensities.csv','start_on_seq':53,'datacolumn':'Intensity','name':'Experiment A'},
	{'file':'gg_601_TS_b_fitted_intensities.csv','start_on_seq':53,'datacolumn':'Intensity','name':'Experiment B'},
	]},
	{
	'name':'Experimental cleavage frequencies from different gel lanes, gg-601 nucleosomes, Bottom Strand',
	'pngfileout':'exp_compar_BS.png',
	'zero_at':74, #move zero on the plot to this position
	'seq':BS_seq, #common sequence for all profiles
	#to align data from different files we set
	#start_on_seq - parameter, 1-based numbering position
	# on 'seq' determines where porfile from corresponding file starts (including any NaNs)
	'rescale':'every',
	'profiles':[
	{'file':'gg_601_BS_a_fitted_intensities.csv','start_on_seq':33,'datacolumn':'Intensity','name':'Experiment A'},
	{'file':'gg_601_BS_b_fitted_intensities.csv','start_on_seq':34,'datacolumn':'Intensity','name':'Experiment B'},
	]},
	{
	'name':'Experimental cleavage frequencies from top and bottom strands are symmetric, gg-601 nucleosomes',
	'pngfileout':'exp_compar_BS_TS.png',
	'zero_at':74, #move zero on the plot to this position
	'seq':TS_seq, #common sequence for all profiles
	#to align data from different files we set
	#start_on_seq - parameter, 1-based numbering position
	# on 'seq' determines where porfile from corresponding file starts (including any NaNs)
	'rescale':'fit',
	'profiles':[
	{'file':'gg_601_TS_a_fitted_intensities.csv','start_on_seq':53,'datacolumn':'Intensity','name':'Top strand'},
	{'file':'gg_601_BS_a_fitted_intensities.csv','start_on_seq':33,'datacolumn':'Intensity','name':'Bottom strand'},
	]}
]

#Common code:
#Combines profiles from different files into one file/dataframe and sends for plotting
#last line plot_options in plot_prof_on_seq can be used to customize plot appearance
###################################
for sets in prof_datasets:
	data={}
	columns=[]
	#Here we combine datasets to one dataframe for plotting,
	#by simply stacking columns from files together
	#then provide column names to plotting functions as array
	for p in sets['profiles']:
		data[p['name']]=pd.read_csv(os.path.join(out_path,p['file']),comment='#')
		data[p['name']].index=range(0+p.get('start_on_seq',1),p.get('start_on_seq',1)+len(data[p['name']]))
		data[p['name']].rename(columns={p['datacolumn']:p['name'],p.get('sitecolumn','Site'):'Site_'+p['name']}, inplace=True)
		columns.append(p['name'])
	site_column_df=pd.DataFrame({'Site':['%d%s'%(n,l) for n,l in zip(range(1,1+len(sets['seq'])),sets['seq'])]},index=range(1,1+len(sets['seq'])))
	result=pd.concat(data.values()+[site_column_df],axis=1)
	# print result
	temp = tempfile.TemporaryFile()
	result.to_csv(temp)
	temp.seek(0)
	plot_prof_on_seq(temp,DNAseq=sets['seq'],\
		graphshow=True,pngfileout=os.path.join(out_path,sets['pngfileout']),title=sets['name'],\
		prof_columns=columns,seq_column="Site",rescale=sets['rescale'],zero_at=sets['zero_at'],\
		ylab=sets.get('ylab',None),
		colorb={'A':'#0b0','T':'#00b','G':'#fff','C':'#fff'},colorf={'A':'#fafafa','T':'#fafafa','G':'#000','C':'#000'},\
		plot_options={'linewidth':1.0,'markersize':8.0,'figsize':(12,3),'fontsize':None,'legendloc':'upper right'})
	temp.close()