#!/usr/bin/env python
"""
HYDROID (HYDroxyl-Radical fOotprinting Interpretation for DNA)
Example 1, centromeric nucleosome of yeast reconstituted on a well-positioning sequence, Shaytan et al. 2017

Comparison stage:
Comparing experimental profiles between top and bottom strand profiles due
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

#Set DNA sequence
#Full experimental sequence
TS_seq=Seq('TCGGGCTGGAGATACCCGGTGCTAAGGCCGCTTAATTGGTCGTAGCAAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCTACCGCGTTTTAACCGCCAATAGGATTACTTACTAGTCTCTAGGCACGTGTAAGATATATACATCCTGTCTGCT')

BS_seq=TS_seq.reverse_complement()

prof_datasets=[
	{
	'name':'Experimental cleavage frequencies from top and bottom strands are symmetric, scCSE4-601TA nucleosome',
	'pngfileout':'exp_compar_BS_TS.png',
	'zero_at':79, #move zero on the plot to this position
	'seq':TS_seq, #common sequence for all profiles
	#to align data from different files we set
	#start_on_seq - parameter, 1-based numbering position
	# on 'seq' determines where porfile from corresponding file starts (including any NaNs)
	'rescale':'fit',
	'profiles':[
	{'file':'scCSE4_601TA_TS_fitted_intensities.csv','start_on_seq':49,'datacolumn':'Intensity','name':'Top strand'},
	{'file':'scCSE4_601TA_BS_fitted_intensities.csv','start_on_seq':48,'datacolumn':'Intensity','name':'Bottom strand'},
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

