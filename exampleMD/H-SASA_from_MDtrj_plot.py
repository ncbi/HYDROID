#!/usr/bin/env python
"""
Plotting H-SASA profile calculated from MD trajectory
"""


from Bio.Seq import Seq
import os,sys,tempfile
import pandas as pd

sys.path.insert(0,'../..')

from HYDROID.HYDROIDexp import plot_prof_on_seq

out_path="results"
try:
	os.mkdir(out_path)
except:
	pass

#Set DNA sequence
TS_seq=Seq('GGAGATACCCGGTGCTAAGGCCGCTTAATTGGTCGTAGCAAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCTACCGCGTTTTAACCGCCAATAGGATTACTTACTAGTCTCTAGGCACGTGTAAGATATATACATCCT')
BS_seq=TS_seq.reverse_complement()

#List of plots below, each plot combines serveral profiles from different files
prof_datasets=[
	{
	'name':'H-SASA, scCSE4-601TA nucleosome, Top Strand, Chain I. Average profile from MD trajectory',
	'pngfileout':'scCSE4_601TA_TS_H-sasa_MDavr.png',
	'zero_at':72, #move zero on the plot to this position
	'seq':TS_seq, 
	'profiles':[
	{'file':'scCSE4_601TA_H-SASA_MDavr.csv','start_on_seq':1,'column':'H-SASA_0','name':'H-SASA_MDavr'},
	]},
	{
	'name':'H-SASA, scCSE4-601TA nucleosome, Bottom Strand, Chain J. Average profile from MD trajectory',
	'pngfileout':'scCSE4_601TA_BS_H-sasa_MDavr.png',
	'zero_at':72, #move zero on the plot to this position
	'seq':BS_seq, 
	'profiles':[
	{'file':'scCSE4_601TA_H-SASA_MDavr.csv','start_on_seq':1,'column':'H-SASA_1','name':'H-SASA_MDavr'},
	]},
	]

#Common code:
#Combines profiles from different files into one file/dataframe and sends for plotting
###################################
for sets in prof_datasets:
	data={}
	columns=[]
	#Here we combine datasets to one dataframe for plotting,
	#by simply stacking columns from files together
	#then provide column names to plotting functions as array
	for p in sets['profiles']:
		data[p['file']]=pd.read_csv(os.path.join(out_path,p['file']),comment='#')
		data[p['file']].index=range(1,1+len(data[p['file']]))
		data[p['file']].rename(columns={p['column']:p['name'],'Site':'Site_'+p['name']}, inplace=True)
		columns.append(p['name'])
	site_column_df=pd.DataFrame({'Site':['%d%s'%(n,l) for n,l in zip(range(1,1+len(sets['seq'])),sets['seq'])]},index=range(1,1+len(sets['seq'])))
	result=pd.concat(data.values()+[site_column_df],axis=1)
		
	temp = tempfile.TemporaryFile()
	result.to_csv(temp)
	temp.seek(0)
	plot_prof_on_seq(temp,DNAseq=sets['seq'],\
		graphshow=True,pngfileout=os.path.join(out_path,sets['pngfileout']),title=sets['name'],\
		prof_columns=columns,seq_column="Site",normalize='together',zero_at=sets['zero_at'],\
		colorb={'A':'#0b0','T':'#00b','G':'#fff','C':'#fff'},colorf={'A':'#fafafa','T':'#fafafa','G':'#000','C':'#000'})
	temp.close()





