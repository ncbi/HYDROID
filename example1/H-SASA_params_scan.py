#!/usr/bin/env python
"""
HYDROID (HYDroxyl-Radical fOotprinting Interpretation for DNA)
Example 1, centromeric nucleosome of yeast reconstituted on a well-positioning sequence, Shaytan et al. 2017

HYDROIDpred, Special script:

This example calculates comprehensive set of H-SASA profiles with different sets of atomic radii,
different probe radius and different combinations of deoxyribose hydrogen atom contributions.

The code is modified with respect to other examples to handle averaging over multiple files (not used here),
and scanning through a set of different parameters.

See get_DNA_H_SASA inline help for explanation of parameters.

"""
PARALLEL=True #if Ture - will run CPUs threads.
CPUs=10 #Max number or process to run simultaneously.
#Problems with CPU affinity while using multiprocessing have been reported for some machines
#following command might help to manually assign CPU affinity to all the processes
# "pgrep -u USERNAME  python | while read l; do taskset -p 0xFFFFFFFFFFFFFFFF $l; done"
# To automatically trigger this call in the script, input username below, and uncomment corresponding lines 143, 150
USERNAME="USERNAME"

import os,sys,tempfile
from Bio.Seq import Seq
import sys
import pandas as pd
import numpy as np

if PARALLEL:
	from multiprocessing import Process, Manager

sys.path.insert(0,'../..')

from HYDROID.HYDROIDpred import get_DNA_H_SASA

out_path="results"
str_path="data/structures"
try:
	os.mkdir(out_path)
except:
	pass



#Set DNA sequence
TS_seq=Seq('GGAGATACCCGGTGCTAAGGCCGCTTAATTGGTCGTAGCAAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCTACCGCGTTTTAACCGCCAATAGGATTACTTACTAGTCTCTAGGCACGTGTAAGATATATACATCCT')
BS_seq=TS_seq.reverse_complement()

prof_data=[
	{'prof_name':'nucl_H_Ndist_multparam','pdb_file':['scCSE4_601TA_nucl_H_Ndist.pdb'],'chain':['I','J'],'resids':[range(-71,72),range(-71,72)],'seq':[TS_seq,BS_seq]},
	{'prof_name':'nucl_H_ELdist_multparam','pdb_file':['scCSE4_601TA_nucl_H_ELdist.pdb'],'chain':['I','J'],'resids':[range(-71,72),range(-71,72)],'seq':[TS_seq,BS_seq]},
	{'prof_name':'nucl_H_MD1ns_multparam','pdb_file':['scCSE4_601TA_nucl_H_MD1ns.pdb'],'chain':['I','J'],'resids':[range(-71,72),range(-71,72)],'seq':[TS_seq,BS_seq]},
	{'prof_name':'nucl_H_MD2ns_multparam','pdb_file':['scCSE4_601TA_nucl_H_MD2ns.pdb'],'chain':['I','J'],'resids':[range(-71,72),range(-71,72)],'seq':[TS_seq,BS_seq]},
	]


#Common code:
#will calculate H-SASA profiles from PDB structure with hydrogen atoms.
#And scan the space of parameters. Will average results over different PDB files provided in the list.
#adjust n_threads to the number of CPU cores in your computer for optimal performance
###################################

if not PARALLEL:

	for p in prof_data:
		sum_df=pd.DataFrame() #sum results with different parameters
		vdw_set=['None','charmm36-rmin','amber10-rmin']
		probe_size=[1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0]
		Hcontrib=[[1,0,0,0,0,0,0],[0,1,1,0,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,1,0,0],[0,0,0,0,0,1,1]]
		# hatoms=['H1\'','H2\'','H2\'','H3\'','H4\'','H5\'','H5\'\'']
		for pr_size in probe_size:
			print p['prof_name']," Testing probe size: ",pr_size
			for v_set in vdw_set:
				sum_dff=pd.DataFrame() #we will average over files
				for pdbfile in p['pdb_file']:
					tempout = tempfile.NamedTemporaryFile()
					get_DNA_H_SASA(os.path.join(str_path,pdbfile),tempout,\
					chain=[p['chain'][0]]*5+[p['chain'][1]]*5,resids=[p['resids'][0]]*5+[p['resids'][1]]*5,seq=[p['seq'][0]]*5+[p['seq'][1]]*5,probe_radius=pr_size,slicen=200,vdw_set=v_set,\
					Hcontrib=Hcontrib*2,n_threads=1,verbose=False)
					tempout.seek(0)
					df=pd.read_csv(tempout,comment='#')
					tempout.close()
					df['pdb_file']=pdbfile
					df['probe_size']=pr_size
					df['vdw_set']=v_set
					# df['Hcontrib']='_'.join(map(str,Hcontrib))
					# df['Hcontrib']='full'
					sum_dff=pd.concat([sum_dff,df])
					del df
				df2=sum_dff.groupby(['resid_%i'%i for i in range(10)]+['Site_%i'%i for i in range(10)]+['probe_size','vdw_set'], as_index=False).agg(np.mean)
				del sum_dff
				sum_df=pd.concat([sum_df,df2])
				del df2
		gv=sum_df.groupby(['resid_%i'%i for i in range(10)]+['Site_%i'%i for i in range(10)]+['probe_size','vdw_set'],as_index=False).agg(np.mean)
		gv.to_csv(os.path.join(out_path,p['prof_name']+'_H-SASA.csv,gz'),compression='gzip')

if PARALLEL:

	for p in prof_data:
	
		Hcontrib=[[1,0,0,0,0,0,0],[0,1,1,0,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,1,0,0],[0,0,0,0,0,1,1]]
		vdw_set=['None','charmm36-rmin','amber10-rmin']
		probe_size=[1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0]
		# hatoms=['H1\'','H2\'','H2\'','H3\'','H4\'','H5\'','H5\'\'']

		manager = Manager()
		return_dict = manager.dict()
		jobs = []

		def worker(num,p,pdbfile,pr_size,v_set,return_dict):
			'''worker function'''
			print "Launching calc:",pdbfile," R=",pr_size," ",v_set
			tempout = tempfile.NamedTemporaryFile()
			get_DNA_H_SASA(os.path.join(str_path,pdbfile),tempout,\
				chain=[p['chain'][0]]*5+[p['chain'][1]]*5,resids=[p['resids'][0]]*5+[p['resids'][1]]*5,seq=[p['seq'][0]]*5+[p['seq'][1]]*5,probe_radius=pr_size,slicen=200,vdw_set=v_set,\
				Hcontrib=Hcontrib*2,n_threads=1,verbose=False)
			tempout.seek(0)
			df=pd.read_csv(tempout,comment='#')
			tempout.close()
			df['pdb_file']=pdbfile
			df['probe_size']=pr_size
			df['vdw_set']=v_set
			return_dict[num]=df


		mapd=dict()
		i=0
		for pr_size in probe_size:
			print p['prof_name']," Testing probe size: ",pr_size
			mapd[pr_size]=dict()
			for v_set in vdw_set:
				mapd[pr_size][v_set]=dict()
				for pdbfile in p['pdb_file']:
					mapd[pr_size][v_set][pdbfile]=i
					prc=Process(target=worker, args=(i,p,pdbfile,pr_size,v_set,return_dict))
					jobs.append(prc)
					prc.start()
					# time.sleep(1)
					i=i+1
					if (i%CPUs)==0:
						# os.system("pgrep -u %s  python | while read l; do taskset -p 0xFFFFFFFFFFFFFFFF $l; done"%USERNAME)
						print '!!!!!Paused and waiting for jobs to finish!!!!'
						for proc in jobs:
							print "Waiting for ",proc
							proc.join()
						print "!!!!!!!Staring submitting jobs!!!!!!!"
		
		# os.system("pgrep -u %s  python | while read l; do taskset -p 0xFFFFFFFFFFFFFFFF $l; done"%USERNAME)
		for proc in jobs:
			print "Waiting for ",proc
			proc.join()

		print "Collecting data..."
		sum_df=pd.DataFrame() 
		for pr_size in probe_size:
			print "for probe size ",pr_size
			for v_set in vdw_set:
				sum_dff=pd.DataFrame() #we will average over files
				for pdbfile in p['pdb_file']:
					sum_dff=pd.concat([sum_dff,return_dict[mapd[pr_size][v_set][pdbfile]]])
				df2=sum_dff.groupby(['resid_%i'%i for i in range(10)]+['Site_%i'%i for i in range(10)]+['probe_size','vdw_set'], as_index=False).agg(np.mean)
				del sum_dff
				sum_df=pd.concat([sum_df,df2])
				del df2
		gv=sum_df.groupby(['resid_%i'%i for i in range(10)]+['Site_%i'%i for i in range(10)]+['probe_size','vdw_set'],as_index=False).agg(np.mean)
		gv.to_csv(os.path.join(out_path,p['prof_name']+'_H-SASA.csv.gz'),compression='gzip')
		del return_dict



