#!/usr/bin/env python

import pkgutil
import os
from Bio.Seq import Seq
import tempfile
import os


def test_exp():

	from hydroid.HYDROIDexp import assign_peaks_interactive

	temp = tempfile.NamedTemporaryFile(delete=False)
	temp2 = tempfile.NamedTemporaryFile(delete=False)
	temp.write(pkgutil.get_data('hydroid', 'pkgdata/test_data/lane_profiles.xls'))
	temp2.write(pkgutil.get_data('hydroid', 'pkgdata/test_data/lane_config.csv'))
	temp.seek(0)
	temp2.seek(0)
	lane_profile_file = temp.name
	lane_config_file = temp2.name

	temp.close()
	temp2.close()
	lane_names=['scCSE4_601TA_BS']

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

	os.remove(temp2.name)
	os.remove(temp.name)


def test_pred():

	from hydroid.HYDROIDpred import get_DNA_H_SASA

	out_path=""
	temp = tempfile.NamedTemporaryFile(delete=False)
	temp.write(pkgutil.get_data('hydroid', 'pkgdata/test_data/test.pdb'))
	temp.seek(0)

	try:
		os.mkdir(out_path)
	except:
		pass

	#Set DNA sequence
	TS_seq=Seq('GGAGATACCCGGTGCTAAGGCCGCTTAATTGGTCGTAGCAAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCTACCGCGTTTTAACCGCCAATAGGATTACTTACTAGTCTCTAGGCACGTGTAAGATATATACATCCT')
	BS_seq=TS_seq.reverse_complement()

	prof_data=[
		{'prof_name':'test','pdb_file':temp.name,'chain':'I','resids':range(-71,72),'seq':TS_seq,'vdw_set':'charmm36-rmin'}
		]

	temp.close()

	#Common code:
	#will calculate H-SASA profiles from PDB structure with hydrogen atoms.
	#adjust n_threads to the number of CPU cores in your computer for optimal performance
	###################################
	for p in prof_data:
		get_DNA_H_SASA(p['pdb_file'],os.path.join(out_path,p['prof_name']+'_H-SASA.csv'),\
			chain=p['chain'],resids=p['resids'],seq=p['seq'],probe_radius=1.4,slicen=200,vdw_set=p['vdw_set'],\
			Hcontrib=[1.0]*7,n_threads=6,verbose=False)
	
	os.remove(temp.name)

def HYDROID_get_ex1():

	import urllib, json
	import os
	url="https://api.github.com/repos/ncbi/HYDROID/contents/examples/example1"

	def get_files_from_git(gitapiurl,savefoldername):
		os.mkdir(savefoldername)
		json_url = urllib.urlopen(gitapiurl)
		data = json.loads(json_url.read())
		for d in data:
			if(d['type']=='file'):
				print("Downloading "+os.path.join(savefoldername,d['name']))
				urllib.urlretrieve(d['download_url'],os.path.join(savefoldername,d['name']))
			if(d['type']=='dir'):
				get_files_from_git(d['url'],os.path.join(savefoldername,d['name']))
	get_files_from_git(url,'example1')

def HYDROID_get_ex2():

	import urllib, json
	import os
	url="https://api.github.com/repos/ncbi/HYDROID/contents/examples/example2"

	def get_files_from_git(gitapiurl,savefoldername):
		os.mkdir(savefoldername)
		json_url = urllib.urlopen(gitapiurl)
		data = json.loads(json_url.read())
		for d in data:
			if(d['type']=='file'):
				print("Downloading "+os.path.join(savefoldername,d['name']))
				urllib.urlretrieve(d['download_url'],os.path.join(savefoldername,d['name']))
			if(d['type']=='dir'):
				get_files_from_git(d['url'],os.path.join(savefoldername,d['name']))
	get_files_from_git(url,'example2')

