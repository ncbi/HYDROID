#!/usr/bin/env python
"""
This is a library to predict expected DNA hydroxyl-radical footprinting profiles
from PDB structures through solven accessible surface area calcualtions.

DNA cleavage by hydroxyl-radicals is thought to proceed through abstraction of
deoxyribose hydrogen atoms. The probablilty of abstraction and hence DNA cleavage
is estimated by calculating the solvent accessible area of these hydrogen atoms resulting
in an H-SASA profile along the DNA sequence.

This library depends on FREESASA library http://freesasa.github.io

"""
import pkgutil
import pandas as pd
import numpy as np
import tempfile
import freesasa

def	get_DNA_H_SASA(pdb_file,csvfileout,chain=None,resids=[],seq=None,probe_radius=1.4,slicen=100,vdw_set=None,Hcontrib=[1.0]*7,n_threads=1,verbose=False):
	"""
	Function is a warapper to the FREESASA library to calculate the Surface Accessible Surface Area out
	atoms in pdb_file, then expreacts the SASA deoxiribose hydrogen atoms and sums it up
	for every nucleotide with coefficients Hcontrib.
	chain - name of the DNA chain of interest in pdb_file, if chain has no name leave blank ('')
	resids - a list of resids to calculate H-SASA values.
	seq - seqeunce of the DNA strand, string or biopython Seq object.
	Hcontrib - coefficients for individual SASA of deoxyribose hydrogens for summing them up into H-SASA profile,
		order [H1' H2' H2'' H3' H4' H5' H5'']
	Note: chain, resids, seq, Hcontrib - can be also a list of two or more instances,
			to make calculation for several chains, spans of resids or combinations of Hcontrib at once.
			In this case number of elements in chain, resids, Hcontrib should be the same,
			and the algorithm will iterate through all list simultaneously (i.e. no combination will be tried).
			Chains should be of the same length.
	probe_radius - size of probe to roll.
	slicen - number of slices per atom, controls precision of the calculation.
	vdw_set - seleting the set of VdW radii:
		None - default for FREESASA used
		charmm36-rmin - rmin from charmm36 forcefield
		abmer10-rmin - rmin from AMBER10 forcefield


	Return
	--------
	CSV file csvfileout with columns of H-SASA profiles along the sequence.
	"""
	chain=[chain] if isinstance(chain,basestring) else list(chain)
	if len(chain)>1:
		assert len(chain)==len(resids)
		assert len(chain)==len(seq)
		assert len(chain)==len(Hcontrib)
	else:
		resids=[resids]
		seq=[seq]
		Hcontrib=[Hcontrib]

	if not verbose:
		freesasa.setVerbosity(freesasa.nowarnings)
	hatoms=['H1\'','H2\'','H2\'','H3\'','H4\'','H5\'','H5\'\'']

	if vdw_set=='charmm36-rmin':
		#Open config from package in a tricky way, independent of package installation mode
		temp2 = tempfile.NamedTemporaryFile()
		conffile = pkgutil.get_data('HYDROID', 'pkgdata/charmm36_rmin.config')
		temp2.write(conffile)
		temp2.seek(0)
		classifier = freesasa.Classifier(temp2.name)
		temp2.close()
		####
		structure = freesasa.Structure(pdb_file,classifier, options={'hydrogen' : True,'hetatm' : True})
	elif vdw_set=='amber10-rmin':
		#Open config from package in a tricky way, independent of package installation mode
		temp2 = tempfile.NamedTemporaryFile()
		conffile = pkgutil.get_data('HYDROID', 'pkgdata/amber10_rmin.config')
		temp2.write(conffile)
		temp2.seek(0)
		classifier = freesasa.Classifier(temp2.name)
		temp2.close()
		####
		structure = freesasa.Structure(pdb_file,classifier, options={'hydrogen' : True,'hetatm' : True})
	else:
		structure = freesasa.Structure(pdb_file, options={'hydrogen' : True,'hetatm' : True})
	print "Launching FreeSASA calculation..."
	result = freesasa.calc(structure,freesasa.Parameters({'algorithm' : freesasa.LeeRichards,'n-slices' : slicen,'probe-radius':probe_radius,'n-threads':n_threads}))
	# result = freesasa.calc(structure,freesasa.Parameters({'algorithm' : freesasa.ShrakeRupley,'n-slices' : slicen,'n-threads':n_threads}))
	print "Calculation done"
	
	print "Extracting SASA values ..."
	
	res=dict()
	
	for ch,rids,Hcont,i in zip(chain,resids,Hcontrib,range(len(chain))):
		res[i]=pd.Series()

		if (np.array(Hcont)==1.0).all():
		#simplified procedure, we can do it faster: we need to calculate all H-SASA at once
			sels=[]
			for resid in rids:
				if len(ch)>0:
					sels.append('%d,(chain %s) and (resi %s%d) and (name %s)'%(resid, ch,'\\' if resid<0 else '', resid, '+'.join(hatoms)))
				else:
					sels.append('%d,(resi %s%d) and (name %s)'%(resid,'\\' if resid<0 else '', resid, '+'.join(hatoms)))
			selections = freesasa.selectArea(sels,structure, result)
			res[i]=res[i].add(pd.Series(selections)*1.0,fill_value=0)
		else:
		#regular procedure
			for hat,hcont in zip(hatoms,Hcont):
				sels=[]
				if hcont!=0:
					for resid in rids:
						if len(ch)>0:
							sels.append('%d,(chain %s) and (resi %s%d) and (name %s)'%(resid, ch,'\\' if resid<0 else '', resid, hat))
						else:
							sels.append('%d,(resi %s%d) and (name %s)'%(resid,'\\' if resid<0 else '', resid, hat))
				selections = freesasa.selectArea(sels,structure, result)
				res[i]=res[i].add(pd.Series(selections)*float(hcont),fill_value=0)

	for i in range(len(chain)):
		res[i].index=res[i].index.map(int)
		res[i]=res[i].sort_index()
	if len(chain)==1:
		df=pd.DataFrame({'resid':res[0].index,'Site':['%d%s'%(n,l) for n,l in zip(range(1,1+len(seq[0])),seq[0])],'H-SASA':res[0].values})
	else:
		df=pd.DataFrame()
		for ch,i in zip(chain,range(len(chain))):
			# print res[i]
			# print seq[i]
			ndf=pd.DataFrame({'resid_%d'%i:res[i].index,'Site_%d'%i:['%d%s'%(n,l) for n,l in zip(range(1,1+len(seq[i])),seq[i])],'H-SASA_%d'%i:res[i].values})
			df=pd.concat([df,ndf],axis=1)
	print "Outputting H-SASA profile to %s"%csvfileout
	df.to_csv(csvfileout)

