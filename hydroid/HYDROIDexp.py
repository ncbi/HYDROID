#!/usr/bin/env python
"""
This is a library to analyze lane profile from DNA gel electrophoresis
of OH radical footprinting experiments.

Plot lane profiles,
identify peaks semi automatically,
fit Lorentzians into peaks,
extract total intensities.

"""
import time
import tempfile
import os
import pkgutil
import collections

import numpy as np
import pandas as pd
import peakutils
import statsmodels.api as sm


import matplotlib.pyplot as plt
from matplotlib import rc,gridspec,rcParams
from matplotlib.widgets import Slider, Button, CheckButtons, RadioButtons

from scipy.optimize import minimize 
from scipy.optimize import leastsq # Levenberg-Marquadt Algorithm #
from scipy.ndimage.filters import gaussian_filter as gfilt


from PIL import Image, ImageDraw, ImageFont

#####Fucntions to read the original data files and read/write parameter files ######
#####

def read_profile(lane_profile_file,lane_config_file,lane_name):
	"""
	Reads dataset with lane profiles (lane_profile_data) and configuration file (lane_config_file),
	returns one profile specified by lane_name (lname column value in configuration file) as PANDAS series.
	"""
	confdf=pd.read_csv(lane_config_file,skipinitialspace=False,delimiter=",\s*",engine='python',comment='#').drop_duplicates(subset='lname',keep='last').set_index('lname',drop=False)
	profdf=pd.read_csv(lane_profile_file,delimiter="\t",engine='python')
	col_name=confdf.loc[lane_name,'column']
	return profdf[col_name]

def read_prof_configs(lane_config_file,lane_name):
	"""
	Reads all parameters from lane_config_file for a given lane_name. See an example of such file for the description of parameters.
	"""
	confdf=pd.read_csv(lane_config_file,skipinitialspace=False,delimiter=r",\s*",engine='python',comment='#',header=0,skip_blank_lines=True).drop_duplicates(subset='lname',keep='last').set_index('lname',drop=False)
	params=confdf.loc[lane_name].to_dict()
	params['delpeaks']=map(lambda x: int(float(x)),str(params['delpeaks']).split()) if not pd.isnull(params['delpeaks']) else []
	params['addpeaks']=map(lambda x: int(float(x)),str(params['addpeaks']).split()) if not pd.isnull(params['addpeaks']) else []
	params['seqpos_num']=int(params['seqpos'][0:-1]) if not pd.isnull(params['seqpos']) else np.nan #The nucleotide is attacked and distructed by OH-radical
	return params

def guess_missing_params(params,lane_length):
	"""
	If certain parameters in lane_config_file are not defined this will guess some values for them to start with.
	"""
	params['alignpos']=lane_length-1 if params['alignpos']=='end' or params['alignpos']=='NaN' or pd.isnull(params['alignpos']) else int(float(params['alignpos']))
	params['seqpeak']=lane_length-1 if pd.isnull(params['seqpeak']) else int(params['seqpeak'])
	params['seqpos_num']=1 if pd.isnull(params['seqpos_num']) else params['seqpos_num']
	params['peakthresh']=0.01 if pd.isnull(params['peakthresh']) else params['peakthresh']
	params['min_dist_right']=lane_length/100 if pd.isnull(params['min_dist_right']) else params['min_dist_right']
	params['min_dist_left']=lane_length/100 if pd.isnull(params['min_dist_left']) else params['min_dist_left']
	params['interpolate']=False if pd.isnull(params['interpolate']) else params['interpolate']
	params['leftlim']=0 if pd.isnull(params['leftlim']) else int(params['leftlim'])
	params['rightlim']=lane_length-1 if pd.isnull(params['rightlim']) else int(params['rightlim'])
	params['segments']=2 if pd.isnull(params['segments']) else params['segments']
	params['base']=False if pd.isnull(params['base']) else params['base']

	return params

def modify_prof_configs(lane_config_file,lane_name,params):
	"""
	Modifies one line in lane_config_file for lane_name with new params.
	Comments out the previous line
	Currently only params modified interactively are supported.
	Interpolation, peakthresh,leftlim,rightlim,min_dist_left, min_dist_right
	"""
	confdf=pd.read_csv(lane_config_file,skipinitialspace=False,delimiter=r",\s*",engine='python',comment='#',header=0,skip_blank_lines=True).drop_duplicates(subset='lname',keep='last').set_index('lname',drop=False)
	pparams=params.copy() #convert some arrays to strings
	pparams['delpeaks']=' '.join(map(str,pparams['delpeaks']))
	pparams['addpeaks']=' '.join(map(str,pparams['addpeaks']))
	df=pd.DataFrame(pparams,index=[lane_name])
	confdf.update(df)
	with open(lane_config_file,'a') as file:
		file.write('#Saved updated values for %s from interactive window on %s\n'%(lane_name,time.strftime("%d, %b %Y, %H:%M:%S")))
	confdf.loc[lane_name:lane_name].to_csv(lane_config_file,na_rep='NaN',mode='a',index=False,header=False)
	return 0


	with open(lane_config_file,'r') as file:
		filedata=file.read()
	fdata=filedata.split('\n')
	comments=[]
	hline=0
	for f,i in zip(fdata,range(len(fdata))):
		if[0]=='#':
			comments.append(f)
			hline=i+1
	header=fdata[hline]
	for l,ind in zip(fdata,range(len(fdata))):
		if l.split(',')[1] in lane_name:
			index=ind
			break
	newline=''
	for t,v,p in zip(re.split(r',\s*',ftitle),re.split(r',\s*',l),re.split(r'\S+,',l)[1:]):
		npar=params.get(t,v)
		print t," ",v," ",p
		if t=='addpeaks' or t=='delpeaks':
			npar='['+'_'.join(map(str,npar))+']'
		newline+=('%.3f'%npar if isinstance(npar,float) else str(npar))+','+p
	print l
	print newline
	fdata[index]=newline
	with open(paramfile,'w') as file:
		file.write('\n'.join(fdata))

#####Fucntions to assign the peaks on a profile, compare profiles and call peaks (assign to positions on DNA sequence)######
#####

def assign_peaks_interactive(lane_profile_file,lane_config_file,lane_name,title=None,qthresh=100,pngfileout=None):
	"""
	Opens an interactive window to adjust parameters for automatic peak location identification.
	lane_profile_file,lane_config_file - files with data and configs - see example provided.
	lane_name - name of the lane from lane_config_file
	qthresh - limit data range view to qthresh percentile of the values.
	fileout - is not None an image would be written to file instead of an interactive figure.
	"""
	lane=read_profile(lane_profile_file,lane_config_file,lane_name)
	params=read_prof_configs(lane_config_file,lane_name)
	lane=lane[~np.isnan(lane)].values #remove any missing values, usually found at the end if profiles were of different length and convert to NUMPY ndarray
	params=guess_missing_params(params,len(lane)) #Guess if certain params are NaN, or use "end" keyword for alignpos.
	peaks,laneSM=find_peaks(lane, params['peakthresh'], params['leftlim'], params['rightlim'],params['base'],params['segments'],params['min_dist_left'],params['min_dist_right'],params['addpeaks'],params['delpeaks'],params['interpolate'])

	global scplot,liml,limr,anchp,seqp,seqpeakannot
	#Generate a plot of the profile and gel mockup
	if not title:
		title='Interactive peak assignment window, Lane: %s'%lane_name
	fig1=plt.figure(figsize=(15, 6))
	gs = gridspec.GridSpec(3, 1,height_ratios=[14,1,5])
	ax1 = fig1.add_subplot(gs[0])
	ax2 = fig1.add_subplot(gs[1])
	ax1.plot(range(len(lane)),lane, linewidth=1, label=lane_name)
	scplot=ax1.scatter(peaks, lane[peaks.astype(int)], marker='x', color='g', s=40,linewidth=2)
	ax1.set_ylabel("Intensity")
	ax1.set_title(title)
	ax1.set_ylim((0,np.percentile(lane,qthresh)*1.05))
	ax1.set_xlim((0,len(lane)))
	ax1.legend(loc="upper right")
	# seqpeakannot=ax1.annotate('Known DNA position: %s'%params['seqpos'], xy=(params['seqpeak'], lane[int(params['seqpeak'])]+0.05*max(lane)), xytext=(params['seqpeak'], lane[int(params['seqpeak'])]+0.2*max(lane)), arrowprops=dict(facecolor='black', shrink=0.01,width=2))

	fakegelimg=np.tile(lane,(100,1))
	ax2.imshow(fakegelimg,cmap='Greys',vmin=0,vmax=np.percentile(lane,qthresh),interpolation='nearest', aspect='auto')
	ax2.set_xlim((0,len(lane)))
	ax2.set_xticks([])
	ax2.set_yticks([])
	try: #try to mark data set limits set in the config file and alignpos if provided.
		liml=ax1.axvline(x=params['leftlim'],color='red')
		limr=ax1.axvline(x=params['rightlim'],color='red')
		anchp=ax1.axvline(x=params['alignpos'],color='pink')
		# seqp=ax1.axvline(x=params['seqpeak'],color='green')
	except:
		pass

	#Block below implements interactive controls for changing many parameters,
	#and recalculate peak positions upon their change.
	#The Save button allows to update lane_config_file with these parameters.
	axcolor = 'lightgoldenrodyellow'
	axleftlim = plt.axes([0.35, 0.22, 0.55, 0.025], axisbg="red")
	axrightlim = plt.axes([0.35, 0.185, 0.55, 0.025], axisbg="red")
	axpeakthresh = plt.axes([0.35, 0.15, 0.55, 0.025], axisbg=axcolor)
	axmin_dist_left= plt.axes([0.35, 0.115, 0.55, 0.025], axisbg=axcolor)
	axmin_dist_right = plt.axes([0.35, 0.08, 0.55, 0.025], axisbg=axcolor)
	axalignpos = plt.axes([0.35, 0.045, 0.55, 0.025], axisbg="pink")
	# axseqpeak = plt.axes([0.35, 0.01, 0.55, 0.025], axisbg="green")

	rax = plt.axes([0.052, 0.165, 0.1, 0.08])
	check = CheckButtons(rax, ['Interpolate'], [params['interpolate']])

	rax2 = plt.axes([0.152, 0.165, 0.1, 0.08])
	check2 = CheckButtons(rax2, ['Baseline'], [params['base']])

	rax3 = plt.axes([0.152, 0.01, 0.1, 0.15])
	radio = RadioButtons(rax3, ('2 Segments','3 Segments','4 Segments','5 Segments'),active=params['segments']-2)

	save = plt.axes([0.052, 0.01, 0.1, 0.05])
	button = Button(save, 'Save', color='0.8', hovercolor='0.975')
	
	leftlim = Slider(axleftlim, 'leftlim', 0, len(lane)-1, valinit=params['leftlim'])
	rightlim = Slider(axrightlim, 'rightlim', 0, len(lane)-1, valinit=params['rightlim'])
	peakthresh = Slider(axpeakthresh, 'peakthresh', 0, 0.50, valinit=params['peakthresh'])
	min_dist_left = Slider(axmin_dist_left, 'min_dist_left', 0, 100, valinit=params['min_dist_left'])
	min_dist_right = Slider(axmin_dist_right, 'min_dist_right', 0, 100, valinit=params['min_dist_right'])
	alignpos = Slider(axalignpos, 'alignpos', 0, len(lane)-1, valinit=params['alignpos'])
	# seqpeak = Slider(axseqpeak, 'seqpeak', 0, len(lane)-1, valinit=params['seqpeak'])

	def save(event): #This is where the lane_config_file is updated
		print 'Modifying parameters in ',lane_config_file
		print 'for ',lane_name
		modify_prof_configs(lane_config_file,lane_name,params)
	button.on_clicked(save)

	def update(val):
		if val == 'Interpolate':
			params['interpolate']= not params['interpolate']
		if val == 'Baseline':
			params['base']= not params['base']
		if isinstance(val, basestring):
			if "Segments" in val:
				params['segments']=int(val[0])

		global scplot,liml,limr,anchp,seqp,seqpeakannot
		scplot.remove()
		liml.remove()
		limr.remove()
		anchp.remove()
		# seqp.remove()
		# seqpeakannot.remove()
		params['min_dist_left']=min_dist_left.val
		params['min_dist_right']=min_dist_right.val
		params['peakthresh']=peakthresh.val
		params['leftlim']=int(leftlim.val)
		params['rightlim']=int(rightlim.val)
		params['alignpos']=int(alignpos.val)
		# params['seqpeak']=int(seqpeak.val)


		peaks,laneSM=find_peaks(lane, params['peakthresh'], params['leftlim'], params['rightlim'],params['base'],params['segments'],params['min_dist_left'],params['min_dist_right'],params['addpeaks'],params['delpeaks'],params['interpolate'])
		scplot=ax1.scatter(peaks, lane[peaks.astype(int)], marker='x', color='g', s=40,linewidth=2)
		try:
			liml=ax1.axvline(x=params['leftlim'],color='red')
			limr=ax1.axvline(x=params['rightlim'],color='red')
			anchp=ax1.axvline(x=params['alignpos'],color='pink')
			# seqp=ax1.axvline(x=params['seqpeak'],color='green')
			# seqpeakannot=ax1.annotate('Known DNA position: %s'%params['seqpos'], xy=(params['seqpeak'], lane[int(params['seqpeak'])]+0.05*max(lane)), xytext=(params['seqpeak'], lane[int(params['seqpeak'])]+0.2*max(lane)), arrowprops=dict(facecolor='black', shrink=0.01,width=2))
		except:
			pass

	min_dist_left.on_changed(update)
	min_dist_right.on_changed(update)
	leftlim.on_changed(update)
	rightlim.on_changed(update)
	peakthresh.on_changed(update)
	# seqpeak.on_changed(update)
	alignpos.on_changed(update)
	check.on_clicked(update)
	check2.on_clicked(update)
	radio.on_clicked(update)
	#========
	fig1.tight_layout()
	if(pngfileout):
		plt.savefig(pngfileout,dpi=(600))
	else:
		plt.show()

	return 0

def find_peaks(lane,peakthresh=0.005,leftlim=0,rightlim=4000,base=False,segments=3,min_dist_left=5,min_dist_right=50,addpeaks=[],delpeaks=[],interpolate=False,interp_thresh=1.25):
	"""
	Implements the peak finding algorithm.
	Takes a lane profile as a 1D numpy ndarray.
	Returns indexes of peaks and lane profile with baseline substracted (if requested by base=True)
	"""
	print "=======Invoking peak identification algorithm======"
	peaks=[] #This list will be filled up with peak indexes

	if base: # remove baseline if requested
		baseline = peakutils.baseline(lane, 2)
		lane=lane-baseline

	segm=int(segments) #data range will be split into a number of segments if requested, and min_dist parameter will be interpolated for each segment
	segml=int(len(lane)/segm)
	for i in range(segm):
		if segm>1:
			md=(min_dist_right-min_dist_left)/(segm-1)*(i)+min_dist_left
		else:
			md=(min_dist_left+min_dist_right)/2
		p=peakutils.indexes(lane, thres=peakthresh, min_dist=md)
		for k in p:
			if k>i*segml and k<=segml*(i+1):
				peaks.append(k)

	for p in delpeaks: #Remove peaks that were manually specified for deletion in lane configuration file
		try:
			peaks.remove(int(p))
		except:
			print "Failed delete peak ", p
			

	for p in addpeaks: #Add peaks that were manually specified for addition in lane configuration file
		peaks.append(int(p))
	peaks.sort()

	peaks=np.array(peaks) #Convert to NUMPY ndarray

	#Try to find missed peaks by interpolation if requested.
	if interpolate:
		print "==Interpolation triggered=="
		dointerp=True
		#this is a cycle run iteratively until no new peaks are found by interpolation
		while dointerp:
			peaksL=peaks[(peaks < rightlim)&(peaks > leftlim)]
			increments=peaksL[1:]-peaksL[0:-1]
			try:
				res=np.polyfit(peaksL[0:-1].astype(float),increments.astype(float),1)
				fitvals=peaksL[0:-1]*res[0]+res[1]
			except:
				fitvals=increments.astype(float)
			msd=((increments-fitvals)**2).sum()
			interp_peaks=[]
			for p,pn,v in zip(peaksL[0:-1],peaksL[1:],increments/fitvals):
				if v>=interp_thresh:
					#trigger interpolation attempts for this segment
					fold=int(v+1.0)
					optsegpeaks=[]
					msdinit=msd
					for num_segm in range(2,fold+2):
						#trying to insert num_ins points
						segpeaks=[]
						for k in range(1,num_segm):
							segpeaks.append(int(p+k*(pn-p)/num_segm))
						#test if this yields better residual
						peaksI=np.append(peaksL,segpeaks)
						peaksI.sort()
						incrementsI=peaksI[1:]-peaksI[0:-1]
						fitvalsI=peaksI[0:-1]*res[0]+res[1]
						msdI=((incrementsI-fitvalsI)**2).sum()
						if msdI<msdinit:
							#means successful insertion
							msdinit=msdI
							optsegpeaks=segpeaks
							#try next attempt - inserting more
						else:
							break
					interp_peaks.extend(optsegpeaks)
			peaks=np.append(peaks,interp_peaks)
			peaks.sort()
			# print interp_peaks           
			dointerp=False
		peaks.sort()
	print "Peaks identified: ",len(peaks)
	# print peaks
	print "==============="
	return peaks,lane

def call_peaks_interactive(lane_profile_file,lane_config_file,DNAseq=None,labeled_end='three_prime',lane_name=None,helper_prof_names=None, normalize_profiles=True,title=None,pngfileout=None):
	"""
	Opens an interactive window to overlay footprinting profiles with helper seqeuncing profiles
	and allows to interactively assign any peak to a location on the DNA seqeunce.
	DNAseq - DNA sequence in usual 5'-3' format.
	labeled_end='three_prime' or 'five_prime' - the labeled end of the DNA
	(Digitized profiles are assumed to run from the top of the gel to the bottom, and not vice versa).
	"""
	lane=read_profile(lane_profile_file,lane_config_file,lane_name)
	params=read_prof_configs(lane_config_file,lane_name)
	lane=lane[~np.isnan(lane)].values #remove any missing values, usually found at the end if profiles were of different length and convert to NUMPY ndarray
	params=guess_missing_params(params,len(lane)) #Guess if certain params are NaN, or use "end" keyword for alignpos.
	peaks,laneSM=find_peaks(lane, params['peakthresh'], params['leftlim'], params['rightlim'],params['base'],params['segments'],params['min_dist_left'],params['min_dist_right'],params['addpeaks'],params['delpeaks'],params['interpolate'])

	helper_lanes={}
	helper_lanes_orig_array={}
	helper_aligned_lanes={}
	helper_params={}
	helper_params_changed=set()
	for l in helper_prof_names:
		helper_lanes[l]=read_profile(lane_profile_file,lane_config_file,l)
		helper_lanes_orig_array[l]=helper_lanes[l][~np.isnan(helper_lanes[l])].values
		helper_params[l]=guess_missing_params(read_prof_configs(lane_config_file,l),len(helper_lanes_orig_array[l]))
		if(normalize_profiles):
			helper_lanes[l]=helper_lanes[l]/max(helper_lanes[l])*max(lane)
		#Shifting the positions of the helper lanes using alignpos values.
		delta=params['alignpos']-helper_params[l]['alignpos']
		helper_lanes[l].index=range(0+delta,len(helper_lanes[l])+delta)
		helper_aligned_lanes[l]=helper_lanes[l].reindex(range(0,len(lane))).values

	global anchp,seqp,seqpeakannot,hplots
	hplots={}
	#Generate a plot of the footprinting profile together with helper profiles.
	#All profiles will be aligned using alingnpeak value at the very righ side of the plot.
	#Sliders will allow to change alignpos value for every profile and thus to slide profiles independently.
	if not title:
		title='Interactive peak calling window, Footprinting lane: %s'%lane_name
	fig1=plt.figure(figsize=(15, 6))
	gs = gridspec.GridSpec(2, 1,height_ratios=[15,5])
	ax1 = fig1.add_subplot(gs[0])
	ax1.plot(range(len(lane)),lane, linewidth=1, label=lane_name,color='C0')
	ax1.set_ylabel("Intensity of main profile")
	ax1.set_title(title)
	ax1.set_ylim((0,max(lane)*1.05))
	ax1.set_xlim((0,len(lane)))
	scplot=ax1.scatter(peaks, lane[peaks.astype(int)], marker='x', color='g', s=40,linewidth=2)
	seqpeakannot=ax1.annotate('Known DNA position:\n %s'%params['seqpos'], xy=(params['seqpeak'], lane[int(params['seqpeak'])]+0.05*max(lane)), xytext=(params['seqpeak'], max(lane)*0.95), arrowprops=dict(facecolor='black', shrink=0.01,width=2))


	hprofscolors=[]
	for l,i in zip(helper_prof_names,range(len(helper_prof_names))):
		hplots[l],=ax1.plot(range(len(helper_aligned_lanes[l])),helper_aligned_lanes[l], linewidth=1, label=l, color='C%d'%(i+1))
	ax1.legend(loc="upper right")

	try: #try to mark data set limits set in the config file and alignpos if provided.
		anchp=ax1.axvline(x=params['alignpos'],color='pink')
		seqp=ax1.axvline(x=params['seqpeak'],color='green')
		liml=ax1.axvline(x=params['leftlim'],color='red')
		limr=ax1.axvline(x=params['rightlim'],color='red')
	except:
		pass

	#Block below implements interactive controls for changing position of seqeunce peak,
	#and interactive aligning of the profiles (changing alignpos)
	#The Save button allows to update lane_config_file with these parameters.
	axcolor = 'lightgoldenrodyellow'
	axalignpos_helper={}
	alignpos_helper={}
	for l,i in zip(helper_prof_names,range(len(helper_prof_names))):
		axalignpos_helper[l] = plt.axes([0.3, 0.095+i*0.035, 0.6, 0.025], axisbg="pink")
		alignpos_helper[l] = Slider(axalignpos_helper[l], '%s alignpos'%l, 0, len(helper_lanes_orig_array[l])-1, valinit=helper_params[l]['alignpos'])


	axseqpeak = plt.axes([0.3, 0.06, 0.6, 0.025], axisbg="green")
	axseq = plt.axes([0.3, 0.01, 0.6, 0.025], axisbg="green")


	save = plt.axes([0.052, 0.01, 0.1, 0.05])
	button = Button(save, 'Save', color='0.8', hovercolor='0.975')
	
	seqpeak = Slider(axseqpeak, 'seqpeak', 0, len(lane)-1, valinit=params['seqpeak'])
	seqpos_num = Slider(axseq, 'sequence', 1, len(DNAseq), valinit=params['seqpos_num'],valfmt='%d')
	if labeled_end=='three_prime':
		seq_text=plt.figtext(0.01,0.22,"DNA(5'>3'): "+str(DNAseq)[0:params['seqpos_num']-1].lower()+str(DNAseq)[params['seqpos_num']-1]+str(DNAseq)[params['seqpos_num']:].lower(),wrap=True)
	elif labeled_end=='five_prime':
		seq_text=plt.figtext(0.01,0.22,"DNA(3'>5'): "+(str(DNAseq)[0:params['seqpos_num']-1].lower()+str(DNAseq)[params['seqpos_num']-1]+str(DNAseq)[params['seqpos_num']:].lower())[::-1],wrap=True)
	else:
		raise Exception("labeled_end should be three_prime of five_prime")


	def save(event): #This is where the lane_config_file is updated
		print 'Modifying parameters in ',lane_config_file
		print 'for ',lane_name
		modify_prof_configs(lane_config_file,lane_name,params)
		for l in helper_params_changed:
			print 'Modifying parameters in ',lane_config_file
			print 'for ',l
			modify_prof_configs(lane_config_file,l,helper_params[l])
	button.on_clicked(save)

	def update(val):
		global seqp,seqpeakannot,hplots
	# 	anchp.remove()
		seqp.remove()
		seqpeakannot.remove()
		for l,i in zip(helper_prof_names,range(len(helper_prof_names))):
			hplots[l].remove()
			if helper_params[l]['alignpos']!=int(alignpos_helper[l].val):
				helper_params[l]['alignpos']=int(alignpos_helper[l].val)
				helper_params_changed.add(l)
			delta=params['alignpos']-helper_params[l]['alignpos']
			helper_lanes[l].index=range(0+delta,len(helper_lanes[l])+delta)
			helper_aligned_lanes[l]=helper_lanes[l].reindex(range(0,len(lane))).values
			hplots[l],=ax1.plot(range(len(helper_aligned_lanes[l])),helper_aligned_lanes[l], linewidth=1, label=l, color='C%d'%(i+1))

	# 	params['alignpos']=int(alignpos.val)
		params['seqpeak']=int(seqpeak.val)
		params['seqpos_num']=int(seqpos_num.val)
		params['seqpos']='%d%s'%(params['seqpos_num'],DNAseq[params['seqpos_num']-1])
		if labeled_end=='three_prime':
			seq_text.set_text("DNA(5'>3'): "+str(DNAseq)[0:params['seqpos_num']-1].lower()+str(DNAseq)[params['seqpos_num']-1]+str(DNAseq)[params['seqpos_num']:].lower())
		elif labeled_end=='five_prime':
			seq_text.set_text("DNA(3'>5'): "+(str(DNAseq)[0:params['seqpos_num']-1].lower()+str(DNAseq)[params['seqpos_num']-1]+str(DNAseq)[params['seqpos_num']:].lower())[::-1])
		else:
			raise Exception("labeled_end should be three_prime of five_prime")


		try:
	# 		anchp=ax1.axvline(x=params['alignpos'],color='pink')
			seqp=ax1.axvline(x=params['seqpeak'],color='green')
			seqpeakannot=ax1.annotate('Known DNA position:\n %s'%params['seqpos'], xy=(params['seqpeak'], lane[int(params['seqpeak'])]+0.05*max(lane)), xytext=(params['seqpeak'], max(lane)*0.95), arrowprops=dict(facecolor='black', shrink=0.01,width=2))
		except:
			pass

	seqpeak.on_changed(update)
	seqpos_num.on_changed(update)
	for l in helper_prof_names:
		alignpos_helper[l].on_changed(update)
	# alignpos.on_changed(update)
	#========
	fig1.tight_layout()
	if(pngfileout):
		plt.savefig(pngfileout,dpi=(600))
	else:
		plt.show()

	return 0

#####Fucntions to model fitting into footprinting profiles ######
#####

##Elementary functions (Gaussian/Lorenzians) and arrays of these functions

def gaussian(x,sigma,x_0):
	"""
	G(s,sigma,x_0)=1/(2*PI*sigma**2)**0.5*exp(-(x-x_0)**2/2/sigma)
	"""
	k=-((x-x_0)/sigma)**2/2
	y = np.exp(k)/(np.abs(sigma)*(3.14159*2)**0.5)

	return y

def lorentzian(x,g_2,x_0): 
	"""
	L(x,Gamma/2,x_0)=(1/PI)*(Gamma/2/((x-x_0)**2+(Gamma/2)**2))
	when g_2 == sigma. Gaussian and Lorentzian have similar shapes
	g_2 is Gamma/2
	"""
	numerator =  g_2
	denominator = ( x - x_0 )**2 + g_2**2
	y = (numerator/denominator)/3.14159
	return y

def sigma_to_alpha(sigma_vector):
	"""
	Converting a vector of with parameters (sigma or gamma_2) to vector of alpha
	This variable substitution allows to constrain width of the bands to increase from left to right during fitting.
	(Sigma)_k=SUM_{i=1}^k alpha_i**2
	alpha_k=((Sigma)_{k}-(Sigma)_{k-1})**0.5
	Sigma_0=0 (technical hack)
	"""
	return (sigma_vector-np.concatenate((np.array([0]),sigma_vector[0:-1])))**0.5

def alpha_to_sigma(alpha_vector):
	"""
	Inverse of sigma_to_alpha
	(Sigma)_k=SUM_{i=1}^k alpha_i**2
	"""
	return np.dot(np.tril(np.ones((len(alpha_vector),len(alpha_vector))),0),(alpha_vector**2))

def pwsigma2alpha(p,M0):
	return np.array([p[0::3],sigma_to_alpha(p[1::3]),p[2::3]]).T.reshape(-1)

def pwalpha2sigma(p,M0):
	return np.array([p[0::3],alpha_to_sigma(p[1::3]),p[2::3]]).T.reshape(-1)

def sum_peakfuncs(x,p,M0,func,pwa2sigmafunc=lambda x,M0: x, Nauxvirtpeaks=0):
	"""
	For a given vector x, return a vector of values of function described by sum of Gaussinan/Lorentzians
	centered at provided positions.
	func - reference to any bell shaped function taking (x, width, position) as input.
	p - vector of function parameters [H,width,position]*N, N - number of functions, H multipier for each function.
	width can be substitited by a vector of other parameters if a conversion function
	p2sigmafunc is provided. For example, to do (Sigma)_k=SUM_{i=1}^k alpha_i**2 substitution.
	Nauxvirtpeaks - adds virtual peaks to the left and to the right of the data range as copies of the last peaks in the data range
	"""
	
	# sigmas=pwa2sigmafunc(p.reshape(-1,3)[:,1],M0)
	pn=pwa2sigmafunc(p,M0)
	if Nauxvirtpeaks>0:
		lh=pn[-3]
		lw=pn[-2]
		lp=pn[-1]
		ld=pn[-1]-pn[-4]
		for i in range(Nauxvirtpeaks+1):
			pn=np.append(pn,[lh,lw,lp+ld*(i+1)])
			pn=np.append(pn,[pn[0],pn[1],pn[2]-(pn[5]-pn[2])*(i+1)])
	return (func(x.reshape(-1,1).astype(float),pn[1::3].astype(float),pn[2::3].astype(float))*np.abs(pn[0::3]).astype(float)).sum(axis=1)

def pwsigma2alpha_wpm(p,M0):
	"""
	p - same for combatibility, but only first and last width param matter.
	sigma_i=(l-f)*(x_0-x_0_f)/(x_0_l-x_0_f)+f
	This function is the same a pwalpha2sigma_wpm
	"""
	l=p[1::3][-1]
	f=p[1::3][0]
	x_0=p[2::3]
	x_0_f=x_0[0]
	x_0_l=x_0[-1]
	sigmas=(l-f)*(x_0-x_0_f)/(x_0_l-x_0_f)+f
	sigmas=sigmas.reshape(-1)
	return np.array([p[0::3],sigmas,p[2::3]]).T.reshape(-1)

def pwalpha2sigma_wpm(p,M0):
	"""
	Identical to pwsigma2alpha_wpm, we are not doing any substitutions here, just for compatibility.
	"""
	return pwsigma2alpha_wpm(p,M0)

def pwsigma2alpha_logquad(p,M0):
	"""
	p - same for combatibility, but only first, last and middle width param matter (f,l,m)
	log(sigma)=P2(log(M)); sigma=exp(a*log(M)**2+b*log(M)+c)
	we need to know M0 - DNA length of the first peak.
	sigma_i=exp(a*log(M0+i)**2+b*log(M)+c)
	"""

	logsl=np.log(p[1::3][-1])
	logsf=np.log(p[1::3][0])
	mi=int(len(p[1::3])/2)
	pl=len(p[1::3])
	logsm=np.log(p[1::3][mi])

	logmf=np.log(M0)
	logml=np.log(M0-pl+1)
	logmm=np.log(M0-mi)

	a,b,c=np.linalg.solve([[logmf**2,logmf,1],[logmm**2,logmm,1],[logml**2,logml,1]],[logsf,logsm,logsl])

	M=M0-np.arange(pl)
	sigmas=np.exp(a*np.log(M)**2+b*np.log(M)+c)
	sigmas=sigmas.reshape(-1)
	return np.array([p[0::3],sigmas,p[2::3]]).T.reshape(-1)

def pwalpha2sigma_logquad(p,M0):
	"""
	it is the same since we do not do any substitutions
	"""
	return pwsigma2alpha_logquad(p,M0)

def pwsigma2alpha_safa(p,M0):
	"""
	p - same for combatibility, but only first and last width param matter (f,l)
	sigma=sigma_0+k*(x_{i+1}-x_{i-1})/2
	"""

	l=p[1::3][-1]
	f=p[1::3][0]
	x_0=p[2::3]

	diffs=x_0[2:]-x_0[:-2]
	diffs=np.insert(diffs,0,(x_0[1]-x_0[0])*2)
	diffs=np.append(diffs,(x_0[-1]-x_0[-2])*2)

	#f=sigma_0+diffs[0]*k
	#l=sigma_0+diffs[-1]*k
	k=(l-f)/(diffs[-1]-diffs[0])
	sigma_0=f-diffs[0]*k

	sigmas=sigma_0+k*diffs
	sigmas=sigmas.reshape(-1)
	return np.array([p[0::3],sigmas,p[2::3]]).T.reshape(-1)

def pwalpha2sigma_safa(p,M0):
	"""
	it is the same since we do not do any substitutions
	"""
	return pwsigma2alpha_safa(p,M0)

def residuals(p,y,x,M0,func,pwa2sigmafunc=lambda x,M0: x,Nauxvirtpeaks=0):
	err = y - sum_peakfuncs(x,p,M0,func,pwa2sigmafunc,Nauxvirtpeaks)
	#We still need penalty for deviation of peaks outside boundaries and
	#penalty if there is an unlikely situation, when peaks change positions.
	x_0=p[2::3]
	diffs=x_0[1:]-x_0[:-1]
	penalty=1+diffs[diffs<0].sum()*(-1)
	if np.min(x_0)<np.min(x):
		penalty+=10*(np.min(x)-np.min(x_0))
	if np.max(x_0)>np.max(x):
		penalty+=10*(np.max(x_0)-np.max(x))

	return err*penalty

def residuals_with_penalty(p,y,x,peaks,M0,func,pwa2sigmafunc=lambda x,M0: x, Nauxvirtpeaks=0):
	"""
	Here we penalize solutions where width of individual Gaussians/Lorentzians
	becomes too wide. This prevents overfitting of the data.
	Contstraint is that width cannot be greater that the distance between peaks *2.
	peaks - array of peak positions. One per band on the gel.
	p2sigmafunc - function for converting width parameters to width if such substitution was made.
	i.e. sigmaconvfunc=alpha_to_sigma
	"""
	err = y - sum_peakfuncs(x,p,M0,func,pwa2sigmafunc,Nauxvirtpeaks)

	sigmas=pwa2sigmafunc(p,M0).reshape(-1,3)[:,1]
	it=peaks[1:]-peaks[0:-1]
	sigmathresh=np.insert(it,0,it[0])*2
	diff=sigmathresh-sigmas
	penalty=1+diff[diff<0].sum()*(-1)

	x_0=p[2::3]
	diffs=x_0[1:]-x_0[:-1]
	penalty+=diffs[diffs<0].sum()*(-1)
	if np.min(x_0)<np.min(x):
		penalty+=10*(np.min(x)-np.min(x_0))
	if np.max(x_0)>np.max(x):
		penalty+=10*(np.max(x_0)-np.max(x))

	return err*penalty

###Function to estimate initial parameters for curve fitting

def getinit_func_param(peakvalues,peaks,func,first_sigma=1,last_sigma=10):
	"""
	Gets initial parameters for Gaussian/Lorentzian fitting.
	Requires that width increases from left to right.
	sigma=k*x_0+b
	last_sigma=k*peaks[-1]+b
	first_sigma=k*peaks[0]+b
	first_sigma and last_simga - initial estimates for width of first and last peaks.
	func - elementary Gaussian or Lorentzian function.
	Returns vector [H,width,position]*N
	"""
	s=[]
	x0=[]
	h=[]
	maxi=len(peaks)
	k=np.abs(float(np.abs(last_sigma)-np.abs(first_sigma))/float(peaks[-1]-peaks[0]))
	b=float(np.abs(last_sigma))-(float(k)*float(peaks[-1]))

	for i in range(maxi):
		x_0=float(peaks[i])
		sigma=k*peaks[i]+b
		s.append(sigma)
		if i==0:
			d=peaks[i+1]-peaks[i]
			hi=(peakvalues[i]/(func(0,sigma,0)+func(d,sigma,0)+func(2*d,sigma,0)+func(3*d,sigma,0)))
		elif i==(maxi-1):
			d=peaks[i]-peaks[i-1]
			hi=(peakvalues[i]/(func(0,sigma,0)+func(d,sigma,0)+func(2*d,sigma,0)+func(3*d,sigma,0)))
		else:
			d1=peaks[i+1]-peaks[i]
			d2=peaks[i]-peaks[i-1]
			hi=(peakvalues[i]/(func(0,sigma,0)+func(d1,sigma,0)+func(2*d1,sigma,0)+func(3*d1,sigma,0)+func(d2,sigma,0)+func(2*d2,sigma,0)+func(3*d2,sigma,0)))
		h.append(hi)
		x0.append(x_0)
	return np.array([h,s,x0]).T.reshape(-1)

###Main fitting routine

def fit_peaks(lane_profile_file,lane_config_file,DNAseq,lane_name,labeled_end='three_prime',\
	out_path='',peaktype='Gaussian',fitting_constraint='dSIGMA>=0',maxfev=10000,\
	Nauxvirtpeaks=0,\
	csvfileout=None,pngfileout=None,graphshow=False,\
	plotcontrib=False,ploterr=False,plotinitialguess=False, title=None,qthresh=100,\
	csvfilein=None):
	"""
	Wrapper for fit_peaks2footprint that loads data from files,
	and triggers automatic peak identification.
	Input
	-----
	lane_profile_file,lane_config_file - paths to files with lane data and configs
	DNAseq - DNA sequence in str of biopython format
	lane_name - name of lane to be analyzed as specified in lane_config_file
	labeled_end='three_prime' or 'five_prime' - the labeled end of the DNA
	out_path - path to directory where to output file with fitted data
	peaktype - 'Gaussian' or 'Lorentzian'
	----- Constrained fit parameters -----
	fitting_constraint equals one of the following:
		None - unconstrained fit.
		'dSIGMA>=0' - peak width (sigma) are constrained to increase monotonically from left to right. This dramatically improves stability of fitting.
		'SIGMA<2*dD' - Peak width (sigma) should not exceed twice the distance between the neighboring peaks (dD). Imlies automatically 'dSIGMA>=0'.
		'SIGMA=k*D+b' - peak width (sigma) will be linearly related to the position of the peak (D), which is in turn proportional to its mobility.
		'log(SIGMA)=P2(log(M))' - logarithm of peak width (sigma) will be related to number of base bairs in DNA (M) via a second degree polynomical.
		'SAFA' - SAFA type contraints, sigma=sigma_0+k*(D_{i+1}-D_{i-1})/2
	----------------------------------------
	Nauxvirtpeaks - number of auxillary virtual Gaussians/Lorenzians to add to the fitting model to the left and
		to the right of the fitted data range. This is needed if contribution from bands outside of the
		data range is expected to influence the profile values near the ends of the data range.
		If not 0, Nauxvirtpeaks will be added to the fitting model at each end of the data range.
		These peaks will be of the same width and height, there positions will be fixed so that the distance to the last peak
		in the data range will be a multiple of the distance between the two last peaks in the data range. 
	csvfileout - name of file to output data relative to output_path, if None the lane_name will we used.
	pngfileout - name of file to output image of the fitting results, relative to output_path
	graphshow - whether to show an interactive plot for exploration after the fitting is done
	ploterr - whether to plot graph of errors between data and fitted curve
	plotcontrib - plot individial Gaussian/Lorentzian functions fitted to each band
	plotinitialguess - whether to plot initial guess curve used to start optimization
	title - custom title for the plot
	qthresh - limit data range view to qthresh percentile of the values.
	csvfilein - readin fitting results produced earlier, fitting will not be performed - only plotting.
	"""
	lane=read_profile(lane_profile_file,lane_config_file,lane_name)
	lane=lane[~np.isnan(lane)].values #remove any missing values, usually found at the end if profiles were of different length and convert to NUMPY ndarray
	params=guess_missing_params(read_prof_configs(lane_config_file,lane_name),len(lane)) #Guess if certain params are NaN, or use "end" keyword for alignpos.
	peaks,laneSM=find_peaks(lane, params['peakthresh'], params['leftlim'], params['rightlim'],params['base'],params['segments'],params['min_dist_left'],params['min_dist_right'],params['addpeaks'],params['delpeaks'],params['interpolate'])

	title = "%s: Fitting experimental profile by modeling each band/peak as a %s."%(lane_name,peaktype)
	fit_peaks2footprint(lane,peaks,DNAseq,params=params,\
		csvfileout=os.path.join(out_path,csvfileout if csvfileout else lane_name+'_fitted_intensities.csv'),\
		peaktype=peaktype,name=lane_name,pngfileout=os.path.join(out_path,pngfileout if pngfileout else None),\
		labeled_end=labeled_end,Nauxvirtpeaks=Nauxvirtpeaks,\
		graphshow=graphshow,plotcontrib=plotcontrib,ploterr=ploterr,plotinitialguess=plotinitialguess,\
		title=title,qthresh=qthresh,fitting_constraint=fitting_constraint,maxfev=maxfev,csvfilein=os.path.join(output_path,csvfilein) if csvfilein else None)

def fit_peaks2footprint(lane,peaks,seq,params={},csvfileout='fit.csv',labeled_end='three_prime',peaktype='Gaussian',fitting_constraint=None,Nauxvirtpeaks=0,name='Profile',pngfileout=None,graphshow=True,plotcontrib=False,ploterr=False,plotinitialguess=False, title='',qthresh=100,maxfev=10000,csvfilein=None):
	"""
	This function fits a sum of Gaussians/Lorentzians functions to the hydroxyl-radical footptining profile
	and outputs a table of estimated cleavage intensities for every position in sequence.
	Input parameters:
	lane - 1D array of datavalues defining the experimental profile.
	peaks - an array of starting positions of peaks where Gaussian/Lorentzians will be fitted, data points are numbered starting from 0.
	seq - DNA sequence of the strand. string or biopython Seq object. 5'-3' always
	params - a dictionary of additional important parameters
		'leftlim' - left border of the region in the data that will be used for fitting.
		'rightlim' - right border of the region in the data that will be used for fitting.
		'seqpeak' - approximate position of the peak, corresponding to nucleotide whose position in DNA sequence will we provided.
		'seqpos' - position of this nucleotide in following format: "1A", or "45G"
	csvfileout - path to file where table with intensities will be written.
	labeled_end='three_prime' or 'five_prime' - the labeled end of the DNA
	peaktype - 'Gaussian' or 'Lorentzian'
	----- Constrained fit parameters -----
	fitting_constraint equals one of the following:
		None - unconstrained fit.
		'dSIGMA>=0' - peak width (sigma) are constrained to increase monotonically from left to right. This dramatically improves stability of fitting.
		'SIGMA<2*dD' - Peak width (sigma) should not exceed twice the distance between the neighboring peaks (dD). Imlies automatically 'dSIGMA>=0'.
		'SIGMA=k*D+b' - peak width (sigma) will be linearly related to the position of the peak (D), which is in turn proportional to its mobility.
		'log(SIGMA)=P2(log(M))' - logarithm of peak width (sigma) will be related to number of base bairs in DNA (M) via a second degree polynomical.
		'SAFA' - SAFA type contraints, sigma=sigma_0+k*(D_{i+1}-D_{i-1})/2
	----------------------------------------
	Nauxvirtpeaks - number of auxillary virtual Gaussians/Lorenzians to add to the fitting model to the left and
		to the right of the fitted data range. This is needed if contribution from bands outside of the
		data range is expected to influence the profile values near the ends of the data range.
		If not 0, Nauxvirtpeaks will be added to the fitting model at each end of the data range.
		These peaks will be of the same width and height, there positions will be fixed so that the distance to the last peak
		in the data range will be a multiple of the distance between the two last peaks in the data range. 
	pngfileout - path to file to output technical plot or None.
	graphshow - True if interactive plot desired.
	plotcontrib - plot individial Gaussian/Lorentzian functions fitted to each band
	ploterr - True if plot error plot in graph.
	plotinitialguess - True if initial guess plot desired
	title - title
	maxfev - maximum number of function evaluations, default 10000
	csvfilein - readin fitting results produced earlier, fitting will not be performed - only plotting.
	"""
	#Get experimental profile
	lanefull=lane

	#Assign parameters
	leftlim=params.get('leftlim',0)
	rightlim=params.get('rightlim',len(lanefull)-1)
	seqpeak=params.get('seqpeak',0)
	seqpos=params.get('seqpos','1%s'%seq[0])
	seqpos_num=params.get('seqpos_num',int(seqpos[0:-1]))

	#Let's chop data to specified limit range
	drange=np.arange(leftlim,rightlim+1)
	lane=np.array(lanefull[leftlim:rightlim+1])
	peaks=np.array(peaks)
	peaks=peaks[(peaks>=leftlim)&(peaks<=rightlim)]

	peakvalues=lane[np.searchsorted(drange,peaks)]

	#we need to know the position of known seq peak.
	known_pos = (np.abs(peaks-int(seqpeak))).argmin()
	print "Known peak is ",known_pos," ",peaks[known_pos]
	M0=(len(seq)-(seqpos_num-known_pos)) # size of DNA corresponding to leftmost peak (with index 0)

	#Let's define helper functions, so that Lorentz and Gauss will be the same:
	if(peaktype=='Lorentzian'):
		func=lorentzian
	elif(peaktype=='Gaussian'):
		func=gaussian
	else:
		print "Peaktype not supported!"
		raise(Exception)

	##### Now it's time to do fitting or read precalculated parameters#####
	if not csvfilein:
		#Define the points that will be used for fittig.
		#10 point between every peak
		dp=np.searchsorted(drange,peaks)
		dpn=dp
		for i,st,dd in zip(range(len(dp[0:-1])),dp[0:-1],(dp[1:]-dp[0:-1])):
			dpn=np.insert(dpn,i+1,[int(st+k*dd/10) for k in range(1,10) ])
		dpn=np.insert(dpn,0,np.where(drange<peaks[0])[0])
		dpn=np.insert(dpn,-1,np.where(drange>peaks[-1])[0])
		dp=np.unique(dpn)
		print len(dp)," datapoints will be used for fitting"

		if(fitting_constraint=='dSIGMA>=0'):
			pwa2sigmafunc=pwalpha2sigma
			pws2alphafunc=pwsigma2alpha
		elif(fitting_constraint=='SIGMA<2*dD'):
			pwa2sigmafunc=pwalpha2sigma
			pws2alphafunc=pwsigma2alpha
		elif(fitting_constraint=='SIGMA=k*D+b'):
			pwa2sigmafunc=pwsigma2alpha_wpm
			pws2alphafunc=pwalpha2sigma_wpm
		elif(fitting_constraint=='log(SIGMA)=P2(log(M))'):
			pwa2sigmafunc=pwsigma2alpha_logquad
			pws2alphafunc=pwalpha2sigma_logquad
		elif(fitting_constraint=='SAFA'):
			pwa2sigmafunc=pwsigma2alpha_safa
			pws2alphafunc=pwalpha2sigma_safa
		else:
			pwa2sigmafunc=lambda x,M0: x
			pws2alphafunc=lambda x,M0: x


		print "=====Starting empirical optimization of initial parameters ====="
		
		def sum_constr_residuals_for_init_fit(p):
			parinit=getinit_func_param(peakvalues,peaks,func,p[0],p[1])
			sresid=(residuals(parinit,lane[dp],drange[dp],M0,func,Nauxvirtpeaks=Nauxvirtpeaks)**2).sum()
			penalty=(p[0]-p[1]+1) if p[1]<p[0] else 1
			return sresid*penalty

		#Fit a line to dependency of peak width with respect to position (should be approximately linear)
		resS=np.polyfit(peaks[0:-1], peaks[1:]-peaks[0:-1],1)
		#Make a conservative estimate of width of first and last peak below the actual values
		# to initiate search by further minimization.
		fs=(peaks[0]*resS[0]+resS[1])/5
		ls=(peaks[-2]*resS[0]+resS[1])/5

		#Get initial parameters before optimization
		pinit_b_o=getinit_func_param(peakvalues,peaks,func,first_sigma=fs,last_sigma=ls)
		#bounds set to distance between peaks/2 - you can not see peaks if distance is equal sigma. even 0.8 sigma.
		#But to start optimiztion we need sharp peaks
		res=minimize(sum_constr_residuals_for_init_fit,[fs,ls],bounds=[(fs,fs*5/2),(ls,ls*5/2)],options={'maxiter':5000, 'disp': True})
		pinit=getinit_func_param(peakvalues,peaks,func,first_sigma=res['x'][0],last_sigma=res['x'][1])

		print "RMSD of initial guess=",((residuals(pinit_b_o,lane,drange,M0,func)**2).sum()/len(drange))**0.5
		print "RMSD of optimized initial guess=",((residuals(pinit,lane,drange,M0,func)**2).sum()/len(drange))**0.5
		print "Params before optimization:", fs,ls
		print "Params after optimization:", res['x']

		print "======================"

		print "=====Starting least square fit ====="


		pinit_subs=pws2alphafunc(pinit,M0)

		if(fitting_constraint=='SIGMA<2*dD'):
			pbest = leastsq(residuals_with_penalty,pinit_subs,args=(lane[dp],drange[dp],peaks,M0,func,pwa2sigmafunc,Nauxvirtpeaks),full_output=1,maxfev=maxfev,xtol=1.49012e-08,ftol=1.49012e-08,factor=1)
		else:
			pbest = leastsq(residuals,pinit_subs,args=(lane[dp],drange[dp],M0,func,pwa2sigmafunc,Nauxvirtpeaks),full_output=1,maxfev=maxfev,xtol=1.49012e-08,ftol=1.49012e-08,factor=1)

		pres = pbest[0] #This is the vector of optimized parameters
		print "====Finished least square fit"
		print "Final RMSD LS=",((pbest[2]['fvec']**2).sum()/len(dp))**0.5
		best_parameters=pwa2sigmafunc(pres,M0)

	else: #We have precalculated parameters in csvfilein	
		df_fp_int=pd.read_csv(csvfilein,comment='#')
		if('gamma_2' in df_fp_int.columns):
			peaktypeF='Lorentzian'
			best_parameters = np.array([[df_fp_int.iloc[i]['Intensity'],df_fp_int.iloc[i]['gamma_2'],df_fp_int.iloc[i]['x_0']] for i in range(len(df_fp_int))]).reshape(1,-1)[0]

		if('sigma' in df_fp_int.columns):
			peaktypeF='Gaussian'
			best_parameters = np.array([[df_fp_int.iloc[i]['Intensity'],df_fp_int.iloc[i]['sigma'],df_fp_int.iloc[i]['x_0']] for i in range(len(df_fp_int))]).reshape(1,-1)[0]
		if peaktype!=peaktypeF:
			raise Exception("peaktype in %s does not match the one specified in paramters of the function"%csvfilein)

	finrmsd=((residuals(best_parameters,lane,drange,M0,func,Nauxvirtpeaks=Nauxvirtpeaks)**2).sum()/len(drange))**0.5
	print "Final RMSD=",finrmsd
	reldev=np.abs(residuals(best_parameters,lane,drange,M0,func,Nauxvirtpeaks=Nauxvirtpeaks)/lane*100)
	relrmsd=100*(((residuals(best_parameters,lane,drange,M0,func,Nauxvirtpeaks=Nauxvirtpeaks)/lane)**2).sum()/len(drange))**0.5
	print "Final relative RMSD=",relrmsd,"%"

	#RMSD of peak area
	fit = sum_peakfuncs(drange,best_parameters,M0,func,Nauxvirtpeaks=Nauxvirtpeaks)
	area_dev=[]
	for i in range(len(peaks)):
		lb=(peaks[i-1]+peaks[i])/2 if i>0 else leftlim
		rb=peaks[i]+peaks[i+1]/2 if i<(len(peaks)-1) else rightlim
		farea=fit[(drange>=lb)&(drange<rb)].sum()
		rarea=lane[(drange>=lb)&(drange<rb)].sum()
		area_dev.append((farea-rarea)/rarea)
	area_dev=np.array(area_dev)
	area_rrmsd=((area_dev**2).sum()/len(area_dev))**0.5
	area_mdev=max(np.abs(area_dev))
	print "Relative area RMSD=%f %% Max rel peak area error=%f %%"%(100*area_rrmsd,100*area_mdev)

	fit = sum_peakfuncs(drange,best_parameters,M0,func,Nauxvirtpeaks=Nauxvirtpeaks)
	if not csvfileout:
		initguess = sum_peakfuncs(drange,pinit_b_o,M0,func,Nauxvirtpeaks=Nauxvirtpeaks)
		initguess_opt = sum_peakfuncs(drange,pinit,M0,func,Nauxvirtpeaks=Nauxvirtpeaks)
	elif plotinitialguess:
		raise Exception("plotinitialguess=True cannot be combined with reading precalculated results from file %s"%csvfilein)

	########Extracting intensities for sequence and output
	#############
	#Footprint intensities are here
	h=np.abs(best_parameters[0::3])

	#The caluculation of peak positions below and outout

	#Test that order of x_0 follows the order of peaks
	if np.array_equal(np.sort(best_parameters[2::3]),best_parameters[2::3]) and not csvfilein:
		if labeled_end=='three_prime':
			print "Outputing data to ",csvfileout
			print("Outputting for 3 prime labeled DNA")
			if(peaktype=='Lorentzian'):
				df_fp_int=pd.DataFrame({'Intensity':h,
				'Site':[ "%d%s"%((i-known_pos+seqpos_num),seq[(i-known_pos+seqpos_num)-1]) for i in range(len(h))],
				'gamma_2':best_parameters[1::3],
				'x_0':best_parameters[2::3],
				'area_dev_percent':area_dev*100})
			if(peaktype=='Gaussian'):
				df_fp_int=pd.DataFrame({'Intensity':h,
				'Site':[ "%d%s"%((i-known_pos+seqpos_num),seq[(i-known_pos+seqpos_num)-1]) for i in range(len(h))],
				'sigma':best_parameters[1::3],
				'x_0':best_parameters[2::3],
				'area_dev_percent':area_dev*100})

		elif labeled_end=='five_prime':
			print "Outputing data to ",csvfileout
			print("Outputting for 5 prime labeled DNA")
			if(peaktype=='Lorentzian'):
				df_fp_int=pd.DataFrame({'Intensity':h[::-1],
				'Site':[ "%d%s"%((known_pos-i+seqpos_num),seq[(known_pos-i+seqpos_num)-1]) for i in (range(len(h))[::-1])],
				'gamma_2':best_parameters[1::3][::-1],
				'x_0':best_parameters[2::3][::-1],
				'area_dev_percent':(area_dev*100)[::-1]})
			if(peaktype=='Gaussian'):
				df_fp_int=pd.DataFrame({'Intensity':h[::-1],
				'Site':[ "%d%s"%((known_pos-i+seqpos_num),seq[(known_pos-i+seqpos_num)-1]) for i in (range(len(h))[::-1])],
				'sigma':best_parameters[1::3][::-1],
				'x_0':best_parameters[2::3][::-1],
				'area_dev_percent':(area_dev*100)[::-1]})

		else:
			raise Exception("labeled_end should be 'three_prime' or 'five_prive'")


		with open(csvfileout,'w') as file:
			file.write('#Model fitted to experimental profiles on %s\n#Profile name: %s\n#Fitted model type: %s\n#Fitting constraints:%s\n'%(time.strftime("%d, %b %Y, %H:%M:%S"),name,peaktype,fitting_constraint))
			file.write('#Fitting errors:\n#RMSD exp vs fit=%f\n#Relative RMSD exp vs fit=%f %%\n#Relative area RMSD=%f %%\n#Max relative peak area error=%f %%\n'%(finrmsd,relrmsd,100*area_rrmsd,100*area_mdev))
			file.write('#Data column description\n#Intensity - total deconvoluted intensity of each band (coefficient before Gaussian or Lorentzian function,\
area under the individual Gaussian/Lorentzian curve)\n#Site - DNA site corresponding to the band (cleaved site)\n\
#area_dev_percent - % difference in area under each peak between fitted curve and original data (the segment used to calculate are for each peak starts and ends at the midpoints between the neighboring peaks)\n\
#sigma or gamma_2 - the width parameter of the functions modeling each band, sgima if Gaussain ,gamma/2 if Lorentzian\n\
#x_0 - position of the Gaussian/Lorentzian modeling each band on the original profile\n')

		df_fp_int.to_csv(csvfileout,mode='a') #Here we output everything!!!
	else:
		if not csvfilein:
			print "Order of the peaks in fit does not follow order in initial data!!!\n Something bad happend, data will not be output!!!"
			print map(str,np.sort(best_parameters[2::3])),'\n',map(str,best_parameters[2::3]),'\n',map(str,peaks)
		else:
			pass
	########Technical plots##############
	#####################################
	fig=plt.figure(figsize=(12,4))

	gs = gridspec.GridSpec(2, 1,height_ratios=[18,2])

	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1])
	ax1.plot(drange,lane, linewidth=1, label=name+", experimental")
	ax1.scatter(peaks, peakvalues, marker='x', color='g', s=40,linewidth=2,label="Band peak locations")
	ax1.scatter(peaks[known_pos], peakvalues[known_pos]*1.05, marker='o', color='orange', s=40,linewidth=2,label=seqpos)

	if(plotinitialguess):
		ax1.plot(drange,initguess,linewidth=1,color='cyan',label='Initial guess')
		ax1.plot(drange,initguess_opt,linewidth=1,color='brown',label='Initial optimized guess')
	ax1.plot(drange,fit,linewidth=1,color='magenta',label='Fitted model curve')
	if(ploterr):
		ax3 = ax1.twinx()
		ax3.plot(drange,reldev,linewidth=1,color='orange',label='Relative deviation, %')
		ax3.set_ylabel('Relative deviation,%')
		ax3.set_ylim((0,100))
		ax3.set_yticks(np.arange(0, 100,10))
		ax3.legend(loc='right', shadow=False)

	if(plotcontrib):
		labflag=True
		for h, g_2, x_0 in zip(*(iter(best_parameters),) * 3):
			if peaktype=='Lorentzian':
				y=lorentzian(drange,g_2,x_0)*h
			if peaktype=='Gaussian':
				y=gaussian(drange,g_2,x_0)*h
			if labflag:
				ax1.plot(drange,y,linewidth=1,color='green',label="Individual "+peaktype+"s")
				labflag=False
			else:
				ax1.plot(drange,y,linewidth=1,color='green',label=None)


	ax1.annotate('RMSD=%.2f\nRelative RMSD=%.2f %% \nRelative peak area RMSD=%.2f %%  \nMax rel peak area err=%.2f %%'%(finrmsd,relrmsd,100*area_rrmsd,100*area_mdev), xy=(0.5, 0.7), xycoords='axes fraction')

	ax1.set_ylabel("Intensity")
	ax1.set_title(title)
	ax1.set_ylim((0,np.percentile(lane,qthresh)))
	ax1.set_xlim((leftlim-1,rightlim+1))

	ax1.legend(loc="upper right")


	fakegelimg=np.tile(lanefull,(100,1))
	ax2.imshow(fakegelimg,cmap='Greys',vmin=0,vmax=np.percentile(lane,qthresh),interpolation='nearest', aspect='auto')
	# ax2.set_xlabel('Coordinate')
	ax2.set_xlim((leftlim-1,rightlim+1))

	ax2.set_yticks([])
	ax2.set_xticks([])

	try:
		ax1.axvline(x=leftlim,color='red')
		ax1.axvline(x=rightlim,color='red')
	except:
		pass

	fig.tight_layout()

	if(pngfileout):
		plt.savefig(os.path.splitext(pngfileout)[0]+'.png',dpi=(300))
		# os.system('convert %s -resize 25%% %s'%(os.path.splitext(pngfileout)[0]+'.png',os.path.splitext(pngfileout)[0]+'_small.png'))

	if(graphshow):
		plt.show()
	plt.close()



##Functions to plot profiles on sequence.
#######################

def plot_prof_on_seq(data_file,DNAseq,prof_names=None,\
			graphshow=False,pngfileout=None,title='',prof_columns="Intensity",seq_column="Site",rescale=None,zero_at=0,ylab=None,colors=None,\
			colorb={'A':'#0b0','T':'#0b0','G':'#00b','C':'#00b'},colorf={'A':'#fafafa','T':'#fafafa','G':'#fafafa','C':'#fafafa'},\
			plot_options={'linewidth':1.0,'markersize':8.0,'figsize':(12,3),'fontsize':None,'legendloc':'upper right'}):
	"""
	Function makes a plot of some values along the DNA sequence.
	data_file - csv file with a dataframe.
	prof_columns - list of column names to plot.
	rescale - None - no normalization,
		'every' - every profile separtely devide by its max value, (Y/Ymax),
		'every_min_max' - (Y-Ymin)/(Ymax-Ymin) for every profile,
		'together' - divide all the profiles by highest among max values,
		'together_min_max' - (Y-Ymin)/(Ymax-Ymin), where min and max are global for all profiles,
		'fit' - a linear fit between profiles without intersect will be made,
			and then they will be rescaled by the highest max value (as in 'together' mode).
		'fit_min_max' - a normal linear fit between profiles (with intersect) will be made
			and then they will be rescaled to [0,1] segement from global min to global max (same as together_min_max). 

	colors - list of colors for plottin several datasets.
	colorb - colors for DNA bases
	"""
	if plot_options['fontsize']:
		rcParams.update({'font.size': plot_options['fontsize']})

	#Let's make an image of a sequence
	temp = tempfile.TemporaryFile()
	im=Image.new('RGB',(24*len(DNAseq),40))
	aspect=(24.0*len(DNAseq))/40.0
	draw = ImageDraw.Draw(im)
	
	#Open font from package in a tricky way, independent of package installation mode
	temp2 = tempfile.TemporaryFile()
	fontfile = pkgutil.get_data('hydroid', 'pkgdata/cnrb.otf')
	temp2.write(fontfile)
	temp2.seek(0)
	font = ImageFont.truetype(temp2, 40)
	####

	for i,l in zip(range(len(DNAseq)),DNAseq):
		draw.rectangle(((24*i,0),(24*(i+1),40)),fill=colorb[l],outline=colorb[l])
		draw.text((24*i+1, 0), l,fill=colorf[l],font=font)
	im.save(temp, "PNG")
	# im.save('temp.png', "PNG")
	temp.seek(0)

	prof_columns=[prof_columns] if isinstance(prof_columns,basestring) else list(prof_columns)
	prof_names=([prof_names] if isinstance(prof_names,basestring) else list(prof_names)) if prof_names else None
	colors=list(colors) if colors else None


	data=pd.read_csv(data_file,comment='#')
	# print data
	yvalues=data[prof_columns].values

	if rescale=='together':
		yvalues=yvalues/np.nanmax(yvalues)
	elif rescale=='together_min_max':
		yvalues=yvalues-np.nanmin(yvalues)
		yvalues=yvalues/np.nanmax(yvalues)
	elif rescale=='every':
		yvalues=yvalues/np.nanmax(yvalues,axis=0)
	elif rescale=='every_min_max':
		yvalues=yvalues-np.nanmin(yvalues,axis=0)
		yvalues=yvalues/np.nanmax(yvalues,axis=0)
	elif rescale=='fit':
		coefs=[]
		for c,i in zip(prof_columns,range(len(prof_columns))):
			fitres = sm.OLS(yvalues[:,i], yvalues[:,0],missing='drop').fit()
			coefs.append(fitres.params[0])
		yvalues=yvalues/np.array(coefs)
		yvalues=yvalues/np.nanmax(yvalues)
	elif rescale=='fit_min_max':
		coefs=[]
		const=[]
		for c,i in zip(prof_columns,range(len(prof_columns))):
			fitres = sm.OLS(yvalues[:,i],sm.add_constant(yvalues[:,0], prepend=False) ,missing='drop').fit()
			coefs.append(fitres.params[0])
			const.append(fitres.params[1])
		yvalues=(yvalues-np.array(const))/np.array(coefs)
		yvalues=yvalues-np.nanmin(yvalues)
		yvalues=yvalues/np.nanmax(yvalues)

	seqimgheight=12/3/aspect*np.nanmax(yvalues)*1.05
	vlinestart=seqimgheight/(np.nanmax(yvalues)*1.05+seqimgheight)

	fig=plt.figure(figsize=plot_options['figsize'])
	ax1 = plt.subplot()

	ind=np.array(map(lambda x:int(x[0:-1]),data[seq_column].values))-zero_at

	ax1.yaxis.grid(linewidth=0.5)
	#Let's draw minor grid ticks
	ticklabs=[]
	for ls in range(0-zero_at,len(DNAseq)-zero_at):
		if ls%5==0:
			ax1.axvline(ls,linewidth=1 if ls%10==0 else 0.5,color='black',ymin=vlinestart)
		if ls%10==0:
			ticklabs.append(ls)

	for c,i in zip(prof_columns,range(len(prof_columns))):
		ax1.plot(ind, yvalues[:,i], '.b-',color= colors[i] if colors else None,label=prof_names[i] if prof_names else c,linewidth=plot_options['linewidth'],markersize=plot_options['markersize'])
	
	ax1.set_xticks(ticklabs)
	ax1.set_ylabel((ylab if ylab else 'Value')+(", Normalized (%s)"%rescale if rescale else ""))
	ax1.set_title(title)
	ax1.set_ylim((-seqimgheight,np.nanmax(yvalues)*1.05))
	ax1.set_xlim((0.5-zero_at,len(DNAseq)+0.5-zero_at))
	ax1.set_yticks(ax1.get_yticks()[np.array(ax1.get_yticks())>=0])
	# legend=ax1.legend(loc="upper right",shadow=False)
	legend=ax1.legend(loc=plot_options['legendloc'],shadow=False)
	# legend=ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	legend.get_frame().set_facecolor('#FFFFFF')
	legend.get_frame().set_alpha(1)

	seqimg = Image.open(temp)
	# im = OffsetImage(seqimg, zoom=0.26, resample=True)
	# ab = AnnotationBbox(im, (1,1), xycoords='data',\
	# xybox=(0.5,0),box_alignment=(0., 0.),frameon=False)
	# ax1.add_artist(ab)
	# ax1.imshow(seqimg,interpolation='nearest', aspect=aspect,extent=(0.5,len(DNAseq)+0.5,0.0,0.2))
	# ax1.imshow(seqimg,interpolation='nearest', aspect=2000)
	ax1.imshow(seqimg,interpolation='bilinear',aspect='auto',extent=(0.5-zero_at,len(DNAseq)+0.5-zero_at,-seqimgheight,0.0))
	


	fig.tight_layout()

	if(pngfileout):
		plt.savefig(os.path.splitext(pngfileout)[0]+'.png',dpi=(300))
	if(graphshow):
		plt.show()
	plt.close(fig)
	temp.close()
	temp2.close()


### Functions for plotting theoretical gels and estimated data
def ogston_pos(M,options={'Mu_0':3193.01,'D_0':350.64,'B':0.00,'C':0.14}):
    """
    This function implements calculations of DNA mobility (traveled distance) 
    in gel electrophoresis based on Ogston mobility model.
    ln(Mu) ~ -M, where Mu is mobility and M is molecular weight.
    
    This function returns positions of peaks on the gel lane based on their 
    molecular weight, expressed as a number of nucleotides.
    M - a numpy array with weights (by default molecular weight is measured in 
    the number of nucleotides)
    The formula to calculate the mobility is as follows
        Mu = Mu_0*exp(-(B+C*M^0.5)^2)-D_0
    options - dict with parameters.
        Default paramerters {'Mu_0':3193.01,'D_0':350.64,'B':0.00,'C':0.14}
        are derived from a sample experimental dataset. 
               	
    """
    return options['Mu_0']*np.exp(-(options['B']+options['C']*M**0.5)**2)-options['D_0']

def make_lane(positions,CleavageInt,V_sigmas,H_sigma,width,ident=None):
    """
    This is a helper function that simulates a gel lane image in the form of 2D
    numpy array based on band positions, widths and intensities (Gaussian approximation is used)
    positions - peak positions in pixel units
    CleavageInt - peak intensities 
    V_sigmas - sigmas for vertical core of gaussian filter for each peak (as they tend to grow)
    !positions,CleavageInt,V_sigmas  must have same size!
    H_sigma - sigma for horisontal core of gaussian filter 
    width - width of the resulting lane
    ident - identation from the last peak in pixels
    """
    if ident==None:
        ident=V_sigmas.max()*2
    H_sigma  =int(H_sigma)
    ident=int(ident)
    lane=np.zeros((int(positions.max()+ident*2),int(width)))
    for pos,intencity,sigma in zip(positions,CleavageInt,V_sigmas):
        temp=np.zeros((int(positions.max()+ident*2),int(width)))
        temp[int(pos),H_sigma:-H_sigma]=intencity #np.repeat(intencity,width)
        lane+=gfilt(temp,(sigma,H_sigma))

    return lane


def simulate_gel(intensities, offset, model='ogston',
				graphshow=False,pngfileout=None,title='',
                options={'Mu_0':3193.01,'D_0':350.64,'B':0.00,'C':0.14},
                Dc=0.00581843,
                v_blur={'A':-0.05, 'B':0.3, 'C':3.81},
                H_sigma=5,width=55):
        '''
        This function simulates the gel lane profile and image of DNA gel
        electrophoresis experiments based on known band intensities and DNA
        molecular weight.
        intensities - list of gel band intensities (DNA cleavage frequencies)
                ordered by increasing distance from the labeled end.
        offset  - distance from the DNA labeled end to the first datapoint
                in the list of intensities
        model   - model of DNA mobility in gel electrophoresis.
         'linear' - peak run distance depends from nucleotide number as B*(Nmax-Ni) + C
         'ogston' - peak run distance depends from nucleotide number as 
                Mu_0*exp(-(B+C*M^0.5)^2)-D_0
        options - dict for model parameters 
                example options={'Mu_0':3193.01,'D_0':350.64,'B':0.00,'C':0.14}
        Dc - dispersion coeffiscient constant, which relates movility and band 
            D = Dc * Mu  
        H_sigma - horisontal blur strength in pixels
        width   - lane width
        '''
        lengths_in_bases=np.arange(offset,offset+intensities.size)
        if model == 'ogston':
            runs=ogston_pos(lengths_in_bases,options)
        elif model == 'linear':
            runs = options['C'] + options['B']*(lengths_in_bases.max()- lengths_in_bases)
        
        V_sigmas=Dc*(runs+options['D_0'])
               
        basic_lane=make_lane(runs[runs>0],
                     intensities[runs>0],
                     V_sigmas[runs>0],
                     H_sigma=H_sigma,width=width)

        plt.rcParams["figure.figsize"] =[12,3]
        figGel, ax = plt.subplots(nrows=2)
        
        ax[0].imshow(basic_lane.T,cmap='Greys')
        
        ax[1].plot(basic_lane[:,50],label='Simulated lane')
        ax[1].set_ylabel('Intensity')
        ax[1].set_xlabel('Travel distance, px')
        ax[1].set_xlim(ax[0].get_xlim())
        
        margins = dict(hspace=0.2, wspace=0.2, top=1, bottom=0.17, left=0.06, right=0.94)
        figGel.subplots_adjust(**margins) 
        ax[0].set_title('%s simulated gel image' % title)
        ax[1].set_title('Lane intensity profile')       
        if pngfileout != None:
	        plt.savefig(pngfileout,dpi=300)
        if (graphshow):
            plt.show()
        plt.close(figGel)


###Functions under development
def mob_plt(filein,intcsv,seq,paramfile,graphout=None,graphshow=True,ploterr=True,title='',qthresh=100,brfcut=-1,D0=False):
	params=read_paramfile(paramfile,filein)
	mobility_plot(filein,intcsv,seq,params,graphout=graphout,graphshow=graphshow,ploterr=ploterr,title=title,qthresh=qthresh,brfcut=brfcut,D0=D0)

def mobility_plot(filein,intcsv,seq,params={},graphout=None,graphshow=True,ploterr=True,title='',qthresh=100,brfcut=-1,D0=False):
	"""
	Plot mobility and dispersion coefficients information as well as fits.
	"""

	lanefull=read_datafile(filein)

	#Assign parameters
	leftlim=params.get('leftlim',0)
	rightlim=params.get('rightlim',len(lanefull)-1)
	seqpeak=params.get('seqpeak',0)
	seqpos=params.get('seqpos','1X')
	seqpos_num=params.get('seqpos_num',int(seqpos[0:-1]))
	drange=np.arange(leftlim,rightlim+1)
	lane=np.array(lanefull[leftlim:rightlim+1])
	try:
		peaks,laneS=find_peaks(lanefull, params['peakthresh'], params['leftlim'], params['rightlim'],params['base'],params['segments'],params['min_dist_left'],params['min_dist_right'],params['addpeaks'],params['delpeaks'],params['interpolate'])
		peaks=peaks[(peaks>=leftlim)&(peaks<=rightlim)]
	except:
		pass
	df_fp_int=pd.read_csv(intcsv)
	if('gamma_2' in df_fp_int.columns):
		peaktype='Lorentzian'
		best_parameters = np.array([[df_fp_int.iloc[i]['Intensity'],df_fp_int.iloc[i]['gamma_2'],df_fp_int.iloc[i]['x_0']] for i in range(len(df_fp_int))]).reshape(1,-1)[0]

	if('sigma' in df_fp_int.columns):
		peaktype='Gaussian'
		best_parameters = np.array([[df_fp_int.iloc[i]['Intensity'],df_fp_int.iloc[i]['sigma'],df_fp_int.iloc[i]['x_0']] for i in range(len(df_fp_int))]).reshape(1,-1)[0]

	known_pos = (np.abs(best_parameters[2::3]-int(seqpeak))).argmin()
	known_pos_x=int(best_parameters[known_pos*3+2])

	distance=best_parameters[2::3]
	bp_size=[(len(seq)-(seqpos_num+i-known_pos))  for i in range(len(best_parameters[2::3]))]

	#Assume that distance is not measured from actual zero, so mobility is (distance + something)/Field.
	#Ogston model we have Log (Dist+D_0) = Log(Mu_0)-(B+C*R_g)**2, where R_g=M**0.5
	#Repatation BRF model we have Dist+D_0=Mu_0*(1/(3*M)+B)

	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	plt.rcParams['ps.useafm'] = True
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	plt.rcParams['pdf.fonttype'] = 42

	#Fitting Ogston model
	fig=plt.figure(1,figsize=(18,6))
	ax1 = fig.add_subplot(1, 1, 1)

	ax1.scatter(bp_size,distance, linewidth=1, label=filein.split('/')[-1])
	ax1.set_title(title+": Peak positions vs Base pair")
	#Fitting Ogston model
	def resid_cur(p,func,mbp,dist):
		return dist-func(p,mbp)
	def func_ogston(p,mbp):
		x=np.array(mbp).astype(float)
		Mu_0=np.abs(p[0])
		B=np.abs(p[1])
		C=np.abs(p[2])
		if D0:
			D_0=D0
		else:
			D_0=np.abs(p[3])
		R_g=np.array(mbp)**0.5
		return Mu_0*np.exp(-(B+C*R_g)**2)-D_0
		# return Mu_0*np.exp(-B*np.array(mbp))-D_0
	Mu_00=distance[-1]*1.5
	if D0:
		D_00=D0
	else:
		D_00=distance[-1]*0.5
	B0=0
	C0=(-np.log((distance[0]-distance[-1])/Mu_00+1)/bp_size[0])**0.5

	print "Starting values"
	print [Mu_00,B0,C0,D_00]
	pbest = leastsq(resid_cur,[Mu_00,B0,C0,D_00],args=(func_ogston,bp_size,distance),full_output=1,maxfev=10000,xtol=1.49012e-08,ftol=1.49012e-08)
	print "Finished fitting Ogston model"
	print pbest[3]
	pres=pbest[0]
	print pres
	finrmsd=((pbest[2]['fvec']**2).sum()/len(bp_size))**0.5
	print "Final Ogston model RMSD LS=",finrmsd
	ax1.plot(bp_size,func_ogston(pres,bp_size), linewidth=1,label='Ogston model fit',color='g')
	# ax1.plot(bp_size,func_ogston([Mu_00,B0,C0,D_00],bp_size), linewidth=1,label='Ogston model initial',color='r')
	ax1.annotate('Ogston model: Mu_0*exp(-(B+C*M^0.5)^2)-D_0\nOgston RMSD=%.2f\nMu_0=%.2f\nD_0=%.2f\nB=%.2f\nC=%.2f'%(finrmsd,pres[0],pres[-1],pres[1],pres[2]), xy=(0.5, 0.7), xycoords='axes fraction')


	#Let's fit BRF model
	#Dist+D_0=Mu_0*(1/(3*M)+B)
	def func_brf(p,mbp):
		x=np.array(mbp).astype(float)
		Mu_0=np.abs(p[0])
		B=np.abs(p[1])
		if D0:
			D_0=D0
		else:
			D_0=np.abs(p[2])
		return Mu_0*(1.0/(3.0*np.array(mbp))+B)-D_0

	Mu_00=distance[-1]*1.5
	if D0:
		D_00=D0
	else:
		D_00=distance[-1]*0.5
	B0=1.0


	print "Starting values"
	print [Mu_00,B0,D_00]
	pbest = leastsq(resid_cur,[Mu_00,B0,D_00],args=(func_brf,bp_size[0:brfcut],distance[0:brfcut]),full_output=1,maxfev=10000,xtol=1.49012e-08,ftol=1.49012e-08)
	print "Finished fitting BRF model"
	print pbest[3]
	pres=pbest[0]
	print pres
	finrmsd=((pbest[2]['fvec']**2).sum()/len(bp_size))**0.5
	print "Final BRF model RMSD LS=",finrmsd
	ax1.plot(bp_size[0:brfcut],func_brf(pres,bp_size[0:brfcut]), linewidth=1,label='BRF model fit',color='r')
	ax1.annotate('BRF model: Mu_0*(1/(3*M)+B)-D_0\nBRF RMSD=%.2f\nMu_0=%.2f\nB=%.2f\nD_0=%.2f'%(finrmsd,pres[0],pres[1],pres[-1]), xy=(0.5, 0.5), xycoords='axes fraction')
	# ax1.plot(bp_size,func_brf([Mu_00,B0,D_00],bp_size), linewidth=1,label='BRF model initial',color='b')

	#Let's fit Empirical model
	#Dist=(Mu_0*B)/(M+B)-D_0
	def func_emp(p,mbp):
		x=np.array(mbp).astype(float)
		Mu_0=np.abs(p[0])
		B=np.abs(p[1])
		if D0:
			D_0=D0
		else:
			D_0=np.abs(p[2])
		return Mu_0*B/(np.array(mbp)+B)-D_0

	Mu_00=distance[-1]*1.5
	if D0:
		D_00=D0
	else:
		D_00=distance[-1]*0.5
	B0=1.0


	print "Starting values"
	print [Mu_00,B0,D_00]
	pbest = leastsq(resid_cur,[Mu_00,B0,D_00],args=(func_emp,bp_size,distance),full_output=1,maxfev=10000,xtol=1.49012e-08,ftol=1.49012e-08)
	print "Finished fitting EMP model"
	print pbest[3]
	pres=pbest[0]
	print pres
	finrmsd=((pbest[2]['fvec']**2).sum()/len(bp_size))**0.5
	print "Final EMP model RMSD LS=",finrmsd
	ax1.plot(bp_size,func_emp(pres,bp_size), linewidth=1,label='EMP model fit',color='c')
	ax1.annotate('EMP model: Mu_0*B(M+B)-D_0\nEMP RMSD=%.2f\nMu_0=%.2f\nB=%.2f\nD_0=%.2f'%(finrmsd,pres[0],pres[1],pres[-1]), xy=(0.8, 0.5), xycoords='axes fraction')
	# ax1.plot(bp_size,func_brf([Mu_00,B0,D_00],bp_size), linewidth=1,label='BRF model initial',color='b')
	
	
	ax1.legend(loc="upper right")
	ax1.set_ylabel("Peak coordinate (distance)")
	ax1.set_xlabel("DNA size (M), bp")

	if(graphout):
		fig.savefig(graphout+'_mob'+'.png',dpi=(300))
		os.system('convert %s -resize 25%% %s'%(graphout+'_mob'+'.png',graphout+'_mob'+'_small.png'))

	#Now let's look at peak width, sigma, or dispersion.
	sigma=best_parameters[1::3]
	fig2=plt.figure(2,figsize=(18,6))
	ax2 = fig2.add_subplot(1, 1, 1)
	ax2.scatter(bp_size,sigma,linewidth=1,color='g',label=filein.split('/')[-1])
	ax2.set_title(title+": Peak width vs Base pair")

	def func_expquad(p,x):
		x=np.array(x).astype(float)
		l=np.log(x)
		a=p[0]
		b=p[1]
		c=p[2]
		return np.exp(a*l*l+b*l+c)

	pbest = leastsq(resid_cur,[-1.0,-1.0,0],args=(func_expquad,bp_size,sigma),full_output=1,maxfev=10000,xtol=1.49012e-08,ftol=1.49012e-08)
	print "Finished fitting expquad model"
	print pbest[3]
	pres=pbest[0]
	print pres
	finrmsd=((pbest[2]['fvec']**2).sum()/len(bp_size))**0.5
	print "Final expquad model RMSD LS=",finrmsd
	ax2.plot(bp_size,func_expquad(pres,bp_size), linewidth=1,label='Expquad fit',color='r')
	# ax1.plot(bp_size,func_ogston([Mu_00,B0,C0,D_00],bp_size), linewidth=1,label='Ogston model initial',color='r')
	ax2.annotate('Exp quad fit: exp(a*l^2+b*l+c)\nRMSD=%.3f\na=%.2f\nb=%.2f\nc=%.2f'%(finrmsd,pres[0],pres[1],pres[-1]), xy=(0.7, 0.7), xycoords='axes fraction')

	ax2.set_ylabel("Width (sigma fitted)")
	ax2.set_xlabel("Base numer (M)")
	ax2.legend(loc='upper right', shadow=False)

	if(graphout):
		fig2.savefig(graphout+'_width'+'.png',dpi=(300))
		os.system('convert %s -resize 25%% %s'%(graphout+'_width'+'.png',graphout+'_width'+'_small.png'))


	fig3=plt.figure(3,figsize=(18,6))
	ax3 = fig3.add_subplot(1, 1, 1)
	ax3.scatter(np.log(bp_size),np.log(sigma),linewidth=1,color='g',label=filein.split('/')[-1])
	ax3.set_title(title+": Peak width vs Base pair, Log-log plot")


	def func_lin(p,x):
		x=np.array(x).astype(float)
		k=p[0]
		b=p[1]
		return k*x+b

	pbest = leastsq(resid_cur,[-1.0,0],args=(func_lin,np.log(bp_size),np.log(sigma)),full_output=1,maxfev=10000,xtol=1.49012e-08,ftol=1.49012e-08)
	print "Finished fitting lin model"
	print pbest[3]
	pres=pbest[0]
	print pres
	finrmsd=((pbest[2]['fvec']**2).sum()/len(bp_size))**0.5
	print "Final lin model RMSD LS=",finrmsd
	ax3.plot(np.log(bp_size),func_lin(pres,np.log(bp_size)), linewidth=1,label='Lin fit',color='r')
	# ax1.plot(bp_size,func_ogston([Mu_00,B0,C0,D_00],bp_size), linewidth=1,label='Ogston model initial',color='r')
	ax3.annotate('Lin fit: k*x+b\nRMSD=%.3f\nk=%.2f\nb=%.2f'%(finrmsd,pres[0],pres[-1]), xy=(0.8, 0.7), xycoords='axes fraction')

	def func_quad(p,x):
		x=np.array(x).astype(float)
		a=p[0]
		b=p[1]
		c=p[2]
		return a*x*x+b*x+c

	pbest = leastsq(resid_cur,[-1.0,-1.0,0.0],args=(func_quad,np.log(bp_size),np.log(sigma)),full_output=1,maxfev=10000,xtol=1.49012e-08,ftol=1.49012e-08)
	print "Finished fitting quad model"
	print pbest[3]
	pres=pbest[0]
	print pres
	finrmsd=((pbest[2]['fvec']**2).sum()/len(bp_size))**0.5
	print "Final quad model RMSD LS=",finrmsd
	ax3.plot(np.log(bp_size),func_quad(pres,np.log(bp_size)), linewidth=1,label='Lin fit',color='b')
	# ax1.plot(bp_size,func_ogston([Mu_00,B0,C0,D_00],bp_size), linewidth=1,label='Ogston model initial',color='r')
	ax3.annotate('Quad fit: a*x*x+b*x+c\nRMSD=%.3f\na=%.2f\nb=%.2f\nc=%.2f'%(finrmsd,pres[0],pres[1],pres[-1]), xy=(0.8, 0.5), xycoords='axes fraction')


	ax3.set_ylabel("Log Width (sigma fitted)")
	ax3.set_xlabel("Log Base numer (M)")
	ax3.legend(loc='upper right', shadow=False)

	if(graphout):
		fig3.savefig(graphout+'_widthlog'+'.png',dpi=(300))
		os.system('convert %s -resize 25%% %s'%(graphout+'_widthlog'+'.png',graphout+'_widthlog'+'_small.png'))


	if(graphshow):
		plt.show()
	plt.close("all")



