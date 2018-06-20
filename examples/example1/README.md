# A stage-by-stage tutorial on analyzing HRF data of a protein-DNA complex and comparing it to a PDB structure.(Example 1)
System: *S. cerevisiae* centromeric nucleosome reconstituted on a well-positioning 601TA DNA sequence, DNA is radioactively labeled on 3' end. Maxam-Gilbert sequencing reactions products were run side-by-side with the HRF gel lanes to assign HRF peaks to the sites on DNA sequence. Data set is taken from [Shaytan et al., NAR (2017)](https://www.ncbi.nlm.nih.gov/pubmed/28934480).

This is a generic example that outlines (1) HRF data quantification from a gel image with the 3' labeled DNA, (2) prediction of theoretical profiles from PDB strucutres, (3) data comparison.

Python files implementing every stage are provided in this directory.

The contents of this directory can be conveniently downloaded by the following command once HYDROID is installed:
```
HYDROID_get_ex1
```

Document [H-SASA_params.md](H-SASA_params.md) further outlines details of theoretical cleavage profiles calculation and their dependence on various parameters.

Video tutorial outlining the key stages is available [here](https://www.youtube.com/playlist?list=PL_GHGdsPyn0nVSvrRnyvuvkRCrNBjqeuC).

## Quantification of HRF gel electrophoresis images using [HYDROIDexp.py](../hydroid/HYDROIDexp.py).
### A stage-by-stage usage example:
- Stage 1: Extract lane profiles from gel images via ImageJ - see instructions in [exp_s1_extract_lp.md](exp_s1_extract_lp.md).
- Stage 2: Assign band peak locations - see details in [exp_s2_assign_peaks.py](exp_s2_assign_peaks.py) and run it:
```
python exp_s2_assign_peaks.py
```
- Stage 3: Assign peaks to DNA sequence - see details in [exp_s3_call_peaks.py](exp_s3_call_peaks.py) and run it:
```
python exp_s3_call_peaks.py
```
- Stage 4: Quantify cleavage frequencies by fitting a model to the data  - see details in  [exp_s4_fit_model.py](exp_s4_fit_model.py) and run it:
```
python exp_s4_fit_model.py
```
See [fitting results](results/scCSE4_601TA_BS_fitted_intensities.png).
- Stage 5: Plot cleavage frequencies along DNA sequence  - see details in [exp_s5_plot_cl_freq.py](exp_s5_plot_cl_freq.py) and run it:
```
python exp_s5_plot_cl_freq.py
```
See the [resulting plots](results/scCSE4_601TA_BS_cl_freq_profile.png).

## Prediction of cleavage intensities from PDB-structures using [HYDROIDpred.py](HYDROIDpred.py).
### A stage-by-stage example:
- Stage 1: Prepare PDB file with hydrogen atoms - see instructions in [pred_s1_prep_pdb.md](pred_s1_prep_pdb.md).
- Stage 2: Estimate theoretical hydroxyl-footprinting cleavage profiles through calculating H-SASA profile - see details in  [pred_s2_calc_H-SASA.py](pred_s2_calc_H-SASA.py) and run it:
```
python pred_s2_calc_H-SASA.py
```
- Stage 3: Plot H-SASA profiles - see details in [pred_s3_plot_H-SASA.py](pred_s3_plot_H-SASA.py) and run it:
```
python pred_s3_plot_H-SASA.py
```
See resulting plots [here](results/scCSE4_601TA_TS_H-SASA.png) and [here](results/scCSE4_601TA_BS_H-SASA.png). Also [simulated gel lane profiles](results/scCSE4_601TA_TS_H-SASA_simulated.png) from theoretical profiles may be plotted.

## Compare profiles
- Experimental vs H-SASA: run [comp_plot_exp_vs_pred.py](comp_plot_exp_vs_pred.py)
```
python comp_plot_exp_vs_pred.py
```
to get comparative plots for [top](results/exp_vs_H-SASA_TS.png) and [bottom](results/exp_vs_H-SASA_BS.png) strands.

- Top and bottom strands of experimental profiles: run [comp_plot_exp.py](comp_plot_exp.py)
```
python comp_plot_exp.py
```
to get compartive [plots](results/exp_compar_BS_TS.png) between top and bottom DNA strands to assess pseudosymmetry of nucleosome with respect to the two-fold dyad symmetry axis.
