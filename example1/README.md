# Step-by-step tutorial on analyzing HRF data of a protein-DNA complex and comparing it to a PDB structure. (Example 1)
System: *S. cerevisiae* centromeric nucleosome reconstituted on a well-positioning 601TA DNA sequence, DNA is radioactively labeled on 3' end. Data set taken from [Shaytan et al., unpublished](https://www.ncbi.nlm.nih.gov/pubmed/)

This is a generic example that outlines (1) HRF data quantification from a gel image with 3' labeled DNA, (2) prediction of theoretical profiles from PDB strucutres, (3) data comparison.

Python files implementing every step are provided in this directory.

Document [H-SASA_params.md](H-SASA_params.md) further outlines details of theoretical cleavage profiles calculation and their dependence on various parameters.


## Quantification of HRF gel electrophoresis images using [HYDROIDexp.py](../HYDROIDexp.py).
###Step-by-step usage example:
- Step 1: Extract lane profiles from gel images via ImageJ - see instructions in [exp_s1_extract_lp.md](exp_s1_extract_lp.md).
- Step 2: Assign band peak locations - run [exp_s2_assign_peaks.py](exp_s2_assign_peaks.py)
- Step 3: Assign peaks to DNA sequence - run [exp_s3_call_peaks.py](exp_s3_call_peaks.py)
- Step 4: Quantify cleavage frequencies by fitting a model to the data  - run [exp_s4_fit_model.py](exp_s4_fit_model.py), see [fitting results](results/scCSE4_601TA_BS_fitted_intensities.png).
- Step 5: Plot cleavage frequencies along DNA sequence  - run [exp_s5_plot_cl_freq.py](exp_s5_plot_cl_freq.py), get the [resulting plots](results/scCSE4_601TA_BS_cl_freq_profile.png).

## Prediction of cleavage intensities from PDB-structures using [HYDROIDpred.py](HYDROIDpred.py).
###Step-by-step example:
- Step 1: Prepare PDB file with hydrogen atoms - see instructions in [pred_s1_prep_pdb.md](pred_s1_prep_pdb.md).
- Step 2: Estimate theoretical hydroxyl-footprinting cleavage profiles through calculating H-SASA profile - run [pred_s2_calc_H-SASA.py](pred_s2_calc_H-SASA.py)
- Step 3: Plot H-SASA profiles - run [pred_s3_plot_H-SASA.py](pred_s3_plot_H-SASA.py), resulting plots [here](results/scCSE4_601TA_TS_H-sasa.png) and [here](results/scCSE4_601TA_TS_H-sasa_MD.png). Plot [simulated gel lane profiles](results/scCSE4_601TA_TS_H-sasa_simulated.png) from theoretical profiles.

## Compare profiles
- Experimental vs H-SASA: run [comp_plot_exp_vs_pred.py](comp_plot_exp_vs_pred.py), get comparative plots for [top](results/exp_vs_H-SASA_TS.png) and [bottom](results/exp_vs_H-SASA_BS.png) strands.

- Top and bottom strands of experimental profiles: run [comp_plot_exp.py](comp_plot_exp.py). Get compartive [plots](results/exp_compar_BS_TS.png) between top and bottom DNA strands to assess pseudosymmetry of nucleosome with respect to the two-fold dyad symmetry axis.
