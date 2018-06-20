# A stage-by-stage tutorial on analyzing HRF data. (Example 2)

Data set: Chicken (Gallus gallus) nucleosomes reconstituted on a well-positioning 601 DNA sequence, DNA is radioactively labeled on 5' end. Products of PCR with dideoxynucleotide triphosphates were run side-by-side the HRF gel lanes to assign HRF peaks to DNA sequence. Data set taken from [Morozov et al., NAR 2009](https://www.ncbi.nlm.nih.gov/pubmed/?term=19509309).

This is a generic example that outlines HRF data quantification from a gel image with 5' labeled DNA. It illustrates algorithm robustness with respect to quantification of independent experiments, noise tolerance and special fitting features such as constraints on the peak width and acconting for the peaks outside of the studied data range.

Python files implementing every stage are provided in this directory.

The contents of this directory can be conveniently downloaded by the following command once HYDROID is installed:
```
HYDROID_get_ex2
```

Video tutorial outlining the key stages is available [here](https://www.youtube.com/playlist?list=PL_GHGdsPyn0nVSvrRnyvuvkRCrNBjqeuC).

## Quantification of HRF gel electrophoresis images using [HYDROIDexp.py](../HYDROIDexp.py).
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
See [fitting results](results/gg_601_BS_a_fitted_intensities.png).
- Stage 5: Plot cleavage frequencies along DNA sequence  - see details in [exp_s5_plot_cl_freq.py](exp_s5_plot_cl_freq.py) and run it:
```
python exp_s5_plot_cl_freq.py
```
See the [resulting plots](results/gg_601_BS_a_cl_freq_profile.png).


## Compare different experimental profiles
-  Run [compare_exp.py](compare_exp.py).
```
python compare_exp.py
```
Get comparative plots for [top](results/exp_compar_TS.png) and [bottom](results/exp_compar_BS.png) strands from different gel lanes. Get compartive [plots](results/exp_compar_BS_TS.png) between top and bottom DNA strands to assess pseudosymmetry of nucleosome with respect to the two-fold dyad symmetry axis.
