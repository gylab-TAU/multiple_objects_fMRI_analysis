# multiple_objects_fMRI_analysis

This code was used to analyze the data for the paper:
Kliger, L., & Yovel, G. (2019). The functional organization of high-level visual cortex determines the representation of complex visual stimuli. (bioRxiv). https://doi.org/10.1101/2019.12.22.852293
Please cite if using this code.

The code is divided to four parts:
1.	fMRI analysis: including preprocessing and standard GLM analysis.
2.	Arranging data: loading the fMRI GLM .nii results files, including betas, percent signal change, contrast t-maps and ROI (region of interest) masks, into a combined .mat file.
3.	Perform searchlight analysis and predicting the response to multiple objects based on the response to single objects.
4.	Perform statistics on the results and plot figures.

The first three parts are written in Matlab and the fourth part is written in R.

The code uses spm12 (https://www.fil.ion.ucl.ac.uk/spm/doc/ ), marsbar (http://marsbar.sourceforge.net/), and stan (https://mc-stan.org/).


