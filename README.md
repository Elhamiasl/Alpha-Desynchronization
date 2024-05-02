# Alpha Desynchronization across Infancy
This repository contains MATLAB and R code utilized for the visualization and analysis of data regarding the development of alpha desynchronization across infancy.

## MATLAB Code
After preprocessing the EEG data and obtaining Fast Fourier Transform (FFT) and desynchronization values per electrode and frequency for each subject, the MATLAB code (alpha_desycnh_ME_2023.m) in this repository is used for the following purposes:

- Generating spectra plots to visualize the pre-processed data.
- Creating headplots to illustrate the distribution of alpha desynchronization across the brain.
- Conducting one-sample t-tests against the null hypothesis value of zero (0) to examine significant desynchronization across each electrode per group and per condition.


## R Code
The R code (alpha_desycnh_crosssectional_longitudinal_ME_2023.Rmd) provided here is used for the following analyses:

- Performing ANOVA analysis to examine between-age and between-condition differences in alpha desynchronization, as well as their interaction effects.
- Conducting follow-up t-tests to further explore significant effects.
- Visualizing main effects and interactions using appropriate plots.
