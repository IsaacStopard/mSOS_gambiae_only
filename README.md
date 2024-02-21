# **Data and code for: Estimating of the effects of temperature on transmission of the human malaria parasite, Plasmodium falciparum**
## **Authors: Eunho Suh, Isaac J. Stopard, Ben Lambert, Jessica L. Waite, Nina L. Dennington, Thomas Churcher, and Matthew B. Thomas**

This GitHub repository provides all code necessary to run the analyses for this paper. We ran the model using R version 4.2.1 and Stan version 2.21.0.

**Stan models**

:one: mSOS_multi_8.stan - implements the pooled model in Stan. Note to reduce computational time the likelihood for each unique combination of datapoint values is only calculated once.

:two: mSOS_multi - implements the independent model fits.

**Running**

:one: run_model_multi_temp_only.R - provides R code to run the Stan models.

:two: vectorial_capacity_cluster.R - runs the vectorial capacity calculations. Will only work on the MRC Global Infectious Disease Analysis cluster.

:three: biting_rate.R - fits the Brière function to the gonotrophic cycle length data and plots the results.

**Helper functions**

:one: plot_model_temp_only.R - plot the model fits and run some analyses.

:two: functions_temp_only.R - consistent functions to load in and wrangle the data.

:three: plotting_functions_temp_only.R - functions to help plotting.

:four: wrangle_Eunho_new_data.R - wrangle the data into a format for R.

:five: read_data.R - reads data in the correct format for model fitting and plotting.

:six: VC_functions.R - functions required to calculate the expected number of infectious bites per infected mosquito.

Notes

⚠️ Please note that the Stan outputs are not saved on the Github repository, so run_model_multi_temp_only.R will need to be run before attempting to visualise the results.

⚠️ This code is released with no support and we cannot endorse any outputs from the model other than those we generate.

⚠️ This code is under development, so code may change without notice.

⚠️ No liability is accepted by the authors.




