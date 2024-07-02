# sero-spatial-test
code associated with the manuscript entitled "Geospatial joint modeling of vector and parasite serology to microstratify malaria transmission"

In order to comply with the ethical regulations of our study, we are unable to share the true data collected in our study.
Instead, we have provided test datasets that use the true geolocations of the villages with the sample size and counts of SG6, CSP and PCR generated as random numbers.
Therefore, please note that the models, outputs and validations are not the same as reported in the manuscript.

We have provided the Shapefiles, Raw data (rasters) and test datasets to run the scripts in the following order:

. 00) Reference_grid - establishes the reference grid (Bago (East), Kayin and Kayah states)
. 01) Clean_static_datat_5sd - cleans, transforms and standardises our satellite-derived covariates (rasters)
. 02) SG6 spatial - model selection, run, prediction and validation for SG6 only
. 02) SG6 spatial seasons - SG6 spatial - run and prediction for SG6 per season only
. 02) Joint SG6 CSP PCR - model selection, run and prediction for joint spatial model 
. VIF_func - iteratively calculates VIF and removes highest 
. 03) joint_validation_gsg6_na - validation for joint model that sets just gsg6 to NA in the test datasets
. 03) joint_validation_all_na - validation for joint model that sets all outcomes (gsg6, csp, pcr) to NA in the test datasets
