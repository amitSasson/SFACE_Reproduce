#  The SF-ACE reproducibility code 
Simulation and data code for reproducibility of results of the paper "The Subtype-Free Average Causal Effect for Heterogeneous Disease Etiology". 

All analyses use the “TheSFACE” R package. 

## Simulations 

The folder **Simulations** contains the following subfolders, each holding code and data of a different simulation study: 

- Study I - The sample size varied. 
- Study II -The unmeasured covariate U had on Y2 varied.
- Study III - The models were misspecified. 

In each subfolder, there is 
- Code to rerun this analysis both for the RR and the difference scale. 
- A  results (sub)subfolder, containing the results in RData files. 
- A combine file, to combine all the results to a TeX table or a plot. 


## Data

- In this paper, we used data from two large US cohorts, the  Nurses Health Study (NHS) and the Health Professionals Follow-up Study (HPFS). 

The folder **Data** contains the following files: 
- calc_weights.R -  script  for calculating weights to adjust for missing MSI values. 
- data_analysis.R -  script for analysing the data (under S-Monotonicity for both subtypes), both for the main scenario and for the two additional scenario (see section C of the SM). 
- data_analysis_sensitivity.R -  script for analysing the data using the suggested sensitivity methods (only for the main scenario).
- pseudo_dataset_for_weights_SFACE.csv - a pseudo dataset with a similar design with the same sample size and with the same number of covariates, same type of covariates, and same marginal distribution for the covariates as in our dataset. 
- table_covariates.R - script for creating descriptive tables, showing the distribution of the covariates and the models estimators (see section C in the SM)
- a results folder with the results of each scenario as an RData file. 

All the scripts can be run using the pseudo dataset. Of course, the results will be different than in the paper. 
## Note  

Please note that the code was run separately for each effect scale and subtype. To use this code for all scales and subtypes effectively, one might want to adjust the code. 