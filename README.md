#  The SF-ACE reproducibility 
Simulation and data code for reproducibility of results from the paper "The Subtype-Free Average Causal Effect for Heterogeneous Disease Etiology". 

All analyses use the “TheSFACE” R package. 

## Simulations 

The folder **Simulations** contains the following folders, each holding code and data from a difference simulation study: 

- Study I - The sample size varied between 1,000 and 50,000.
- Study II -The unmeasured covariate U had on Y2 varied
- Study III - The models were misspecified. 

In each file, there is code available to rerun this analysis both for the RR and the difference scale. 
There is also a results folder with the results of each scenario as an RData file. 


## Data

- In this paper, we used data from two large US cohorts, the  Nurses Health Study (NHS) and the Health Professionals Follow-up Study (HPFS). 

The folder **Data** contains the following folders: 
-  Main scenario - keep all
- Missing data - has outcome data
- Missing data - weights 
along with a pre process file which is relevant for all three scenarios. 

The first folder holds the code and results relevant for the main scenario analysed in the paper. The other two folders hold code and results for the two other scenarios referenced in the SM. 

In each scenario, the following steps were made: 
1. pre process 
2. calc weights (for missing subtype and missing outcome, if needed) 
3. data analysis 
4. sensitivity analysis (only for the main scenario). 

In each file, there is code available to rerun this analysis both for the RR and the difference scale. 
There is also a results folder with the results of each scenatro as an RData file. 

