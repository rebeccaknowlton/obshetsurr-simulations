This repository contains all code to reproduce the simulation results in Knowlton & Parast, L. (2025). Assessing Surrogate Heterogeneity in Real World Data Using Meta-Learners.

To reproduce results, proceed through the following steps. Download all files in this repository and place them in a single folder. These include function files, helper files, and two main files. 
The two main files are: sims_obs_master_011725.R and sims_obs_readin.R. The files called sims_runX.R (Setting 1), sims_runX_2.R (Setting 2), sims_runX_3.R (Setting 3), and sims_runX_3.R (Setting 4) 
where X is a number, each run 50 iterations of the sims_obs_master_011725.R file. For Setting 1 for example, you will need to run all 20 files: sims_runX.R (designed to run in parallel); 
this will write results to various files in the folder. Then use the sims_obs_readin.R file to read in all the results and create the table; make sure to set the setting and directory correctly 
at the top of the sims_obs_readin.R file.
