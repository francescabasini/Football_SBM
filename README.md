# Football_SBM (2021)

This repository contains the data and code associated with the article/preprint by *F. Basini, V. Tsouli, I Ntzoufras, N. Friel*, 
*"Assessing competitive balance in the English Premier League for over forty seasons using a stochastic block model" (2021)*, currently deposited in [arXiv?].

## Repository structure:  

* ```README.md``` - you're reading it
* [/Data_Premier](https://github.com/basins95/Football_SBM/tree/master/Data_Premier) - contains the tables of results (match grids) of the Premier League championship from season 1977/78 to 2019/20 in csv format. e.g. *Result_Premier_0102.csv* for season 2001/02.
* [```LABEL_CORRECTION_AND_ANALYSIS.R```](https://github.com/basins95/Football_SBM/blob/master/LABEL_CORRECTION_AND_ANALYSIS.R): code to apply the label switching algorithm ([collpcm]) and carry out post-hoc analysis of the chain. 
* [```SBM_FUNCTIONS.R```](https://github.com/basins95/Football_SBM/blob/master/SBM_FUNCTIONS.R) contains all functions used in the MCMC algorithm. e.g. *get_loglik* returns the collapsed logllikelihood. 
* [```READ_TABLE_RESULTS.R```](https://github.com/basins95/Football_SBM/blob/master/READ_TABLE_RESULTS.R) code to load the result table and extract the relational pattern **y**.
* [```MCMC_main.R```](https://github.com/basins95/Football_SBM/blob/master/MCMC_main.R) the heart of the whole code which calls the other sources and runs the MCMC algorithm.
* [/Inference_results](https://github.com/basins95/Football_SBM/tree/master/Inference_results) - now empty, folders for each season analysed will be created inside this folder one the code is run, e.g. [/Inference_results/mcmc_Premier_Season_0102] for season 2001/02.

## Usage

* Clone repository.
* Open [```MCMC_main.R```] in RStudio and set Football_SBM as you working directory.
* In [line 28] set [season] to the season you want to analyse, e.g. 01/02 for 2001/02. (Provided that it is between 1977/78 and 2019/20)
* Run it all.

## Results
In the associated folder [/Inference_results/mcmc_Premier_Season_*season*] the following items will be saved:
* []
* [```WS_```]
*
*

