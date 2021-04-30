# Football_SBM (2021)

This repository contains the data and code associated with the article/preprint by *F. Basini, V. Tsouli, I Ntzoufras, N. Friel*, 
*"Assessing competitive balance in the English Premier League for over forty seasons using a stochastic block model" (2021)*, currently deposited in [arXiv?].

## Repository structure:  

* ```README.md``` - you're reading it
* [/Data_Premier](https://github.com/basins95/Football_SBM/tree/master/Data_Premier) - contains the tables of results (match grids) of the Premier League championship from season 1978/79 to 2019/20 in csv format. e.g. *Result_Premier_0102.csv* for season 2001/02.
	* Note on the data: each csv file contains the results table where entries of the main diagonal are blank and the entry of match results is written as | 4~3 | for a match won 4 to 3 by row team against col team, where row team is playing home.
* [```LABEL_CORRECTION_AND_ANALYSIS.R```](https://github.com/basins95/Football_SBM/blob/master/LABEL_CORRECTION_AND_ANALYSIS.R) - code to apply the label switching algorithm (```collpcm```) and carry out post-hoc analysis of the chain. 
* [```SBM_FUNCTIONS.R```](https://github.com/basins95/Football_SBM/blob/master/SBM_FUNCTIONS.R) - contains all functions used in the MCMC algorithm. e.g. *get_loglik* returns the collapsed loglikelihood. 
* [```READ_TABLE_RESULTS.R```](https://github.com/basins95/Football_SBM/blob/master/READ_TABLE_RESULTS.R) - code to load the result table and extract the relational pattern **y**.
* [```MCMC_main.R```](https://github.com/basins95/Football_SBM/blob/master/MCMC_main.R) - the heart (:heartpulse:) of the whole code which calls the other source files and runs the MCMC algorithm.
* [/Inference_results](https://github.com/basins95/Football_SBM/tree/master/Inference_results) - now empty, folders for each season analysed will be created inside this folder once the code is run, e.g. [/Inference_results/mcmc_Premier_Season_0102] for season 2001/02.


## Usage

* Clone repository.
* Open ```MCMC_main.R``` in RStudio and set ```Football_SBM``` as your working directory.
* In ```line 28``` set ```season``` to the season you want to analyse, e.g. 01/02 for 2001/02. (Provided that it is between 1978/79 and 2019/20)
* Run it all.

## Results
In the associated folder ```/Inference_results/mcmc_Premier_Season_*season*``` that will be created, the following items will be available:
* ```Heatmap_Season_*season*.pdf``` the match grid of the season with results categorised by colour (See fig. 2 (a) in paper).
* ```Heatmap_Estimated_Season_*season*.pdf``` the permuted match grid of the season, from *a posteriori* analysis, with teams listed according to the most likely block membership, results categorised by colour (See fig. 2 (b) in paper). 
	Plotted and saved only when the MAP for K is larger than 1.
* ```StackedPlot_Premier_Season_*season*_*iter*k_seed_*seed*.pdf``` is the stacked plot of posterior allocations of the teams in the league ordered by the final ranking.
	* Note: The subtitle will print ```K = *K*```, the biggest number K of non empty clusters found by the algorithm search.
* ```Summary_table_Premier_Season_season_iterk_seed.txt``` is the table of posterior allocations of the teams in the league, in latex table format.
* ```FinalTable_Premier_Season_*season*_*iter*k_seed_*seed*.txt``` is the final league table for the season, in latex format.
* ```Ktrue_NF_Premier_Season_*season*_*iter*k_seed_*seed*.pdf``` is the jittered traceplot of K as plotted in Nobile and Fearnside (2007).
* ```WS_Premier_Season_*season*_*iter*k_seed_*seed*.RData``` the whole workspace is saved at the end of the analysis.

## Dependencies:
All packages used are available on CRAN.
* #### Data wrangling: 
```
install.packages("plyr")
install.packages("stringi")
install.packages("seqinr")
```
* #### Plotting: 
```
install.packages("RColorBrewer")
install.packages("lattice")
```
* #### Printing tables in Latex format: 
```
install.packages("xtable")
```
* #### Algorithm to adjust for label switching: 
```
install.packages("collpcm")
```
Jason Wyse and Caitriona Ryan (2019). collpcm: Collapsed Latent Position Cluster Model
  for Social Networks. R package version 1.1. <https://CRAN.R-project.org/package=collpcm>

