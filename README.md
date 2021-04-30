# Football_SBM (2021)

This repository contains the data and code associated with the article/preprint by *F. Basini, V. Tsouli, I Ntzoufras, N. Friel*, 
*"Assessing competitive balance in the English Premier League for over forty seasons using a stochastic block model" (2021)*, currently deposited in [arXiv?].

## Dependencies:
All packages used are available on CRAN:
* [Data wrangling]: 
```
install.packages("plyr")
install.packages("stringi")
install.packages("seqinr")
```
* [Plotting]: 
```
install.packages("RColorBrewer")
install.packages("igraph")
```
* [Printing tables in Latex format]: 
```
install.packages("xtable")
```
* [Algorithm to adjust for label switching]: 
```
install.packages("collpcm")
```
Jason Wyse and Caitriona Ryan (2019). collpcm: Collapsed Latent Position Cluster Model
  for Social Networks. R package version 1.1. <https://CRAN.R-project.org/package=collpcm>

## Repository structure:  

* ```README.md``` - you're reading it
* [/Data_Premier](https://github.com/basins95/Football_SBM/tree/master/Data_Premier) - contains the tables of results (match grids) of the Premier League championship from season 1977/78 to 2019/20 in csv format. e.g. *Result_Premier_0102.csv* for season 2001/02.
* [```LABEL_CORRECTION_AND_ANALYSIS.R```](https://github.com/basins95/Football_SBM/blob/master/LABEL_CORRECTION_AND_ANALYSIS.R): code to apply the label switching algorithm ([collpcm]) and carry out post-hoc analysis of the chain. 
* [```SBM_FUNCTIONS.R```](https://github.com/basins95/Football_SBM/blob/master/SBM_FUNCTIONS.R) contains all functions used in the MCMC algorithm. e.g. *get_loglik* returns the collapsed loglikelihood. 
* [```READ_TABLE_RESULTS.R```](https://github.com/basins95/Football_SBM/blob/master/READ_TABLE_RESULTS.R) code to load the result table and extract the relational pattern **y**.
* [```MCMC_main.R```](https://github.com/basins95/Football_SBM/blob/master/MCMC_main.R) the heart of the whole code which calls the other sources and runs the MCMC algorithm.
* [/Inference_results](https://github.com/basins95/Football_SBM/tree/master/Inference_results) - now empty, folders for each season analysed will be created inside this folder one the code is run, e.g. [/Inference_results/mcmc_Premier_Season_0102] for season 2001/02.

## Usage

* Clone repository.
* Open ```MCMC_main.R``` in RStudio and set Football_SBM as your working directory.
* In ```line 28``` set ```season``` to the season you want to analyse, e.g. 01/02 for 2001/02. (Provided that it is between 1977/78 and 2019/20)
* Run it all.

## Results
In the associated folder [/Inference_results/mcmc_Premier_Season_*season*] that will be created, the following items will be available:
* [```Heatmap_Season_season.pdf```] the match grid of the season with results categorised by colour (See fig. 2 (b) in paper).
* [```StackedPlot_Premier_Season_season_iterk_seed.pdf```] is the stacked plot of posterior allocations of the teams in the league ordered by the final ranking.
* [```Summary_table_Premier_Season_season_iterk_seed.txt```] is the table of posterior allocations of the teams in the league, in latex table format.
* [```Tabellone_Premier_Season_season_iterk_seed.txt```] is the final ranking table or scoreboard for the season, in latex format.
* [```Ktrue_NF_Premier_Season_season_iterk_seed.pdf```] is the traceplot of K as plotted in Nobile and Fearnside (2007).
* [```Loglik_Premier_Season_season_iterk_seed.pdf```] is the traceplot of the collapsed loglikelihood.
* [```PostDensities_Premier_Season_season_iterk_seed.pdf```] contains plots for the kernel densities for posterior allocation in each cluster for each team in the league.
* [```TrueK_density_Premier_Season_season_iterk_seed.pdf```] shows the kernel density of teh posterior of K. 
* [```WS_Premier_Season_season_iterk_seed.RData```] the whole workspace is saved at the end of the analysis.

