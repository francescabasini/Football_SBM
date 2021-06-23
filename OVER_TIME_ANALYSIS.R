# OVER TIME ANALYSIS - OVER ALL SEASONS

# IMPORTANT: Set this as current working directory

# RUN THE ANALYSIS FOR ALL SEASONS IN THE DATA FOLDER
####################################################################################
years = c(78:99, 0:20)
season_labels = c(0)
season_numbers = c(0)
for(yy in 1:(length(years)-1)){
  year1<-years[yy]
  year2<-years[yy+1]
  if(nchar(year1)==1){
    year1<-as.character(paste0("0",year1))
  }
  if(nchar(year2)==1){
    year2= as.character(paste0("0",year2))
  }
  season<-paste0(year1, year2) 
  season_lab = paste0(year1, "/",year2) 
  
  season_numbers = append(season_numbers, season)
  season_labels = append(season_labels, season_lab)
}

season_labels = season_labels[-1]
season_numbers = season_numbers[-1]

posterior_list = vector(mode = "list", length = length(season_numbers))

for(yy in 1:length(season_numbers)){
  season = season_numbers[yy]
  source("MCMC_main.R")
  posterior_list[[yy]] = posterior_k
  rm(list=setdiff(ls(), c("season_numbers","season_labels", "posterior_list")))
  cat("\014")
}


#new folder path
folder_path_OVER = "Inference_results//mcmc_Premier_OVER_TIME_ANALYSIS"
# Create OVER_TIME_ANALYSIS directory insinde "Inference_results" folder
dir.create(folder_path_OVER)
# Saving the posterior of K on a separate WorkSpace inside "OVER_TIME_ANALYSIS" folder
save.image(paste0(folder_path_OVER, "//mcmc_Premier_OVER_TIME.RData"))
#load(paste0("Inference_results//mcmc_Premier_OVER_TIME_ANALYSIS//mcmc_Premier_OVER_TIME.RData"))


  # Print Posterior of K table
####################################################################################
library(xtable)

# Using what is written below
max_k_over = max(unlist(lapply(posterior_list, FUN = function(x){length(x)}),
                        recursive = TRUE, use.names = TRUE))

my_fun = function(x, max_k_over){
  ll = length(x)
  if (ll<max_k_over){
    clust_diff = max_k_over-ll
    for (jj in 1:clust_diff){
      x = cbind(x, 0.0)
    }
    colnames(x) = 1:max_k_over
    x
  }else{
    x
  }
}

posterior_list_updates = lapply(posterior_list, FUN = my_fun, max_k_over = max_k_over)

posterior_k_matrix <- matrix(unlist(posterior_list_updates), ncol = max_k_over, byrow = TRUE)
colnames(posterior_k_matrix) = 1:max_k_over
row.names(posterior_k_matrix) = season_labels

# into x table
posterior_k_matrix_table<-xtable(posterior_k_matrix, digits = 2)
print.xtable(posterior_k_matrix_table, type="latex", file=paste0(folder_path_OVER,
                                                        "//Posterior_K_table_OVER_TIME.txt"))
####################################################################################


# Computing the probability of belonging to the top block
####################################################################################
# P(top block we need posterior of K and the )
# P(k=1) + P(K=2)*P(top|k=2) + ...

#Useful Info
S = 250000
burn_in_level<-50000
my_seed<-1909

for (yy in 1:length(season_numbers)){
  season = season_numbers[yy]
  folder_path<-paste0("Inference_results//mcmc_Premier_Season_",season,"//")
  after_object<-paste0("Premier_Season_",season,"_", (S-burn_in_level)/1000 ,"k_seed",my_seed)
  #load
  load(paste0("Inference_results//mcmc_Premier_Season_",season,
              "//WS_Premier_Season_", season, "_", (S-burn_in_level)/1000,
              "k_seed_",my_seed,".RData"))
  #------------------
  # Computing P(top)
  max_k_season = length(posterior_k)
  acronym_winner<-colnames(O)[pos_winner<-which(rownames(Results)==rownames(Ordered_Tabellone)[1])]
  
  P_top = matrix(NA, N, max_k_season)
  P_top[,1] = rep(posterior_k[1], N)*1
  for (kk in 2:max_k_season){
    current_clust_model = get(paste0("Cluster_Percentages_Model", kk))
    # find the label of top block
    pos_winner_labs = which(colnames(current_clust_model)==acronym_winner)
    winner_lab_model = which.max(current_clust_model[,pos_winner_labs])
    p_top_block_k = current_clust_model[winner_lab_model,]*(posterior_k[kk]/100)
    P_top[,kk] = p_top_block_k
  }
 # sum over
  p_top_over_k = t(as.matrix(rowSums(P_top)))
  colnames(p_top_over_k) = colnames(O)
  row.names(p_top_over_k) = season
  assign(paste0("P_top_block_season", season), p_top_over_k)
}
####################################################################################


# Probability of Top block for each Championship Over time (Figure 6)
####################################################################################
balance_or_not = apply(posterior_k_matrix, MARGIN = 1, FUN = function(x) which.max(x))
colour_bean = ifelse(balance_or_not==1, "navajowhite2", "#B0E0E6")

pdf(paste0(folder_path_OVER, "//TopBlock_prob_datapoints_JITTERED.pdf"), width = 20, height=10)

season_points = c(1:length(season_numbers))
plot(season_points,rep(c(0,1),length(season_points)/2), "n", xact = NULL, xlab = "Season", ylab = "P(top block)",
     main = "Posterior probabilities of belonging to the top block over time", 
     cex.main = 1.9,xaxt="n", las = 1, cex.axis = 1.2, cex.lab = 1.4)
axis(1, at=1:length(season_numbers), labels=season_labels, cex.axis = 1.1)
for (yy in 1:length(season_numbers)){
  season = season_numbers[yy]
  current_p_top = get(paste0("P_top_block_season", season))/100

  x_vals = rep(yy, length(current_p_top))
  x_jitt = rnorm(length(x_vals), 0, 0.1)
  y_jitt= rnorm(length(current_p_top), 0, 0.01)
  points(x_vals+x_jitt, current_p_top+y_jitt, cex = 1.35, pch = 21, col = "black", bg = colour_bean[yy])
  # adding mean
  #lines(c(seas-0.4, seas +0.4), rep(mean(probs_top),2), lwd = 2.5)
}
dev.off()


# Table for number of teams in top block
####################################################################################
how_many_on_top = rep(NA, length(season_numbers))
for (yy in 1:length(season_numbers)){
  season = season_numbers[yy]
  current_p_top = get(paste0("P_top_block_season", season))/100
  how_many_top_this_season = sum(current_p_top>=0.5)
  how_many_on_top[yy] = how_many_top_this_season
}

# Print
pdf(paste0(folder_path_OVER, "//TopBlock_Size_barplot.pdf"), width = 20, height=10)
bp  =barplot(how_many_on_top, main= "Size of top block over seasons", xlab = "Season",
             ylab= "Number of teams in top block", col = colour_bean, ylim = c(0, max(how_many_on_top)+1.5))
grid(nx=NA, ny=NULL)
axis(1, at=bp, labels=season_labels, cex.axis = 1.1)
# Add numbers
#text(x = bp, y = how_many_on_top, label = how_many_on_top, pos = 3, cex = 0.8, col = "black")
dev.off()
####################################################################################


# Print table
####################################################################################
top_block_size_table = cbind(how_many_on_top)
rownames(top_block_size_table) = season_labels
colnames(top_block_size_table) = "Size of Top Block"
top_block_size_table_table<-xtable(top_block_size_table, digits = 0)
print.xtable(top_block_size_table_table, type="latex", file=paste0(folder_path_OVER,
                                                                 "//Size_top_block_table.txt"))

####################################################################################


# Top block prob over time for KNOWN TEAMS    
####################################################################################

# list of all teams labels
All_teams_labs = c()
for (yy in 1:length(season_numbers)){
  season = season_numbers[yy]
  current_p_top = get(paste0("P_top_block_season", season))/100
  All_teams_labs = c(All_teams_labs, colnames(current_p_top))
}
All_teams_labs = unique(All_teams_labs)

# Extracting p(top block) for each team in all leagues over time
# absent are assigned NA
Team_prob_over_time = matrix(NA, length(All_teams_labs), length(season_numbers))
for (yy in 1:length(season_numbers)){
  season = season_numbers[yy]
  current_p_top = get(paste0("P_top_block_season", season))/100
  current_colnames = colnames(current_p_top)
  for(tt in 1:length(All_teams_labs)){
    current_team_pos = which(current_colnames == All_teams_labs[tt])
    #Not in that league
    if (identical(current_team_pos, integer(0))){
      Team_prob_over_time[tt, yy] = NA
    }else{
      Team_prob_over_time[tt, yy] = current_p_top[current_team_pos]
    }
  }
}

colnames(Team_prob_over_time) = season_labels
row.names(Team_prob_over_time) = All_teams_labs


# PLOTTING P(TOP BLOCK) OVER TIME FOR EACH TEAM
####################################################################################
xleft_team<-c(1:length(season_labels)-0.5)
ybottom_team<-rep(0, length(xleft_team))
xright_team<-c(2:(length(season_labels)+1)-0.5)
ytop_BALANCE<-ifelse(balance_or_not==1, 1.05, 0)
ytop_UNBALANCE<-ifelse(balance_or_not!=1, 1.05, 0)

# Creating a separate directory inside "OVER TIME" for all the figures for each team
folder_path_OVER_for_TEAMS = paste0(folder_path_OVER, "//Over_time_Teams")
dir.create(folder_path_OVER_for_TEAMS)

for(tt in 1:length(All_teams_labs)){
  current_team = row.names(Team_prob_over_time)[tt]
  pdf(paste0(folder_path_OVER_for_TEAMS,"//",current_team ,"_Top_block_over_time.pdf"),
      width = 15, height=8)
  # p(top) for team tt
  pp_top_team_over_time = Team_prob_over_time[tt, ]
  plot(pp_top_team_over_time, pch=19,ylim=c(0,1), ylab="Pr(strong cluster)",
       main=current_team, xaxt = "n",
       xlab=c("Season"), col=ifelse(pp_top_team_over_time>=0.5,"purple", "chocolate1"), yaxt = "n", cex=1)
  axis(1, at=seq(2,length(season_labels),by=3), labels=season_labels[seq(2, length(season_labels), by=3)])
  axis(2, at=seq(0,1, by=0.2), labels=seq(0,1, by=0.2))
  # Adding balance bars
  rect(xleft=xleft_team, ybottom=ybottom_team, xright=xright_team,
       ytop=ytop_BALANCE, col="navajowhite2", border=NA,density=400)
  rect(xleft=xleft_team, ybottom=ybottom_team, xright=xright_team,
       ytop=ytop_UNBALANCE, col="#B0E0E6", border=NA,density=400)
  #
  for(yy in 1:length(season_labels)){
    if(!(is.na(pp_top_team_over_time[yy])))
      segments(x0=yy, y0=0, x1=yy,
               y1=pp_top_team_over_time[yy], lwd = 3.5) 
  }
  points(pp_top_team_over_time, col = "black",bg=ifelse(pp_top_team_over_time>=0.5,"purple", "chocolate1"),
         cex = 4.5, pch =21)
  dev.off()
}


# Save all
save.image(paste0(folder_path_OVER, "//mcmc_Premier_OVER_TIME.RData"))


# Print P(Top block) as xtable for season 18/19 (Paper)
####################################################################################
# print.xtable(xtable(t(P_top_block_season1819)), type="latex", file=paste0("Inference_results",
#                                                                  "//P_top_block_1819.txt"))

