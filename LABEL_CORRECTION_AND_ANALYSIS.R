#### SCRIPT FOR LABEL SWITCHING CORRECTION AND ANALYSIS

# paths for elements
folder_path<-paste0("Inference_results//mcmc_Premier_Season_",season,"//")
after_object<-paste0("Premier_Season_",season,"_", (S-burn_in_level)/1000 ,"k_seed",my_seed)

# In case one wants to load it separately

# burn_in_level<-50000
# S<-250000
# 
# MCMC<-"AE_GS_MK_corrected"
# Championship<-"Premiere"
# my_seed<-1909
# 
# load(paste0("Inference_results//mcmc_Premier_Season_",season,
# "//WS_Premier_Season_", season, "_", (S-burn_in_level)/1000,
# "k_seed_",my_seed,".RData"))

####################################################################################
## MODELS SEPARATION AFTER LABEL SWITCH

# parameters to feed into label correction algorithm
SS<-dim(Z_seq_burned)[1]
Knumb<-rep(max(K_seq), SS)

# timing the agorithm
start_time_undo <- Sys.time()
# object output
Undone_Z_results<-collpcm::collpcm.undo.label.switching(Z=Z_seq_burned,                                                      Gsamp=Knumb)
end_time_undo <- Sys.time()
end_time_undo-start_time_undo

# Relabelled clusters
Undone_Z_seq<-Undone_Z_results$relab

# Record no. clusters for each draw
how_many_clust_labels = apply(Undone_Z_seq, MARGIN = 1, FUN = function(x){
  length(unique(x))
})
## Max no. of clusters obtained
K_max_seq = max(how_many_clust_labels)

# joining
Z_seq_burned_label_undone = cbind(Undone_Z_seq, how_many_clust_labels)

# SEPARATION
how_many_draws_k = rep(NA, K_max_seq)
for (kk in 1:K_max_seq){
  selected_Z_lab = Z_seq_burned_label_undone[Z_seq_burned_label_undone[,ncol(Z_seq_burned_label_undone)]==kk,]
  if (is.matrix(selected_Z_lab)){
    assign(paste0("Z_seq_burned_K", kk, "undone"), selected_Z_lab[,1:(ncol(selected_Z_lab)-1)])
    how_many_draws_k[kk] = dim(selected_Z_lab)[1]
  }else{
    assign(paste0("Z_seq_burned_K", kk, "undone"), selected_Z_lab[-length(selected_Z_lab)])
    how_many_draws_k[kk] = 1
  }
}


####################################################################################
## ANALYSIS AND PLOTS


# Plotting Traceplot for K as in Nobile and Fearnside
##########################################################
Noisy_True_K<-True_K_seq_burned
for(ii in 1:SS){
  Noisy_True_K[ii]<-Noisy_True_K[ii]+runif(1,min=0.1, max=0.9)
}
K_maxmax<-max(True_K_seq_burned)
K_minmin<-min(True_K_seq_burned)
#Plot directly the points
pdf(paste0(folder_path, "Ktrue_NF_",after_object,".pdf"),
    width = 11, height=8.5, paper="USr")
plot(Noisy_True_K, col="#00000033",ylab="K", xlab=" ", cex.lab=1.7,
     ylim=c(1,K_maxmax+1), yaxt="n")

at_points<-seq(1, K_maxmax)+.5
axis(side=2, at=at_points, labels=1:K_maxmax)
dev.off() 
##########################################################

## POSTERIOR OF K
##########################################################
posterior_k = rbind(how_many_draws_k/(S-burn_in_level)*100)
colnames(posterior_k) = 1:K_max_seq
print("Model percentages")
posterior_k
##########################################################


# Extracting allocation probabilities: P(Z|K)
# Write Summary table of posterior allocations for Teams
##########################################################
library(xtable)

for (kk in 1:K_max_seq){
  selected_Z_lab = get(paste0("Z_seq_burned_K", kk, "undone"))
  if (is.matrix(selected_Z_lab)){
    not_null = dim(selected_Z_lab)[1]!=0
    n_teams = ncol(selected_Z_lab)
  }else{
    not_null = is.vector(selected_Z_lab)
    n_teams = length(selected_Z_lab)
  }
  if (not_null){
    s_k = how_many_draws_k[kk]
    my_levs_current = unique(as.vector(selected_Z_lab))
    if (length(my_levs_current)!=kk){
    # at times, the clusterings might differ
      max_kk = max(length(my_levs_current), kk)
      current_alloc_matrix =  matrix(NA, n_teams, max_kk)
    }else{
      current_alloc_matrix =  matrix(NA, n_teams, kk)
    }
    for (tt in 1:n_teams){
      if (is.vector(selected_Z_lab)){
        current_alloc_matrix[tt,] = table(factor(selected_Z_lab[tt], levels = my_levs_current))/s_k*100
      }else{
        current_alloc_matrix[tt,] = table(factor(selected_Z_lab[,tt], levels = my_levs_current))/s_k*100
      }
    }
    assign(paste0("Allocation_prob_Model", kk),current_alloc_matrix)
    current_cluster_percentages = t(current_alloc_matrix)
    colnames(current_cluster_percentages)<-colnames(O)
    rownames(current_cluster_percentages)<- paste0("Cluster ", 1:length(my_levs_current))
    # Save matrix for each cluster model
    assign(paste0("Cluster_Percentages_Model", kk),current_cluster_percentages)
    
    summary_table<-xtable(current_cluster_percentages,
                          caption = paste0("Season ", substr(season, start = 1, stop = 2), "/",
                                           substr(season, start = 3, stop = 4), " for model with K = ", kk))
    # Save Latex Table
    print.xtable(summary_table, type="latex", caption.placement = 'top',
                 file=paste0(folder_path,"Summary_table_K",kk,"_",after_object,".txt"))
  }
}
##########################################################

# Print Results Table
##########################################################
Point_Table<-mapvalues(O, from=c(1,2,3), to=c(3,1,0))
Reversed_Table<-mapvalues(Point_Table, from=c(3,0), to=c(0,3))
# Italian tip: "Tabellone" means scoreboard or finale ranking table :) 
Tabellone<-as.matrix(rowSums(Point_Table, na.rm = TRUE)+
                       colSums(Reversed_Table, na.rm = TRUE))

Ordered_Tabellone<-cbind(sort(Tabellone, decreasing = TRUE))
rownames(Ordered_Tabellone)<-rownames(Tabellone)[order(Tabellone, decreasing = TRUE)]
Tabellone_table<-xtable(Ordered_Tabellone, digits = 0)
print.xtable(Tabellone_table, type="latex", file=paste0(folder_path,
                                                        "FinalTable_",after_object,".txt"))
##########################################################


## Stacked Plot: of posterior allocations for each team in each cluster
##########################################################
library(RColorBrewer)

coul = brewer.pal(K_max_seq, "Set1") 

for (kk in 2:K_max_seq){
  current_clust_model = get(paste0("Cluster_Percentages_Model", kk))
  dim_k = dim(current_clust_model)[1]
  rev_clust<-matrix(NA,dim_k, N)
  for(i in 1:dim_k){
    rev_clust[i,]<-rev(current_clust_model[i,])*100
  }
  table_clusters<-as.table(rev_clust)
  colnames(table_clusters)<-rev(colnames(O))
  
  # order table according to final table (Tabellone)
  Team_Names = cbind(colnames(Results), rownames(Results))
  rownames(Team_Names) = Team_Names[,2]
  # ordered_team_names (according to Tabellone)
  Team_Names_ordered = Team_Names[rownames(Ordered_Tabellone),]
  
  # ordered stacked table according to tabellone
  ordered_stacked_table = table_clusters[,rev(Team_Names_ordered[,1])]
  
  pdf(paste0(folder_path,"StackedPlot_K",kk,after_object,".pdf"),width = 7, height=20, paper="a4") 
  coul_kk = coul[1:dim_k]
  barplot(ordered_stacked_table,col=coul_kk, horiz = TRUE, border="white", ylab="Teams", cex.names=0.8)
  mtext(side=3,"Posterior allocation probabilities",line=2,cex=1.8)
  mtext(side=3,paste("K =", kk),line=0,cex=1.5)
  dev.off() 
}
##########################################################


# Mode, needed for estimated Heatmap
##########################################################
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

K_estimated = Mode(True_K_seq_burned)

