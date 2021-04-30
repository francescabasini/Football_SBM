#### SCRIPT FOR LABEL SWITCHING CORRECTION AND ANALYSIS

# paths for elements
folder_path<-paste0("Inference_results//mcmc_Premier_Season_",season,"//")
after_object<-paste0("Premier_Season_",season,"_", (S-burn_in_level)/1000 ,"k_seed",my_seed)

# checks
is.matrix(Z_seq_burned)
is.vector(K_seq_burned)

dim(Z_seq_burned)

####################################################################################
## LABEL CORRECTION

# parameters to feed into label correction algorithm
SS<-dim(Z_seq_burned)[1]
Knumb<-rep(max(K_seq), SS)

# timing the agorithm
start_time_undo <- Sys.time()
#object output
Undone_Z_results<-collpcm::collpcm.undo.label.switching(Z=Z_seq_burned,
                                                        Gsamp=Knumb)
end_time_undo <- Sys.time()
end_time_undo-start_time_undo

# Used all "Undone" elements which are label corrected

Undone_Z_seq<-Undone_Z_results$relab
K_est<-Undone_Z_results$numcomponents
# access allocation probabilities from object output of collpcm
Allocation_Probs<-Undone_Z_results$label.probs
Cluster_percentages<-t(Allocation_Probs[[1]])

is.matrix(Cluster_percentages)

# Names of teams
colnames(Cluster_percentages)<-colnames(O)
rownames(Cluster_percentages)<-seq(1,K_est)
for(ii in 1:K_est){
  rownames(Cluster_percentages)[ii]<-paste0("Cluster ",ii)
}

####################################################################################
## ANALYSIS AND PLOTS

## Plotting posterior densities, although continuous here
############################################################
pdf(paste0(folder_path, "PostDensities_",after_object,".pdf"),
      width = 11, height=8.5, paper="USr")
# 2. Create a plot
par(mfrow=c(2,3))
for(i in 1:N){
  plot(density(Undone_Z_seq[,i]), col=i,
       main=paste("Posterior allocation density of team", i, colnames(O)[i]))
}
# Close the pdf file
dev.off() 
##########################################################


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



# Write Summary table of posterior allocations for Teams
##########################################################
library(xtable)
colnames(Cluster_percentages)<-colnames(O)
summary_table<-xtable(Cluster_percentages)
print.xtable(summary_table, type="latex", 
             file=paste0(folder_path,"Summary_table_",after_object,".txt"))
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
                                                        "Scoreboard_",after_object,".txt"))
##########################################################


## Stacked Plot: of posterior allocations for each team in each cluster
##########################################################
rev_clust<-matrix(NA,K_est, N)
for(i in 1:K_est){
  rev_clust[i,]<-rev(Cluster_percentages[i,])*100
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

library(RColorBrewer)
pdf(paste0(folder_path,"StackedPlot_",after_object,".pdf"),width = 7, height=20, paper="a4") 
coul = brewer.pal(K_est, "Set1") 
barplot(ordered_stacked_table,col=coul, horiz = TRUE, border="white", ylab="Teams", cex.names=0.8)
mtext(side=3,"Posterior allocation probabilities",line=2,cex=1.8)
mtext(side=3,paste("K =", K_est),line=0,cex=1.5)
dev.off() 
##########################################################


## Plotting loglikelihood 
##########################################################
pdf(paste0(folder_path, "Loglik_traceplot_",after_object,".pdf"),
    width = 11, height=8.5, paper="USr")
plot(log_lik_seq_burned, type="l", main = "Loglikelihood traceplot")
dev.off()
##########################################################


## kernel density of True K
##########################################################
True_K_from_Z_burned<-BURN_BABY_BURN(seq_to_burn = True_K_from_Z,
                                     burn_in_level = burn_in_level,
                                     maxS=S)
pdf(paste0(folder_path, "TrueK_density_",after_object,".pdf"),
    width = 11, height=8.5, paper="USr")
plot(density(True_K_from_Z_burned), type="l", ylab="K")
dev.off()
##########################################################


