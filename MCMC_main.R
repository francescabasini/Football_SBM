######################################
####       FOOTBALL SBM     ##########
######################################

# Written by Francesca Basini (Francesca.Basini@warwick.ac.uk), for the preprint
#
#     F. Basini, V. Tsouli, I Ntzoufras, N. Friel
#     "Assessing competitive balance in the English Premier League
#     for over forty seasons using a stochastic block model" (2021)
#

# The MCMC algorithm is made of the following moves, choosen uniformly at random:
# (GS move) - Metropolis-within-Gibbs move, reallocates one node at time until full sweep,
#             changes Z only
# (MK move) - Metropolis move, to increase or decrease K, leaves Z unchanged.
# (AE move) - Metropolis-Hastings move to eject or absorb a cluster,
#             changes both Z and K. 

####################################################################################

#### IMPORTANT: Set THIS as working directory
#### Uncomment line29 if you want to run it for only one season

## Clean Environment
#rm(list = ls())
#cat("\014")

## Choose Premier Season
#season<-"2021"
  
####################################################################################

#Create directory inside Inference_results for the selected season
dir.create(paste0("Inference_results//mcmc_Premier_Season_", season))

## set seed
my_seed<-1909   

## READ DATA FROM RESULTS TABLE
source("READ_TABLE_RESULTS.R")

###   VISUALIZE AND SAVE HEATMAP FOR RESULTS TABLE
####################################################################################
library("lattice")
palf <-colorRampPalette(c("green3", "yellow", "red1"))

# option to print season with 4 figures
year1 = substr(season, 1, 2)
year2 = substr(season, 3, 4)
if((year1)>20){
  year1<-as.character(paste0("19",year1))
}else{
  year1<-as.character(paste0("20",year1))
}
if((year2)==0){
  year2
  year2<-as.character(paste0("200",year2))
}
season_lab = paste0(year1, "/",year2) 

print(levelplot(t(O[nrow(O):1,]),
                col.regions=palf(100), xlab = NULL, ylab = NULL, colorkey = FALSE,
                main =  list(label=paste0("Results table season: ", season_lab), cex=2.4),
                scales = list(list(alternating=1), x=list(cex=1),y=list(cex=0.8)),
))


# save it in folder
pdf(paste0("Inference_results//mcmc_Premier_Season_",season,
           "//Heatmap_Season_", season,".pdf"),width = 10, height=10)

print(levelplot(t(O[nrow(O):1,]),
                col.regions=palf(100), xlab = NULL, ylab = NULL, colorkey = FALSE,
                main =  list(label=paste0("Results table season: ", season_lab), cex=2.4),
                scales = list(list(alternating=1), x=list(cex=1),y=list(cex=0.8)),
                ))

dev.off()
####################################################################################



## Load functions used in the MCMC algorithm
source("SBM_FUNCTIONS.R")

library(seqinr)   #necessary for the swap function
library(plyr)     #to use mapvalues

####################################################################################
## Parameters setting

#hyperparameters
gamma_0=1
beta_0=1

#Fixing K max 
K_max=5

set.seed(1605)

#random initialization of K and Z
K_current<-sample.int(n=K_max, size=1)
Z_current<-sample.int(n=K_current, size=N, replace=TRUE)

# set number of steps to discard
burn_in_level<-50000
# Number of steps to retain + to discard
S = 200000+burn_in_level

####################################################################################

# Initialize structures for chain output
K_seq=numeric(S)
K_seq[1]<-K_current

Z_seq=matrix(nrow=S, ncol=N)
Z_seq[1,]=Z_current

labels_available<-seq(1:K_current)
# K_from labels: number K of non empty clusters
K_from_labels<-K_seq
K_from_labels[1]<-length(labels_available)

# Number of non empty clusters based on allocations in Z
True_K_from_Z<-K_from_labels
True_K_from_Z[1]<-length(table(Z_seq[1,]))

# Start of algorithm
N_current<-get_N(Lab_vec=labels_available, y_array = y, Z_vector = Z_current)
n_current<-get_n(Z_vector = Z_current,Lab_vec=labels_available)
A_current<-get_loglik(K_clusters=K_current, N_array = N_current, beta_value = beta_0)

# We keep track of "True K", that is K without the empty clusters
# equivalent to True_K_from_Z but based on n_k (to double check)
True_K_seq<-K_seq
True_K_seq[1]<-K_current-sum(n_current==0)

# we take track of the sequence of the loglik which is just
# A + log((gamma(3*beta_0)/((gamma(beta_0)^3)))
log_lik_seq<-numeric(S)
log_lik_seq[1]=A_current

####################################################################################

## MCMC ALGORITHM START  

# first step is the initialization
s=1
count_MK<-0
count_GS<-0
count_AE<-0

#couting the time the algorithm takes
start_time <- Sys.time()

####################################################################################

while(s<S){
  chosen_move= sample(3,size=1)
  MK_chosen<- chosen_move==1
  GS_chosen<- chosen_move==2
  AE_chosen<- chosen_move==3
  ##  MK move
  if (MK_chosen){
    count_MK<-count_MK+1
    add_or_remove=sample.int(2, size=1)
    add_attempt= add_or_remove==1
    # we ADD only if we don't exceed K_max
    not_more_than_max= K_current<K_max
    # Add accepted
    add_accepted<- K_current/((N+K_current)*(K_current+1)) >= runif(1)
    # we REMOVE only if not less than 1
    at_least_two= K_current>=2
    # remove if there are empty clusters
    empty_cluster= any(n_current==0)
    
    ## insert attempt 
    if (add_attempt & not_more_than_max & add_accepted){
      # update available labels
      labels_available<-labels_add(LABELS=labels_available)$NEW_LABEL_VECTOR
      # update K --- true update
      K_seq[s+1]<-K_current+1
      K_current=K_current+1
      # compute also
      N_current<-get_N(Lab_vec=labels_available, y_array = y, Z_vector = Z_current)
      n_current<-get_n(Lab_vec = labels_available, Z_vector = Z_current)
    ## delete attempt
    } else if (!add_attempt & at_least_two & empty_cluster){
      lab_to_delete<-which(n_current==0)[1]
      labels_available<-labels_available[-lab_to_delete]
      # update K --- true update
      K_seq[s+1]<-K_current-1
      K_current=K_current-1
      # compute also
      N_current<-get_N(Lab_vec=labels_available, y_array = y, Z_vector = Z_current)
      n_current<-get_n(Lab_vec = labels_available, Z_vector = Z_current)
    } else {
      # labels_available are the same
      K_seq[s+1]<-K_current
      # N and n current stay the same
    }
    # no changes on Z
    Z_seq[s+1,]=Z_current
    # A_current is always the same
    log_lik_seq[s+1]=A_current
    # Remains same as before
    True_K_seq[s+1]<-True_K_seq[s]
    # Other measures to double check True K
    True_K_from_Z[s+1]<-length(table(Z_seq[s+1,]))
    K_from_labels[s+1]<-length(labels_available)
    # new step
    s=s+1
    # print every 10k steps
    multiple_10k<-s%%10000 == 0
    if(multiple_10k){
      print(paste("Iteration number", s, "Time:", Sys.time()))
    }
  }
  ## AE move
  if(AE_chosen){
    count_AE<-count_AE+1
    K_one= K_current==1
    K_between= K_current>1 & K_current<K_max
    K_is_max= K_current==K_max
    ejection_attempt= sample(x=2, size=1) ==1
    ejection_case_2= K_between & ejection_attempt
    absorption= (K_between & !ejection_attempt) || K_is_max
    p_E<-runif(1)
    # ejection move with K = 1
    if(K_one){
      K_prime=2
      Z_prime=Z_current
      # nodes to eject
      nodes_to_eject1=as.logical(rbinom(n=N, size=1, prob=p_E))
      # check if it is eject all nodes or none
      all_TRUEorFALSE= all(nodes_to_eject1) || all(!nodes_to_eject1)
      if(!all_TRUEorFALSE){
        labels_available_prime<-labels_add(LABELS=labels_available)$NEW_LABEL_VECTOR
        lab_j2<-labels_add(LABELS=labels_available)$LABEL_ADDED
        Z_prime[nodes_to_eject1]<-lab_j2
        N_prime=get_N(Lab_vec=labels_available_prime, y_array=y, Z_vector=Z_prime)
        A_prime=get_loglik(K_clusters=K_prime, N_array = N_prime, beta_value = beta_0)
        n_prime=get_n(Lab_vec=labels_available_prime, Z_vector = Z_prime)
        B_prime=get_B(n_vector=n_prime)
        B_current=get_B(n_vector=n_current)
        log_R_eject=(A_prime-A_current+B_prime-B_current)+log(K_current/(K_current+1)/(N+K_current))+log(N+1)
        eject_condition_1= min(log_R_eject,0)>=log(runif(1))
        if(eject_condition_1){
          #update Z
          Z_seq[s+1,]<-Z_prime
          Z_current<-Z_prime
          #update K
          K_seq[s+1]<-K_prime
          K_current<-K_prime
          #update related quantities
          N_current<-N_prime
          n_current<-n_prime
          A_current<-A_prime
          #update loglik sequence
          log_lik_seq[s+1]<-A_current
          #labels
          labels_available<-labels_available_prime
          #update true k
          True_K_seq[s+1]<-K_current-sum(n_current==0)
        }else{
          # ejection rejected
          Z_seq[s+1,]<-Z_current
          K_seq[s+1]<-K_current
          log_lik_seq[s+1]<-A_current
          True_K_seq[s+1]<-True_K_seq[s]
        }
      } else{
        # ejection rejected
        Z_seq[s+1,]<-Z_current
        K_seq[s+1]<-K_current
        log_lik_seq[s+1]<-A_current
        True_K_seq[s+1]<-True_K_seq[s]
      }
    }
    # ejection attempt with K>1
    if(ejection_case_2){
      K_prime= K_current+1
      j1= sample(labels_available, size=1)
      n_j1= n_current[which(labels_available==j1)]
      ejecting_j1_empty= n_j1==0
      # ejecting a non empty cluster
      if(!ejecting_j1_empty){
        nodes_to_eject2=as.logical(rbinom(n=n_j1, size=1, prob=p_E))
        all_FALSE_or_TRUE_2= all(nodes_to_eject2) || all(!nodes_to_eject2)
        if(!all_FALSE_or_TRUE_2){
          labels_available_prime<-labels_add(labels_available)$NEW_LABEL_VECTOR
          j2_added=labels_add(labels_available)$LABEL_ADDED
          Z_prime=Z_current
          Z_prime[Z_prime==j1][nodes_to_eject2]<-j2_added
          # swap if new lab is different
          new_j2= sample(labels_available_prime, size=1)
          diff_new_lab= new_j2 != j2_added
          if(diff_new_lab){
            Z_prime=mapvalues(Z_prime, from=c(new_j2, j2_added),
                              to=c(j2_added,new_j2),warn_missing = FALSE)
          }
          N_prime=get_N(Lab_vec=labels_available_prime, y_array=y, Z_vector=Z_prime)
          A_prime=get_loglik(K_clusters=K_prime, N_array = N_prime, beta_value = beta_0)
          n_prime=get_n(Lab_vec=labels_available_prime, Z_vector = Z_prime)
          B_prime=get_B(n_vector=n_prime)
          B_current=get_B(n_vector=n_current)
          log_prob_ratio=A_prime-A_current+B_prime-B_current+log(K_current/(K_current+1)/(N+K_current))
          log_R_eject=log_prob_ratio+log(n_j1+1)
          eject_condition_2= min(log_R_eject,0)>=log(runif(1))
          ## ejection accepted
          if(eject_condition_2){
            #update Z
            Z_seq[s+1,]<-Z_prime
            Z_current<-Z_prime
            #update K
            K_seq[s+1]<-K_prime
            K_current<-K_prime
            #update related quantities
            N_current<-N_prime
            n_current<-n_prime
            A_current<-A_prime
            #update loglik sequence
            log_lik_seq[s+1]<-A_current
            #labels
            labels_available<-labels_available_prime
            #update true k
            True_K_seq[s+1]<-K_current-sum(n_current==0)
          }else{
            # ejection rejected
            Z_seq[s+1,]<-Z_current
            K_seq[s+1]<-K_current
            log_lik_seq[s+1]<-A_current
            True_K_seq[s+1]<-True_K_seq[s]
          }
        } else{
          Z_seq[s+1,]<-Z_current
          K_seq[s+1]<-K_current
          log_lik_seq[s+1]<-A_current
          True_K_seq[s+1]<-True_K_seq[s]
        }
      } else{
        Z_seq[s+1,]<-Z_current
        K_seq[s+1]<-K_current
        log_lik_seq[s+1]<-A_current
        True_K_seq[s+1]<-True_K_seq[s]
      }
    }
    ## absorption attempt
    if(absorption){
      K_prime= K_current-1
      j1_absorbing=sample(x=labels_available,size=1)
      K_prime_one= K_prime==1
      if(K_prime_one){
        j2_absorbed= labels_available[-which(labels_available==j1_absorbing)]
      } else{
        j2_absorbed=sample(x=labels_available[-which(labels_available==j1_absorbing)],
                           size=1)
      }
      Z_prime=Z_current
      Z_prime[Z_prime==j2_absorbed]=j1_absorbing
      labels_available_prime=labels_available[-which(labels_available==j2_absorbed)]
      j1_prime_position=which(labels_available_prime==j1_absorbing)
      N_prime=get_N(Lab_vec=labels_available_prime, y_array=y, Z_vector=Z_prime)
      A_prime=get_loglik(K_clusters=K_prime, N_array = N_prime, beta_value = beta_0)
      n_prime=get_n(Lab_vec=labels_available_prime, Z_vector = Z_prime)
      B_prime=get_B(n_vector=n_prime)
      B_current=get_B(n_vector=n_current)
      log_prob_ratio=A_prime-A_current+B_prime-B_current+log((K_current*(N+K_current-1))/(K_current-1))
      log_R_abs=log_prob_ratio-log(n_prime[j1_prime_position]+1)
      abs_condition=min(log_R_abs,0)>= log(runif(1))
      if(abs_condition){
        #update Z
        Z_seq[s+1,]=Z_prime
        Z_current=Z_prime
        #update K
        K_seq[s+1]<-K_prime
        K_current=K_prime
        #update labels
        labels_available=labels_available_prime
        #update n
        n_current=n_prime
        #update True value of K
        True_K_seq[s+1]<-K_current-sum(n_current==0)
        #update M
        N_current=N_prime
        #update loglik
        log_lik_seq[s+1]=A_prime
        A_current<-A_prime
      } else{
        #keep Z
        Z_seq[s+1,]=Z_current
        #keep K
        K_seq[s+1]<-K_current
        #keep True K as before
        True_K_seq[s+1]<-True_K_seq[s]
        #keep loglik
        log_lik_seq[s+1]=A_current
      }
    }
    # should already be incremented but double check
    True_K_from_Z[s+1]<-length(table(Z_seq[s+1,]))
    K_from_labels[s+1]<-length(labels_available)
    s=s+1
    multiple_10k<-s%%10000 == 0
    if(multiple_10k){
      print(paste("Iteration number", s, "Time:", Sys.time()))
    }
  }
  ## GS move
  if (GS_chosen) {
    count_GS<-count_GS+1
    # we reallocate if K_current is more than 1
    one_cluster= K_current ==1
    if (!one_cluster) {
      Z_prime=Z_current
      # full sweep
      for(ii in 1:N){
        Z_scanning=Z_prime
        #save current label of z_ii
        k_cur<-Z_current[ii]
        k_cur_position=which(labels_available==k_cur)
        # PROPOSE A Z_PRIME for z_prime
        # new clustering setup Z prime with only one node changed
        # reassign Z_[ii] randomly from a discrete uniform
        k_prime=sample(x=labels_available, size=1)
        k_prime_position= which(labels_available==k_prime)
        Z_scanning[ii]=k_prime
        #compute new N and A
        # K_current unchanged 'cause we are in GS context
        N_prime<-get_N(Lab_vec=labels_available, y_array=y, Z_vector = Z_scanning)
        A_prime<-get_loglik(K_clusters=K_current, N_array = N_prime, beta_value = beta_0)
        log_alpha=A_prime-A_current+log((n_current[k_prime_position]+gamma_0)/
                                          (n_current[k_cur_position]+gamma_0-1))
        #create statements that check conditiond to accept move
        GS_condition= min(log_alpha,0)>=log(runif(1))
        if(GS_condition){
          Z_prime<-Z_scanning
          N_current=N_prime
          n_current=get_n(Lab_vec=labels_available, Z_vector = Z_prime)
          A_current<-A_prime
        }
        #labels_available are the same
        #else, Z_prime[ii] stays equal to Z_current[ii]
      }
      #in any case, at the end of scanning
      Z_seq[s+1,]=Z_prime
      Z_current<-Z_prime
      True_K_seq[s+1]<-K_current-sum(n_current==0)
    } else {
      Z_seq[s+1,]=Z_current
      True_K_seq[s+1]<-True_K_seq[s]
    }
    K_seq[s+1]<-K_current
    log_lik_seq[s+1]=A_current
    #should already be incremented
    True_K_from_Z[s+1]<-length(table(Z_seq[s+1,]))
    K_from_labels[s+1]<-length(labels_available)
    s=s+1
    multiple_10k<-s%%10000 == 0
    if(multiple_10k){
      print(paste("Iteration number", s, "Time:", Sys.time()))
    }
  }
}

end_time <- Sys.time()
end_time - start_time
print(paste("For", S, "iterations"))

####################################################################################

## BURN IN

#check dimensions
dim(Z_seq)
length(K_seq)
length(True_K_seq)
length(log_lik_seq)

## BURN IN
Z_seq_burned<-BURN_BABY_BURN(seq_to_burn =Z_seq,maxS=S,
                             burn_in_level = burn_in_level,
                             is_matrix=TRUE)

K_seq_burned<-BURN_BABY_BURN(seq_to_burn =K_seq,maxS=S,
                             burn_in_level = burn_in_level,
                             is_matrix=FALSE)
True_K_seq_burned<-BURN_BABY_BURN(seq_to_burn =True_K_seq,maxS=S,
                                  burn_in_level = burn_in_level,
                                  is_matrix=FALSE)
K_from_labels_burned<-BURN_BABY_BURN(seq_to_burn =K_from_labels,maxS=S,
                                     burn_in_level = burn_in_level,
                                     is_matrix=FALSE)
True_K_from_Z_burned<-BURN_BABY_BURN(seq_to_burn =True_K_from_Z,maxS=S,
                                     burn_in_level = burn_in_level,
                                     is_matrix=FALSE)

log_lik_seq_burned<-BURN_BABY_BURN(seq_to_burn =log_lik_seq,maxS=S,
                                   burn_in_level = burn_in_level,
                                   is_matrix=FALSE)

#This two should coincide
table(K_seq_burned)
table(K_from_labels_burned)
#This two should coincide
table(True_K_seq_burned)
table(True_K_from_Z_burned)


####################################################################################
## Label correction and analysis

library(collpcm) # used for label switching algorithm
source("LABEL_CORRECTION_AND_ANALYSIS.R")


# VISUALIZE AND SAVE PERMUTED MATCH GRID AFTER ESTIMATING K AND Z
####################################################################################
if (K_estimated>1){
  Winner_label = Team_Names[rownames(Ordered_Tabellone)[1],][1]
  # select cluster percentages
  Cluster_percentages = get(paste0("Cluster_Percentages_Model", 
                                                         K_estimated))
  
  Top_block = as.numeric(which.max(Cluster_percentages[,Winner_label]))
  
  # Team is in topblock if the posterior allocation is >=0.5
  Top_block_Teams = as.numeric(which(Cluster_percentages[Top_block,]>50))
  how_many_top = length(Top_block_Teams)
  
  all = 1:(dim(O)[1])
  the_others <- all[!all %in% Top_block_Teams]
  
  new_block_order = c(Top_block_Teams, the_others)
  Reordered_O = O[new_block_order,]
  Final_O = Reordered_O[,new_block_order]
  
  print(levelplot(t(Final_O[nrow(Final_O):1,]),
                  col.regions=palf(100), xlab = NULL, ylab = NULL, colorkey = FALSE,
                  main = list(paste0("Results table ordered by block membership for season: ",
                                     season), cex= 2), 
                  scales = list(list(alternating=1), x=list(cex=1),y=list(cex=0.8)),
                  panel = function(...){
                    panel.levelplot(...)
                    panel.abline(h = (20-how_many_top)+0.5, lw = 2.5)
                    panel.abline(v = how_many_top+0.5, lw =2.5)
                  }))
  
  # save it in folder
  pdf(paste0("Inference_results//mcmc_Premier_Season_",season,
             "//Heatmap_Estimated_Season_", season,".pdf"),width = 10, height=10)
  
  print(levelplot(t(Final_O[nrow(Final_O):1,]),
                  col.regions=palf(100), xlab = NULL, ylab = NULL, colorkey = FALSE,
                  main = list(paste0("Results table ordered by block membership for season: ",
                                     season_lab), cex= 1.9), 
                  scales = list(list(alternating=1), x=list(cex=1),y=list(cex=0.8)),
                  panel = function(...){
                    panel.levelplot(...)
                    panel.abline(h = (20-how_many_top)+0.5, lw = 2.5)
                    panel.abline(v = how_many_top+0.5, lw =2.5)
                  }))
  dev.off()
}

## Save workspace
save.image(paste0("Inference_results//mcmc_Premier_Season_",season,
                  "//WS_Premier_Season_", season, "_", (S-burn_in_level)/1000,
                  "k_seed_",my_seed,".RData"))
