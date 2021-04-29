##### SCRIPT OF ALL SBM FUNCTIONS #######

#Function to have the array N_kl^omega
get_N<-function(Lab_vec, y_array, Z_vector){
  lab_values_ok<-sum(unique(Z_vector) %in% Lab_vec) == length(unique(Z_vector))
  Kclus<-length(Lab_vec)
  if(lab_values_ok){
    N_kl_omega=array(NA,dim=c(Kclus,Kclus,3))
    for(k in 1:Kclus){
      where_k=Z_vector==Lab_vec[k]
      for(l in 1:Kclus){
        where_l=Z_vector==Lab_vec[l]
        for(omega in 1:3){
          N_kl_omega[k,l,omega]=sum(y_array[where_k,where_l,omega], na.rm=TRUE)
        }
      }
    }
    return(N_kl_omega)
  }else{
    message("[Error in get_N]: Some nodes are assigned to a non existing label:")
    which_not<-unique(Z_vector)[!(unique(Z_vector) %in% Lab_vec)]
    return(which_not)
  }
}

# Function for n_k vector 
get_n<-function(Lab_vec, Z_vector){
  repeated_lab= any(diff(Lab_vec)==0)
  nodes_in_lab<- is.element(Z_vector, Lab_vec)
  each_node_labeled<- all(nodes_in_lab)
  if(repeated_lab){
    message("[Error in get_n]: One or more labels in Lab_vec are repeated")
    return(NULL)
  }else if(!each_node_labeled){
    message("[Error in get_n]: One or more nodes are assigned to a non existing label:")
    nodes_unassigned<-which(nodes_in_lab==FALSE)
    return(list(NODE_POSITION=nodes_unassigned, NODE_ASSIGMENT=Z_vector[nodes_unassigned]))
  } else {
    K_clusters<-length(Lab_vec)
    n<-numeric(K_clusters)
    for(i in 1:K_clusters){
      n[i]<-sum(Z_vector==Lab_vec[i])
    }
    return(n)
  }
}

# Function to compute the log of approximate likelihood (integrated out)
get_loglik<-function(N_array, beta_value, K_clusters){
  A<-sum(lgamma(N_array+beta_value))-sum(lgamma(rowSums(N_array, dims=2)+3*beta_value))
  C<-K_clusters^2*lgamma(3*beta_value)-K_clusters^2*log(gamma(beta_value)^3)
  log_lik<-A+C
  return(log_lik)
}


# compute sum of logGamma (nk +1) over all k
get_B<-function(n_vector){
  sum(lgamma(n_vector+1))
}

#Function to have the vector of labels available to assign
labels_add<-function(LABELS){
  CLUSTERS<-length(LABELS)
  #first scenario is the first is missing
  first_one_missing<- LABELS[1]!=1
  if(first_one_missing){
    #scenario where the first is missing
    NEW_LABELS_first_missing<-c(1,LABELS)
    return(list(NEW_LABEL_VECTOR=NEW_LABELS_first_missing, LABEL_ADDED=1))
  }else{
    diff_lab<-diff(LABELS)
    all_ones<- all(diff_lab==1)
    if(all_ones){
      #I need to add the last label
      NEW_LABELS_last_missing<-c(LABELS, CLUSTERS+1)
      return(list(NEW_LABEL_VECTOR=NEW_LABELS_last_missing, LABEL_ADDED=CLUSTERS+1))
    }else{
      label_to_insert<-which(diff_lab!=1)[1]+1
      NEW_LABELS_middle<-append(x=LABELS,values=label_to_insert, after=label_to_insert-1 )
      return(list(NEW_LABEL_VECTOR=NEW_LABELS_middle, LABEL_ADDED=label_to_insert))
    }
  }
}

# function to compute burn in of chain results
BURN_BABY_BURN<-function(seq_to_burn, maxS, burn_in_level, is_matrix=FALSE){
  if(is_matrix){
    BURNED_SEQ<-seq_to_burn[(burn_in_level+1):maxS,]
  }else{
    BURNED_SEQ<-seq_to_burn[(burn_in_level+1):maxS]
  }
  return(BURNED_SEQ)
}
  