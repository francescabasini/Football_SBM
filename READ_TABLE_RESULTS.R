####### SCRIPT FOR LOADING AND CLEANING DATA ###  

# Read results table
Results<-read.csv(paste0("Data_Premier//Result_Premier_",season, ".csv"), header = TRUE, sep=",")


rownames(Results)=Results$Home...Away
Results=Results[,-1]


library(stringi)

N<-dim(Results)[1]

O<-matrix(NA,N,N)
for(i in 1:N){
  for(j in 1:N){
    if(i==j){
      O[i,j]=NA
    }else{
      b=as.character(Results[i,j])
      a=stri_split_fixed(b,"~")
      if(as.numeric(a[[1]][1])>as.numeric(a[[1]][2])){
        O[i,j]=1
      }
      if(as.numeric(a[[1]][1])==as.numeric(a[[1]][2])){
        O[i,j]=2
      }
      if(as.numeric(a[[1]][1])<as.numeric(a[[1]][2])){
        O[i,j]=3
      }
    }
  }
}

rownames(O)=rownames(Results)
colnames(O)=colnames(Results)


##### BUILDING THE ADJACENCY MATRIX from OUTCOME MATRIX 
####################################################################################

#Function for the adjacency matrix from O
to_adjacency<-function(O_matrix){
  N_nodes=dim(O_matrix)[1]
  adja_y=array(data=NA, dim=c(N,N,3))
  for(i in 1:N_nodes){
    for(j in 1:N_nodes){
      #do not consider self loops
      if(i!=j){
        #value 1 for W
        if(O_matrix[i,j]==1)
          adja_y[i,j,]=c(1,0,0)
        #Draw
        if(O_matrix[i,j]==2)
          adja_y[i,j,]=c(0,1,0)
        #Loss
        if(O_matrix[i,j]==3)
          adja_y[i,j,]=c(0,0,1)
      }
    }
  }
  return(adja_y)
}

y<-to_adjacency(O)

#number of nodes
N= dim(O)[1]

#do not need to build a generalised adjacency matrix as it would be a waste of
#memory and more complicate to work with
O_value=O
O_value<-ifelse(O==1, "W", ifelse(O==2, "D", "L"))

O_value
