library(rjson)
library(igraph)
library(foreach)
library(doParallel)


#############
# functions #
#############

getLargeAdjMatrix <- function(Lambda, Gamma){
  n_obs <- nrow(Lambda)
  n_lat <- ncol(Gamma)
  n = n_obs+n_lat
  L = matrix(0,n,n)
  L[1:n_obs, 1:n_obs] <- Lambda
  L[(n_obs+1):n, 1:n_obs] <- t(Gamma)
  return(L)
}

getAllGraphs <- function(Gamma, max_obs_edges, acyclic=TRUE, cores=20){
  
  if (!acyclic){
    stop("Cyclic graphs is not implemented yet!")
  }
  
  n = nrow(Gamma)
  if (n==1){
    return(list("0"=list(matrix(0,1,1))))
  }
  
  nOffDiag <- n * (n-1)
  adjMatrix <- matrix(0,n,n)
  storage <- list(adjMatrix)
  
  # Initialize the cluster
  cores = min(cores, max_obs_edges) 
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl)
  
  results <- foreach(m = 1:max_obs_edges, 
                     .combine="c", 
                     .multicombine=TRUE, 
                     .packages=c("igraph"), .export=c("getLargeAdjMatrix")) %dopar%{
    
    print(paste(m, choose(nOffDiag,m)))
    res = NULL
    
    edgeIndices <- t(combn(nOffDiag,m))
    if (is.vector(edgeIndices)){
      edgeIndices <- matrix(edgeIndices, nrow=1)
    }
    
    for (i in 1:nrow(edgeIndices)){
      if (i%%1000==0){print(i)}
      Lambda <- matrix(0,n,n)
      Lambda[upper.tri(adjMatrix)| lower.tri(adjMatrix)][edgeIndices[i,]] <- 1
      if (is_dag(graph_from_adjacency_matrix(Lambda, mode="directed"))){
        largeAdjMatrix = getLargeAdjMatrix(Lambda, Gamma)
        g = graph_from_adjacency_matrix(largeAdjMatrix, mode="directed")
        # Test pairwise isomorphism
        is_new = TRUE
        if (!is.null(res)){
          for (j in 1:length(res)){
            if (isomorphic(res[[j]], g)){
              is_new = FALSE
              break
            }
          }
        }
        if (is_new){
          res = c(res, list(g))
        }
      }
    }
    res
  }
  stopCluster(cl)
  
  # Graph to adjacency matrix
  for (i in 1:length(results)){
    storage[[i+1]] = as_adj(results[[i]], sparse=FALSE)[1:n,1:n]
  }
return(storage)
}

######################################
# Load latent adjacency matrix Gamma #
######################################

exp_name = "DAG6Nodes2Latent/"
prefix = paste(getwd(), "/experiments/", sep="")
filename = paste(prefix, exp_name, "meta_data.json", sep="")
meta_data = fromJSON(file=filename)
Gamma = do.call(rbind, meta_data$Gamma)
max_obs_edges = 6 #meta_data$codim_Omega
cores = 4



############################
# Get all graphs for Gamma #
############################

res <- getAllGraphs(Gamma, max_obs_edges, acyclic=meta_data$acyclic, cores=cores)

# Change format of list
jsonList = list()
for (k in 1:length(res)){
  obj <- list(list("edges" = sum(res[[k]]), "adjMatrix" = c(t(res[[k]]))))
  names(obj) <- k
  jsonList = c(jsonList, obj)
}
print(length(res))

# Save as json
name =  paste(prefix, exp_name, "results.json", sep="")
jsonData = toJSON(jsonList)
write(jsonData, name)
