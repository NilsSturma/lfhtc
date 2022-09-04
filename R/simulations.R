#library(devtools)
#install_github("NilsSturma/SEMID", ref = "lfhtc")

library(rjson)
library(SEMID)
library(foreach)
library(doParallel)
source(paste(getwd(), "/R/utils.R", sep=""))


#############
# Read data #
#############

# Set location of experiment
exp_name = "DAG6Nodes2Latent/"
prefix = paste(getwd(), "/experiments/", sep="")


# Read metadata
meta_data = fromJSON(file=paste(prefix, exp_name, "meta_data.json", sep=""))
Gamma = do.call(rbind, meta_data$Gamma)
max_obs_edges = meta_data$codim_Omega
n_obs = nrow(Gamma)
n_lat = ncol(Gamma)
n = n_obs + n_lat
print(Gamma)

# Read graphs
graphs = fromJSON(file=paste(prefix, exp_name, "results.json", sep=""))
print(length(graphs))

# Check identifiability
for (i in 1:length(graphs)){
  if ((i%%100)==0){print(i)}
  LambdaVec <- graphs[[as.character(i)]]$adjMatrix
  Lambda <- matrix(LambdaVec, nrow=sqrt(length(LambdaVec)), byrow=TRUE)
  g <- getLatentDigraph(Lambda, Gamma)
  graphs[[as.character(i)]]$maxObsParents <- max(sapply(lapply(g$observedNodes(), g$observedParents), length))
  graphs[[as.character(i)]]$lfhtcID <- isIdentifiable(g)
  
  gMixed <- g$getMixedGraph()
  resHTC <- htcID(gMixed, tianDecompose = FALSE)
  graphs[[as.character(i)]]$htcID <- (sum(sapply(resHTC$unsolvedParents, length)) == 0)
  
  resEdgewiseHTC <- generalGenericID(gMixed, list(edgewiseIdentifyStep), tianDecompose = FALSE)
  graphs[[as.character(i)]]$edgewiseID <- (sum(sapply(resEdgewiseHTC$unsolvedParents, length)) == 0)
}


# How many are identifiable by LF-HTC?
countLFHTC = 0
countHTC = 0
countEdgewise = 0
for (i in 1:length(graphs)){
  if (graphs[[as.character(i)]]$lfhtcID==TRUE){
    countLFHTC = countLFHTC + 1
  }
  if (graphs[[as.character(i)]]$htcID==TRUE){
    countHTC = countHTC + 1
  }
  if (graphs[[as.character(i)]]$edgewiseID==TRUE){
    countEdgewise = countEdgewise + 1
  }
}
print(countLFHTC)
print(countHTC)
print(countEdgewise)

# Save results
jsonData = toJSON(graphs)
write(jsonData, paste(prefix, exp_name, "results.json", sep=""))
