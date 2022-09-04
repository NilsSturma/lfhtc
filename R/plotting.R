library(rjson)
library(SEMID)
source(paste(getwd(), "/R/utils.R", sep=""))

# Set location of experiment
exp_name = "DAG6Nodes2Latent/"
prefix = paste(getwd(), "/experiments/", sep="")

# Read metadata
meta_data = fromJSON(file=paste(prefix, exp_name, "meta_data.json", sep=""))
Gamma = do.call(rbind, meta_data$Gamma)

# Read graphs
graphs = fromJSON(file=paste(prefix, exp_name, "results.json", sep=""))

# Define function
plotGraph <- function(g, Gamma){
  print(paste("Edges: ", g$edges, sep=""))
  print(paste("Identifiable: ", g$identifiable, sep=""))
  print(paste("LF-HTC identifiable: ", g$lfhtcID, sep=""))
  
  LambdaVec <- g$adjMatrix
  Lambda <- matrix(LambdaVec, nrow=sqrt(length(LambdaVec)), byrow=TRUE)
  g <- getLatentDigraph(Lambda, Gamma)
  plot(g)
}

# Plot specific graph
nr = 5
plotGraph(graphs[[as.character(nr)]], Gamma)
