library(SEMID)

getLatentDigraph <- function(Lambda, Gamma){
  n_obs <- nrow(Lambda)
  n_lat <- ncol(Gamma)
  n = n_obs+n_lat
  L = matrix(0,n,n)
  L[1:n_obs, 1:n_obs] <- Lambda
  L[(n_obs+1):n, 1:n_obs] <- t(Gamma)
  g <- LatentDigraph(L, 1:n_obs, (n_obs+1):(n_obs+n_lat))
  return(g)
}

isIdentifiable <- function(graph){
  res = lfhtcID(graph)
  if (sum(sapply(res$unsolvedParents, length)) == 0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
