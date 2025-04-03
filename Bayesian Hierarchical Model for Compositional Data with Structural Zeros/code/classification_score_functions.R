CohenKappa <- function(tab)
{
  # Cohen's kappa, see Agresti, A. (1990) Categorical Data Analysis, p. 366
  subtab <- substitute(tab)
  if(!is.matrix(tab)) stop(paste(subtab, "is not a square matrix"))
  d <- dim(tab)
  if(d[1] != d[2]) stop(paste(subtab, "is not a square matrix"))
  tot <- sum(tab)
  if(tot == 0) stop(paste("all cells in", subtab, "are 0"))
  prob <- tab / tot
  rowprob <- rowSums(prob)
  colprob <- colSums(prob)
  Pi0 <- sum(diag(prob))
  Pie <- sum(rowprob * colprob)
  kappa <- (Pi0 - Pie) / (1 - Pie)
  return(kappa)
}

GKtau <- function(tab)
{
  # Goodman and Kruskal's tau,
  #   see Agresti, A. (1990) Categorical Data Analysis, p. 24
  subtab <- substitute(tab)
  if(!is.matrix(tab)) stop(paste(subtab, "is not a matrix"))
  tot <- sum(tab)
  if(tot == 0) stop(paste("all cells in", subtab, "are 0"))
  prob <- tab / tot
  rowprob <- rowSums(prob)
  colprob <- colSums(prob)
  if(sum(colprob > 0) <= 1)
    stop(paste("Less than 2 columns in", subtab, "have positive counts"))
  A <- sum( rowSums(prob*prob) / rowprob )
  B <- sum(colprob * colprob)
  tau <- (A - B) / (1 - B)
  return(tau)
}

TheilU <- function(tab)
{
  # Theil's U, see Agresti, A. (1990) Categorical Data Analysis, p. 25
  subtab <- substitute(tab)
  if(!is.matrix(tab)) stop(paste(subtab, "is not a matrix"))
  tot <- sum(tab)
  if(tot == 0) stop(paste("all cells in", subtab, "are 0"))
  prob <- tab / tot
  rowprob <- rowSums(prob)
  colprob <- colSums(prob)
  if(sum(colprob > 0) <= 1)
    stop(paste("Less than 2 columns in", subtab, "have positive counts"))
  NumMat <- prob * log(prob / outer(rowprob, colprob))
  NumMat[is.nan(NumMat)] <- 0
  Num <- sum(NumMat)
  DenVec <- colprob * log(colprob)
  DenVec[is.nan(DenVec)] <- 0
  Den <- sum(DenVec)
  U <- - Num / Den
  return(U)
}
