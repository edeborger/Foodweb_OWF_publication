### Plot ranges wrapper

# To debug
# range <- Flows.330[284:317, ]

PRanges <- function(range, mtitle, labcex = 0.5, textcex = 0.6) {
  
  ord   <- order(range[,2], decreasing = F)
  xlabz <- expression(paste("Flow value ") *
                         group("(", "mmol C " * m^{-2} * d^{-1}, ")"))
  l    <- length(ord)
  rmin <- range[, 2] - range[, 3]
  rmin <- min(rmin[rmin >= 0])
  rmax <- max(range[, 2] + range[, 3])
  unks <- range$flow
  
  if((l %% 2) == 0) {
    even   <- seq(2, l  , 2)
    uneven <- seq(1, l-1, 2)
  } else {
    even   <- seq(2, l-1, 2)
    uneven <- seq(1, l  , 2)    
  }
  
  par(mar=c(4, 4, 4, 6))
  
  Plotranges(min     = range[, 2][ord] - range[, 3][ord], 
             max     = range[, 2][ord] + range[, 3][ord], 
             value   = range[, 2][ord]     ,
             log     = "x"                 ,
             main    = mtitle              ,
             xlab    = xlabz               ,
             labels  = rep("", length(ord)),
             xlim    = c(rmin, rmax)       ,
             lab.cex = labcex,
             pch     = 20    ,
             pch.col = c("black", "grey"))
  mtext((unks[ord][even]  ), 
        side = 4, line = 0.5, at = even  , adj = 0, las = 2, cex = textcex)
  mtext((unks[ord][uneven]), 
        side = 2, line = 2.5, at = uneven, adj = 0, las = 2, cex = textcex)
}



### Extract from multicore runs
# Extracts x-values from multi-core runs.
# It drops the first result from each run since this is always the same initial 
# solution.

mcExtract <- function(lst) {
  
  x <- NULL
  
    for(i in (1:length(lst))){
      newx <- lst[[i]]$X[-1, ]
      x    <- rbind(x, newx)
    }
  
  return(x)
}


# getFlowMatrix

# @ DANIELLE DE JONGE
# This function creates a flowmatrix with flow values.
# It produces the same as Flowmatrix() - which is a build in function -, but
# takes into the account the possibility that multiple flows may occur between compartments.
# For every flow it checks whether there are other flows with the same source and sink
# compartments. The values of all flows with the same source and sink are summed to
# obtain the total flow between those two compartments. The total flow value is stored
# at the location of the original flownumber. So if flow1 and flow2 both occur between
# compartment A and B, the total flow value at index 1 is flow1 + flow2 and the total flow value
# at index 2 is also flow1 + flow2.
# It then replaces all flownumbers in the original flowmatrix with the corresponding total flow
# value. So, if flow1 is not included in the flowmatrix, but flow2 is, then the total flux value at
# index 2 is taken, which already also includes the flow value of flow1.

getFlowMatrix <- function (originalFM = NULL, flows = NULL, answers = NULL) {
  # Check function input
  if(is.null(originalFM)) {stop("Supply flowmatrix with flownumbers from LIM$Flowmatrix")}
  if(is.null(flows)) {stop("Supply the flows dataframe from readLIM$flows")}
  if(is.null(answers)) {stop("Supply a vector with all flow values")}
  
  # Create a vector which will store all total flow values
  totalflows <- numeric(length(flows[,1]))
  # Create a vector with all flownumbers
  flownrs <- c(1:length(flows[,1]))
  
  for(flownr in flownrs) {
    same <- findParallelFlows(flows = flows, flownr = flownr)
    totalflow <- sum(answers[same], na.rm = TRUE) # sum the flow values of all flows with the same source and sink.
    totalflows[flownr] <- totalflow # store the value at the index of the original flow
  }
  
  # Find the locations in the matrix with flownumbers
  indices <- match(flownrs, originalFM)
  indices <- indices[!is.na(indices)] # remove NA's (no matches) from the list.
  # Replace the flownumbers in the original matrix with the corresponding
  # total flow values (so if the flownr is 3 in the matrix, but it should be flow
  # 2 and 3, it will insert the totalflow of flow 2 and 3).
  FM <- replace(originalFM, indices, totalflows[originalFM[indices]])
  # Return a flowmatrix with flow values.
  return(FM)
}


# findParallelFlows

# @DANIELLE DE JONGE
# Input flows = readLIM$flows and a flow number.
# Returns vector with row numbers in readLIM$flows of all flows
# with the same source and sink as the flow defined by flownr
# (so all parallel flows).
findParallelFlows <- function(flows = NULL, flownr = NULL) {
  sc1 <- which(flows[,1] == flows[flownr, 1]) # flows with same source compartment
  sc2 <- which(flows[,2] == flows[flownr, 2]) # flows with same sink compartment
  sc3 <- c(sc1, sc2) # all flows ...
  same <- sc3[duplicated(sc3)] # if a flownr occurs twice it has the same source and sink
  return(same)
}


# stdfromstd()

# @ DANIELLE DE JONGE
# Function that calculates the new standard deviation when means
# with their own standard deviations are summed. To sum standard
# deviations one uses the following formula:
# Square root(std1^2 + std2^2 + ... + stdn^2)
stdfromstd <- function(stdevs = NULL){
  vars <- stdevs^2
  sum <- sum(vars, na.rm = TRUE)
  std <- sqrt(sum)
  return(std)
}


# GetTrophInd

# @ DANIELLE DE JONGE
# is a copy of the TrophInd function from the NetInd package. In the original 
# function the TL is the average of the trophic levels of the input compartments
# plus 1. We want the carrion compartment to have a trophic level equal to the average
# trophic levels of all input compartments, but without the addition of 1, because
# the carrion pool does not really 'feed' on anything. We also don't just want to
# set it to trophic level 1.
getTrophInd <- function (Flow = NULL, Tij = t(Flow), Import = NULL, Export = NULL, 
                         Dead = NULL) 
{
  if (is.character(Dead)) {
    dead <- which(rownames(Tij) %in% Dead)
  }else {dead <- Dead}
  if (length(dead) != length(Dead)) {
    stop("Dead not recognized")
  }
  N <- NetIndices:::InternalNetwork(Tij, Import, Export)
  p <- NetIndices:::Diet(N$Tint, dead, N$iN)
  ncomp <- ncol(N$Tint)
  A <- -p
  diag(A) <- 1
  B <- rep(1, ncomp)
  B[9] <- 0 # Carcass has index 9, should not add 1 but add 0.
  TL <- ginv(A) %*% B
  OI <- vector(length = ncomp)
  for (i in 1:ncomp) OI[i] <- sum((TL - (TL[i] - 1))^2 * p[i, 
  ])
  return(data.frame(TL, OI, row.names = rownames(N$Tint)))
}


