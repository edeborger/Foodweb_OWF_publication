# Solving the Coarse foodweb.
# Emil De Borger
# July- 2023

# Loading required packages
require(LIM)
require(splus2R)
require(NetIndices)
require(MASS)
require(parallel)
require(doParallel)

# Load input file
file_coarse <- "./Foodweb_models/final_versions/coarse.input"

# Reading LIM
readLIM_coarse <- LIM::Read(file = file_coarse, checkLinear = TRUE)

#Make initial solution and ranges
LIM_coarse <- Setup(readLIM_coarse)
Solldei_coarse <- Ldei(LIM_coarse, tol = 1e-6)
Ranges_coarse <- Xranges(LIM_coarse, central = T, ispos = T)

# Run looped batches
numCores <- 90

doParallel::registerDoParallel(numCores)

mcmccoarse <- foreach(i = 1:50) %dopar% {
  
    xsnew <- xsample(E    = LIM_coarse$A,
                     F    = LIM_coarse$B,
                     G    = LIM_coarse$G,
                     H    = LIM_coarse$H,
                     jmp  = (Ranges_coarse[,2] - Ranges_coarse[,1])/100,
                     x0   = Solldei_coarse$X,
                     iter = 200)
}
    
doParallel::stopImplicitCluster()


save(mcmccoarse     , file = "./Foodweb_models/final_versions/MCMCcoarse.rda")