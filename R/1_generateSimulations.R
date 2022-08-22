# Filename: 1_generateSimulations.R
# Purpose: to generate the simulated datasets through
# the functions present in "./R_modules/wellMixed.R"

# +---------------------+
# | IMPORTING LIBRARIES | ===================================
# +---------------------+

source("./R_modules/wellMixed.R")

# +-----------------------------+
# | DEFINING INITIAL CONDITIONS | ===================================
# +-----------------------------+
seed = 29343064
set.seed(seed = seed)
totalPopVec <- c(75,100,125,150)
fitnessMutantVal <- c(1.01, 1.1, 1.3, 1.5, 1.7)
mutantPopVal <- c(1, 10, 30, 50, 70)
numTrajec <- 1000
outputDir <- "./data/raw"


# iterate through all the fitnessMutantVal values
# and save them in ./data/raw
if(!dir.exists(outputDir)) {dir.create(outputDir, recursive = TRUE)}

lapply(X = fitnessMutantVal,
       FUN = execHoranForFitnessValue,
       dirOutput = outputDir,
       simFunction = oneIterMoran,
       totalPopVec = totalPopVec,
       numTrajec = numTrajec,
       mutantPopVal = mutantPopVal)
