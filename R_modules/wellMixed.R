# Filename: wellMixed.R
# Purpose: to hold functions related to Moran process modeling
# in a well-mixed space

# +---------------------+
# | IMPORTING LIBRARIES | ===================================
# +---------------------+

if(!(any(rownames(installed.packages())=="stringr"))) {
  install.packages("stringr")
}
library(stringr)

# +--------------------+
# | DEFINING FUNCTIONS | ===================================
# +--------------------+

probMoran <- function(totalPop, mutantPop, fitnessMutant) {
  # Calculate the probability of increment, decrement or staying
  # for random drift with constant selection.
  # The formula considers that the fitness difference modifies
  # the probability of reproduction, and not death.
  # Assumes that the fitness of the non-mutant is 1,
  # and the fitness of mutant is fitnessMutant.

  # probability that the mutant is chosen to reproduce
  # and that the non-mutant is chosen to die
  probIncrement <- ((fitnessMutant*mutantPop) / (fitnessMutant*mutantPop + totalPop - mutantPop)) * ((totalPop - mutantPop)/(totalPop))

  # probability that the non-mutant is chosen to reproduce
  # and that the mutant is chosen to die
  probDecrement <- ((totalPop -  mutantPop) / (fitnessMutant*mutantPop + totalPop - mutantPop)) * ((mutantPop)/(totalPop))

  # probability of staying in the same state
  # (i.e., that a same species -  mutant or non-mutant is chosen to reproduce
  # and the same species is chosen to die)
  probStaying <- 1 - probIncrement - probDecrement

  probList <- c(probIncrement, probDecrement, probStaying)

  return(probList)
}

#testing one function at a time
#everything seems right
#probMoran(10, 5, 1.1)
#probMoran(10, 10, 1.1)
#probMoran(10, 0, 1.1)


oneIterMoran <- function(totalPop,
                         mutantPop,
                         fitnessMutant) {
  # Calculates the next iteration of population.
  # The iteration assumes that the total population is constant
  # and for each iteration,
  # one random individual is replaced by another random individual,
  # given their abundances and their fitness.

  #reestimate probabilities for each iteration
  probList <- probMoran(totalPop, mutantPop, fitnessMutant)

  #probIncrement <- probList[1]
  #probDecrement <- probList[2]

  randomNumber <- runif(1)
  if(randomNumber < probList[1]) {mutantPop <- mutantPop + 1}
  else if(randomNumber < probList[1] + probList[2]) {mutantPop <- mutantPop - 1}
  # else is the case of not incrementing nothing,
  # so no alteration in the population number is necessary

  iterMoran <- c(totalPop, mutantPop, fitnessMutant)
  return(iterMoran)
}

# testing one function at a time
# everything seems fine
#oneIterMoran(totalPop = 10, mutantPop = 5, fitnessMutant = 1.1)


oneIterMoranWithoutStay <- function(totalPop,
                                    mutantPop,
                                    fitnessMutant) {
  # Similar to oneIterMoran,
  # but it ignores the probability of staying in the same state.
  # Inappropriate to answer questions related to time until fixation

  #reestimate probabilities for each iteration
  probList <- probMoran(totalPop, mutantPop, fitnessMutant)

  #probIncrement <- probList[1]
  #probDecrement <- probList[2]

  randomNumber <- runif(1)
  #message(paste("randomNumber = ", randomNumber))
  #message(paste("probList[1] = ", probList[1]))
  #message(paste("probList[2] = ", probList[2]))
  #message(paste("probList[3] = ", probList[3]))

  if(randomNumber < (probList[1] / (probList[1] + probList[2]))) {
    mutantPop <- mutantPop + 1
  }
  else {mutantPop <- mutantPop - 1}

  iterMoran <- c(totalPop, mutantPop, fitnessMutant)
  return(iterMoran)
}

# testing one function at a time
#ola <- oneIterMoranWithoutStay(totalPop = 10, mutantPop = 5, fitnessMutant = 1.1)





genDriftTrajectory <- function(simFunction = oneIterMoranWithoutStay,
                               totalPop = 100,
                               mutantPop = 1,
                               minLimit = 0,
                               maxLimit = totalPop,
                               fitnessMutant) {
  # generate a random drift trajectory
  # until the value of mutantPop is zero or totalPop.


  # appending numbers is very slow in R
  # to avoid this, it is better to allocate a temporary vector
  # and then, when it is completely filled,
  # transfer its values to the permanent value
  mutantPopVec <- integer(100)
  i = 1 # counter

  # save initial value in the first position of the iteration
  mutantPopVec[i] <- mutantPop

  while(mutantPop != minLimit && mutantPop != maxLimit) {

    i <- i + 1

    #message(paste("ANT:totalPop = ", totalPop))
    #message(paste("ANT:mutantPop = ", mutantPop))
    #message(paste("ANT:fitnessMutant = ", fitnessMutant))

    iterMoran <- simFunction(totalPop, mutantPop, fitnessMutant)
    totalPop <- iterMoran[1]
    mutantPop <- iterMoran[2]
    fitnessMutant <- iterMoran[3]

    #message(paste("DEP:totalPop = ", totalPop))
    #message(paste("DEP:mutantPop = ", mutantPop))
    #message(paste("DEP:fitnessMutant = ", fitnessMutant))

    if(i > length(mutantPopVec)) {
      mutantPopVec <- c(mutantPopVec, integer(100))
    }

    mutantPopVec[i] <- mutantPop

    #increment variable
  }

  #shorten the vector to not contain neutral values
  mutantPopVec <- mutantPopVec[1:i]
  return(mutantPopVec)
}

repeatTrajec <- function(simFunction = oneIterMoranWithoutStay,
                         totalPop,
                         mutantPop,
                         fitnessMutant,
                         numRepet) {
  # execute the function genDriftTrajectory
  # numRepet times,
  # saving in distinct vectors
  # (1) the iteration in which each simulation finished, and
  # (2) if the iteration ended with the number of mutants at maximum or zero.


  endIter <- integer(numRepet)
  upOrDown <- logical(numRepet)

  for(i in 1:numRepet) {
    trajec <- genDriftTrajectory(simFunction = oneIterMoranWithoutStay,
                                 totalPop = totalPop,
                                 mutantPop = mutantPop,
                                 fitnessMutant = fitnessMutant)
    trajecLen <- length(trajec)

    endIter[i] <- trajecLen
    #message(paste("endIter[i]=", endIter[i], sep=""))
    if(trajec[trajecLen] == totalPop) {upOrDown[i] <- TRUE}
    #message(paste("upOrDown[i]=", upOrDown[i], sep=""))
  }

  endIterAndUp <- list(endIter, upOrDown)

  return(endIterAndUp)
}


nameForInitMutantPopCol <- function(mutantPop) {
  # generate automatically a name for
  # the column in the dataframes "endIter"  and "upOrDown"

  mutantPopStr <- as.character(mutantPop)
  colName <- paste("initMutPop_", mutantPopStr, sep="")
  return(colName)
}


nameForSimOutput <- function(prefixColumn,
                             fitnessMutant,
                             totalPop) {
  # return the names for output files,
  # following the layout
  # <prefixColumn><fitnessMutantVal without "." char>.csv

  fitnessMutantStr <- as.character(fitnessMutant)
  procFitnessMutantStr <- stringr::str_replace(fitnessMutantStr,
                                               "\\.",
                                               "_")
  fileName <- paste(prefixColumn,
                    procFitnessMutantStr,
                    "_totalPop_",
                    totalPop,
                    ".csv",
                    sep = "")
  return(fileName)
}


simForMutantPopVal <- function(simFunction = oneIterMoranWithoutStay,
                               totalPop,
                               fitnessMutant,
                               numTrajec,
                               mutantPopVal) {
  # generate the dataframes in which each column
  # contains the results of all simulations.

  endIter <- data.frame(obs_ID = 1:numTrajec)
  upOrDown <- data.frame(obs_ID = 1:numTrajec)

  #message(paste(mutantPopVal, sep = ""))

  for(j in 1:length(mutantPopVal)) {

    #message(paste("mutantPopVal[j] = ", mutantPopVal[j], sep = ""))

    # generate the trajectories
    # repeatTrajec returns two vectors:
    # (1) the iteration in which each simulation finished, and
    # (2) if the iteration ended with
    # the number of mutants at maximum or zero.
    endIterAndUp <- repeatTrajec(simFunction = oneIterMoranWithoutStay,
                                 totalPop = totalPop,
                                 mutantPop = mutantPopVal[j],
                                 fitnessMutant = fitnessMutant,
                                 numRepet = numTrajec)

    #message(paste(endIterAndUp[1][[1]]))
    endIterVec <- endIterAndUp[1][[1]]
    upOrDownVec <- endIterAndUp[2][[1]]

    # generate automatic name for the dataframes' column
    # for a specific value of mutantPopVal
    colName <- nameForInitMutantPopCol(mutantPop = mutantPopVal[j])
    #message(paste("colName=", colName, sep=""))

    # record the information of the trajectories
    # in the dataframes
    endIter[colName] <- endIterVec
    upOrDown[colName] <- as.integer(upOrDownVec)
  }

  #message(paste("colnames=", colnames(endIter), sep=" "))

  resultDataframes <- list(endIter, upOrDown)
  #message(paste(resultDataframes[1][[1]], sep=" "))
  return(resultDataframes)
}




execHoranForFitnessValue <- function(fitnessMutant,
                                     dirOutput,
                                     simFunction,
                                     totalPopVec,
                                     numTrajec,
                                     mutantPopVal) {
  # Execute and save the simulation of a well-mixed Moran model
  # for a single value of mutant relative fitness.
  # It assumes that the absolute fitness of the nonMutant is 1,
  # that the absolute fitness of the mutant is fitnessMutant,
  # and that the relative fitness is fitnessMutant/1.
  # The functions accepted for FUN are
  # "oneIterMoranWithoutStay" and "oneIterMoran".
  # The variables "totalPop" and "numTrajec" are
  # numbers higher than 0.
  # "mutantPopVal" is the vector of values of
  # initial mutant populations that will be simulated for
  # a value of fitness.
  # The returned files are:
  # (1) "endIter_<fitnessMutant>.csv",
  # which records the number of iterations until fixation
  # for each one of the <numTrajec> simulations;
  # and (2) "upOrDown_<fitnessMutant>.csv"
  # which records the end result of the iterations
  # ("0" when the mutantPop is completely eliminated,
  # "1" when the nonMutant pop is completely eliminated).
  # Each column in these files refer to simulations
  # executed with a specific value of
  # initial mutant population quantity
  # (i.e., a element of the vector "mutantPopVal").


  for(i in totalPopVec) {
    resultDataframes <- simForMutantPopVal(simFunction = oneIterMoranWithoutStay,
                                           totalPop = i,
                                           fitnessMutant = fitnessMutant,
                                           numTrajec = numTrajec,
                                           mutantPopVal = mutantPopVal)

    endIter <- resultDataframes[1][[1]]
    upOrDown <- resultDataframes[2][[1]]

    filenameEndIter <- nameForSimOutput("endIter_",
                                        fitnessMutant,
                                        totalPop = i)
    pathEndIter <- file.path(dirOutput, filenameEndIter)

    filenameUpOrDown <- nameForSimOutput("upOrDown_",
                                         fitnessMutant,
                                         totalPop = i)
    pathUpOrDown <- file.path(dirOutput, filenameUpOrDown)

    # generate a file for each value of fitnessMutantVal
    write.csv(x = endIter,
              file = pathEndIter,
              row.names = FALSE,
              quote = FALSE)
    write.csv(x = upOrDown,
              file = pathUpOrDown,
              row.names = FALSE,
              quote = FALSE)

    # output message to confirm the success to the user
    message(paste(filenameEndIter,
                  " and ",
                  filenameUpOrDown,
                  " were written in ",
                  dirOutput,
                  sep = ""))
  }
}
