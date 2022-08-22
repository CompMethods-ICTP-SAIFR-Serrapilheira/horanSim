errorSymGraph <- function(numRow,
                          numCol,
                          fitnessMutant,
                          initMutants) {
  # verifica se os parametros passados para symGraph
  # possuem valores adequados

  messageInit <- ""
  messageLowValue <- ""
  messageFitness <- ""

  # messageInit
  if (initMutants > numCol*numRow) {
    messageInit <- paste("initMutants (=",
                         initMutants,
                         ") has value higher than numRow*numCol (=",
                         numRow,
                         "*",
                         numCol,
                         ")",
                         sep = "")
  }

  #messageLowValue
  if (numRow < 2 || numCol < 2) {
    messageLowValue <- paste("numCol or numRow is less than 2.")
  }

  #messageFitness
  if(fitnessMutant < 0) {
    messageFitness <- paste("fitnessMutant is less than one")
  }


  fullMessage <- paste(messageInit,
                       messageLowValue,
                       messageFitness,
                       collapse = "",
                       sep = "\n")

  return(fullMessage)
}



symGraph <- function(numRow = 5,
                     numCol = 5,
                     fitnessMutant = 1,
                     initMutants = 10) {
  # This function returns the necessary vectors
  # for the Moran iteration.
  # Using the strategy described here,
  # we translate the 2D lattice to a vector,
  # through a neighbor vector and a weighProb vector.

  #verifica se há erros nos valores passados
  errorMessage <- errorSymGraph(numRow = numRow,
                                numCol = numCol,
                                fitnessMutant = fitnessMutant,
                                initMutants = initMutants)
  if(errorMessage != "\n\n") {stop(message)}


  # choose random positions for the mutants
  initPosMutant <- sample.int(n = numRow*numCol,
                              size = initMutants,
                              replace = FALSE)
  #record in popVec the positions of mutants
  popVec <- rep(x = c(FALSE), times = numRow*numCol)
  popVec[initPosMutant] <- TRUE

  # produce the weighted probabilities array
  weighProb <- rep(x = 1/4, times = numRow*numCol)
  weighProb[initPosMutant] <- fitnessMutant/4
  weighProb <- rep(x = weighProb, each = 4)


  #produce the up, down, left and right positions vectors
  totalLen = numRow*numCol
  # UP VECTOR
  upPos <- array(data = c(1), dim = c(numRow, numCol))
  upPos[,1] <- c(numRow, 2:numRow-1)
  for(i in 2:numCol) {upPos[,i] <- upPos[,i-1]+numRow}

  # DOWN VECTOR
  downPos <- array(data = c(1), dim = c(numRow, numCol))
  downPos[,1] <- c(2:numRow, 1)
  for(i in 2:numCol) {downPos[,i] <- downPos[,i-1]+numRow}

  # LEFT VECTOR
  leftPos <- array(data = c(1), dim = c(numRow, numCol))
  leftPos[,1] <- (totalLen-numRow+1):totalLen
  leftPos[,2] <- 1:numRow
  if(numCol>2) {
    for(i in 3:numCol) {leftPos[,i] <- leftPos[,i-1]+numRow}
  }

  # RIGHT VECTOR
  rightPos <- array(data = c(1), dim = c(numRow, numCol))
  rightPos[,numCol] <- 1:numRow
  rightPos[,1] <- rightPos[,numCol] + numRow
  if(numCol>2) {
    for(i in 2:(numCol-1)) {rightPos[,i] <- rightPos[,i-1]+numRow}
  }

  # finally, produce a single vector,
  # from the same size of weighProb,
  # that contain the position of the selected neighbors
  neighVector <- integer(numRow*numCol*4)
  for(i in 0:(numRow*numCol-1)) {
    neighVector[(i*4)+1] <- upPos[i+1]
    neighVector[(i*4)+2] <- downPos[i+1]
    neighVector[(i*4)+3] <- leftPos[i+1]
    neighVector[(i*4)+4] <- rightPos[i+1]
  }

  neighAndWeigh <- list(popVec, neighVector, weighProb)
  return(neighAndWeigh)
}



symTrajec <- function(numRow = 5,
                     numCol = 5,
                     fitnessMutant = 1,
                     initMutants = 10,
                     lowLim = 0,
                     upLim = numRow*numCol) {
  # this function generates the results
  # for a single Moran lattice execution.

  #1. First, generate the initial vectors
  neighAndWeigh <- symGraph(numRow = numRow,
                            numCol = numCol,
                            fitnessMutant = fitnessMutant,
                            initMutants = initMutants)
  popVec <- neighAndWeigh[1][[1]]
  neighVec <- neighAndWeigh[2][[1]]
  weighVec <- neighAndWeigh[3][[1]]

  # initialize variables for the iteration
  quantMut <- initMutants
  # mutantPopVec vai armazenar a quantidade de mutantes
  # na iteracao i
  mutantPopVec <- integer(20000)
  mutantPopVec[1] <- quantMut
  i = 1 # counter
  # replacementType eh uma matriz que codifica
  # o tipo de transicao que esta sendo feito
  # caso 0: FALSE substituindo um FALSE
  # caso 1: FALSE substituindo um TRUE
  # caso 2: TRUE substituindo um FALSE
  # caso 3: TRUE substituindo um TRUE
  replacementType <- matrix(data = c(0,2,1,3),
                            nrow = 2,
                            ncol = 2)

  while(quantMut < upLim && quantMut > lowLim){

    i <- i + 1
    # produce the iteration
    # we use the value of the cell in weighVec to say if
    # the cell is a mutant or a nonMutant

    # indexToReplace holds the index of the spatial position of the cell
    # that will be replaced
    indexReplaced <- sample(x = neighVec, size = 1, prob = weighVec)
    # replacedIndividual shows the cell that will be replaced
    replacedIndividual  <- neighVec[indexReplaced]
    # individualToRep holds the index of the cell that will reproduce
    reprIndividual <- ceiling(indexReplaced/4)

    isReprMutant <- popVec[reprIndividual]
    isReplacedMutant <- popVec[replacedIndividual]

    # identifique o tipo de transicao que esta acontecendo
    # no caso, se esta mudando o valor
    transitionType <- replacementType[isReprMutant+1, isReplacedMutant+1]
    #message(paste("transitionType:",transitionType,sep=""))

    # execute the replacement
    #1. First, check if the replacement produced
    # a increment or a decrement in the number of mutants
    if(transitionType == 1) { #FALSE replacing TRUE
      # decrement mutant number
      quantMut <- quantMut - 1
      # change the value of the replaced cell in popVec
      popVec[replacedIndividual] <- FALSE
      # change the value of the weighted probabilities
      weighVec[(replacedIndividual*4-3):(replacedIndividual*4)] <- rep(x = 1/4,
                                                                       times = 4)
    }
    else if(transitionType == 2) { #TRUE replacing FALSE
      # increment mutant number
      quantMut <- quantMut + 1
      # change the value of the replaced cell in popVec
      popVec[replacedIndividual] <- TRUE
      # change the value of the weighted probabilities
      weighVec[(replacedIndividual*4-3):(replacedIndividual*4)] <- rep(x = fitnessMutant/4,
                                                                       times = 4)
    }

    # record the quantity of mutants in the matrix
    # makes the vector larger only when it fills the vector entirely
    if(i > length(mutantPopVec)) {mutantPopVec <- c(mutantPopVec, integer(20000))}
    mutantPopVec[i] <- quantMut
    if(i%%10000 == 0) {
      message(paste("i = ", i, ", quantMut = ", quantMut, sep=""))
      message(paste("popVec"))
      message(paste(popVec))
      message(paste("weighVec"))
      message(paste(weighVec))
      }
  }

  #shorten the vector to not contain neutral values
  mutantPopVec <- mutantPopVec[1:i]
  return(mutantPopVec)
}



