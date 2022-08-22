# filename: 2_processRawData.R
# Purpose: generate the tables for plotting from
# the datasets contained in "./output/raw"


# +--------------------+
# | IMPORTING PACKAGES | ===================================
# +--------------------+


# +--------------------------+
# | ABSORBING endIter TABLES | ===================================
# +--------------------------+

# endIter_1_01_initMutPop_100.csv
endIterFiles <- list.files(path = "data/raw",
                           pattern = "endIter_.+\\.csv$",
                           full.names = TRUE)

# endIter_1_01_initMutPop_100
endIterFilesWithoutExt <- gsub(".csv",
                               "",
                               basename(endIterFiles),
                               fixed = TRUE)

# _1_01_initMutPop_100
endIterValsWithoutPrefix <- gsub("endIter",
                    "",
                    endIterFilesWithoutExt,
                    fixed = TRUE)

# _1_01
endIterFitnessValsNoTotalPopInfo <- sub("_totalPop_[[:digit:]]+",
                                 "",
                                 endIterValsWithoutPrefix)
# 1_01
endIterFitnessValsNoUnder <- sub("_",
                                 "",
                                 endIterFitnessValsNoTotalPopInfo,
                                 fixed = TRUE)
# 1.01 (numeric)
endIterFitnessValsDec <- as.numeric(gsub("_",
                                         ".",
                                         endIterFitnessValsNoUnder,
                                         fixed = TRUE))

# _100 (string)
endIterTotalPop <- sub(".+_totalPop",
                      "",
                      endIterValsWithoutPrefix)

endIterTotalPopDec <- as.integer(sub("_",
                                     "",
                                     endIterTotalPop))


endIterDataList <- lapply(endIterFiles, read.csv)

# names of variables in R cannot begin with numbers
# we'll have to process this values before plotting
names(endIterDataList) <- paste(endIterFitnessValsNoTotalPopInfo,
                                endIterTotalPop,
                                sep = "")


# +----------------------------+
# | ABSORBING upOrDown TABLES | ===================================
# +----------------------------+

# upOrDown_1_01_initMutPop_100.csv
upOrDownFiles <- list.files(path = "data/raw",
                           pattern = "upOrDown_.+\\.csv$",
                           full.names = TRUE)

# upOrDown_1_01_initMutPop_100
upOrDownFilesWithoutExt <- gsub(".csv",
                               "",
                               basename(upOrDownFiles),
                               fixed = TRUE)

# _1_01_initMutPop_100
upOrDownValsWithoutPrefix <- gsub("upOrDown",
                                  "",
                                  upOrDownFilesWithoutExt,
                                  fixed = TRUE)

# _1_01
upOrDownFitnessValsNoTotalPopInfo <- sub("_totalPop_[[:digit:]]+",
                                         "",
                                         upOrDownValsWithoutPrefix)
# 1_01
upOrDownFitnessValsNoUnder <- sub("_",
                                 "",
                                 upOrDownFitnessValsNoTotalPopInfo,
                                 fixed = TRUE)
# 1.01 (numeric)
upOrDownFitnessValsDec <- as.numeric(gsub("_",
                                         ".",
                                         upOrDownFitnessValsNoUnder,
                                         fixed = TRUE))

# _100 (string)
upOrDownTotalPop <- sub(".+_totalPop",
                       "",
                       upOrDownValsWithoutPrefix)

upOrDownDataList <- lapply(upOrDownFiles, read.csv)

# names of variables in R cannot begin with numbers
# we'll have to process this values before plotting
names(upOrDownDataList) <- paste(upOrDownFitnessValsNoTotalPopInfo,
                                upOrDownTotalPop,
                                sep = "")


# +--------------------------------+
# | CHECKING IF EVERY endIter FILE | ======================
# | HAS ITS upOrDown PAIR          |
# +--------------------------------+

#First, check first the fitnessVals
if(sum(endIterValsWithoutPrefix == upOrDownValsWithoutPrefix) != length(upOrDownValsWithoutPrefix)) {
  message(paste("Not every endIter file in \\'./data\\'",
                "has a upOrDown pair file. ",
                "Probably an error will be generated."))
}

# +-------------------------------------------+
# | TABLE FOR GRAPH1: PROBABILITY OF FIXATION | ======================
# |      (initial mutant population = 1)      |
# +-------------------------------------------+

# obtain the possible values for fitnessMutant
# and totalPop to build the dataframe.
# PROBLEM: I don't know if unique()
# will order the vectors in the same way,
# but it seems that I'm right
testedFitnessVals <- c(unique(upOrDownFitnessValsDec))
testedFitnessValsStr <- c(unique(upOrDownFitnessValsNoTotalPopInfo))
testedTotalPopVals <- c(unique(upOrDownTotalPop))


# build data.frame to receive the values
probInitPop1 <- data.frame(fitnessValue = testedFitnessVals)
# name the rows using the same
rownames(probInitPop1) <- testedFitnessValsStr
for(i in testedTotalPopVals) {
  probInitPop1[i] <- numeric(length = length(testedFitnessVals))
}

# obtain automatically the name of column for initialPopulation 1
# PROBLEM to solve in the future: variable being defined statically
initialPop1ColName <- "initMutPop_1"

# second: extract probabilities from the matrix
dataframeNames <- names(upOrDownDataList)
for(i in 1:length(dataframeNames)) {

    dataframeName <- dataframeNames[i]

    # extract vector of events
    desiredDataframe <- upOrDownDataList[[dataframeName]]
    upOrDownVec <- desiredDataframe[,initialPop1ColName]
    numSucesses <- sum(upOrDownVec)
    numAttempts <- length(upOrDownVec)

    # discover in what cell of probDataframe you have to record them
    fitnessValueStr <- upOrDownFitnessValsNoTotalPopInfo[i]
    totalPopStr <- upOrDownTotalPop[i]

    probInitPop1[fitnessValueStr,totalPopStr] <- numSucesses/numAttempts
}


# +-----------------------+
# | SAVE TABLE FOR GRAPH1 | ===============================
# +-----------------------+

# modify the column names to make them acceptable for loading
colnames(probInitPop1) <- sub(pattern="_",
                              replacement="totalPop",
                              colnames(probInitPop1))

# generate path for the file
outputFilePath = "./data/processed/"
outputFileName = "probInitPop1.csv"
if(!dir.exists(outputFilePath)) {dir.create(outputFilePath, recursive = TRUE)}
fullPath = paste(outputFilePath, outputFileName, sep="")

write.csv(x = probInitPop1,
          file = fullPath,
          row.names = FALSE,
          quote = FALSE)
