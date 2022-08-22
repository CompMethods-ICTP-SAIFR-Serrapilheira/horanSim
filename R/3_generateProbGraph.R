# Filename: 3_generateProbGraph.R
# Purpose: generate the probability graphs
# from the table "./data/processed/probInitPop1.csv"
# (which is produced at the end of "2_generateProcessedTable.R")

# +--------------------+
# | IMPORTING PACKAGES | =========================================
# +--------------------+

packages <- c("ggplot2","reshape2", "dplyr")

# method recommended by Antoine Soetewey
# <https://statsandr.com/blog/an-efficient-way-to-install-and-load-r-packages/#what-is-a-r-package-and-how-to-use-it>
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Loading package
invisible(lapply(packages, library, character.only = TRUE))


# +------------------+
# | CREATE GRAPH DIR | =========================================
# +------------------+

outputDir <- "./graph"
# iterate through all the fitnessMutantVal values
# and save them in ./data/raw
if(!dir.exists(outputDir)) {dir.create(outputDir, recursive = TRUE)}


# +------------------------------------------+
# | READING AND PREPARING TABLE FOR PLOTTING | =========================================
# +------------------------------------------+

probInitPop1 <- read.csv(file = "./data/processed/probInitPop1.csv",
                row.names = NULL)

meltedProbInitPop1 <- reshape2::melt(data = probInitPop1,
                                     id = "fitnessValue",
                                     variable_name = "totalPop")

names(meltedProbInitPop1) <- c("fitnessValue","totalPop","prob")

# turn the totalPop column in "numeric characters"
meltedProbInitPop1$totalPop <- as.factor(meltedProbInitPop1$totalPop)
levelsTotalPop <- levels(meltedProbInitPop1$totalPop)
levelsTotalPop <- sub("totalPop",
                      "",
                      levelsTotalPop,
                      fixed=TRUE)
levels(meltedProbInitPop1$totalPop) <- levelsTotalPop

# reorder by totalPop in ascending order
ascendingTotalPop <- as.character(sort(as.numeric(levels(meltedProbInitPop1$totalPop))))
meltedProbInitPop1$totalPop <- factor(meltedProbInitPop1$totalPop, levels=ascendingTotalPop)
meltedProbInitPop1 <- arrange(meltedProbInitPop1, totalPop)



# +-------------------------+
# | PLOTTING THE TABLE DATA | =========================================
# +-------------------------+

graphFilename <- "fixationProb.png"
xAxisTitle <- "Relative fitness"
yAxisTitle <- "Fixation probability"
legendTitle <- "Total population"
lineSize <- 1
marginRight <- 1
textSize <- 12

plotProb <- ggplot(data = meltedProbInitPop1, aes(x = fitnessValue, y = prob)) +
            geom_line(aes(colour = totalPop), na.rm = TRUE, size = lineSize*1.5) +
            labs(x = xAxisTitle, y = yAxisTitle, colour = legendTitle) +
            scale_x_continuous(labels = scales::comma_format(big.mark=",", decimal.mark="."),
                               expand = c(0.0,0.01),
                                limits = c(1, 1.8)) +
            scale_y_continuous(labels = scales::comma_format(big.mark=",", decimal.mark="."),
                               expand = c(0.0,0.01),
                               limits = c(0, 0.5)) +
            theme(panel.background = element_blank()) +
            theme(panel.spacing = unit(1, "lines")) +
            theme(panel.grid.major = element_line(colour = "grey90")) +
            theme(plot.margin = margin(8, marginRight, 8, 8)) +
            theme(strip.text = element_text(size=textSize*1.2)) +
            theme(strip.background = element_rect(fill = "white", colour = NULL)) +
            theme(axis.line = element_line(colour = "black", size = lineSize)) +
            theme(axis.text = element_text(size = textSize*1.2)) +
            theme(axis.title.x = element_text(size = textSize*1.2, margin = margin(t=10, r=0, b=0, l=0))) +
            theme(axis.title.y = element_text(size = textSize*1.2, margin = margin(t=0, r=10, b=0, l=0))) +
            theme(legend.key.size = unit(2, "lines")) +
            theme(legend.text = element_text(size = textSize)) +
            theme(legend.title = element_text(size = textSize, face = "bold"))
            #facet_wrap(~sample,  ncol = numOfColsGraph, nrow = numOfRowsGraph, scales = 'free')

fullPath = file.path(outputDir, graphFilename)

ggsave(filename = fullPath,
       plot = plotProb,
       dpi = 600,
       width = 8,
       height = 4)
