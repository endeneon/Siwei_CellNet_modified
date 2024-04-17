# Siwei 13 Sept 2023
# CellNet pre-configure file

# install packages ####
install.packages("devtools")
library(devtools)
install_github("pcahan1/CellNet",
               ref = "master")
# ERROR: dependencies ‘GO.db’, ‘org.Hs.eg.db’
# are not available for package ‘cancerCellNet’
install_github("pcahan1/cancerCellNet@v0.1.1",
               ref = "master")


source("pacnet_utils.R")

# load required R packages ####
{
  library(CellNet)
  library(cancerCellNet)
  library(plyr)
  library(ggplot2)
  library(RColorBrewer)
  library(pheatmap)
  library(plotly)
  library(igraph)
}
