# Siwei 15 Oct 2023
# Use MiNND batch 1 data set
# Check if homogeneity found between different samples
# (issue: no difference found between risk/wt samples)

# Load libraries
{
  library(CellNet)
  library(cancerCellNet)
  library(plyr)
  library(ggplot2)
  library(RColorBrewer)
  library(pheatmap)
  library(plotly)
  library(igraph)
  source("pacnet_utils.R")

  library(readr)
  library(stringr)

  library(edgeR)
}

# source("pacnet_utils.R")

# init #####
set.seed(42)
# load ENSG Geneid-Gene_Symbol ref list
load("~/backuped_space/Siwei_misc_R_projects/CellNet/ENSG_gene_index.RData")

# load data #####
df_raw <-
  read_delim("MiNND/MiNND_06Oct2023.txt",
             delim = "\t", escape_double = FALSE,
             trim_ws = TRUE)

# convert gene id to gene symbol for intersection
df_4_DGE <- df_raw
df_4_DGE$Geneid <-
  str_split(string = df_4_DGE$Geneid,
            pattern = '\\.',
            simplify = T)[, 1]

df_4_DGE <-
  merge(x = df_4_DGE,
        y = ENSG_anno_gene_indexed,
        by = "Geneid")
df_4_DGE <-
  df_4_DGE[!duplicated(df_4_DGE$Gene_Symbol), ]
# Keep gene name as a separate vector
gene_name_list <-
  df_4_DGE$Gene_Symbol

# remove all non-numerical columns
df_4_DGE$Geneid <- NULL
df_4_DGE$Gene_Symbol <- NULL

# remove all genes with zero counts
gene_name_list <-
  gene_name_list[rowSums(df_4_DGE, na.rm = T) != 0]
df_4_DGE <-
  df_4_DGE[rowSums(df_4_DGE, na.rm = T) != 0, ]
# reassign rownames!
rownames(df_4_DGE) <-
  gene_name_list
colnames(df_4_DGE) <-
  str_remove_all(string = colnames(df_4_DGE),
                 pattern = "_")

# build DGE object to calculate cpm
df_DGE <-
  DGEList(counts = as.matrix(df_4_DGE),
          samples = colnames(df_4_DGE),
          genes = rownames(df_4_DGE),
          remove.zeros = T)

df_DGE <-
  calcNormFactors(df_DGE)
df_DGE <-
  estimateDisp(df_DGE)

## calc cpm, requires log transform #####
df_queryExpDat <-
  cpm(df_DGE,
      normalized.lib.sizes = T,
      log = T)

## construct sample metadata table #####
df_querySampTab <-
  data.frame(sample_name = colnames(df_queryExpDat),
             pe_se = "PAIRED",
             description1 = str_replace_all(string = colnames(df_queryExpDat),
                                            pattern = "V",
                                            replacement = "V00087"),
             description2 = colnames(df_queryExpDat),
             stringsAsFactors = F)
rownames(df_querySampTab) <-
  df_querySampTab$sample_name

# make sure the count matrix and metadata have the same row/col orders
df_queryExpDat <-
  df_queryExpDat[, match(df_querySampTab$sample_name,
                         colnames(df_queryExpDat))]

# remove Sum gene counts == 0 (should be already removed, just to confirm...)
gene_name_list <- rownames(df_queryExpDat)
df_queryExpDat_raw <- df_queryExpDat
df_queryExpDat <-
  df_queryExpDat[!rowSums(df_queryExpDat, na.rm = T) == 0, ]
rownames(df_queryExpDat) <-
  gene_name_list[!rowSums(df_queryExpDat_raw, na.rm = T) == 0]
min(rowSums(df_queryExpDat))
rm(df_queryExpDat_raw)

# Assign variables ####
study_name <-
  "MiNND_16Oct2023"
queryExpDat <- df_queryExpDat
querySampTab <- df_querySampTab

# Load training data: ####
expTrain <-
  utils_loadObject("Hs_expTrain_Jun-20-2017.rda")
stTrain <-
  utils_loadObject("Hs_stTrain_Jun-20-2017.rda")

iGenes <-
  intersect(rownames(expTrain),
            rownames(queryExpDat))

expTrain <- expTrain[iGenes, ]

# Split the training data into a training subset and a validation subset: ####
stList <-
  splitCommon_proportion(sampTab = stTrain,
                         proportion = 0.66,
                         dLevel = "description1") # Use 2/3 of training data for training and 1/3 for validation
stTrainSubset <- stList$trainingSet
expTrainSubset <-
  expTrain[, rownames(stTrainSubset)]

# See number of samples of each unique type in description1 in training subset
table(stTrainSubset$description1)
min(table(stTrainSubset$description1))

stValSubset <-
  stList$validationSet
expValSubset <-
  expTrain[, rownames(stValSubset)]
# See number of samples of each unique type in description1 in validation subset
table(stValSubset$description1)

# Train the random forest classifier, takes 3-10 minutes depending on memory availability:
system.time(expr = (my_classifier =
                       broadClass_train(stTrain =
                                           stTrainSubset,
                                        expTrain = expTrainSubset,
                                        colName_cat = "description1",
                                        # colName_samp = "sra_id",
                                        nRand = 70,
                                        nTopGenes = 100,
                                        nTopGenePairs = 100,
                                        nTrees = 2000,
                                        stratify = TRUE,
                                        sampsize = min(table(stTrainSubset$description1)) - 1, # Must be less than the smallest n in table(stTrainSubset$description1)
                                        quickPairs = TRUE))) # Increasing the number of top genes and top gene pairs increases the resolution of the classifier but increases the computing time
saveRDS(my_classifier,
        file = "MiNND/cellnet_classifier_100topGenes_100genePairs.rda")

### section #######
# Classifier Validation

stValSubsetOrdered <-
  stValSubset[order(stValSubset$description1), ] #order samples by classification name
expValSubset <-
  expValSubset[, rownames(stValSubsetOrdered)]
cnProc <-
  my_classifier$cnProc #select the cnProc from the earlier class training

classMatrix <-
  broadClass_predict(cnProc,
                     expValSubset,
                     nrand = 60)
stValRand <-
  addRandToSampTab(classRes = classMatrix,
                   sampTab = stValSubsetOrdered,
                   desc = "description1",
                   id = "sra_id")

grps <-
  as.vector(stValRand$description1)
names(grps) <-
  rownames(stValRand)

# Gene pair validation
genePairs <- cnProc$xpairs
# Get gene to gene comparison of each gene pair in the expression table
expTransform <-
  query_transform(expTrainSubset,
                  genePairs)
avgGenePair_train <-
  avgGeneCat(expDat = expTransform,
             sampTab = stTrainSubset,
             dLevel = "description1",
             sampID = "sra_id")

genePairs_val <-
  query_transform(expValSubset,
                  genePairs)
geneCompareMatrix <-
  makeGeneCompareTab(queryExpTab = genePairs_val,
                     avgGeneTab = avgGenePair_train,
                     geneSamples = genePairs)
val_grps <-
  stValSubset[ ,"description1"]
val_grps <- c(val_grps,
              colnames(avgGenePair_train))
names(val_grps) <-
  c(rownames(stValSubset),
    colnames(avgGenePair_train))

# Create and save xpairs_list object
xpairs_list <- vector("list", 14)
for (pair in rownames(avgGenePair_train)) {
  for (j in 1:ncol(avgGenePair_train)) {
    if (avgGenePair_train[pair,j] >= 0.5) {
      if (is.null(xpairs_list[[j]])) {
        xpairs_list[[j]] <- c(pair)
      } else {
        xpairs_list[[j]] <- c(xpairs_list[[j]], pair)
      }
    }
  }
}

xpair_names <- colnames(avgGenePair_train)
xpair_names <- sub(pattern = "_Avg",
                   replacement = "",
                   x = xpair_names)
names(xpairs_list) <- xpair_names

for (type in names(xpairs_list)) {
  names(xpairs_list[[type]]) <- xpairs_list[[type]]
}

classMatrixQuery <-
  broadClass_predict(cnProc = cnProc,
                     expDat = queryExpDat,
                     nrand = 3)
grp_names <-
  c(as.character(querySampTab$description1),
    rep("random", 3))
names(grp_names) <-
  c(as.character(rownames(querySampTab)),
    paste0("rand_",
           c(1:3)))
# Re-order classMatrixQuery to match order of rows in querySampTab
classMatrixQuery <-
  classMatrixQuery[, names(grp_names)]
# save(classMatrixQuery,
#      file = "MiNND/classificationMatrix.rda")

# classMatrixQuery <-
#   load("MiNND/classificationMatrix.rda")
#Plot classification heatmap:
# png(file = "example_outputs/query_classification_heatmap.png",
#     height = 4,
#     width = 8,
#     units = "in",
#     res = 300)
# This function can be found in pacnet_utils.R
acn_queryClassHm(classMatrixQuery,
                 main = paste0("Classification Heatmap"),
                 grps = grp_names,
                 fontsize_row = 10,
                 fontsize_col = 10,
                 isBig = F,
                 cCol = T, cRow = F,
                 scale = F)
dev.off()
write.table(classMatrixQuery,
            file = "MiNND_iPS_classificationScore_16Oct2023.tsv",
            sep = "\t")
########


