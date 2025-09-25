#Useful tutorials:

#WGCNA Colonoids human
rm(list=ls(all=TRUE))

setwd("") # Set wrorking directory to folder where your input files are present
getwd()

#Loading required libraries
library(dplyr)
library(WGCNA)
library(dplyr)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(ggplot2)

#Preparing metadata and counts data
list.files()

SampleInfoFile <- "Organoids_meta_human.txt"

traitData <- read.delim(SampleInfoFile, stringsAsFactors=FALSE, sep="\t", header=TRUE, row.names=1)
head(traitData)
NORM_DATA_DOXO_COLON <-"Organoids_NORM_human.txt"

norm_counts_R_ODAF <- read.delim(NORM_DATA_DOXO_COLON, sep='\t', stringsAsFactors=FALSE, header=TRUE,row.names = 1)
norm_counts_R_ODAF[1:3,1:3]
norm_counts_R_ODAF[1:50,]
ncol(norm_counts_R_ODAF)
nrow(traitData)

traitData
traitData <- subset(traitData, FactorValuedose != 30.0 & FactorValuedose != 60.0 & FactorValuedose != 0.0 )
traitData <- traitData[ c(1,2,3,5) ]

colnames(traitData)<-c( "Treatment", "severity", "Dose", "time" )

traitData$Treatment[traitData$Treatment == 'dimethyl sulfoxide'] <- 'Controls'

traitData$Treatment[traitData$Treatment == 'doxorubicin'] <- 'Case'

nrow(traitData)
ncol(norm_counts_R_ODAF)

# Identify the matching identifiers between the row names of traitData and the column names of norm_counts_R_ODAF
matching_ids <- intersect(row.names(traitData), colnames(norm_counts_R_ODAF))
# Subset traitData based on the matching identifiers
traitData <- traitData[matching_ids, ]
# Subset norm_counts_R_ODAF based on the matching identifiers
norm_counts_R_ODAF <- norm_counts_R_ODAF[, matching_ids]
ncol(norm_counts_R_ODAF)
nrow(traitData)


#Parent genes file 
library(readxl)

excelfile <- read_xlsx("Final_Database_genes_Ensembl.xlsx") #Provided in supplementary files

# Extract the Mouse gene stable IDs from the second column of excelfile
human_gene_ids <- excelfile$`Human Gene stable ID`

# Filter norm_counts_R_ODAF based on matching Mouse gene stable IDs
nrow(norm_counts_R_ODAF)
filtered_norm_counts <- norm_counts_R_ODAF[row.names(norm_counts_R_ODAF) %in% human_gene_ids, ]
# Print the filtered data frame
print(filtered_norm_counts)
colnames(excelfile)
excelfile_subset <- excelfile[, c("Human Gene stable ID", "Human Gene name")]
# Merge the Human Gene name column into filtered_norm_counts based on Mouse gene stable ID
merged_data <- merge(filtered_norm_counts, excelfile_subset, by.x = "row.names", by.y = "Human Gene stable ID", all.x = TRUE)

nrow(merged_data)
nrow(norm_counts_R_ODAF)

# Print the merged data frame
print(merged_data)

merged_data <- merged_data[, c(ncol(merged_data), 1:(ncol(merged_data)-1))]

merged_data <- merged_data[, -2]

norm_counts_R_ODAF<-merged_data

duplicates <- duplicated(norm_counts_R_ODAF$`Human Gene name`)

# Check for duplicate values in the column
duplicated_names <- duplicated(norm_counts_R_ODAF$`Human Gene name`)

# Identify the duplicate values
duplicate_values <- norm_counts_R_ODAF$`Human Gene name`[duplicated_names]

nrow(norm_counts_R_ODAF)
# Remove duplicate values from the data frame or handle them in an appropriate way
norm_counts_R_ODAF <- norm_counts_R_ODAF[!duplicated_names, ]

nrow(norm_counts_R_ODAF)
# Assign the unique column as row names
rownames(norm_counts_R_ODAF) <- norm_counts_R_ODAF$`Human Gene name`
nrow(norm_counts_R_ODAF)

#Removing column with Human gene names as the are already rownames.
norm_counts_R_ODAFX <- norm_counts_R_ODAF[, -1]

#Traitdata numbering
traitData$time[traitData$time == '24'] <- '1'
traitData$time[traitData$time == '48'] <- '2'
traitData$time[traitData$time == '72'] <- '3'

traitData
traitData = traitData[order(rownames(traitData)),]
colnames(traitData)
traitData$Dose<-as.numeric(as.factor(traitData$Dose))

# Transform severity to have levels in the order "control", "low", "high"
traitData$severity <- factor(traitData$severity, levels = c("control", "low", "high"))
# Convert the ordered factor to numeric for further analysis if needed
traitData$severity <- as.numeric(traitData$severity)

# Display the transformed data
traitData

traits <- traitData %>% 
  mutate(case_controls = ifelse(grepl('Case', Treatment), 1, 0)) %>% 
  select(5)

traitData$Treatment <- factor(traitData$Treatment, levels = c("Controls", "Case"))

traitData$Treatment<-as.numeric(as.factor(traitData$Treatment))

traitData<- traitData[ , c(1,2,4)]
colnames(traitData)<-c("Dose+Time","Dose","Time")

traits <- cbind(traits, traitData)
nrow(traits)

colnames(traits)

traits<- traits[ , c(2,3,4)]

data<-norm_counts_R_ODAFX

#Cleaning data
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

traitData
# Outlier removal
traitData <- traitData [c(grep("C47",rownames(traitData), invert=TRUE)),]

traitData <- traitData [c(grep("C48",rownames(traitData), invert=TRUE)),]

traitData <- traitData [c(grep("C30",rownames(traitData), invert=TRUE)),]
norm_counts_R_ODAFX <- norm_counts_R_ODAFX[, rownames(traitData) ]
NCOL(norm_counts_R_ODAFX)
NROW(traitData)

NCOL(norm_counts_R_ODAFX)
NROW(traitData)

data<-norm_counts_R_ODAFX

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")

ncol(data)

nrow(traitData)

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


norm_counts_R_ODAF_TRANSPOSE<- t(norm_counts_R_ODAFX)
#Renaming
norm.counts <- norm_counts_R_ODAF_TRANSPOSE
rownames(norm.counts)

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

# visualization to pick power

sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

library(gridExtra)

grid.arrange(a2, a1, nrow = 2)

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 25
temp_cor <- cor
cor <- WGCNA::cor

getwd()


getwd()
#This takes time, grab a coffee!
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs

nrow(module_eigengenes)

getwd()

# Print out a preview
head(module_eigengenes)
# get number of genes for each module

SAVETABLE_colonoids_human<-table(bwnet$colors)

SAVETABLE_colonoids_human <-as.data.frame(SAVETABLE_colonoids_human)

colnames(SAVETABLE_colonoids_human)<- c("Colors","Genes")

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

nrow(module_eigengenes)
nrow(traits)

nrow(module_eigengenes)
nrow(traits)

# Print the modified 'traits' data frame
print(traits)

library(dplyr)
str(traits)
nrow(module_eigengenes)
nrow(traits)

# Check for common sample names
common_samples <- intersect(rownames(module_eigengenes), rownames(traits))

# Subset both data frames to include only common samples
module_eigengenes_aligned <- module_eigengenes[common_samples, , drop = FALSE]
traits_aligned <- traits[common_samples, , drop = FALSE]

# Compute correlation
module.trait.corr <- cor(module_eigengenes_aligned, traits_aligned, use = 'p')

# Optional: display the correlation matrix
print(module.trait.corr)
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)
nrow(module_eigengenes)

# visualize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

colnames(module_eigengenes)

colnames(heatmap.data)
# Original module names
original_names <- colnames(module_eigengenes) 

# New module names (optional)
new_names <- paste("Module", 1:length(original_names))

# Rename columns
colnames(module_eigengenes) <- new_names    # here

colnames(module_eigengenes)

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

# Rename only the first 12 columns
colnames(heatmap.data)[1:12] <- paste0("Module ", 1:12)

colnames(heatmap.data)

library(dplyr)

X<-CorLevelPlot(heatmap.data,
                x = names(heatmap.data)[13:15],
                y = names(heatmap.data)[1:12],
                col = c("blue1", "skyblue", "white", "pink", "red"),
                titleY = "", rotTitleY = 100,cexLabY = 2,
                rotLabY = 0, fontLabY = 3,cexCorval = 2,cexTitleY = 4,titleX = "Module-Trait Relationship", main =                                         "Mouse Colonoid",cexTitleX = 2, 
                signifSymbols = c("***", "**", "*", ""),
                signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),cexLabX = 2,rotTitleX = 0)



#Extracting genes for pathway analysis

module.gene.mapping <- as.data.frame(bwnet$colors)

#Examples on saving genes for significant module for downstream analysis based on color name

MEred<-module.gene.mapping %>% 
  filter(`bwnet$colors` == 'red') %>% 
  rownames()

ME_turq<-module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()

write.csv(ME_turq, file = "ME_turq.csv", row.names = FALSE)
write.csv(MEred, file = "MEred.csv", row.names = FALSE)

# End of script