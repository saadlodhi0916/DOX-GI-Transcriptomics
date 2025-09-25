# DEG analysis example

  ## **Applied tool:**
  
  ######################################################
#DESeq2 for normalization and DE identification 
######################################################

#DESeq2 was used according to R-ODAF selection criteria for differentially expressed genes
#For more information, refer to the paper: https://www.sciencedirect.com/science/article/pii/S0273230022000307
### **PARAMETERS**                            		 

RunID<-"VIVO_DOX_TREATMENT_Colon_Basic_DEGanalysis"  # You can use a name tailored to your data
species<-"human"
Platform <- "RNA-Seq"

# Specify which groups need to be compared 
Samp4compare<- c("case")
Cont4compare<- c("control")

# Applied normalization strategy
NormType<-paste0(Samp4compare, "_Subset-Matched_Control")
DESIGN<- "treatment" #Design

# Set thresholds for Differential Expression (based on R-ODAF)
minCoverage <- 5000000
MinCount<- 1
pAdjValue<- 0.01 

# Names of files to load
SampleDataFile <- "SAMPLEDATA_COLON.txt"

SampleInfoFile <- "vivoDox_colon.txt"

##### **Normalization strategy**
#All samples where normalized together: Lower amount of DEGs, but DEGs more robust and effects can be compared between DIs.

################################
# Load R libraries					  #
################################ 

library(knitr)
library(DT)
library(data.table)
library(purrr)
library(GenomicFeatures)
library(DESeq2)
require("DESeq2")
require("edgeR")
require("pheatmap")
library("ggplot2")
library("lattice")
library(AnnotationHub)
library(AnnotationFilter)
library(ensembldb)
library(multtest)
library(tibble)
library(dbplyr)
library(ggplot2)
library(ComplexUpset)
library(dbplyr)
library(corrplot)
library("Hmisc")
options(timeout=30000)

# Set directories 

Glob.Data.Dir <- "C:/Users/...../...../..../..../Colon/R-ODAFoutput/"

Glob.output.dir<- "C:/Users/...../...../..../..../Colon/Global_output_directory/"


if (!file.exists(Glob.output.dir)) {dir.create(Glob.output.dir)}
knitr::opts_knit$set(root.dir = Glob.output.dir)

DEGoutputdir<-"C:/Users/..../..../..../...../..../Colon/DEGoutputdir/"


if (!file.exists(DEGoutputdir)) {dir.create(DEGoutputdir)}

knitr::opts_knit$set(root.dir = DEGoutputdir)

DEGoutputdir_NORM<-paste0(DEGoutputdir, NormType, "/")

if (!file.exists(DEGoutputdir_NORM)) {dir.create(DEGoutputdir_NORM)}

# Set file locations
sampledir <- Glob.output.dir
outputdir <- DEGoutputdir

# Load input files 
sampleData <- read.delim(SampleDataFile, sep='\t', stringsAsFactors=FALSE, header=TRUE,row.names = 1)

DESeqDesign <- read.delim(SampleInfoFile, stringsAsFactors=FALSE, sep="\t", header=TRUE, row.names=1)
x=1 

condition1<- Cont4compare[x]	    		
condition2<- Samp4compare[x]

sapply(DESeqDesign, class)                                    # Print classes of all colums

DE_Design <- matrix(data=NA, ncol=2)

DE_Design <- DESeqDesign [c(grep(condition1,DESeqDesign[,DESIGN[x]]), grep(condition2,DESeqDesign[,DESIGN[x]])),]
sampleData <- sampleData[, rownames(DESeqDesign) ]
ncol(sampleData)
NROW(DESeqDesign)  

samples <- sampleData[, rownames(DE_Design) ]

print(paste("Comparison: ", condition2, " vs ", condition1, ":", NormType))	
outputdir<-DEGoutputdir_NORM

##### *Data clean-up: replace NA with 0*
ZeroDetected<-count(samples[ is.na(samples) ])
samples[ is.na(samples) ] <- 0 
print(paste0(ZeroDetected, "x replacement of NA values"))

##### *Remove samples with total readcount < threshold (5M)*
Keep<-ncol(samples[,(colSums(samples)> minCoverage)])
Remove<-ncol(samples[,(colSums(samples)< minCoverage)])
print(paste0("From the total of ", ncol(samples), " samples, ", Remove, " samples had to be removed due to low sequencing depth (<5M) -> ", Keep, " samples remaining"))
samples<- samples[,(colSums(samples)> minCoverage)]
DE_Design <- DE_Design[rownames(DE_Design) %in% colnames(samples),] #Pheno file

##### *Remove outliers (>20% variance between samples)*
# FUNCTION: CHECK FOR OUTLIERS
check_outliers_vst = function(input, RunID, DESIGN, OUTPUT.DIR) {  
  OutNames<-NULL
  vsd = vst(input, blind = T) %>% assay()
  
  p = PCAtools::pca(vsd, 
                    metadata = DE_Design, 
                    removeVar = 0.1)
  
  scree_plot = PCAtools::screeplot(p, axisLabSize = 10, titleLabSize = 22)
  PCA = PCAtools::biplot(p, x = 'PC1', y = 'PC2', lab = rownames(p$metadata),
                         drawConnectors = F, colby = DESIGN, 
                         legendPosition = "bottom", legendLabSize = 13, legendIconSize = 6.0) + coord_fixed(ratio = 1)
  PCA
  #ggrepel::geom_text_repel()
  ggsave(paste0(OUTPUT.DIR, '/', RunID, '_pca.png'))
  
  cowplot::plot_grid(scree_plot, PCA, ncol = 2, nrow = 1)
  ggsave(paste0(OUTPUT.DIR, '/', RunID, '_scree-pca.png'))
  
  # get replicates group
  #  levels_compare = p$metadata$Dose %>% droplevels() %>% unique() #original code Martha
  levels_compare = unique(c(condition1, condition2))
  
  
  # find PC with var >= 20
  pc_check = p$variance[which(p$variance >= 20)]
  
  # if there are no PC >= 20%, print message and stop
  if (length(pc_check) < 1) {
    
    print('No PC contributes to 20% or more of total variance')
    
  } else {
    
    # check outliers in all PC that contribute for at least 20% of total variance
    for (x in seq(1, length(pc_check))) {
      
      print(pc_check[x])
      # max distance along PC
      max_pc = abs(max(p$rotated[x]) - min(p$rotated[x]))
      
      for (i in seq(1, length(levels_compare))) {
        # name of all replicates
        all_reps = p$metadata %>% dplyr::filter(DI %in% levels_compare[i]) %>% rownames()
        
        # PC values for all replicates
        all_reps_value = p$rotated %>% 
          .[rownames(.) %in% all_reps, x, drop = F] %>% 
          as.data.frame() %>% 
          .[order(.[, 1], decreasing = T), , drop = F]
        
        # distance on PC for all replicates. Pairwise and between closest
        all_reps_dist = all_reps_value %>% 
          dplyr::select(1) %>% 
          dist %>%  # get distance
          as.matrix() %>% 
          .[row(.) == col(.) + 1]  # get diagonal
        
        # check if any distance is >= 20% of `max_pc`
        # max_pc : var(PC) = all_reps_dist : x%
        # TRUE = outlier 
        outliers = ( (as.numeric(p$variance[x]) * all_reps_dist) / max_pc ) > 20
        
        # Obtaining names of outliers
        outlierNames1<-NULL
        check_first_half<-sum(outliers[1:(length(outliers)/2)])
        if(check_first_half>0){
          TRUEpos<-NULL
          for(pos in 1:(length(outliers)/2)) {
            if (outliers[pos]==TRUE){TRUEpos<-c(TRUEpos, pos)}
          }
          outlierNames1<-row.names(all_reps_value)[1:max(TRUEpos)]  
        }
        
        outlierNames2<-NULL
        check_second_half<-sum(outliers[(length(outliers)/2):length(outliers)])
        if (check_second_half>0){
          TRUEpos<-NULL
          for(pos in (length(outliers)/2):length(outliers)) {
            if (outliers[pos]==TRUE){TRUEpos<-c(TRUEpos, pos)}
          }
          outlierNames2<-row.names(all_reps_value)[(min(TRUEpos)+1):(length(outliers)+1)]  
        }
        outlierNames<-c(outlierNames1, outlierNames2)
        
        print(paste0(levels_compare[i], " contained ", length(outlierNames), " outliers: "))
        print(outlierNames)
        OutNames<-c(OutNames, outlierNames)
      }
    }
  }
  OutNames<<-OutNames
  All_OutNames<<-c(All_OutNames, OutNames)
}

# FUNCTION: CHECK FOR OUTLIERS
check_outliers_norm_data = function(input, RunID, DESIGN, OUTPUT.DIR) {  
  OutNames<-NULL
  p = PCAtools::pca(input, 
                    metadata = DE_Design, 
                    removeVar = 0.1)
  
  scree_plot = PCAtools::screeplot(p, axisLabSize = 10, titleLabSize = 22)
  PCA = PCAtools::biplot(p, x = 'PC1', y = 'PC2', lab = rownames(p$metadata),
                         drawConnectors = F, colby = DESIGN, 
                         legendPosition = "bottom", legendLabSize = 13, legendIconSize = 6.0) + coord_fixed(ratio = 1)
  PCA
  ggsave(paste0(OUTPUT.DIR, '/', RunID, '_pca.png'))
  
  cowplot::plot_grid(scree_plot, PCA, ncol = 2, nrow = 1)
  ggsave(paste0(OUTPUT.DIR, '/', RunID, '_scree-pca.png'))
  
  # get replicates group
  levels_compare = unique(c(condition1, condition2))
  
  
  # find PC with var >= 20
  pc_check = p$variance[which(p$variance >= 20)]
  
  # if there are no PC >= 20%, print message and stop
  if (length(pc_check) < 1) {
    
    print('No PC contributes to 20% or more of total variance')
    
  } else {
    
    # check outliers in all PC that contribute for at least 20% of total variance
    
    for (x in seq(1, length(pc_check))) {
      
      print(pc_check[x])
      # max distance along PC
      max_pc = abs(max(p$rotated[x]) - min(p$rotated[x]))
      
      for (i in seq(1, length(levels_compare))) {
        # name of all replicates
        all_reps = p$metadata %>% dplyr::filter(DI %in% levels_compare[i]) %>% rownames()
        
        # PC values for all replicates
        all_reps_value = p$rotated %>% 
          .[rownames(.) %in% all_reps, x, drop = F] %>% 
          as.data.frame() %>% 
          .[order(.[, 1], decreasing = T), , drop = F]
        
        # distance on PC for all replicates. Pairwise and between closest
        all_reps_dist = all_reps_value %>% 
          dplyr::select(1) %>% 
          dist %>%  # get distance
          as.matrix() %>% 
          .[row(.) == col(.) + 1]  # get diagonal
        
        # check if any distance is >= 20% of `max_pc`
        # max_pc : var(PC) = all_reps_dist : x%
        # TRUE = outlier 
        outliers = ( (as.numeric(p$variance[x]) * all_reps_dist) / max_pc ) > 20
        
        # Obtaining names of outliers
        outlierNames1<-NULL
        check_first_half<-sum(outliers[1:(length(outliers)/2)])
        if(check_first_half>0){
          TRUEpos<-NULL
          for(pos in 1:(length(outliers)/2)) {
            if (outliers[pos]==TRUE){TRUEpos<-c(TRUEpos, pos)}
          }
          outlierNames1<-row.names(all_reps_value)[1:max(TRUEpos)]  
        }
        
        outlierNames2<-NULL
        check_second_half<-sum(outliers[(length(outliers)/2):length(outliers)])
        if (check_second_half>0){
          TRUEpos<-NULL
          for(pos in (length(outliers)/2):length(outliers)) {
            if (outliers[pos]==TRUE){TRUEpos<-c(TRUEpos, pos)}
          }
          outlierNames2<-row.names(all_reps_value)[(min(TRUEpos)+1):(length(outliers)+1)]  
        }
        outlierNames<-c(outlierNames1, outlierNames2)
        
        print(paste0(levels_compare[i], " contained ", length(outlierNames), " outliers: "))
        print(outlierNames)
        OutNames<-c(OutNames, outlierNames)
      }
      OutNames<<-OutNames
      All_OutNames<<-c(All_OutNames, OutNames)
    }
  }
  
}

#setwd(outputdir)

############################# Normalize data; valuable for both with and without R-ODAF analysis  ##############################

dds <- DESeqDataSetFromMatrix(countData = round(samples), colData = as.data.frame(DE_Design), design = as.formula(paste0("~", DESIGN)))

dds <- DESeq(dds, quiet=TRUE) #"Get some coffee, next step will take a while"


#counts(dds) == norm_data

norm_data <<- counts(dds,normalized=TRUE)

#View(norm_data)

#View(counts(dds))

vsd = vst(dds, blind = T) %>% assay()

#View(vsd) # Variance stabilized transformation data after normalization counts

p = PCAtools::pca(vsd, 
                  metadata = DE_Design, 
                  removeVar = 0.1)
PCA<-PCAtools::biplot(p, x = 'PC1', y = 'PC2', lab = rownames(p$metadata),
                      drawConnectors = F, colby = DESIGN, 
                      legendPosition = "bottom", legendLabSize = 13, legendIconSize = 6.0) + coord_fixed(ratio = 1)

PCA

# first outlier check
All_OutNames<-NULL

##The function below was prepared in the beginning of this  tutorial, this one runs vst
check_outliers_vst(dds, RunID, DESIGN, outputdir)

for (name in OutNames) {
  Remove_Pcode<-paste0(unlist(strsplit(name, "_"))[1], "_")
  DE_Design <- DE_Design [c(grep(Remove_Pcode,rownames(DE_Design), invert=TRUE)),]
  samples <- samples[, rownames(DE_Design) ]
}
print("These outliers (and matched samples) have been removed from the dataset. Check for more outliers")

setwd(outputdir)
# Rerun normalization & outlier check until no outliers are left
while (length(OutNames)>0) {
  # Normalize data 
  dds <- DESeqDataSetFromMatrix(countData = round(samples), colData = as.data.frame(DE_Design), design = as.formula(paste0("~", DESIGN)))
  dds <- DESeq(dds, quiet=TRUE) #"Get some coffee, next step will take a while"
  norm_data <<- counts(dds,normalized=TRUE)
  
  ##was _vst
  check_outliers_norm_data(dds, RunID, DESIGN, outputdir)
  
  for (name in OutNames) {
    Remove_Pcode<-paste0(unlist(strsplit(name, "_"))[1], "_")
    DE_Design <- DE_Design [c(grep(Remove_Pcode,rownames(DE_Design), invert=TRUE)),]
    samples <- samples[, rownames(DE_Design) ]
  }
  print("These outliers (and matched samples) have been removed from the dataset. Check for more outliers")
}
print(paste0("In total, ", length(All_OutNames), " have been excluded:"))
print(All_OutNames)

##### **PCA plot after exclusion**
vsd = vst(dds, blind = T) %>% assay()

p = PCAtools::pca(vsd, 
                  metadata = DE_Design, 
                  removeVar = 0.1)

scree_plot = PCAtools::screeplot(p, axisLabSize = 10, titleLabSize = 22)
PCA = PCAtools::biplot(p, x = 'PC1', y = 'PC2', lab = rownames(p$metadata),
                       drawConnectors = F, colby = DESIGN, 
                       legendPosition = "bottom", legendLabSize = 13, legendIconSize = 6.0) + coord_fixed(ratio = 1)

#Plot scree_plot and PCA
scree_plot 
PCA 

############################################
## Differential expression analysis: DESeq2 
############################################

datatable(DE_Design)
setwd(DEGoutputdir_NORM)

# Filtering genes with low readcounts: 
# 75% of at least 1 group need to be above MinCount CPM

SampPerGroup<-table(DE_Design[,DESIGN])
kable(SampPerGroup, caption = "Amount of samples per group")

Counts<-counts(dds, normalized=TRUE)
CPMdds<-cpm(counts(dds, normalized=TRUE))

Filter <- matrix(data=NA, ncol=3, nrow= nrow(Counts))
rownames(Filter) <- rownames(Counts)
colnames(Filter) <- c("Low readcounts","quantile","spike")

for (gene in 1:nrow(dds)) {
  
  CountsPass<-NULL
  for (group in 1:length(SampPerGroup)) { 
    sampleCols<-grep(dimnames(SampPerGroup)[[1]][group],DE_Design[,DESIGN])
    Check<-sum(CPMdds[gene,sampleCols] >= MinCount)>= 0.75*SampPerGroup[group]
    CountsPass<-c(CountsPass, Check)
  }
  
  if ( sum(CountsPass) > 0 ) {Filter[gene,1] <- 1 }	else { Filter[gene,1] <- 0 }
  
}	

compte <- Counts[Filter[,1] == 1,]
Filter <- Filter[rownames(Filter) %in% rownames(compte),]

##### *Filtering genes with low readcounts: 75% of at least 1 group need to be above "MinCount" CPM*
print(paste("low readcount filtering removed ",nrow(dds)- nrow(Filter)," genes from the ",nrow(dds)," assessed. ", nrow(Filter)," genes remaining",sep=""))

##### *Obtaining the DESeq2 results*
# compute the DEGs on the genes passing the Relevance condition

res <- results(dds[rownames(compte),], contrast=c(DESIGN, condition2, condition1), pAdjustMethod= 'fdr')

setwd(outputdir)
FileName<-paste(NormType, condition2,"vs",condition1, "FDR", pAdjValue, sep="_")

#Save output tables		
norm_data <<- counts(dds[rownames(compte)],normalized=TRUE) 
DEsamples <<- subset(res,res$padj < pAdjValue)	

if (nrow(DEsamples)>0) {
  
  print(paste("A total of ", nrow(DEsamples), " DEGs were selected (before filtering)"))
  
  if (nrow(DEsamples)>2) {
    DECounts <- compte[rownames(compte) %in% rownames(DEsamples),]
    Filter <- Filter[rownames(Filter) %in% rownames(DECounts),]
  } else { 
    DECounts <- t(data.frame(compte[rownames(compte) %in% rownames(DEsamples),]))
    row.names(DECounts)<-row.names(DEsamples) 
    Filter <- t(data.frame(Filter[rownames(Filter) %in% rownames(DECounts),]))
    row.names(Filter)<-row.names(DEsamples)
  }
  
  ####################################################################################################
  print("Filtering for relevance (Check median against third quantile) & Filtering spurious spikes")
  
  for (gene in 1:nrow(DECounts)) {
    # Check the median against third quantile
    quantilePass <-NULL
    sampleColsg1 <- grep(dimnames(SampPerGroup)[[1]][1],DE_Design[,DESIGN])
    sampleColsg2 <- grep(dimnames(SampPerGroup)[[1]][2],DE_Design[,DESIGN])
    
    Check <- median(DECounts[gene,sampleColsg1]) > quantile(DECounts[gene,sampleColsg2], 0.75)[[1]]
    quantilePass <-c(quantilePass, Check)
    Check <- median(DECounts[gene,sampleColsg2]) > quantile(DECounts[gene,sampleColsg1], 0.75)[[1]]
    quantilePass <-c(quantilePass, Check)
    
    if ( sum(quantilePass) > 0 ) {Filter[gene,2] <- 1 }	else { Filter[gene,2] <- 0 }
    
    # Check for spike 
    spikePass <- NULL
    for (group in 1:length(SampPerGroup)) { 
      sampleCols<-grep(dimnames(SampPerGroup)[[1]][group],DE_Design[,DESIGN])
      if (max(DECounts[gene,sampleCols]) ==0) {Check <- FALSE} else {
        Check <- (max(DECounts[gene,sampleCols])/sum(DECounts[gene,sampleCols])) >= 1.4*(SampPerGroup[group])^(-0.66)
        spikePass<-c(spikePass, Check)
      }
    }
    if ( sum(spikePass) > 1 ) {Filter[gene,3] <- 0 }	else { Filter[gene,3] <- 1 }
  }		
  
  print("extract the final list of DEGs")
  
  DECounts_real <- DEsamples[rowSums(Filter) == 3 ,]
  DECounts_no_quant <- DEsamples[Filter[,2] == 0 ,]
  DECounts_spike <- DEsamples[Filter[,3] == 0 ,]
  
  #DECounts_real2<-annotate_DEGs(DECounts_real)      
  
  print(paste0(condition2, " vs ", condition1, ": A total of ",nrow(DECounts_real), " DEGs were selected, after ",nrow(DECounts_no_quant)," genes(s) removed by the quantile rule and ", nrow(DECounts_spike)," gene(s) with a spike --> saved in: ", paste0(RunID, "_DEGtable_", NormType, "_", condition2, "_vs_", condition1, "_WITHrodafCriteria.csv")))
  write.table(DECounts_real, paste0(RunID, "_DEGtable_", NormType, "_", condition2, "_vs_", condition1, "_WITHrodafCriteria.csv"), sep="\t")
  write.table(norm_data,file=paste0(paste0(RunID, "_NormData_", NormType, "_", condition2, "_vs_", condition1, "_WITHrodafCriteria.csv")), sep="\t", quote=FALSE)
  
}else{
  print(paste0("Applying R-ODAF filtering criteria resulted in NO DEGs identified for ", condition2, " vs ", condition1))
}

#DESeq2 identification of DEGs finished. 

# **DESeq2 Results**

### **Heatmap of all genes**

cormM = cor(norm_data, method = "pearson")

mydata.cor = cor(norm_data, method = c("pearson"))

corrplot(mydata.cor)

palette = colorRampPalette(c("green", "white", "red")) (20)

heatmap(x = mydata.cor, col = palette, symm = TRUE)

### **Heatmap of DEGs**

if (nrow(DECounts_real)>2) {
  DEG_counts<-norm_data[c(row.names(DECounts_real)),]
  heatmap(DEG_counts)
} else { 
  print("unable to generate a heatmap with <2 DEGs")
}

### **Table of identified DEGs (FDR<0.01)**
DECounts_real <-as.data.frame(DECounts_real)

datatable(DECounts_real)
NROW(DECounts_real) #FINAL NUMBER OF DEGs.

# Saving DEGs and related statistics for downstream analysis
write.csv(
  DECounts_real,
  file = "DECounts_real.csv",
  row.names = TRUE
)
### **End of DEGanalysis**

