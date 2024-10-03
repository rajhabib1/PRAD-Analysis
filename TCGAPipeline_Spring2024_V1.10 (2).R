opt <- list() #Instantiate opt
opt$project    <- "TCGA-PRAD"       


opt$wd <- path.expand(file.path("~","bio321g","DiffExp"))
opt$pckDir  <- c(
    file.path(opt$wd,"R_gdc"),              
    "/stor/scratch/Bio321G_NB_Spring2024/R" 
)
opt$gdcPath <- "/stor/scratch/Bio321G_NB_Spring2024/GDCdata" 

opt$biocPackages <- c(
    "ggplot2","ggrepel","TCGAbiolinks","DESeq2"
)



if(!dir.exists(opt$wd)){dir.create(opt$wd,recursive = T,showWarnings = T)}
setwd(opt$wd)


opt$backupPath    <- file.path(opt$wd,"R_gdc")
opt$backupGdcPath <- file.path(opt$wd,"GDCdata")


if(!dir.exists(opt$pckDir[1])){
    message("Creating home package directory:\n",opt$backupPath)
    dir.create(opt$backupPath,recursive = T,showWarnings = T)
}
if(!all(dir.exists(opt$pckDir))){
    message("Changing package directory to:\n",opt$backupPath)
    opt$pckDir <- opt$backupPath
}


.libPaths(opt$pckDir)
# Print the values: order is the priority locations when looking for packages
message(paste0(.libPaths(),sep = "\n"))


if(!all(opt$biocPackages%in%installed.packages())){

    if(!"BiocManager"%in%installed.packages()){
        install.packages("BiocManager", lib = opt$pckDir)
    }
    

    update.packages(instlib = opt$pckDir,ask = F)
    

    BiocManager::install(
        opt$biocPackages, lib = opt$pckDir,
        ask = F,update = T
    )

}


if(!all(opt$pckDir%in%.libPaths())){
    print("Setting lib paths...")
    .libPaths(opt$pckDir)
}
library(ggplot2)
library(ggrepel)
install.packages("tidyselect")
install.packages("readr")
library(TCGAbiolinks) #Takes a bit to load
library(DESeq2) #Should produces lots of "warnings", but no errors. 
library(readr)



query1 <- GDCquery(
    project = opt$project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    access = "open"
)


samDf <- query1$results[[1]]

table(samDf$sample_type)


group1 <- read.delim("HighMutationGdc.tsv",fill=T,header = 1)
group2 <- read.delim("NoMutationGcd.tsv",fill=T,header = 1)

set.seed(6481)  # Setting seed for reproducibility
if (nrow(group2) > 50) {
  group2 <- group2[sample(nrow(group2), 50), ]
}

kk <- samDf[samDf$sample.submitter_id %in% group1$Sample.ID,]
ak <- samDf[samDf$sample.submitter_id %in% group2$Sample.ID,]


samDf2 <- rbind(kk,ak)



write.csv(samDf2, file = "project1.csv")
samDf2


sum(group1$Sample.ID%in%group2$Sample.ID) #Should be 0
sum(group2$Sample.ID%in%group1$Sample.ID) #Should be 0




nrow(samDf2)             #Should be >0 
length(group1$Sample.ID) #Should be >0
length(group2$Sample.ID) #Should be >0


sum(samDf2$sample.submitter_id%in%group1$Sample.ID) #Should be >=6 and not >>>50
sum(samDf2$sample.submitter_id%in%group2$Sample.ID) #Should be >=6 and not >>>50
length(unique(samDf2$cases)) 


query2 <- GDCquery(
    project = opt$project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    access = "open",
    barcode = samDf2$cases
)


if(!file.exists(opt$gdcPath)){
    warning(
        "GDC download directory (",opt$gdcPath,") not found...\n",
        "Creating directory to store GDC data in the wd:\n",
        opt$backupGdcPath
    )
    opt$gdcPath <- opt$backupGdcPath
    dir.create(opt$gdcPath,recursive = T,showWarnings = T)
}


GDCdownload(
    query = query2,
    method = "api",
    files.per.chunk = 10,
    directory = opt$gdcPath
)



dds <- GDCprepare(query = query2,directory = opt$gdcPath)



class(dds) #Should be RangedSummarizedExperiment and SummarizedExperiment

## Check the clinical data
dim(colData(dds))   #Should have a number of rows equal to the number of samples
dim(as.data.frame(rowRanges(dds))) #Should have ~60660 rows (number of characterized loci [~genes])

## Check the sample data section
assays(dds) # First of the downloaded data sets should be "unstranded" (what DESeq2 uses)
dim(assays(dds)[[1]]) #Should be the number of human genes by the number of samples

## Check for multiple observations / undesired tissues, etc
# The clinical data can also be directly accessed with $
if(sum(duplicated(dds$sample_submitter_id))>0){
    warning("Some IDs present more than once. Consider subsetting if this was unintentional.")
}
if(length(unique(dds$sample_type))>1){
    warning("More than one type of sample tissue type present. Consider subsetting if this was unintentional.")
}



query_mutations <- GDCquery(
    project = opt$project,
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
)
GDCdownload(
    query = query_mutations,
    method = "api",
    files.per.chunk = 10,
    directory = opt$gdcPath
)
dds2 <- GDCprepare(
    query = query_mutations, 
    directory = opt$gdcPath
)


##### Extract patient IDs ####
dds2$patientID <- gsub(
    pattern = "^(.*?)-(.*?)-(.*?)-.*", #Capture the first three identifiers delimited with hyphens
    "\\1-\\2-\\3",
    dds2$Tumor_Sample_Barcode
)


##### Subset to rows of interest to analysis #####
dds2     <- dds2[dds2$patientID%in%dds$patient,]
dds2_APC <- dds2[dds2$Hugo_Symbol=="APC",]


##### Compare to analysis object #####
dds$patient%in%dds2_APC$patientID #My cases are correctly ordered because of how I subset samDf.




####.......................................####
#---------------------------------------------#
#### Filter the samples and loci analyzed #####
#---------------------------------------------#

###### Notes on exploring the summarized experiment object ####
## Reading:
# https://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html

## Visualize **some** of the clinical data
# View(as.data.frame(colData(dds)))


###### Identify the groupings of samples ######
## Create new "columns" in the sample data based on grouping
dds$group1 <- dds$sample_submitter_id%in%group1$Sample.ID
dds$group2 <- dds$sample_submitter_id%in%group2$Sample.ID

## These should be the opposite of each other if barcode based subsetting was successful
if(!all(dds$group1==!dds$group2)){
    stop("Your groupings are not mutually exclusive")
}


dds$comp <- factor(dds$group2,levels = c("FALSE","TRUE"))



dds <- DESeqDataSet(dds, design = ~ comp)


dds <- estimateSizeFactors(dds)


# Step 1: Remove loci with 90% or more zero raw reads in either group
remove_loci_with_zero_reads <- function(data, threshold=0.9) {
  zero_reads_proportion <- rowMeans(data == 0)
  to_remove <- zero_reads_proportion >= threshold
  return(data[!to_remove, ])
}

# Step 2: Remove samples with extremely different raw read depths
remove_samples_with_different_depths <- function(data, threshold=3) {
  numeric_data <- data[, sapply(data, is.numeric)]  # Filter out non-numeric columns
  depths <- colSums(numeric_data) # Total reads per sample
  mean_depth <- mean(depths)
  std_depth <- sd(depths)
  cat("Mean depth:", mean_depth, "\n")
  cat("Standard deviation of depth:", std_depth, "\n")
  to_remove <- abs(depths - mean_depth) >= threshold * std_depth
  print(to_remove) # Print which samples are marked for removal
  return(data[, !to_remove])
}



# Step 3: Remove samples with extremely different patterns of normalized read depth
remove_samples_with_different_patterns <- function(data, threshold=0.95, num_components=2) {
  variance_per_locus <- apply(data, 1, var)
  most_variable_loci_indices <- tail(order(variance_per_locus), round(threshold * length(variance_per_locus)))
  pca_data <- prcomp(data[most_variable_loci_indices, ])
  component_scores <- pca_data$x[, 1:num_components]
  distances <- apply(component_scores, 1, function(x) sqrt(sum((x - mean(x))^2)))
  to_remove <- distances >= threshold * mean(distances)
  return(data[, !to_remove])
}


samDf9 <- remove_loci_with_zero_reads(group1)
samDf8 <- remove_samples_with_different_depths(group1)
samDf7 <- remove_samples_with_different_patterns(group1)


install.packages("ggplot2")
install.packages("prcomp")

library(ggplot2)
library(prcomp)


pca_result <- prcomp()



# My hidden version: assumes dds exists
source("~/bio321g/DiffExp/nolanHiddenFilteringScript.R")


## Load my filtered dds object to continue forward with example:
download.file(destfile = "dds.zip","https://utexas.box.com/shared/static/7bwh3xdp6qjqphhg0ufhb6aho2qqsqxy.zip")
unzip("dds.zip")
dds_Nolans <- readRDS("dds.RDS")



dds <- DESeqDataSet(dds, design = ~ comp)


dds <- estimateSizeFactors(dds)

## Run the Deseq2 package analysis
dds <- DESeq(dds)
res <- results(dds) # Organizes the results

##### Accumulate data into a data.frame #####
## Add in gene data from the rowRanges section of the SummarizedExperiment object

colsAdded <- !colnames(as.data.frame(rowRanges(dds)))%in%colnames(as.data.frame(res))
resOutput <- cbind(
    as.data.frame(res),
    as.data.frame(rowRanges(dds))[,colsAdded]
)


source("~/bio321g/DiffExp/nolanHiddenDESeq2VisScript.R")


