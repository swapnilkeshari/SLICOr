setwd("/ix/djishnu/Swapnil/CellOracle/primaryBCell")
rm(list=ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"))

install.packages("devtools",,repos = "http://cran.us.r-project.org")
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
library(monocle3)
library(cicero)

# print("Loaded the libraries")
# You can substitute the data path below to your scATAC-seq data.
data_folder <- "filtered_peak_bc_matrix"

# # Create a folder to save results
output_folder <- "cicero_output"
dir.create(output_folder)

# Read in matrix data using the Matrix package
indata <- Matrix::readMM(("peakPercell_matrix.mtx")) 
# Binarize the matrix
indata@x[indata@x > 0] <- 1
indata<- t(indata)

# Format cell info
cellinfo <- read.table("/ix/djishnu/Swapnil/CellOracle/primaryBCell/filtered_feature_bc_matrix/barcodes.tsv")
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

# Format peak info
peakinfo <- read.table(paste0("atac_peaks.bed"))
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# Make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
cell_metadata = cellinfo,
gene_metadata = peakinfo))


input_cds <- monocle3::detect_genes(input_cds)

print("Made CDS")
#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

# Filter cells by peak_count
# Please set an appropriate threshold values according to your data 
max_count <-  15000
min_count <- 2000
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) >= min_count] 
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) <= max_count] 

set.seed(2017)

input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

# Dimensional reduction with umap
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")
plot_cells(input_cds)
umap_coords <- reducedDims(input_cds)$UMAP
print("Ran UMAP")

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# Save cds object if you want
saveRDS(cicero_cds, paste0(output_folder, "/cicero_cds.Rds"))
print("Saved cds object")  

cicero_cds <- readRDS("/ix/djishnu/Swapnil/CellOracle/primaryBCell/cicero_output/cicero_cds.Rds")
chromosome_length <- read.table("./hg38_chromosome_length.txt")

# Run the main function
conns <- run_cicero(cicero_cds, chromosome_length) # Takes a few minutes to run

write.csv(x = conns, file ="cicero_connections.csv")

# Save results if you want
saveRDS(conns, paste0(output_folder, "/cicero_connections.Rds"))

# Check results
print(head(conns))
print("Ran cicero")
all_peaks <- row.names(exprs(input_cds))
write.csv(x = all_peaks, file = paste0(output_folder, "/all_peaks.csv"))
write.csv(x = all_peaks, file = "all_peaks.csv")
write.csv(x = conns, file = paste0(output_folder, "/cicero_connections.csv"))

print("Wrote results")