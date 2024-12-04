#!/usr/bin/env Rscript

# Get the library path from command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the library path is provided
if (length(args) < 1) {
  stop("Please provide the R library path as the first argument.")
}

R_LIB_PATH <- args[1]  # First argument is the R library path

# Set the library path
.libPaths(R_LIB_PATH)
cat("Current R library paths:", .libPaths(), "\n")

# Install required packages
install.packages("Matrix", lib=R_LIB_PATH)
install.packages("SeuratObject", lib=R_LIB_PATH)
install.packages("Seurat", lib=R_LIB_PATH)
install.packages("ggplot2", lib=R_LIB_PATH)
install.packages("dplyr", lib=R_LIB_PATH)
install.packages("foreach", lib=R_LIB_PATH)
install.packages("doParallel", lib=R_LIB_PATH)
install.packages("ggpubr", lib=R_LIB_PATH)
install.packages("R.utils", lib=R_LIB_PATH)
