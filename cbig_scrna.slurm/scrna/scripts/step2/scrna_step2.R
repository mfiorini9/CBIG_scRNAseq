#!/usr/bin/env Rscript

## This code is used to to generate a unified Seurat object for all samples in the study

## Initiate progress tracker
stepp0="Generation of unified Seurat object"
cat("##########################################################################\n")
start_time0 <- Sys.time()
cat(stepp0,"has commenced.\n")
cat("##########################################################################\n")

################################################################################
# Load R libraries and job parameters
################################################################################

## Initiate progress tracker 
stepp="Loading R libraries and configuring job parameters"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

## Load job parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1] 
#output_dir="/home/fiorini9/scratch/CBIG_smajic/SOP"

r_lib_path=args[2] 
#r_lib_path="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R"

source(paste(output_dir,'/job_info/parameters/parameters.txt',sep="")) 

## Load R libraries
.libPaths(r_lib_path)
packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','Matrix', 'ggpubr', 'scCustomize', 'R.utils')
invisible(lapply(packages, library, character.only = TRUE))

## Detect number of available cores
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl) 

## Report step progress
cat(stepp,"has completed. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

################################################################################
# Load author-provided metadata to Seurat object
################################################################################

## Import author-provided metadata
meta <- read.delim(par_meta_data, sep = ",")

## Rename select columns according to standard CBIG nomenclature
colnames(meta)[colnames(meta) == par_barcodes] <- "barcodes_CBIG"
colnames(meta)[colnames(meta) == par_sample_id] <- "sample_id_CBIG"
colnames(meta)[colnames(meta) == par_subject_id] <- "subect_id_CBIG"
colnames(meta)[colnames(meta) == par_cell_type_main] <- "cell_type_main_CBIG"
colnames(meta)[colnames(meta) == par_disease_status_main] <- "disease_status_main_CBIG"

## Reorder column names
cbig_columns <- c("barcodes_CBIG", "sample_id_CBIG", "subect_id_CBIG", 
                  "cell_type_main_CBIG", "disease_status_main_CBIG")
remaining_columns <- setdiff(colnames(meta), cbig_columns)
meta <- meta[, c(cbig_columns, remaining_columns)]

################################################################################
# Generate individual Seurat objects for each sample in the study
################################################################################

## Initiate progress tracker 
stepp="Generation of individual Seurat objects"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

## Identify available feature-barcode expression matrices produced by CellRanger
list<-dir(path = paste(output_dir, "/step1",sep=""),full.names = TRUE)
sample_name<-dir(path = paste(output_dir, "/step1",sep=""))


## Create Seurat object for each sample
seu_list<-foreach (i=1:length(sample_name)) %do% {    
    datadirs <- file.path(list[i],   "ouput_folder","outs","raw_feature_bc_matrix")
    names(datadirs)=sample_name[i]
    sparse_matrix <- Seurat::Read10X(data.dir = datadirs)
    seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=0, min.features= 1, project = sample_name[i])
    nam <- paste("seurat_object", sample_name[i], sep = ".")
    assign(nam, seurat_object)
    seu<-get(nam)
    seu <- subset(seu, cells = meta$barcodes_CBIG)
    seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
    seu <- ScaleData(seu)
}

## Report step progress
cat(stepp,"has completed. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

################################################################################
# Generate a unified Seurat object for all samples in the study
################################################################################

## Initiate progress tracker 
stepp="Generation of unified Seurat object"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

## Create unified Seurat object 
seu_int <- Merge_Seurat_List(
    list_seurat = seu_list,
    add.cell.ids = NULL,
    merge.data = TRUE,
    project = "MergeSeurat"
)
seurat_list <- SplitObject(seu_int, split.by = "orig.ident")
first_obj <- seurat_list[[1]]  
remaining_objs <- seurat_list[-1]  
seu_int <- merge(x = first_obj, y = remaining_objs)
seu_int[["RNA"]] <- JoinLayers(seu_int[["RNA"]]) 

## Report step progress
cat(stepp,"has completed. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

################################################################################
# Add author-provided metadata to Seurat object
################################################################################

## Initiate progress tracker 
stepp="Addition of author-provided metadata"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

## Remove possible column duplicates
overlap_columns <- intersect(colnames(seu_int@meta.data), colnames(meta))
meta <- meta[, !colnames(meta) %in% overlap_columns]

## Rename column names in metadata object
colnames(seu_int@meta.data) <- paste0(colnames(seu_int@meta.data), "_CBIG")
head(seu_int@meta.data)

## Make sure barcodes match between author-provided metadata and Seurat object
rownames(seu_int@meta.data)
head(meta)
meta$barcodes_CBIG <- paste0("_", meta$barcodes_CBIG)

## check if meta barcodes are in the Seurat object.
meta$barcodes_CBIG %in% rownames(seu_int@meta.data)

## Set rownames of metdata as barcodes
rownames(meta) <- meta$barcodes_CBIG
temp_meta <- meta[rownames(seu_int@meta.data),]

## Add metadata to Seurat
seu_int <- AddMetaData(seu_int, temp_meta)
head(seu_int@meta.data)

## Report step progress
cat(stepp,"has completed. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

################################################################################
# Export final data files
################################################################################

## Initiate progress tracker 
stepp="Exporting final data objects"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

## Create output directory
dir.create(paste0(output_dir,'/final_outs'), recursive = FALSE, showWarnings = TRUE)

## Extract metadata and export
meta_new <- data.frame(seu_int@meta.data)
write.csv(meta_new, file= paste0(output_dir,'/final_outs/metadata.csv'))

## Save Seurat RDS file
saveRDS(seu_int, paste0(output_dir,'/final_outs/seurat.RDS'),compress=TRUE)

## Export barcodes.tsv file
barcodes <- colnames(seu_int)
write.table(barcodes, file = paste0(output_dir,'/final_outs/barcodes.tsv'), quote = FALSE, row.names = FALSE, col.names = FALSE)

## Export features.tsv file
features <- rownames(seu_int)
features <- data.frame(
  Ensemble = features,
  Feature = features,
  Expression_Type = rep("Gene Expression", length(features))
)
write.table(features, file = paste0(output_dir,'/final_outs/features.tsv'), quote = FALSE, row.names = FALSE, col.names = FALSE)

## Export matrix.mtx
matrix_data <- GetAssayData(object = seu_int, assay = "RNA", layer = "counts")
writeMM(matrix_data, paste0(output_dir,'/final_outs/matrix.mtx'))

## Compress files
gzip(paste0(output_dir,'/final_outs/barcodes.tsv'), remove = TRUE)  
gzip(paste0(output_dir,'/final_outs/features.tsv'), remove = TRUE)
gzip(paste0(output_dir,'/final_outs/matrix.mtx'), remove = TRUE)

## Test if the output files can be used to creat new seurat object
sparse_matrix <- Seurat::Read10X(data.dir = paste0(output_dir,'/final_outs'))
seu_int <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=0, min.features= 1, project = "merge")

## Report step progress
cat(stepp,"has completed. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

################################################################################
# End
################################################################################

## Print final progress report
cat("##########################################################################\n")
cat(stepp0,"successfully completed. Total time:",as.numeric (Sys.time() - start_time0, units = "mins"),"minutes\n")
cat("##########################################################################\n")
