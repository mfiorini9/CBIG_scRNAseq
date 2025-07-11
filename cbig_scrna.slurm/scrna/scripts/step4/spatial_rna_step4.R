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
#output_dir="/home/fiorini9/scratch/CBIG_smajic/pipeline_test/out_dir_spaceranger"

r_lib_path=args[2] 
#r_lib_path="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R"

source(paste(output_dir,'/job_info/parameters/parameters.txt',sep="")) 

## Load R libraries
#.libPaths(r_lib_path)
#packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','Matrix', 'ggpubr', 'scCustomize', 'R.utils', 'tiff')
#invisible(lapply(packages, library, character.only = TRUE))

library(SeuratDisk, lib=r_lib_path)
library(Seurat, lib=r_lib_path)
library(dplyr, lib=r_lib_path)
library(foreach, lib=r_lib_path)
library(doParallel, lib=r_lib_path)
library(Matrix, lib=r_lib_path)
library(scCustomize, lib=r_lib_path)
library(R.utils, lib=r_lib_path)
library(tiff, lib=r_lib_path)

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

## add condition to run only if par is true.
#par_include_metadata_spatial = "no"         ## TEMP
#par_meta_data_spatial   = "/home/fiorini9/scratch/CBIG_smajic/pipeline_test/out_dir_spaceranger/metadata_temp.csv"   ## TEMP

if (tolower(par_include_metadata_spatial)=="yes") {
## Import author-provided metadata
meta <- read.delim(par_meta_data_spatial, sep = ",")

## Rename select columns according to standard CBIG nomenclature
colnames(meta)[colnames(meta) == par_barcodes_spatial] <- "barcodes_CBIG"
colnames(meta)[colnames(meta) == par_sample_id_spatial] <- "sample_id_CBIG"
colnames(meta)[colnames(meta) == par_subject_id_spatial] <- "subect_id_CBIG"
colnames(meta)[colnames(meta) == par_cell_type_main_spatial] <- "cell_type_main_CBIG"
colnames(meta)[colnames(meta) == par_disease_status_main_spatial] <- "disease_status_main_CBIG"
colnames(meta)[colnames(meta) == par_tissue_spatial] <- "tissue_CBIG"
colnames(meta)[colnames(meta) == par_tissue_region_spatial] <- "tissue_region_CBIG"

## Reorder column names
cbig_columns <- c("barcodes_CBIG", "sample_id_CBIG", "subect_id_CBIG", "tissue_CBIG", "tissue_region_CBIG",
                  "cell_type_main_CBIG", "disease_status_main_CBIG")
remaining_columns <- setdiff(colnames(meta), cbig_columns)
meta <- meta[, c(cbig_columns, remaining_columns)]
}

################################################################################
# Generate individual Seurat objects for each sample in the study
################################################################################

## Initiate progress tracker 
stepp="Generation of individual Seurat objects"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

## Identify available feature-barcode expression matrices produced by CellRanger
list<-dir(path = paste(output_dir, "/step3",sep=""),full.names = TRUE)
sample_name<-dir(path = paste(output_dir, "/step3",sep=""))

## Create output directory
dir.create(paste0(output_dir,'/final_outs'), recursive = FALSE, showWarnings = TRUE)


## Create Seurat object for each sample
seu_list<-foreach (i=1:length(sample_name)) %do% { 
    datadirs <- file.path(list[i],   "ouput_folder","outs","filtered_feature_bc_matrix") ## BEGIN WITH FILETRED BUT EVENTUALLY SWITCH TO RAW: raw_feature_bc_matrix 
    names(datadirs)=sample_name[i]
    sparse_matrix = Seurat::Read10X(data.dir = datadirs) 
    seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=0, min.features= 1, assay = 'Spatial', project = sample_name[i]) 
    imgpath = file.path(list[i],   "ouput_folder","outs","spatial")

    ## Establish destination directory for spatial data
    dest_dir <- file.path(paste0(output_dir,'/final_outs/spatial.', sample_name[i] ))

    # Create destination directory if it doesn't exist
    if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
    }
    
    ## copy contents of the source spatial directory
    files_to_copy <- list.files(imgpath, full.names = TRUE, recursive = TRUE)

    # Copy each file preserving relative path
    for (j in seq_along(files_to_copy)) {
    file.copy(files_to_copy[j], file.path(dest_dir), overwrite = TRUE)
    }
    
    img = Seurat::Read10X_Image(image.dir = imgpath) 
    Seurat::DefaultAssay(object = img) <- 'Spatial' 
    img@boundaries$centroids@cells <- paste0(sample_name[i], '_',img@boundaries$centroids@cells )
    img@key <- paste0(sample_name[i], '_image' )
    img = img[colnames(x = seurat_object)] 
    seurat_object[['image']] = img  
    
    nam <- paste("seurat_object", sample_name[i], sep = ".")
    assign(nam, seurat_object)
    seu<-get(nam)

    if (tolower(par_include_metadata_spatial)=="yes") {
    seu <- subset(seu, cells = meta$barcodes_CBIG)
    }
    
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
    add.cell.ids = paste0("sample", seq_along(seu_list)),
    merge.data = TRUE,
    project = "MergeSeurat"
)

unique(seu_int@meta.data$orig.ident)

seurat_list <- SplitObject(seu_int, split.by = "orig.ident")
first_obj <- seurat_list[[1]]  
remaining_objs <- seurat_list[-1]  
seu_int <- merge(x = first_obj, y = remaining_objs)
seu_int[["Spatial"]] <- JoinLayers(seu_int[["Spatial"]]) 

## Report step progress
cat(stepp,"has completed. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")

#metadata_temp <- seu_int@meta.data
#write.csv(metadata_temp, '/home/fiorini9/scratch/CBIG_smajic/pipeline_test/out_dir_spaceranger/metadata_temp.csv' )

################################################################################
# Add author-provided metadata to Seurat object
################################################################################

## Initiate progress tracker 
stepp="Addition of author-provided metadata"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

if (tolower(par_include_metadata_spatial)=="yes") {

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

    ## maybe remove all barcodes that are present in metadata
    ## FILL THIS IN. 

    ## Report step progress
    cat(stepp,"has completed. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
    cat("#####################################\n")

}

################################################################################
# Export final data files
################################################################################

## Initiate progress tracker 
stepp="Exporting final data objects"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

## Extract metadata and export
meta_new <- data.frame(seu_int@meta.data)
write.csv(meta_new, file= paste0(output_dir,'/final_outs/metadata.csv'))

## rename the image files
current_images <- names(seu_int@images)
orig_idents <- unique(seu_int@meta.data$orig.ident)
if (length(current_images) != length(orig_idents)) {
  stop("Number of images does not match number of unique orig.ident values.")
}
names(seu_int@images) <- paste0('image.',orig_idents)

## Save Seurat RDS file
saveRDS(seu_int, paste0(output_dir,'/final_outs/seurat.RDS'),compress=TRUE)

## Export barcodes.tsv file
barcodes <- colnames(seu_int)
write.table(barcodes, file = paste0(output_dir,'/final_outs/barcodes.tsv'), quote = FALSE, row.names = FALSE, col.names = FALSE)

## Export features.tsv file
features <- rownames(seu_int)
features <- data.frame(
  Feature = features,
  Expression_Type = rep("Gene Expression", length(features))
)

ensembl <- read.delim(paste0(list[1],   "/ouput_folder/outs/filtered_feature_bc_matrix/features.tsv.gz"), sep = '\t', header = FALSE)
ensembl <- ensembl[,1:2]
colnames(ensembl) <- c("Ensembl", "Feature")
ensembl$Feature <- rownames(seu_int)

features2 <- merge(features, ensembl, by = "Feature")
features2_ordered <- features2[match(features$Feature, features2$Feature), ]
features2_ordered <- features2_ordered %>% dplyr::select(Ensembl, Feature, Expression_Type)

write.table(features2_ordered, file = paste0(output_dir,'/final_outs/features.tsv'), quote = FALSE, row.names = FALSE, col.names = FALSE)

## Export matrix.mtx
matrix_data <- GetAssayData(object = seu_int, assay = "Spatial", layer = "counts")
writeMM(matrix_data, paste0(output_dir,'/final_outs/matrix.mtx'))

## Compress files
gzip(paste0(output_dir,'/final_outs/barcodes.tsv'), remove = TRUE)  
gzip(paste0(output_dir,'/final_outs/features.tsv'), remove = TRUE)
gzip(paste0(output_dir,'/final_outs/matrix.mtx'), remove = TRUE)

### Print the .tiff files
# Loop over all images in the Seurat object
for (img_name in names(seu_int@images)) {
    print(img_name)
    # Extract the image data
    img_data <- seu_int@images[[img_name]]@image
    # Create a temporary variable
    img_name_short <- sub("^image\\.", "", img_name)
    # Construct output filename
    file_name <- paste0(output_dir,'/final_outs/spatial.',img_name_short, '/', img_name, ".tif")
    # Write to file
    writeTIFF(img_data, file_name)
    
    cat("Saved:", file_name, "\n")
}


## Test if the output files can be used to creat new seurat object
#sparse_matrix <- Seurat::Read10X(data.dir = paste0(output_dir,'/final_outs'))
#seu_int_test <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=0, min.features= 1, project = "merge")

## Test if we can load in the tif thing. 
#imgpath = "/home/fiorini9/scratch/CBIG_smajic/pipeline_test/out_dir_spaceranger/final_outs/spatial.Test"

#img = Seurat::Read10X_Image(image.dir = imgpath) 
#Seurat::DefaultAssay(object = img) <- 'Spatial' 
#img@boundaries$centroids@cells <- paste0('sample2_', sample_name[i], '_',img@boundaries$centroids@cells )
#img@key <- paste0(sample_name[i], '_image' )

#colnames(seu_int_test)
#rownames(img)

#img = img[colnames(x = seu_int_test)] 
#seu_int_test[['image']] = img  

## Report step progress
cat(stepp,"has completed. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")


################################################################################
# Convert to Scanpy test
################################################################################
## Initiate progress tracker 
stepp="Create AnnData objects"
cat("#####################################\n")
cat(stepp, "started\n")
start_time <- Sys.time()

############ TEMP TEST
#seu_int <- readRDS('/home/fiorini9/scratch/CBIG_smajic/pipeline_test/out_dir_spaceranger/final_outs/seurat.RDS')
#str(seu_int)

# Extract counts, data, and scale.data matrices from the original assay
counts_mat <- as.matrix(seu_int@assays$Spatial@layers$counts)
#data_mat <- as.matrix(seu_int@assays$Spatial@layers$data)
scale_data_mat <- as.matrix(seu_int@assays$Spatial@layers$scale.data)

# Set row and column names for these matrices
rownames(counts_mat) <- rownames(seu_int)
colnames(counts_mat) <- colnames(seu_int)

#rownames(data_mat) <- rownames(seu_int)
#colnames(data_mat) <- colnames(seu_int)

rownames(scale_data_mat) <- rownames(seu_int)
colnames(scale_data_mat) <- colnames(seu_int)

# Create new assay with counts, then manually add data and scale.data
assay.v5 <- CreateAssayObject(counts = counts_mat)
#assay.v5@data <- data_mat
assay.v5@scale.data <- scale_data_mat

# Create new Seurat object using this assay
seu_int_scanpy <- CreateSeuratObject(assay.v5)

# Rename assay and remove default RNA assay
seu_int_scanpy[["Spatial"]] <- seu_int_scanpy[["RNA"]]
DefaultAssay(seu_int_scanpy) <- "Spatial"
seu_int_scanpy[["RNA"]] <- NULL

# Add metadata
df_meta <- data.frame(seu_int@meta.data)
seu_int_scanpy <- AddMetaData(seu_int_scanpy, df_meta)

########################################## 

########################################## Original code -- works
#counts_mat <- as.matrix(seu_int@assays$Spatial@layers$counts)
#colnames(counts_mat) <- colnames(seu_int)
#rownames(counts_mat) <- rownames(seu_int)
#assay.v5 <- CreateAssayObject(counts = counts_mat, )

#seu_int_scanpy <- CreateSeuratObject(assay.v5)

#seu_int_scanpy[["Spatial"]] <- seu_int_scanpy[["RNA"]]  
#DefaultAssay(seu_int_scanpy) <- "Spatial"
#seu_int_scanpy[["RNA"]] <- NULL                         

#df_meta <-data.frame(seu_int@meta.data)
#seu_int_scanpy <- AddMetaData(seu_int_scanpy, df_meta)
##########################################


## Need to add the images -- turns out we need the tissue_lowres_image.png file ************

## create an info dataframe
idents <- unique(seu_int_scanpy@meta.data$orig.ident)
samples <- as.character(unique(seu_int_scanpy@active.ident))
info_df <- data.frame(idents_col = idents, samples_col = samples)

for (row in 1:nrow(info_df)){
    print(row)
    imgpath = paste0(output_dir,'/final_outs/spatial.', info_df[row,1])
    img = Seurat::Read10X_Image(image.dir = imgpath)
    Seurat::DefaultAssay(object = img) <- 'Spatial'    
    img@boundaries$centroids@cells <- paste0(info_df[row,2], '_',info_df[row,1], '_',img@boundaries$centroids@cells )
    img@key <- paste0(info_df[row,1], '_image' )
    img = img[colnames(x = seu_int_scanpy)]
    seu_int_scanpy[[paste0('image.',info_df[row,1])]]  = img
}

## Test
#for (img in names(seu_int_scanpy@images)) {
#  bbox <- seu_int_scanpy@images[[img]]@boundaries[[1]]@bbox
#  # Add dimnames if missing or incorrect
#  if (is.null(dimnames(bbox)) || !identical(dimnames(bbox), list(c("x", "y"), c("min", "max")))) {
#   dimnames(bbox) <- list(c("x", "y"), c("min", "max"))
#    seu_int_scanpy@images[[img]]@boundaries[[1]]@bbox <- bbox
#  }
#}

## TEMP
#SaveH5Seurat(seu_int_scanpy, filename = paste0(output_dir, "/final_outs/seurat.h5Seurat"))
#seu_test <- LoadH5Seurat(paste0(output_dir, "/final_outs/seurat.h5Seurat"))

SaveH5Seurat(seu_int_scanpy, 
        filename = paste0(output_dir,'/final_outs/seurat.h5Seurat'),
        overwrite = T,
        verbose = T
        )

## lets read in the H5 object to see what is stored in it. 
#seu_test <- LoadH5Seurat(paste0(output_dir, '/final_outs/seurat.h5Seurat'))

Convert(paste0(output_dir,'/final_outs/seurat.h5Seurat'), dest = "h5ad", assay = "Spatial")

## Report step progress
cat(stepp,"has completed. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
cat("#####################################\n")


####### TEST
#for (img in names(seu_int_scanpy@images)) {
#  bbox <- seu_int_scanpy@images[[img]]@boundaries[[1]]@bbox
#  
#  # Force correct dimnames
#  dimnames(bbox) <- list(c("x", "y"), c("min", "max"))
  
#  # Reassign to object
#  seu_int_scanpy@images[[img]]@boundaries[[1]]@bbox <- bbox
#}

#SaveH5Seurat(seu_int_scanpy, filename = paste0(output_dir, "/final_outs/seurat2.h5Seurat"))
#seu_test <- LoadH5Seurat(paste0(output_dir, "/final_outs/seurat2.h5Seurat"))

#Convert(paste0(output_dir,'/final_outs/seurat2.h5Seurat'), dest = "h5ad", assay = "Spatial")



################################################################################
# End
################################################################################

## Print final progress report
cat("##########################################################################\n")
cat(stepp0,"successfully completed. Total time:",as.numeric (Sys.time() - start_time0, units = "mins"),"minutes\n")
cat("##########################################################################\n")
