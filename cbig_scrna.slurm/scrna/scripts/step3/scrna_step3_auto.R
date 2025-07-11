## Load parameters
args = commandArgs(trailingOnly=TRUE)
output_dir=args[1]
r_lib_path=args[2]

## Load libraries
.libPaths(r_lib_path)

## Load parameters
source(paste(output_dir,'/job_info/parameters/parameters.txt',sep=""))

## Automated library prep and do not rename samples
if (tolower(par_automated_library_prep_spaceranger)=='yes' & (tolower(par_rename_samples_spaceranger)=='no')) {
    dir.create(paste(output_dir, "/samples_info_spaceranger", sep = ""),showWarnings = FALSE)
    new_dir <- paste(output_dir, "/samples_info_spaceranger/", sep = "")

for(i in par_sample_names_spaceranger) {
    #create directory
    dir.create(paste(new_dir,i, sep = ""),showWarnings = FALSE)
    
    #write the csv file 
    fastqs <- par_fastq_directory_spaceranger
    sample <- i
    library_type <- "Gene Expression"
    df <- data.frame(fastqs, sample,library_type)
    write.table(df, file = paste(new_dir,i,"/library.csv", sep=""),sep=",", row.names = FALSE, quote = FALSE)  
}
}

## Automated library prep and rename samples
if (tolower(par_automated_library_prep_spaceranger)=='yes' & (tolower(par_rename_samples_spaceranger)=='yes')) {
    dir.create(paste(output_dir, "/samples_info_spaceranger", sep = ""),showWarnings = FALSE)
    new_dir <- paste(output_dir, "/samples_info_spaceranger/", sep = "")
for(i in par_sample_names_spaceranger) {
    # Parse new sample names
    old_names <- par_sample_names_spaceranger
    new_names <- par_new_sample_names_spaceranger
    names_frame <- data.frame(old_names,new_names )
    new_i <- names_frame$new_names[old_names == i]    
    
    # Create directory
    dir.create(paste(new_dir,new_i, sep = ""),showWarnings = FALSE)
    
    # Write the csv file 
    fastqs <- par_fastq_directory_spaceranger
    sample <- i
    library_type <- "Gene Expression"
    df <- data.frame(fastqs, sample,library_type)
    write.table(df, file = paste(new_dir,new_i,"/library.csv", sep=""),sep=",", row.names = FALSE, quote = FALSE)  
}
}

## Do not perform automated library prep
if (tolower(par_automated_library_prep_spaceranger)=='no') {
print("Skipping automated library prep for CellRanger")  
}
