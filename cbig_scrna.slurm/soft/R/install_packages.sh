# Load R module
module load StdEnv/2023
module load r/4.4.0

# Designated directory for R packages 
R_PATH=/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R
mkdir -p $R_PATH

# Install R packages 
Rscript ./scrnabox.slurm/soft/R/install_packages.R $R_PATH
