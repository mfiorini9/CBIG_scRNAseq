import sys
import time
import os
import ast
import re

################################################################################
# Load Python libraries and job parameters
################################################################################

# Initiate progress tracker
stepp = "Loading Python libraries and configuring job parameters"
print("#####################################")
print(f"{stepp} started")
start_time = time.time()

# Load job parameters from command line
args = sys.argv[1:]

# Assign arguments
output_dir = args[0]
#output_dir = "/home/fiorini9/scratch/CBIG_smajic/pipeline_test/out_dir_spaceranger"

#pythonenv_path = args[1]
#pythonenv_path = "/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA/lib/python3.8/site-packages"
#pythonenv_path = "/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_Py_ENV"

# Load parameters from text file
parameters_path = os.path.join(output_dir, 'job_info', 'parameters', 'parameters.txt')

#with open(parameters_path, 'r') as f:
#    exec(f.read())

def parse_r_parameters(file_path):
    parameters = {}
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue  # Skip comments and blank lines
            
            # Split key=value
            if "=" in line:
                key, value = line.split("=", 1)
                key = key.strip()
                value = value.strip()
      
                # Convert TRUE/FALSE to Python
                value = value.replace("TRUE", "True").replace("FALSE", "False")
      
                # Convert c("a", "b", ...) to list
                if value.startswith("c(") and value.endswith(")"):
                    items = re.findall(r'"(.*?)"', value)
                    value = items  # becomes a Python list
      
                # Remove surrounding quotes from strings
                elif value.startswith('"') and value.endswith('"'):
                    value = value[1:-1]
      
                # Try to evaluate as int/float if applicable
                else:
                    try:
                        value = ast.literal_eval(value)
                    except Exception:
                        pass  # leave as string if not evaluable
      
                parameters[key] = value
      
    return parameters

params = parse_r_parameters(parameters_path)

## Import required libraries
#sys.path.insert(0, pythonenv_path)

import scanpy as sc
import squidpy as sq
from squidpy.im import ImageContainer
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt

# Report step progress
elapsed_time = (time.time() - start_time) / 60  # in minutes
print(f"{stepp} has completed. Total time: {elapsed_time:.2f} minutes")
print("#####################################")


################################################################################
# Load AnnData object and add image assays.
################################################################################
# Initiate progress tracker
stepp = "Loading Python libraries and configuring job parameters"
print("#####################################")
print(f"{stepp} started")
start_time = time.time()

## try again
# Define sample ID to barcode prefix mapping
#sample_to_prefix = {
#    "Test": "sample1_Test",
#    "Test2": "sample2_Test2"
#}

# Load your merged AnnData
adata_path = f"{output_dir}/final_outs/seurat.h5ad"
adata = sc.read_h5ad(adata_path)

unique_samples = set(adata.obs['orig.ident'])
sample_to_prefix = {sample: f"sample{idx+1}_{sample}" for idx, sample in enumerate(sorted(unique_samples))}

print(sample_to_prefix)

# Initialize spatial slot
adata.uns["spatial"] = {}

# Allocate empty spatial coordinate array
spatial_coords = np.zeros((adata.n_obs, 2))

# Loop over each sample
for sample_id, sample_prefix in sample_to_prefix.items():
    img_path = f"{output_dir}/final_outs/spatial.{sample_id}"
    
    # 1. Read spatial image
    img = ImageContainer(os.path.join(img_path, "tissue_hires_image.png"), layer="image")
    
    # 2. Read scalefactors
    with open(os.path.join(img_path, "scalefactors_json.json"), "r") as f:
        scalefactors = json.load(f)
    
    # 3. Read coordinates
    coords = pd.read_csv(os.path.join(img_path, "tissue_positions.csv"), header=0)
    coords.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    coords["barcode"] = sample_prefix + "_" + coords["barcode"].astype(str)
    coords = coords.set_index("barcode")
    
    # 4. Match cells in adata
    cell_mask = adata.obs_names.str.startswith(sample_prefix)
    adata_subset = adata[cell_mask].copy()
    coords_matched = coords.loc[adata_subset.obs_names]
    
    # 5. Assign coordinates to the full matrix
    #spatial_coords[cell_mask, :] = coords_matched[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
    spatial_coords[cell_mask, :] = coords_matched[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values
    
    # 6. Store image and scalefactors
    adata.uns["spatial"][sample_id] = {
        #"images": {"hires": img},
        "images": {"hires": img["image"].values.squeeze()},
        "scalefactors": scalefactors,
        "metadata": {}
    }

# 7. Save final spatial coordinates
adata.obsm["spatial"] = spatial_coords

######################################################
## TESTS
######################################################
## FIND OUT WHY IT PRINTS IT TWICE
#adata_test2 = adata[adata.obs_names.str.startswith("sample1_Test")].copy()
#adata_test2.uns["spatial"] = {
#    "Test": adata_test2.uns["spatial"]["Test"]
#}

# Plot only for Test2
#sq.pl.spatial_scatter(
#    adata_test2,
#   color=["Sox8", "orig.ident"],
#    library_key="Test",
#    use_raw=False
#)

# Save to PDF
#fig = plt.gcf()
#fig.savefig("/home/fiorini9/scratch/CBIG_smajic/pipeline_test/out_dir_spaceranger/Rorb_spatial_plot.pdf", bbox_inches='tight')
#plt.close(fig)

############################################################################################################

#########################
## Save the object test
#########################
adata

raw = adata.raw.to_adata()

if "_index" in raw.var.columns:
    raw.var = raw.var.drop(columns="_index")

adata.raw = raw

## Delete the old file
file_path = f"{output_dir}/final_outs/seurat.h5ad"

if os.path.exists(file_path):
    os.remove(file_path)
    print(f"{file_path} deleted successfully.")
else:
    print(f"{file_path} does not exist.")

## Delete existing file in case it is present
file_path = f"{output_dir}/final_outs/Scanpy.h5ad"

if os.path.exists(file_path):
    os.remove(file_path)
    print(f"{file_path} deleted successfully.")
else:
    print(f"{file_path} does not exist.")

## Save new file
adata.write(f"{output_dir}/final_outs/Scanpy.h5ad")

## TEST
#adata = sc.read_h5ad('/home/fiorini9/scratch/CBIG_smajic/pipeline_test/out_dir_spaceranger/final_outs/output_file.h5ad')

################################################################################################### END HERE. 



#raw_df = pd.DataFrame(
#    adata.raw.X[:5].toarray() if hasattr(adata.raw.X, "toarray") else adata.raw.X[:5],
#    index=adata.raw.obs_names[:5],
#    columns=adata.raw.var_names
#)

#raw_df.head()


#X_df = pd.DataFrame(
#    adata.X[:5].toarray() if hasattr(adata.X, "toarray") else adata.X[:5],
#    index=adata.obs_names[:5],
#    columns=adata.var_names
#)

#X_df.head()

################################################################################################### END HERE.




#if "_index" in adata.raw.var.columns:
#    adata.raw.var = adata.raw.var.drop(columns="_index")

#adata.write("/home/fiorini9/scratch/CBIG_smajic/pipeline_test/out_dir_spaceranger/final_outs/output_file.h5ad")
#seurat.h5ad





## Lets try dropping ray to see if this fixes the duplicate plots ssue. 
#adata.raw = None

#adata.obs

## Plot to test
#sq.pl.spatial_scatter(adata, color=["Sox8","orig.ident"], library_key = 'Test', use_raw = False)
#plt.savefig("/home/fiorini9/scratch/CBIG_smajic/pipeline_test/out_dir_spaceranger/Rorb_spatial_plot.pdf", bbox_inches='tight')  # Save as PDF

# Temp
#adata_test = adata[adata.obs["orig.ident"] == "Test"]
#sq.pl.spatial_scatter(
#    adata_test,
#    color=["Sox8","orig.ident"],            # or ["Sox8", "orig.ident"]
#    library_key="Test",        # matches key in adata.uns["spatial"]
#    use_raw=False
#)
#plt.savefig("/home/fiorini9/scratch/CBIG_smajic/pipeline_test/out_dir_spaceranger/Rorb_spatial_plot.pdf", bbox_inches='tight')  # Save as PDF

#adata_test = adata[adata.obs["orig.ident"] == "Test2"]
#sq.pl.spatial_scatter(
#    adata_test,
#    color=["Sox8","orig.ident"],            # or ["Sox8", "orig.ident"]
#    library_key="Test2",        # matches key in adata.uns["spatial"]
#    use_raw=False,
#   show=False
#)
#fig = plt.gcf()
#fig.savefig(
#    "/home/fiorini9/scratch/CBIG_smajic/pipeline_test/out_dir_spaceranger/Rorb_spatial_plot.pdf",
#    bbox_inches="tight"
#)

## FIND OUT WHY IT PRINTS IT TWICE
#adata_test2 = adata_test[adata_test.obs_names.str.startswith("sample2_Test2")].copy()
#adata_test2.uns["spatial"] = {
#    "Test2": adata_test2.uns["spatial"]["Test2"]
#}

# Plot only for Test2
#sq.pl.spatial_scatter(
#    adata_test2,
#    color=["Sox8", "orig.ident"],
#    library_key="Test2",
#    use_raw=False
#)

# Save to PDF
#fig = plt.gcf()
#fig.savefig("/home/fiorini9/scratch/CBIG_smajic/pipeline_test/out_dir_spaceranger/Rorb_spatial_plot.pdf", bbox_inches='tight')
#plt.close(fig)
