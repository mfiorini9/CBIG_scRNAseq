# Spatial RNAseq pipeline: 10X Genomics Visium HD
- [Installation](#installation)
- [Step 0: Pipeline Set Up](#step-0-pipeline-set-up)
- [Step 1: Run SpaceRanger](#step-1-run-spaceranger)
- [Step 2: Data Preprocessing and Standardization](#step-2-data-preprocessing-and-standardization)
- [Data upload to C-BIG](#step-2-data-preprocessing-and-standardization)


- - - -

This tutorial provides step-by-step instructions for preparing Visium HD spatial RNA sequencing data from 10X Genomics for upload to the C-BIG repository.

- - - -

## Installation
To use the C-BIG RNAseq proessing pipeline, the [CBIG_scRNAseq GitHub repository](https://github.com/mfiorini9/CBIG_scRNAseq) must be cloned onto your high performance computing (HPC) system. Administrators using the Digital Research Alliance of Canada HPC system can skip this step, as the pipeline is already installed on Narval at `/lustre07/scratch/fiorini9/CBIG_pipe_07_11_25/CBIG_scRNAseq/cbig_scrna.slurm`. Note that this is a temporary location. 

To clone the repository run the following code:

```
mkdir /path/to/pipeline
cd /path/to/pipeline
git clone https://github.com/mfiorini9/CBIG_scRNAseq
```

After running the above code, the following files will be available in `/path/to/pipeline`:

```
CBIG_scRNAseq
├── cbig_scrna.slurm
│   ├── launch
│   │   └── launch_cbig_scrna_slurm.sh
│   ├── launch_cbig_scrna.sh
│   ├── scrna
│   │   ├── configs
│   │   │   └── cbig_scrna_config.ini
│   │   ├── pars
│   │   │   └── parameters.txt
│   │   └── scripts
│   │       ├── step1
│   │       │   ├── create_cellranger_scrna_local.sh
│   │       │   ├── create_cellranger_scrna.sh
│   │       │   ├── scrna_step1_auto.R
│   │       │   └── slurm.template
│   │       ├── step2
│   │       │   ├── pipeline_step2.sh
│   │       │   └── scrna_step2.R
│   │       ├── step3
│   │       │   ├── create_spaceranger_scrna.sh
│   │       │   ├── scrna_step3_auto.R
│   │       │   └── slurm.template
│   │       └── step4
│   │           ├── pipeline_step4.sh
│   │           ├── spatial_rna_step4.py
│   │           └── spatial_rna_step4.R
│   ├── soft
│   │   └── R
│   │       ├── install_packages.R
│   │       ├── install_packages.sh
│   │       └── R.library.ini
│   └── tools
│       └── utils.sh
├── LICENSE
└── README.md
```

**Note:** The repository only needs to be cloned onto an HPC system once and can be used to process any number of independent datasets.

- - - -

## Step 0: Pipeline Set Up

After installation, a working directory and a pipeline instance specific to the dataset being processed must be created. Step 0 initiates the pipeline and generates all necessary parameter and configuration files.

To execute Step 0, run the following code:

```
mkdir /path/to/working/directory/Dataset1
cd /path/to/working/directory/Dataset1

PIPELINE_HOME=/path/to/pipeline/CBIG_scRNAseq/cbig_scrna.slurm
PIPELINE_PWD=/path/to/working/directory/Dataset1

cd $PIPELINE_PWD

bash $PIPELINE_HOME/launch_cbig_scrna.sh \
-d ${PIPELINE_PWD} \
--steps 0 
```

After running the Step 0, the following files will be available in `/path/to/working/directory/Dataset1`:

```
Dataset1
└── job_info
    ├── configs
    │   └── cbig_scrna_config.ini
    ├── logs
    │   └── step_3_spacerangerV1_Adult_Mouse_Brain_S5.log
    ├── parameters
    │   └── parameters.txt
    └── summary_report.txt
```

- The `configs/` directory contains a configuration file which allows users to specify their job allocations (memory, threads, and walltime) for each analytical step using the Slurm Workload Manager; 
- The `logs/` directory records the events of each analytical step; 
- The `parameters/` directory contains adjustable, step-specific text files which allow users to define the execution parameters for each analytical step. 

**Note:** The text files can be opened and modified through nano, vim, or a file manager system like cyberduck. 

- - - -

## Step 1: Run SpaceRangerer

In Step 1, spatial gene expression matrices are generated from 10X Genomics FASTQ files provided by the original investigators using the SpaceRanger _counts_ pipeline. 

The following parameters are adjustable for Step 1 of the spatial RNAseq Visium HD track (`/path/to/working/directory/Dataset1/job_info/parameters/parameters.txt`):

|Column name|Default|Description|
|:--|:--|:--|
|**par_automated_library_prep_spaceranger**|Yes|Whether or not to perform automated library prep. Alternatively, you may set this parameter to "no" and manually prepare the libraries.|
|**par_fastq_directory_spaceranger**|NULL|Path to directory containing the FASTQ files. This directory should only contain FASTQ files for the experiment.|
|**par_sample_names_spaceranger**|NULL|The sample names used to name the FASTQ files according to CellRanger nomeclature.|
|**par_rename_sample_spacerangers**|Yes|Whether or not you want to rename your samples. These names will be used to identify cells in the processed data objects.|
|**par_new_sample_names_spaceranger**|NULL|New sample names. Make sure they are defined in the same order as 'par_sample_names'.|
|**par_paired_end_seq_spaceranger**|Yes|Whether or not paired-end sequencing was performed.|
|**par_ref_dir_grch_spaceranger**|/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/10xGenomics/refdata-cellranger-GRCh38-3.0.0|Path to reference genome for FASTQ alignment. 10X Genomics reference genomes are available for download. For more information see the 10X Genomics documentation.|
|**par_fastqs_spaceranger**|NULL|Path to directory containing the FASTQ files. This directory should only contain FASTQ files for the experiment.|
|**par_image_spaceranger**|NULL|Path to the brightfield tissue image (.tif) for the corresponding tissue section.|
|**par_slide_spaceranger**|NULL|Identifier for the Visium slide used. This must match the slide ID printed on the physical slide and encoded in the .gpr file.|
|**par_slidefile_spaceranger**|TEMP|Path to the GenePix results (.gpr) file for the slide, which contains spatial barcode layout information. It's used to map tissue spots to spatial barcodes.|
|**par_area_spaceranger**|NULL|Indicates the capture area (i.e., section) of the slide being processed. Commonly labeled as A1, B1, C1, D1, etc. This must match the area used for imaging and sample preparation.|
|**par_localcores_spaceranger**|16|Number of CPU cores to use for the analysis.|
|**par_localmem_spaceranger=10**|10|Amount of RAM (in GB) to allocate to the job.|

Given that SpaceRanger runs a user interface, it is recommended to run Step 1 in a _screen_, which will allow the the task to keep running if the connection to the HPC is broken. 

To execute Step 1, run the following code:

```
screen -S pipeline

PIPELINE_HOME=/path/to/pipeline/CBIG_scRNAseq/cbig_scrna.slurm
PIPELINE_PWD=/path/to/working/directory/Dataset1

cd $PIPELINE_PWD

bash $PIPELINE_HOME/launch_cbig_scrna.sh \
-d ${PIPELINE_PWD} \
--steps 3

```

The resulting output files are deposited into `/path/to/working/directory/Dataset1/step1`. The filtered expression matrix, features, and barcode files outputed by SpaceRanger are located in `/path/to/working/directory/Dataset1/step3/sample/ouput_folder/outs/filtered_feature_bc_matrix`.

**Note:** If SpaceRanger was successfull, it will display _Pipestance completed successfully!_. If this message is not displayed, you should check the error logs in `/path/to/working/directory/Dataset1/step3/sample/ouput_folder.log` and re-run Step 1.

- - - -

## Step 2: Data Preprocessing and Standardization

In Step 2, the spatial feature-barcode expression matrices are log-transformed and scaled using a standardized approach. Individual cells are filtered to retain only those barcodes specified in the metadata file provided by the original investigators. Metadata column names are standardized and appended to the count matrices. The processed data are saved as both Seurat and Scanpy objects.


The following parameters are adjustable for Step 2 of the spatial RNAseq Visium HD track (`/path/to/working/directory/Dataset1/job_info/parameters/parameters.txt`):

|Column name|Default|Description|
|:--|:--|:--|
|**par_meta_data_spatial**|NULL|Path to investigator-provided metadata .csv file.|
|**par_barcodes_spatial**|NULL|Column name in the investigator-provided metadata file describing the cell barcodes. |
|**par_sample_id_spatial**|NULL|Column name in the investigator-provided metadata file describing the sample IDs.|
|**par_subject_id_spatial**|NULL|Column name in the investigator-provided metadata file describing the Subject IDs.|
|**par_tissue_spatial**|NULL|Column name in the investigator-provided metadata file describing the tissue type. |
|**par_tissue_region_spatial**|NULL|Column name in the investigator-provided metadata file describing the tissue region.|
|**par_cell_type_main_spatial**|NULL|Column name in the investigator-provided metadata file describing the main cell type.|
|**par_disease_status_main_spatial**|NULL|Column name in the investigator-provided metadata file describing the disease status of the subject from which the cells were obtained.|


To execute Step 2, run the following code:

```
PIPELINE_HOME=/path/to/pipeline/CBIG_scRNAseq/cbig_scrna.slurm
PIPELINE_PWD=/path/to/working/directory/Dataset1

cd $PIPELINE_PWD

bash $PIPELINE_HOME/launch_cbig_scrna.sh \
-d ${PIPELINE_PWD} \
--steps 4
```

After running the Step 2, the following files will be available in `/path/to/working/directory/Dataset1/step4`:

```
Dataset1
└── job_info
    ├── configs
    │   └── cbig_scrna_config.ini
    ├── logs
    │   └── step_3_spacerangerV1_Adult_Mouse_Brain_S5.log
    ├── parameters
    │   └── parameters.txt
    └── summary_report.txt
```

- - - -

## Data upload to C-BIG


