# Available data
- [FASTQ files](#fastq-files)
- [Feature-barcode expression matrix](#feature-barcode-expression-matrix)
- [Seurat and Scanpy objects](#seurat-and-scanpy-objects)
- [Metadata](#metadata)

- - - -

For each scRNAseq study hosted on C-BIG, we provide both raw and processed data at various stages of the processing pipeline. Specifically, each study includes the following:

1.	Raw FASTQ files for each sample (.fastq);
2.	A unified feature-barcode expression matrix (barcodes.tsv, features.tsv, matrix.mtx);
3.	A unified Seurat object (.RDS);
4.	A standardized metadata file provided by the authors (.csv)


Ensuring the availability of multiple data types provides C-BIG users with the flexibility to choose the most appropriate data for their specific analytical needs. We further describe each data type in the subsequent sections. 

## FASTQ files

Independent FASTQ files for each sample obtained directly from the original investigators are uploaded to C-BIG. The FASTQ files contain raw sequence reads, which, in the context of scRNAseq, represent transcriptomic data from individual cells. The FASTQ files are published without modification; however, their naming is adjusted according to standard 10xGenomics nomenclature: 

```
[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
```

Where `Read Type` is one of:

- R1: Read 1
- R2: Read 2

For example, in a scRNAseq study involving two distinct samples sequenced using paired-end technology, where each sample represents a different individual—Parkinson1 and Control1—the following FASTQ files will be available:

```
Parkinson1
├── Parkinson1_S1_L001_R1_001.fastq.gz
└── Parkinson1_S1_L001_R2_001.fastq.gz

Control1
├── Control1_S1_L001_R1_001.fastq.gz
└── Control1_S1_L001_R2_001.fastq.gz 
```

Raw FASTQ files are best suited for investigators aiming to replicate the original study’s findings and apply the exact processing methods described in the manuscript.

## Feature-barcode expression matrix
A single feature-barcode expression matrix for each study is uploaded to C-BIG. This matrix consists of three files: barcodes.tsv, features.tsv, and matrix.mtx, which collectively map the read counts of specific genes (features) to individual cells (barcodes). The feature-barcode expression matrix is derived from the unified Seurat object generated using the standardized pipeline described above. Regardless of the number of samples or subjects included in the study, the following files will be available:

- barcodes.tsv.gz
- features.tsv.gz
- matrix.mtx.gz

The feature-barcode expression matrix is best suited for investigators who wish to perform a comprehensive analysis of the data without the need to run the time-consuming and memory-intensive CellRanger counts pipeline. Additionally, the feature-barcode expression matrix can be imported into R for analysis using various single-cell frameworks, such as SingleCellExperiment, or into Scanpy in Python, depending on the user’s preference.

## Seurat and Scanpy objects
A unified Seurat object in R Data Serialization format (RDS) obtained using the pipeline described above is uploaded to C-BIG. Seurat objects store all the data and metadata associated with scRNAseq analysis, including gene expression data and cell-level metadata. To maintain clarity and avoid ambiguity, we intentionally keep the uploaded Seurat object simple, including only the raw, normalized, and scaled RNA assays, as well as the author-provided metadata. This approach ensures that users have the flexibility to perform their own analyses—such as dimensionality reduction, clustering, and differential gene expression—without being constrained by predefined steps. Regardless of the number of samples or subjects included in the study, the following file will be available:

- seurat.RDS
- scanpy.h5ad

The unified Seurat object is best suited for investigators familiar with R and the Seurat framework who wish to perform a comprehensive analysis of the data without the need to run the time-consuming and memory-intensive CellRanger counts pipeline.

## Metadata
An author-provided metadata file for each scRNAseq study is uploaded to C-BIG. The metadata file is standardized to reduce confusion for external investigators, ensuring that essential cell-level information is provided for effective use of the data. Specifically, each metadata file will contain, at the very minimum, the following columns:

|Column name|Description|
|:--|:--|
|**barcodes_CBIG**|Cell barcodes|
|**orig.ident_CBIG**|Sample identifiers matching those in the FASTQ nomenclature|
|**nCount_RNA_CBIG**|Number of RNA counts per cell when processed by the CBIG pipeline|
|**nFeature_RNA_CBIG**|Number of features expressed per cell when processed by the CBIG pipeline|
|**sample_id_CBIG**|Sample identifier|
|**subect_id_CBIG**|Subject identifier|
|**cell_type_main_CBIG**|Primary cell type annotation|
|**disease_status_CBIG**|Subject disease status |
|**Tissue_CBIG**|Tissue from which individuals cells derived|

