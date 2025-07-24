# Standardized processing protocol

A general outline of the RNAseq processing pipelines for data uploaded to C-BIG is illustrated in Figure 1. All analytical steps are carried out using the scRNAbox pipeline (1) and require collaborative input from both the original investigators and the administrators (The NeuroBioinformatics Core). <br />

- - - -

![Image](https://github.com/user-attachments/assets/b53f60b6-f6c3-447d-81d8-f6abad83d767)

 **Figure 1. General workflow of standardized preprocessing framework.** 

- - - -

The pipeline begins with raw sequencing files in FASTQ format, obtained from commercial platforms such as 10xGenomics. Investigators looking to upload their data must send the sample-specific FASTQ files to the administrators. The sample-specific FASTQ files are processed to generate feature-barcode expression matrices with default parameters. Intronic reads are optionally retained depending on the author’s original protocol and the sequences are aligned to the GRCh38 reference genome. During this step, the sample-specific FASTQ files provided by the investigator are renamed to include descriptive sample identifiers, following the standard 10xGenomics nomenclature (for further details, see Available Data).

The feature-barcode expression matrices are used to generate sample-specific Seurat (v5.1.0) objects. Seurat is an R-based package designed for the analysis and visualization of RNAseq data (2). Initially, a separate Seurat object is created for each sample. These individual Seurat objects undergo basic preprocessing steps, including log normalization and scaling. To ensure that the quality control thresholds match that of the original study, we filter the Seurat object to only retain barcodes in the author-provided metadata file. Upon processing the sample-specific Seurat objects, we generate a unified Seurat object that includes all samples in the study by merging the RNA counts layer. Finally, we add the author-provided metadata to the unified Seurat object. This unified Seurat object is made available on C-BIG and used to generate a single feature-barcode expression matrix describing all samples in the study, as well as a Scanpy (Python) object.

The code is available on [GitHub](https://github.com/mfiorini9/CBIG_scRNAseq). 

- - - -
**References**

1. Thomas, Rhalena A., Michael R. Fiorini, Saeid Amiri, Edward A. Fon, and Sali MK Farhan. "ScRNAbox: empowering single-cell RNA sequencing on high performance computing systems." BMC bioinformatics 25, no. 1 (2024): 319.  <br />
2. Hao, Yuhan, Tim Stuart, Madeline H. Kowalski, Saket Choudhary, Paul Hoffman, Austin Hartman, Avi Srivastava et al. "Dictionary learning for integrative, multimodal and scalable single-cell analysis." Nature biotechnology 42, no. 2 (2024): 293-304. <br />
