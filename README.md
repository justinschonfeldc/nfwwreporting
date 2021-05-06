# Wastwater VOC pipeline with automatic report generation

Convenient workflow for wastewater COVID-19 report generation on Variants of Concern (VOCs)
powered by NextFlow, Conda and Singularity

## Requirements
* Nextflow >=20.07.1 (for DSL 2 support)
* [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps)
* Conda
* Git LFS (optional)

## Installation
1. Create nextflow conda environment with latest NextFlow installed
    ```bash
      conda create -n nextflow@latest nextflow && conda activate nextflow@latest
    ```

1. Install `git-lfs` to download Singularity image to `data` folder
    ```
        # if running inside conda environment
        conda install -c conda-forge git-lfs
        # if running locally on Ubuntu
        apt-get install git-lfs
    ```
1. Clone repository with Git and Git LFS
    ```
        git clone https://github.com/justinschonfeldc/nfwwreporting.git
        git lfs pull
    ```
1. Run NextFlow on inputs (see `Running` and `Usage` sections)   

**Note:** Singularity is already installed on the `waffles` server, but NextFlow conda environment is outdated. 
We recommend installating custom NextFlow conda environment 

## Running
1. Activate NextFlow Conda environment if NextFlow run from Conda 
    ```
    conda activate -n nextflow@latest
    
    ```
1. Setup inputs and run NextFlow workflow. Additional run parameters are in the `Usage` section
    ```bash
    nextflow run main.nf -profile conda --resume
    ```


## Usage
```
Pipeline that automates COVID-19 wastewater report generation
    Usage:
    nextflow run --consensus_file consensus.fasta --bam_file consensus.bam  -profile <profile>
    Mandatory arguments:
      --consensus_file [file]                   Path to the consensus COVID-19 genome wastewater assembly
      --bam_file [file]                         Path to the aligned wastewater reads to the MN908947.3 reference
      --vcfparser_batch_file                    Path to VCFParser input batch file with sample name, abs path to TSV/VCF file and BAM file
      -profile                                  Available: conda, singularity, standard
      --resume                                  Resume from the last failed step
      --help                                    Display this help message
```

## Inputs

### Sample Inputs:
Batch File - batch_file - Assign a relative path in nextflow.config.  Default: "${params.indir}/batch.tsv"

### Reference Inputs:
GISAID Metadata - gisaid_msa_file - Assign a relative path in nextflow.config.  Default: gisaid_msa_file = "${params.indir}/gisaid_msa.fasta 

GISAID Multiple Sequence Alignment - gisaid_metadata_file - Assign a relative path in nextflow.config.  Default: gisaid_metadata_file = "${params.indir}/gisaid_metadata.fasta

Variant Specification - variants_file - Assign a relative path in nextflow.config. Default: variants_file = "${params.indir}/cov_lineage_variants.tsv"


### Batch Example:
(tab delimeted)

SampleID    BAMFile ConsensusFile   IvarTSVFile
id1 /path/inputs/sample.bam /path/inputs/consensus.fasta    /path/inputs/ivar.tsv
id2 /path/inputs/sample2.bam /path/inputs/consensus2.fasta    /path/inputs/ivar2.tsv
