# Wastwater VOC pipeline with automatic report generation


## Installation

```

```

## Running

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
Sample BAM File

Sample Consensus FIle

IVAR Output File

### Reference Inputs:
GISAID Metadata

GISAID Multiple Sequence Alignment

Variant Specification
