#!/usr/bin/env nextflow

// Note: Nextflow.enable.dsl = 2 enables the workflow structure and functions
nextflow.enable.dsl = 2


// Helper functions 
def checkFileExists(file_path) {
  f = file(file_path)
  if ( !f.isFile() || !f.exists() ) {
    exit 1, "File '$file_path' does not exist!"
  }
}

if (params.help){
    helpMessage()
    exit 0
} 

def helpMessage() {
    // Add to this help message with new command line parameters
    log.info"""
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
             """
                                                
}


//=============================================================================
// PROCESSES
//=============================================================================

process generate_bai {
    publishDir "$projectDir/inputs/"
    output:
    val 1
    script:
    """
        samtools index $projectDir/${params.bam_file}
    """

}

process collect_subset_for_trees {
    output:
    file 'canada_recent_subset.fasta' 
    file 'canada_oneper_subset.fasta'
    file 'global_oneper_subset.fasta'

    script:
    println "Generating subsampling of GISAID data for tree construction."
    """
    python $projectDir/bin/subsample.py $projectDir/${params.gisaid_metadata_file} $projectDir/${params.gisaid_msa_file}
    """
}

process align_canada_recent {
    input:
    file 'canada_recent_subset.fasta'
    output:
    file 'canada_recent_subset_aligned.fasta'
    script:
    println "Align the query sequence against the subsamples."
    """
     muscle -profile -in1 "canada_recent_subset.fasta" -in2 "$projectDir/${params.consensus_file}" -out "canada_recent_subset_aligned.fasta"
    """
}

process align_canada_oneper {
    input:
    file 'canada_oneper_subset.fasta'
    output:
    file 'canada_oneper_subset_aligned.fasta'
    script:
    println "Align the query sequence against the subsamples."
    """
     muscle -profile -in1 "canada_oneper_subset.fasta" -in2 "$projectDir/${params.consensus_file}" -out "canada_oneper_subset_aligned.fasta"
    """
}

process align_global_oneper {
    input:
    file 'global_oneper_subset.fasta'
    output:
    file 'global_oneper_subset_aligned.fasta'
    script:
    println "Align the query sequence against the subsamples."
    """
     muscle -profile -in1 "global_oneper_subset.fasta" -in2 "$projectDir/${params.consensus_file}" -out "global_oneper_subset_aligned.fasta"
    """
}

process build_global_oneper {
    publishDir "$projectDir/outputs/trees"
    input:
    file 'global_oneper_subset_aligned.fasta'
    output:
    file 'global.treefile'
    script:
    println "Build the global tree."
    """
    iqtree -s global_oneper_subset_aligned.fasta -m GTR -T AUTO -o Wuhan --prefix global_oneper
    """
}

process build_canada_oneper {
    publishDir "$projectDir/outputs/trees"
    input:
    file 'canada_oneper_subset_aligned.fasta'
    output:
    file 'canada_oneper.treefile'
    script:
    println "Build the canadian lineage tree."
    """
    iqtree -s canada_oneper_subset_aligned.fasta -m GTR -T AUTO -o Wuhan --prefix canada_oneper
    """
}

process build_canada_recent {
    publishDir "$projectDir/outputs/trees"
    input:
    file 'canada_recent_subset_aligned.fasta'
    output:
    file 'canada_recent.treefile'
    script:
    println "Build the canadian recent tree."
    """
    iqtree -s canada_recent_subset_aligned.fasta -m GTR -T AUTO -o 'Wuhan|402124' --prefix canada_recent
    """
}


process draw_canada_recent {
    publishDir "$projectDir/outputs/trees"
    input:
    file 'canada_recent.treefile'
    output:
    file 'canada_recent.png'
    script:
    println "Draw the canadian tree of recent sequences."
    """
    python $projectDir/bin/plot_trees.py canada_recent.treefile canada_recent.png
    """
}

process draw_canada_oneper {
    publishDir "$projectDir/outputs/trees"
    input:
    file 'canada_oneper.treefile'
    output:
    file 'canada_oneper.png'
    script:
    println "Draw the canadian tree of lineage sequences."
    """
    python $projectDir/bin/plot_trees.py canada_oneper.treefile canada_oneper.png
    """
}

process draw_global_oneper {
    publishDir "$projectDir/outputs/trees"
    input:
    file 'global_oneper.treefile'
    output:
    file 'global_oneper.png'
    script:
    println "Draw the canadian tree of lineage sequences."
    """
    python $projectDir/bin/plot_trees.py global_oneper.treefile global_oneper.png
    """
}



process plot_coverage { 
    publishDir "$projectDir/outputs/coverage"
    input:
    file "$projectDir/${params.bam_file}.bai"
    output:
    file 'coverage.png' 
 
    script:
    println "Generate a plot of the coverage."
    """
    python $projectDir/bin/plot_coverage.py $projectDir/${params.bam_file} 
    """
}

process plot_heatmaps{
    publishDir "$projectDir/outputs/heatmaps", mode: 'copy'
    input:
    file "$projectDir/${params.bam_file}.bai"
    output:
    file 'heatmap*.png' 
 
    script:
    log.info "Generate heatmap plots for each VOC using VCFParser"
    """
    vcfparser -i $projectDir/${params.ivar_file} -bam $projectDir/${params.bam_file} -voc all --subplots_mode oneplotperfile --annotate
    """
}

process build_report {
    input:
    file 'coverage.png'
    file 'canada_oneper.png'
    file 'canada_recent.png'
    output:
    file 'covid_variants_report.tex'
    script:
    println "Generating the latex version of the report."
    """
    python $projectDir/bin/report.py $projectDir ${params.consensus_file} ${params.bam_file}
    """
}

process convert_report_to_pdf {

    publishDir "$projectDir/outputs/report", mode: 'copy'
    input:
    file 'covid_variants_report.tex' 
    output:
    file 'covid_variants_report.pdf'   
    println "Generate the PDF file."
    """
    pdflatex covid_variants_report.tex
    """
}

//=============================================================================
// WORKFLOW
//=============================================================================
workflow {
    println "Checking for inputs:"
    checkFileExists(params.gisaid_metadata_file)
    checkFileExists(params.gisaid_msa_file)
    checkFileExists(params.consensus_file)
    checkFileExists(params.bam_file)
    checkFileExists(params.variants_file)

    // Check that that .bai file exists
    generate_bai()

    // Build the trees
    // Step 1: Generate an appropriate subset of GISAID sequences for each tree
    collect_subset_for_trees()

    // Step 2: Align the query sequence to each of the subsets    
    align_canada_recent(collect_subset_for_trees.out[0])
    align_canada_oneper(collect_subset_for_trees.out[1])
    // align_global_oneper(collect_subset_for_trees.out[2])
    
    // Step 3: Build the trees using iqtree
    // build_global_oneper(align_global_oneper.out)
    build_canada_oneper(align_canada_oneper.out)
    build_canada_recent(align_canada_recent.out)

    // Step 4: Build 
    // draw_global_oneper(build_global_oneper.out)
    draw_canada_oneper(build_canada_oneper.out)
    draw_canada_recent(build_canada_recent.out)

    // Plot the coverage
    plot_coverage(generate_bai.out)

    // Plot the heatmaps
    plot_heatmaps(generate_bai.out)

    // Build the report
    build_report(plot_coverage.out,draw_canada_oneper.out,draw_canada_recent.out)

    // Convert the report to PDF
    convert_report_to_pdf(build_report.out)
}



//=============================================================================
// INTROSPECTION
//=============================================================================
workflow.onComplete {
    log.info """\n
    Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Results Dir  : ${file(params.outdir)}
    Work Dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """.stripIndent()
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
