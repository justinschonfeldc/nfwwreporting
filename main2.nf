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
    publishDir "$projectDir/inputs/", mode: 'copy'
    input:
        tuple val(sampleID), path(bamPath), path(consensusPath), path(ivarPath)
    output:
        tuple val(sampleID), path(bamPath), path("${bamPath}.bai"), path(ivarPath)
    script:
    """
        echo $sampleID
        echo $bamPath
        echo $consensusPath
        echo $ivarPath
        samtools index ${bamPath}
    """
}

process generate_consensus_fasta {
    input:
        val(sampleList)
    output:
        path(con_aligned_fasta)
        path(sample_names)
    script:
    log.info "Generate Consensus fasta"
    log.info "${sampleList}"
    // transformedList = sampleList.collate(4)

    fileList = ""

    for (entry in sampleList) {
        log.info "${entry[2]}"
        fileList = fileList + " ${entry[2]}"
    }
    log.info "Concatenate consensus fasta files into one file"

    con_fasta = "combined_consensus.fasta"
    con_aligned_fasta = "con_aligned.fasta"
    sample_names = "sample_names.txt"

    """
        cat ${fileList} > ${con_fasta}
        awk '/>/ {print}' ${con_fasta} | sed 's/^.//' > ${sample_names}
        muscle -in ${con_fasta} -out ${con_aligned_fasta} -maxiters 2
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
    path(con_aligned_fasta)
    file 'canada_recent_subset.fasta'
    output:
    file 'canada_recent_subset_aligned.fasta'
    script:
    println "Align the query sequence against the subsamples."
    """
     muscle -profile -in1 "canada_recent_subset.fasta" -in2 ${con_aligned_fasta} -out "canada_recent_subset_aligned.fasta"
    """
}

process align_canada_oneper {
    input:
    path(con_aligned_fasta)
    file 'canada_oneper_subset.fasta'
    output:
    file 'canada_oneper_subset_aligned.fasta'
    script:
    println "Align the query sequence against the subsamples."
    """
     muscle -profile -in1 "canada_oneper_subset.fasta" -in2 "${con_aligned_fasta}" -out "canada_oneper_subset_aligned.fasta"
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
    path(sampleFile)
    output:
    file 'canada_recent.png'
    script:
    println "Draw the canadian tree of recent sequences."
    """
    echo "SAMPLEFILE"
    cat ${sampleFile} > $projectDir/outputs/sample_names.txt
    python $projectDir/bin/plot_trees.py canada_recent.treefile canada_recent.png --samples ${sampleFile}
    """
}

process draw_canada_oneper {
    publishDir "$projectDir/outputs/trees"
    input:
    file 'canada_oneper.treefile'
    path(sampleFile)
    output:
    file 'canada_oneper.png'
    script:
    println "Draw the canadian tree of lineage sequences."
    """
    python $projectDir/bin/plot_trees.py canada_oneper.treefile canada_oneper.png --samples ${sampleFile}
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
    // file "$projectDir/${params.bam_file}.bai"
    // path(bamPath)
    tuple val(sampleID), path(bamPath), path(baiPath), path(ivarPath)
    output:
    file "${sampleID}_coverage.png"
 
    script:
    println "Generate a plot of the coverage."
    println bamPath
    """
    python $projectDir/bin/plot_coverage.py ${bamPath} ${sampleID}_coverage.png
    """
}


process plot_heatmaps {
    publishDir "$projectDir/outputs/heatmaps", mode: 'copy'
    input:
    val(sampleList)
    output:
    file "heatmap*.png" 
    val 1
 
    script:
    // log.info "HERE IT IS"
    // log.info "${sampleList}"
    transformedList = sampleList.collate(4)
    // log.info "${transformedList}"

    // Write batch file
    // Create a file handler
    File file = new File("${projectDir}/inputs/vcfparser_batch.tsv")
    // Initialize the file for writing by overwriting any existing file
    file.write("")
    for ( entry in transformedList) {
        log.info "Entry: ${entry}, ${entry[0]}"
        // Append each entry to the existing file
        file << "${entry[0]}\t${entry[3]}\t${entry[1]}\n"
    }
    log.info "Generate heatmap plots for each VOC using VCFParser"
    """
    vcfparser -f ${projectDir}/inputs/vcfparser_batch.tsv -voc all --subplots_mode oneplotperfile --annotate
    """
}

process build_report {
    input:
    val(coverage)
    val(heatmaps)
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
    checkFileExists(params.batch_file)

    inFile = file(params.batch_file)
    samples = []
    allLines = inFile.readLines()
    for ( line : allLines ) {
        if ( line.split("\t")[0] == "SampleID" ) {
            continue
        }
        if ( line.split("\t").size() != 4 ) {
            println "Incorrect number of parameters.  Skipping line: $line"
            continue
        }
        samples.add(line.split("\t"))
    } 

    // println samples

    // inputChannel = Channel.from(samples).map( it -> tuple(it[0], it[1], it[2], it[3]))

    // generate_bai(inputChannel)

    // baiChannel = Channel.from(generate_bai.out).view()

    // plot_coverage(baiChannel)

    generate_bai(Channel.from(samples).map(it -> tuple(it[0], it[1], it[2], it[3]))) 

    // generate_bai.out.view()

    plot_coverage(generate_bai.out)

    plot_heatmaps(generate_bai.out.collect().view())
    // Channel.from(samples).map(it -> tuple(it[0], it[1], it[2], it[3])) | generate_bai | plot_coverage

    Channel.from(samples).collect() | generate_consensus_fasta

    collect_subset_for_trees()

    align_canada_recent(generate_consensus_fasta.out[0],collect_subset_for_trees.out[0])
    align_canada_oneper(generate_consensus_fasta.out[0],collect_subset_for_trees.out[1])
    
    // Step 3: Build the trees using iqtree
    build_canada_oneper(align_canada_oneper.out)
    build_canada_recent(align_canada_recent.out)

    // Step 4: Build 
    draw_canada_oneper(build_canada_oneper.out,generate_consensus_fasta.out[1])
    draw_canada_recent(build_canada_recent.out,generate_consensus_fasta.out[1])

    // Build the report
    build_report(plot_coverage.out.collect().toList(),plot_heatmaps.out[1],draw_canada_oneper.out,draw_canada_recent.out)

    // Convert the report to PDF
    convert_report_to_pdf(build_report.out)

    // Channel.from(generate_bai.out).view()

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