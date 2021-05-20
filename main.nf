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
    //publishDir "$projectDir/inputs/", mode: 'copy'
    input:
        tuple val(sampleID), path(bamPath), path(consensusPath), path(ivarPath)
    output:
        tuple val(sampleID), path(bamPath), path("${bamPath}.bai"), path(ivarPath)
    script:
    log.info "Generate BAI files for ${sampleID}"
    """
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
    
    
    // transformedList = sampleList.collate(4)

    fileList = ""

    for (entry in sampleList.collate(4)) {
        log.info "Concat consensus fasta file: ${entry[2]}"
        fileList = fileList + " ${entry[2]}"
    }
    

    con_fasta = "combined_consensus.fasta"
    con_aligned_fasta = "con_aligned.fasta"
    sample_names = "sample_names.txt"

    log.info "Generate a single consensus multi-fasta file ${con_aligned_fasta}"
    //log.info "${sampleList}"

    """
        cat ${fileList} > ${con_fasta}
        awk '/>/ {print}' ${con_fasta} | sed 's/^.//' > ${sample_names}
        muscle -in ${con_fasta} -out ${con_aligned_fasta} -maxiters 2
    """
}

process assign_nextclade_2_consensus_seqs {
    publishDir "$projectDir/outputs/nextclade", mode: 'copy'
    
    input:
    path( concat_consensus_multi_fasta )
    output:
    path('nextclade-output-file.tsv')
    
    script:
    """
    nextclade -i $concat_consensus_multi_fasta --output-tsv  nextclade-output-file.tsv
    """
}

process prepate_input_4_nextclade {
    publishDir "$projectDir/temp/", mode: 'copy'
    
    input:
    val(sampleList)
    output:
    path 'concat_consensus_multi.fasta'
    
    script:
    fileFastaPathList = ""
    concat_consensus_fasta="nextclade-output-file.tsv"

    for (entry in sampleList.collate(4)) {
        fileFastaPathList = fileFastaPathList + " ${entry[2]}"
    }

    script:
    """
    fastaFileNames=()
    for sample_path in ${fileFastaPathList};do 
     filename=\$(basename \$sample_path); 
     cp \$sample_path \$filename;
     sed -i "1 s/^.\\+/>\$filename/g" \$filename
     fastaFileNames[\${#fastaFileNames[@]}]=\$filename
    done
    
    cat \${fastaFileNames[@]} > concat_consensus_multi.fasta
    """
}

process collect_subset_for_trees {
    output:
    file 'canada_recent_subset.fasta' 
    file 'canada_oneper_subset.fasta'
    file 'global_oneper_subset.fasta'

    script:
    log.info "Generating subsampling of GISAID data for tree construction."
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
    publishDir "$projectDir/outputs/trees", mode: 'copy'
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
    publishDir "$projectDir/outputs/trees", mode: 'copy'
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
    publishDir "$projectDir/outputs/trees", mode: 'copy'
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
    publishDir "$projectDir/outputs/coverage", mode: 'copy'
    input:
    tuple val(sampleID), path(bamPath), path(baiPath), path(ivarPath)
    
    output:
    file "*_coverage.png"
    
 
    script:
    log.info "Generate a plot of the coverage for ${sampleID}"
    
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
    sampleIDs=[]
    for ( entry in transformedList) {
        //log.info "Entry: ${entry}, ${entry[0]}"
        // Append each entry to the existing file
        file << "${entry[0]}\t${entry[3]}\t${entry[1]}\n"
        sampleIDs.add("${entry[0]}")
    }
    log.info "Generate heatmap plots for each VOC using VCFParser for samples: ${sampleIDs}"
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
    log.info "Generating the latex version of the report. Heatmap files: ${heatmaps}"
    """
    python $projectDir/bin/report.py $projectDir ${params.batch_file}
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

    samples_batch_ch = Channel.from(samples).map(it -> tuple(it[0], it[1], it[2], it[3]))
    generate_bai(samples_batch_ch)

    //Step 0: Assign NextClade to consensus COVID-19 sequences
    assign_nextclade_2_consensus_seqs(prepate_input_4_nextclade(samples_batch_ch.collect()))

    //Step 1: Generate coverage plots and multi-samples heatmaps for each VOC
    plot_coverage(generate_bai.out)
    plot_heatmaps(generate_bai.out.collect())

    //// Channel.from(samples).map(it -> tuple(it[0], it[1], it[2], it[3])) | generate_bai | plot_coverage

    ///Channel.from(samples).collect() | 
    //Step 2: Generate alignment of input samples to GISAID Canada sequences
    generate_consensus_fasta(samples_batch_ch.collect())
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
    build_report(plot_coverage.out.collect().toList(),
        plot_heatmaps.out[0].collect(),
        draw_canada_oneper.out,
        draw_canada_recent.out)

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
