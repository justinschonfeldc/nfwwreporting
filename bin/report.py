import click
import pandas as pd
import pysam
from Bio import SeqIO
import os

def generate_table(df,report_fh):
    """Generates a latex table from a pandas data frame."""
    num_columns = len(df.columns)
    report_fh.write("\\begin{center}\n")
    report_fh.write("\\begin{tabular}{")
    for i in range(num_columns):
        report_fh.write("|l")
    report_fh.write("|}\n")
    report_fh.write("\\hline\n")
    for i in range(num_columns):
        report_fh.write(f"{df.columns[i]}")
        if i == num_columns - 1:
            report_fh.write(" \\\\\n")
        else:
            report_fh.write(" & ")
    report_fh.write("\\hline\n")
    for j in range(len(df)):
        for i in range(num_columns):
            report_fh.write(f"{df.iloc[j,i]}")
            if i == num_columns - 1:
                report_fh.write(" \\\\\n")
            else:
                report_fh.write(" & ")
    report_fh.write("\\hline\n")
    report_fh.write("\\end{tabular}\n")
    report_fh.write("\\end{center}\n")

def check_bam_indexing(batch_path):
    sdata = {}
    inf = open(batch_path,"r")
    # Skip the lines
    inf.readline()
    for line in inf:
        sline = line.strip().split("\t")
        if len(sline) != 4:
            continue
        sdata[sline[0]] = dict()
        sdata[sline[0]]['bam_path'] = sline[1]
        sdata[sline[0]]['con_path'] = sline[2]
        sdata[sline[0]]['ivar_path'] = sline[3] 
    inf.close()

    for sample in sdata:
        bai_path = os.path.abspath(sdata[sample]['bam_path']) + '.bai'
        if not os.path.exists(bai_path):
            print("Creating bam index for {}".format(sdata[sample]['bam_path']))
            pysam.index(sdata[sample]['bam_path'], bai_path)

def generate_document_header(report_fh):
    report_fh.write("\\documentclass{article}\n")
    report_fh.write("\\usepackage[final]{pdfpages}\n")
    report_fh.write("\\usepackage{graphicx}\n")
    report_fh.write("\\usepackage{grffile}\n")
    report_fh.write("\\usepackage{float}\n")

    report_fh.write("\\title{SARS-CoV-2 Wastewater Variant Report}\n")
    report_fh.write("\\date{\\today}\n")
    report_fh.write("\\author{PHAC}\n")

    report_fh.write("\\begin{document}\n")
    report_fh.write("\\begin{titlepage}\n")
    report_fh.write("\\maketitle\n")
    report_fh.write("\\end{titlepage}\n")


def generate_document_footer(report_fh):
    report_fh.write("\\end{document}\n")

def generate_qa_section(report_fh, base_path, batch_path):
    """Generates a quality assurance section."""
    
    sdata = {}
    inf = open(batch_path,"r")
    # Skip the lines
    inf.readline()
    for line in inf:
        sline = line.strip().split("\t")
        if len(sline) != 4:
            continue
        sdata[sline[0]] = dict()
        sdata[sline[0]]['bam_path'] = sline[1]
        sdata[sline[0]]['con_path'] = sline[2]
        sdata[sline[0]]['ivar_path'] = sline[3] 
    inf.close()

    for sample in sdata:
        for seq_record in SeqIO.parse(sdata[sample]['con_path'], "fasta"):
            sdata[sample]['numn'] = seq_record.seq.upper().count("N")
            sdata[sample]['numg'] = seq_record.seq.strip(" -").count("-")
            break

        samfile = pysam.AlignmentFile(sdata[sample]['bam_path'],"rb")
        sdata[sample]['unmapped'] = samfile.unmapped
        sdata[sample]['mapped'] = samfile.mapped 

    results = {}
    results['SampleID'] = list()
    results['No. N'] = list()
    results['No. Gap'] = list()
    results['Unmapped'] = list()
    results['Mapped'] = list()

    for sample in sdata:
        results['SampleID'].append(sample)
        results['No. N'].append(sdata[sample]['numn'])
        results['No. Gap'].append(sdata[sample]['numg'])
        results['Unmapped'].append(sdata[sample]['unmapped'])
        results['Mapped'].append(sdata[sample]['mapped'])

    report_fh.write("\\section{QA Metrics and Visualizations}\n")
    report_fh.write("Summary Statistics:\n")
    df = pd.DataFrame(results)
    generate_table(df,report_fh)
    # Show the coverage plots
    for sample in sdata:
        report_fh.write("\\begin{figure}[H]\n")
        report_fh.write(f"\\includegraphics[width=7in]{{{base_path+'/outputs/coverage/'+sample+'_coverage.png'}}}")
        report_fh.write(f"\\caption{{{sample}}}")
        report_fh.write("\\end{figure}\n")
 
def generate_qa_consensus_tree(report_fh,path):
    """Generates a section of the report displaying the consensus placed in a tree."""
    report_fh.write("\\section{Consensus Sequence Tree}\n")
    report_fh.write("\\begin{figure}[H]\n")
    report_fh.write(f"\\includegraphics[width=6in,height=6in]{{{path+'canada_oneper.png'}}}\n")
    report_fh.write("\\caption{ML-Tree (iqtree 2.0.3) composed of the sample, the Wuhan reference, and 1 randomly sampled sequence per Canadian Pangolin lineage.}")
    report_fh.write("\\end{figure}\n")
    # report_fh.write("\\begin{figure}\n")
    # report_fh.write(f"\\includegraphics[width=6in,height=6in]{{{path+'global_sample_tree.png'}}}\n")
    # report_fh.write("\\caption{ML-Tree (iqtree 2.0.3) composed of the sample, the Wuhan reference, and 1 randomly sampled sequence per global Pangolin lineage.}")
    # report_fh.write("\\end{figure}\n")
    report_fh.write("\\begin{figure}[H]\n")
    report_fh.write(f"\\includegraphics[width=6in,height=6in]{{{path+'canada_recent.png'}}}\n")
    report_fh.write("\\caption{ML-Tree (iqtree 2.0.3) composed of the sample, the Wuhan reference, and the 100 most recent Canadian sequences present in the GISAID MSA.}")
    report_fh.write("\\end{figure}\n")

def generate_variant_summary_section(report_fh, path, variant_results, batch_path):
    """Generates a variant summary section for the report."""
    report_fh.write("\\section{Variant Presence/Absence Summary}\n")

    sdata = {}
    inf = open(batch_path,"r")
    # Skip the lines
    inf.readline()
    for line in inf:
        sline = line.strip().split("\t")
        if len(sline) != 4:
            continue
        sdata[sline[0]] = dict()
        sdata[sline[0]]['bam_path'] = sline[1]
        sdata[sline[0]]['con_path'] = sline[2]
        sdata[sline[0]]['ivar_path'] = sline[3] 
    inf.close()

    for sample in sdata:
        summary = {}
        summary['Sample'] = list()
        summary['Variant'] = list()
        summary['Key Mut Present'] = list()
        summary['Key Mut Total'] = list()
        summary['Sig Mut Present'] = list()
        summary['Sig Mut Total'] = list()
        for vnt in variant_results[sample]:
            summary['Sample'].append(sample)
            summary['Variant'].append(vnt)
            summary['Key Mut Present'].append(variant_results[sample][vnt]['key_cnt'])
            summary['Key Mut Total'].append(variant_results[sample][vnt]['key_tot'])
            summary['Sig Mut Present'].append(variant_results[sample][vnt]['sig_cnt'])
            summary['Sig Mut Total'].append(variant_results[sample][vnt]['sig_tot'])
        df = pd.DataFrame(summary)
        generate_table(df,report_fh)

def generate_variant_specific_sections(report_fh, base_path, batch_path, variant_results):
    """Generate a table for each variant specific section."""
    report_fh.write("\\section{Variant Specific Support}\n")

    sdata = {}
    inf = open(batch_path,"r")
    # Skip the lines
    inf.readline()
    for line in inf:
        sline = line.strip().split("\t")
        if len(sline) != 4:
            continue
        sdata[sline[0]] = dict()
        sdata[sline[0]]['bam_path'] = sline[1]
        sdata[sline[0]]['con_path'] = sline[2]
        sdata[sline[0]]['ivar_path'] = sline[3] 
    inf.close()

    for sample in sdata: 
        variant_results[sample] = dict()
        report_fh.write(f"{sample}\n")
        df = pd.read_csv(base_path+"/inputs/cov_lineage_variants.tsv",sep="\t")
        df['RS']=0
        df['TR']=0
        df['Freq']="0.0"
        df['Present']="No"

        samfile = pysam.AlignmentFile(sdata[sample]['bam_path'],"rb")
        rsl = df.columns.get_loc("RS")
        trl = df.columns.get_loc("TR")
        freql = df.columns.get_loc("Freq")
        presentl = df.columns.get_loc("Present")

        for i in range(len(df)):
            crs = 0
            ctr = 0
            start = df['Position'][i]
            stop = start+df['Length'][i]
            itr = samfile.fetch("MN908947.3",start,stop)
            for read in itr:
                ap = read.get_aligned_pairs()
                vir = list()
                # print(f"start: {start} stop: {stop}\n{read}\n{read.get_aligned_pairs()}")
                for pair in ap:
                    if pair[1] == None:
                        continue
                    if pair[1] >= start and pair[1] < stop:
                        vir.append(pair[0])
                # print(vir)
                if len(vir)<df['Length'][i]:
                    continue
                subseq = ""
                for el in vir:
                    if el == None:
                        subseq+='-'
                    else:
                        subseq+=read.seq[el-1]
                # print(f"subseq: {subseq} alt: {df['Alt'][i]}"
                if df['Type'][i] == "Sub" and subseq == df['Alt'][i]:
                    crs += 1
                elif df['Type'][i] == "Del" and subseq == ("-"*df['Length'][i]):
                    crs += 1
                ctr += 1

            df.iat[i,rsl] = crs
            df.iat[i,trl] = ctr
            if ctr > 0:
                freq_val = crs/ctr
                df.iat[i,freql] = f"{freq_val:.3f}"
            else:
                df.iat[i,freql] = 0
            if crs>20:
                df.iat[i,presentl] = "Yes"

        for i in range(len(df)):
            variant = df['VOC'][i]
            if variant in variant_results[sample]:
                if df['Key'][i] == True:
                    variant_results[sample][variant]['key_tot'] += 1
                    if df['Present'][i] == "Yes":
                        variant_results[sample][variant]['key_cnt'] += 1
                if df['SignatureSNV'][i] == True:
                    variant_results[sample][variant]['sig_tot'] += 1
                    if df['Present'][i] == "Yes":
                        variant_results[sample][variant]['sig_cnt'] += 1        
            else:
                variant_results[sample][variant] = dict()
                variant_results[sample][variant]['key_cnt']=0
                variant_results[sample][variant]['sig_cnt']=0
                variant_results[sample][variant]['key_tot']=0
                variant_results[sample][variant]['sig_tot']=0
                if df['Key'][i] == True:
                    variant_results[sample][variant]['key_tot'] += 1
                    if df['Present'][i] == "Yes":
                        variant_results[sample][variant]['key_cnt'] += 1
                if df['SignatureSNV'][i] == True:
                    variant_results[sample][variant]['sig_tot'] += 1
                    if df['Present'][i] == "Yes":
                        variant_results[sample][variant]['sig_cnt'] += 1
        
        df = df.drop("Key",1)
        df = df.drop("Type",1)
        df = df.drop("Length",1)
        df = df.drop("Ref",1)
        df = df.drop("Alt",1)
        df = df.drop("PangoLineage",1)
        df = df.drop("NextStrainClade",1)
        df = df.drop("SignatureSNV",1)

        vocl = df['VOC'].value_counts().index
        for el in vocl:
            dfs = df[df['VOC']==el]
            generate_table(dfs,report_fh)



    file_list = os.listdir(base_path+"/outputs/heatmaps/")
    png_file_list = [x for x in file_list if x.endswith('.png')]
    
    for file_name in png_file_list:        
        report_fh.write("\\begin{figure}[H]\n")
        report_fh.write(f"\\includegraphics{{{base_path+'/outputs/heatmaps/'+file_name}}}\n")
        report_fh.write("\\end{figure}\n")

def generate_unassociated_variant_section(report_fh, path, bam_file_path):
    """Generate a list of important mutations which are unassociated with a ny particluar variant"""
    report_fh.write("\\section{Significant Mutations Unassociated with Specific Variants}\n")
    df = pd.read_csv(path+"cov_lineage_variants.tsv",sep="\t")
    df = df[df['Key']==True]
    df.reset_index(drop=True,inplace=True)  

    df['RS']=0
    df['TR']=0
    df['Freq']="0.0"
    df['Present']="No"

    samfile = pysam.AlignmentFile(bam_file_path,"rb")
    rsl = df.columns.get_loc("RS")
    trl = df.columns.get_loc("TR")
    freql = df.columns.get_loc("Freq")
    presentl = df.columns.get_loc("Present")

    for i in range(len(df)):
        # if df["Type"][i] != "Sub":
        #     continue
        crs = 0
        ctr = 0
        start = df['Position'][i]
        stop = start+df['Length'][i]
        itr = samfile.fetch("MN908947.3",start,stop)
        for read in itr:
            ap = read.get_aligned_pairs()
            vir = list()
            # print(f"start: {start} stop: {stop}\n{read}\n{read.get_aligned_pairs()}")
            for pair in ap:
                if pair[1] == None:
                    continue
                if pair[1] >= start and pair[1] < stop:
                    vir.append(pair[0])
            # print(vir)
            if len(vir)<df['Length'][i]:
                continue
            subseq = ""
            for el in vir:
                if el == None:
                    subseq+='-'
                else:
                    subseq+=read.seq[el-1]
            if df['Type'][i] == "Sub" and subseq == df['Alt'][i]:
                crs += 1
            elif df['Type'][i] == "Del" and subseq == ("-"*df['Length'][i]):
                crs += 1
            ctr += 1

        df.iat[i,rsl] = crs
        df.iat[i,trl] = ctr
        if ctr > 0:
            freq_val = crs/ctr
            df.iat[i,freql] = f"{freq_val:.3f}"
        else:
            df.iat[i,freql] = 0
        if crs>20:
            df.iat[i,presentl] = "Yes"
    df = df.drop("Key",1)
    df = df.drop("Type",1)
    df = df.drop("Length",1)
    df = df.drop("Ref",1)
    df = df.drop("Alt",1)
    df = df.drop("PangoLineage",1)
    df = df.drop("NextStrainClade",1)
    df = df.drop("SignatureSNV",1)

    generate_table(df,report_fh)

@click.command()
@click.argument('base_path')
@click.argument('batch_path')
def report(base_path, batch_path):
    """Generate a wastewater variant report file."""
    # Open the report file
    report_fh = open("covid_variants_report.tex","w")


    check_bam_indexing(base_path+"/"+batch_path)

    generate_document_header(report_fh)
    
    generate_qa_section(report_fh, base_path, base_path+"/"+batch_path)

    generate_qa_consensus_tree(report_fh,base_path+"/outputs/trees/")


    variant_results = dict()
    generate_variant_specific_sections(report_fh, base_path, base_path+"/"+batch_path, variant_results)

    generate_variant_summary_section(report_fh, base_path+"/inputs/", variant_results, base_path+"/"+batch_path)

    # generate_unassociated_variant_section(report_fh,base_path+"/inputs/", os.path.join(base_path,bam_file))

    generate_document_footer(report_fh)

    # Close the report file handler  
    report_fh.close()




if __name__ == "__main__":
    report()