import click
import pandas as pd
import random
from os import path
from Bio import SeqIO

@click.command()
@click.argument('metadata_fn')
@click.argument('msa_fn')
@click.option('--nums',default=1,help='Number of samples to get from each lineage.')
@click.option('--rseed',default=42,help='Random number seed used in sampling.  Default seed is 42.')
@click.option('--nrecent',default=100,help='Number of recent Canadian sequences in MSA to use in tree.')
def subsample(metadata_fn, msa_fn, nums, rseed, nrecent):
    # Check that the metadata file exists.
    if not (path.exists(metadata_fn) and path.isfile(metadata_fn)) :
        print(f"Error! {metadata_fn} is an invalid path or not a file!")
        return
    
    # Check that the multiple sequence alignment file exists.
    if not (path.exists(msa_fn) and path.isfile(msa_fn)):
        print(f"Error! {msa_fn} is an invalid path or not a file!")
        return
    
    # Set the random number seed
    random.seed(rseed)

    # Read the metadata file into a pandas data frame
    df = pd.read_csv(metadata_fn,sep='\t')

    # Read through the msa and calcuate the number of N's per sequence
    acc_ndict = {}
    for seq_record in SeqIO.parse(msa_fn, "fasta"):
        lid = seq_record.id.split("|")
        if len(lid) < 2:
            continue
        gei = lid[1]
        acc_ndict[gei] = seq_record.seq.upper().count("N")

    # Create two subsets, one with one per lineage and one with one per province
    # Create the one per lineage

    gisaid_id_acc_name_list = ['gisaid_epi_isl', 'Accession ID']
    if df.columns[df.columns.isin(gisaid_id_acc_name_list)].shape[0] == 0:
        print("No valid accession key found "+",".join(gisaid_id_acc_name_list)+".")
        return 1
    else:
        acckey = df.columns[df.columns.isin(gisaid_id_acc_name_list)][0]

   # if 'gisaid_epi_isl' in df.columns:
   #     acckey = 'gisaid_epi_isl'
   # elif 'Accession ID' in df.columns:
   #     acckey = 'Accession ID'
   # else:
   #     print("No valid accession key found (gisai_epi_isl,Accession ID).")
   #     return




    # Subset the data to rows with an accession (acckey)
    df = df[df[acckey].isna()==False]

    gisaid_pango_lin_name_list = ['pangolin_lineage', 'Pango lineage']
    if df.columns[df.columns.isin(gisaid_pango_lin_name_list)].shape[0] == 0:
        print("Pangoline lineage column not found. Tried "+",".join(gisaid_pango_lin_name_list)+"")
        return 1
    else:
        gisaid_pango_lin_name = df.columns[df.columns.isin(gisaid_pango_lin_name_list)][0]


    map_to_pango = {}
    for i in range(len(df)):
        map_to_pango[df[acckey][i]]=df[gisaid_pango_lin_name][i]

    lineage_list = df[gisaid_pango_lin_name].value_counts().index
    lineage_dict = {}
    for lineage in lineage_list:
        dfs = df[df[gisaid_pango_lin_name]==lineage]
        lineage_dict[lineage]=list()
        for accession in dfs[acckey]:
            if accession not in acc_ndict:
                continue
            lineage_dict[lineage].append([accession,acc_ndict[accession]])
        # Sort the set based on the number of Ns
        # Sorted default is ascending so the sequence with the fewest Ns should be first
        lineage_dict[lineage] = sorted(lineage_dict[lineage],key=lambda x:x[1])

    lineage_subset = set()
    for lineage in lineage_dict.keys():
        if len(lineage_dict[lineage]) == 0:
            continue
        lineage_subset.add(lineage_dict[lineage][0][0])
    
    lineage_subset.add("EPI_ISL_402124")

    outf = open("global_oneper_subset.fasta","w")
    for seq_record in SeqIO.parse(msa_fn, "fasta"):
        lid = seq_record.id.split("|")
        if len(lid) < 2:
            continue
        gei = lid[1]
        if gei in lineage_subset:
            if gei == "EPI_ISL_402124":
                outf.write(f">Wuhan\n")
            else:
                outf.write(f">{map_to_pango[gei]}\n")
            outf.write(f"{seq_record.seq}\n")
    outf.close()

    # Create a Canadian subset
    # JSTD: Clean up country field
    # JSTD: Province list is fragile, clean up assumes 'Canada' is the smallest value
    if 'country' in df.columns:
        dfc = df[df['country']=='Canada']
        dfc = dfc[dfc['division'].isna() == False]
        #province_list = dfc['division'].value_counts().index[:-1]
    elif 'Location' in df.columns:
        dfc = df.loc[df['Location'].str.contains("Canada")]
    else:
        print("There is no \'country\' field. Check metadata")
        return  1


    lineage_list = dfc[gisaid_pango_lin_name].value_counts().index
    lineage_dict = {}
    for lineage in lineage_list:
        dfs = dfc[dfc[gisaid_pango_lin_name]==lineage]
        lineage_dict[lineage]=list()
        for accession in dfs[acckey]:
            if accession not in acc_ndict:
                continue
            lineage_dict[lineage].append([accession,acc_ndict[accession]])
        # Sort the set based on the number of Ns
        # Sorted default is ascending so the sequence with the fewest Ns should be first
        lineage_dict[lineage] = sorted(lineage_dict[lineage],key=lambda x:x[1])

    lineage_subset = set()
    for lineage in lineage_dict.keys():
        if len(lineage_dict[lineage]) == 0:
            continue
        lineage_subset.add(lineage_dict[lineage][0][0])
    
    #If the reference is not in the lineage_subset then add it
    lineage_subset.add("EPI_ISL_402124")

    outf = open("canada_oneper_subset.fasta","w")
    for seq_record in SeqIO.parse(msa_fn, "fasta"):
        lid = seq_record.id.split("|")
        if len(lid) < 2:
            continue
        gei = lid[1]
        if gei in lineage_subset:
            if gei == "EPI_ISL_402124":
                outf.write(f">Wuhan\n")
            else:
                outf.write(f">{map_to_pango[gei]}\n")
            outf.write(f"{seq_record.seq}\n")
    outf.close()

    ids_in_msa = set()
    for seq_record in SeqIO.parse(msa_fn, "fasta"):
        lid = seq_record.id.split("|")
        if len(lid) < 2:
            continue
        ids_in_msa.add(lid[1])

    gisaid_data_submitted_name_list = ['Submission date', 'date_submitted']
    if df.columns[df.columns.isin(gisaid_data_submitted_name_list)].shape[0] == 0:
        print("Date submitted column not found. Tried " + ",".join(gisaid_data_submitted_name_list) + "")
        return 1
    else:
        gisaid_data_submitted_name = df.columns[df.columns.isin(gisaid_data_submitted_name_list)][0]

    dfcs = dfc.sort_values(gisaid_data_submitted_name,ascending = False)
    dfcs.reset_index(drop=True,inplace=True)  
    sequence_list = []
    lineage_list = []
    date_list = []
    found_recent = 0
    for i in range(len(dfcs)):
        if dfcs[acckey][i] in ids_in_msa:
            sequence_list.append(dfcs[acckey][i])
            lineage_list.append(dfcs[gisaid_pango_lin_name][i])
            date_list.append(dfcs[gisaid_data_submitted_name][i])
            found_recent += 1
        if found_recent == nrecent:
            break

    if "EPI_ISL_402124" not in sequence_list:
        sequence_list.append("EPI_ISL_402124")

    outf = open("canada_recent_subset.fasta","w")
    for seq_record in SeqIO.parse(msa_fn, "fasta"):
        lid = seq_record.id.split("|")
        if len(lid) < 2:
            continue
        gei = lid[1]
        if gei in sequence_list:
            pos = sequence_list.index(gei)
            if gei == "EPI_ISL_402124":
                outf.write(f">Wuhan|402124\n")
            else:
                outf.write(f">{lineage_list[pos]}|{date_list[pos]}|{sequence_list[pos][8:]}\n")
            outf.write(f"{seq_record.seq}\n")
    outf.close()

if __name__ == "__main__":
    subsample()