import pandas as pd
import os
import ast
from Bio import SeqIO

def write_fasta(path_to_df_species, path_to_orthologs_groups, path_to_dir_GCF, path_to_fasta_dir, path_to_df_dir):
    df_species = pd.read_csv(path_to_df_species,sep='\t')
    df_genes_ortho = pd.read_csv(path_to_orthologs_groups, index_col=0)

    dictionaries_seq={}
    for gen_dir in os.listdir(path_to_dir_GCF):
        if gen_dir.startswith('GCF'):
            path_to_dir = str(path_to_dir_GCF+ gen_dir + '/')
            path_list=[]
            for file in os.listdir(path_to_dir):
                path_list.append(str(path_to_dir_GCF + gen_dir + '/' + file))
            fasta_genome = path_list[0]
            dictionaries_seq[gen_dir] = SeqIO.to_dict(SeqIO.parse(fasta_genome, "fasta"))

    df_genes_ortho['orthologs'] = df_genes_ortho['orthologs'].apply(ast.literal_eval)
    genome_dict = {}
    for index, row in df_genes_ortho.iterrows():
        ortho_fasta = []
        all_genes = [row['main']] + row['orthologs']
        fasta_path = path_to_fasta_dir + row['main'] + ".fasta"
        for ortho in all_genes:
            tax = int(ortho.split('_')[1])
            genome = df_species.loc[df_species['Organism Taxonomic ID']==tax,'Assembly Accession'].iloc[0]
            if genome not in genome_dict:
                path = path_to_df_dir + genome + "/df_exons"
                genome_dict[genome] = pd.read_csv(path, index_col=0)
            df = genome_dict[genome]
            gene_id = f"gene-{ortho.split('_')[0]}"
            gene_df = df[df['gene'] == gene_id]
            for gene_index, gene_row in gene_df.iterrows():
                exons_list = ast.literal_eval(gene_row['exons'])
                header = f">{gene_row['gene'][5:]}_{str(exons_list[0])[5:]}_{tax}"
                chromosome = gene_row['chromosome']
                borders = [int(gene_row['start']) - 1, int(gene_row['end'])]
                species_dictionary = dictionaries_seq[genome]
                record = species_dictionary[chromosome]
                if gene_row['strand'] == '-':
                    exon_seq = record.seq[borders[0]:borders[1]].reverse_complement()
                else:
                    exon_seq = record.seq[borders[0]:borders[1]]
                entry = f"{header}\n{exon_seq}"
                ortho_fasta.append(entry)

        with open(fasta_path, 'w') as f:
                f.write('\n'.join(ortho_fasta))

