import pandas as pd
import os
import numpy as np

def df_species_update (path_to_otiginal_csv, path_to_new_tsv): #Update : new genome accession id is replaced by the old one (currently used)
    df_species = pd.read_csv(path_to_otiginal_csv,sep='\t')
    df_species = df_species.sort_values(by='Assembly Accession')
    df_species['Assembly Accession'] = df_species['Assembly Accession'].replace('GCF_043159975.1', 'GCF_000956065.1')
    pd.DataFrame.to_csv(df_species,path_to_new_tsv,sep='\t', index=False)

def create_table (path_to_df_species,path_to_df_dir, path_to_output) : 
    """Create a table with statistics for all species"""
    list_species = []
    df_species = pd.read_csv(path_to_df_species, sep='\t')
    for index, row in df_species.iterrows(): 
        path_dir= os.path.join(path_to_df_dir,row['Assembly Accession'])
        df = pd.read_csv(os.path.join(path_dir,'df'))
        df_exons = pd.read_csv(os.path.join(path_dir,'df_exons'))
        df_introns = pd.read_csv(os.path.join(path_dir,'df_introns'))
        df_stats = pd.read_csv(os.path.join(path_dir,'df_stats'))
        species = []
        species.append(row['Organism Name'])
        species.append(row['Assembly Accession'])
        genes_mean_length = np.mean(df_stats['length'])
        species.append(genes_mean_length)
        nb_genes = len(df['gene'].unique())
        species.append(nb_genes)        
        nb_trancripts = len(df['transcript'].unique())
        species.append(nb_trancripts)        
        nb_exons = len(df_exons)
        species.append(nb_exons) 
        exons_mean_length = np.mean(df_exons['length'])       
        species.append(exons_mean_length)
        exon_ratio = np.mean(df_stats['ratio_exons'])        
        species.append(exon_ratio) 
        nb_introns = len(df_introns)
        species.append(nb_introns)  
        intron_mean_length = np.mean(df_introns['length']) 
        species.append(intron_mean_length)             
        intron_ratio = np.mean(df_stats['ratio_introns'])
        species.append(intron_ratio)
        list_species.append(species)
    table_species = pd.DataFrame(list_species, 
                                 columns=('species_name','genome_id','genes_mean_length','nb_genes','nb_transcripts','nb_exons','exons_mean_length','exon_ratio','nb_introns','introns_mean_length','intron_ratio'))
    table_species = table_species.sort_values(by='genome_id')
    table_species.to_csv(path_to_output)
    return table_species
create_table('stage/collot/collot/out_stats/df_species.tsv', '/home/collot/stage/collot/collot/out_stats/output_dataframes','stage/collot/collot/out_stats/table_species.csv' )
