import pandas as pd
import argparse

def table (path_to_scpecies_table, path_to_new_species_df):
    row_df=[]
    df = pd.read_csv(path_to_scpecies_table, sep='\t', names=['query', 'subject', 'event', 'query_coverage', 'gene', 'query_species', 'subject_species'])
    df = df.drop_duplicates()
    df['query_coverage'] = pd.to_numeric(df['query_coverage'], errors='coerce')
    for gene in df['gene'].unique():
        df_gene = df[df['gene']==gene]
        row = []
        nb_whole = (pd.isna(df_gene['event']) & (df_gene['query_coverage'] >= 80)).sum()
        nb_other = (pd.isna(df_gene['event']) & (df_gene['query_coverage'] < 80)).sum()
        nb_fusion = (df_gene['event'] == 'fusion').sum()
        nb_fission = (df_gene['event'] == 'fission').sum()
        nb_full_duplication = ((df_gene['event'] == 'duplication') & (df_gene['query_coverage'] >= 80)).sum()
        nb_other += ((df_gene['event'] == 'duplication') & (df_gene['query_coverage'] < 80)).sum()
        row.append(gene)
        row.append(df_gene['query_species'].iloc[0])
        row.append(int(len(df_gene['query'].unique())))
        row.append(int(len(df_gene)))
        row.append(int(nb_whole))
        row.append(int(nb_fusion))
        row.append(int(nb_fission))
        row.append(int(nb_full_duplication))
        row.append(int(nb_other))
        row_df.append(row)
    gene_ali_df = pd.DataFrame(row_df, columns= ['gene', 'species', 'nb_exons','nb_hits', 'no_event', 'fusion', 'fission', 'full_duplication', 'others'])      
    with open(path_to_new_species_df, 'w') as s:
        gene_ali_df.to_csv(s, sep='\t', index=False)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("path_to_species_table")
    parser.add_argument("path_to_new_species_df")
    args = parser.parse_args()
    table(args.path_to_species_table, args.path_to_new_species_df)




    

