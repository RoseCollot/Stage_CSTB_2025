import pandas as pd
import gffutils
import os
import json

def dictionaries (path_to_databases, path_to_dictionaries):
    """Make a json dictionarie for each species containing the protein id and gene (format : gene-genename)"""
    for species in os.listdir(path_to_databases):
        species_dir = str(os.path.join(path_to_databases + species))
        db = str(species[:16])  #db name is db_GCF
        db = gffutils.FeatureDB(species_dir)
        cds_genes = str(species[3:16]) #cds dict name is GCF
        cds_genes = {}
        for cds in db.features_of_type('CDS'): 
            for gene_name in cds.attributes["gene"]:
                if cds.id.startswith('cds-'):
                    key = cds.id[4:]
                else: 
                    key = cds.id 
                cds_genes[key] = gene_name       
        file_path = path_to_dictionaries + species[3:16]+'.json'
        with open(file_path, 'w') as fp:
            json.dump(cds_genes, fp)

def group_dict (path_to_dictionaries) :
    """Group all species dictionaries into one (the key to each species dictionary is the genome_id) """
    dictionaries = {}
    for ortho_dict in os.listdir(path_to_dictionaries):
        dict_path = path_to_dictionaries + ortho_dict
        with open(dict_path) as prot_dict:
            dictionaries[ortho_dict[:-5]] = json.load(prot_dict)
    return dictionaries

def create_orthologous_genes_df (path_to_ortho_proteins, path_to_species_info_df, path_to_new_csv) :   
    df_prot_orth = pd.read_csv(path_to_ortho_proteins, sep='\t', index_col=0)   #dataframe containing orthologous proteins (first columns : protein used as query, second : list of all the orthologs) 
    df_species = pd.read_csv(path_to_species_info_df,sep='\t')  #dataframe with informations on species

    df_genes_ortho = pd.DataFrame(columns=['main', 'orthologs'])    #create an empty dataframe
    for index, row in df_prot_orth.iterrows() : 
        query_prot = df_prot_orth.loc[index,'protein_id']   
        for ortho_dict_name, prot in dictionaries.items(): 
                tax_id = str(df_species.loc[df_species['Assembly Accession'].str[:13]==ortho_dict_name, 'Organism Taxonomic ID'].iloc[0])
                query_gene = prot[query_prot] + '_' + tax_id    #add the taxa id to the gene name
                df_genes_ortho.loc[index,'main'] = query_gene   

        ortho_list = []
        prot_ortho_list=list(str(df_prot_orth.loc[index,'orthologs_protein_id']).split(sep=','))    #correction of the list format
        for i in range(len(prot_ortho_list)): 
            for ortho_dict_name, prot in dictionaries.items(): 
                if prot_ortho_list[i] in prot : 
                    tax_id = str(df_species.loc[df_species['Assembly Accession'].str[:13]==ortho_dict_name, 'Organism Taxonomic ID'].iloc[0])
                    gene_ortho = prot[prot_ortho_list[i]] + '_' + tax_id
                    ortho_list.append(gene_ortho)
        df_genes_ortho.loc[index, 'orthologs']=ortho_list

    pd.DataFrame.to_csv(df_genes_ortho,path_to_new_csv)