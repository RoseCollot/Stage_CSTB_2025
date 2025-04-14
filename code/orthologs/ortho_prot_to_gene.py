import pandas as pd
import gffutils
import os
import json

#make a json dictionaries for each species containing prot_id and gene (gene-genename)
for species in os.listdir('/home/collot/stage/collot/collot/out_stats/databases/'):
    species_dir = str(os.path.join('/home/collot/stage/collot/collot/out_stats/databases/' + species))
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
    file_path = '/home/collot/stage/collot/collot/out_ortho/dictionaries/'+ species[3:16]+'.json'
    with open(file_path, 'w') as fp:
        json.dump(cds_genes, fp)


#group all species dictionaries into one (key to species dictionaries = genome_id) 
dictionaries = {}
for ortho_dict in os.listdir('/home/collot/stage/collot/collot/out_ortho/dictionaries/'):
    dict_path = '/home/collot/stage/collot/collot/out_ortho/dictionaries/' + ortho_dict
    with open(dict_path) as prot_dict:
        dictionaries[ortho_dict[:-5]] = json.load(prot_dict)


#get orthologs dataframe with prot_id   
df_prot_orth = pd.read_csv('/home/collot/stage/collot/collot/out_ortho/prot_oto_homo_sapiens.tsv', sep='\t', index_col=0)
df_prot_orth = df_prot_orth.dropna(ignore_index=True)

df_species = pd.read_csv('/home/collot/stage/collot/collot/out_stats/df_species.tsv',sep='\t')
    
#df with HS gene and orthologs genes 
df_genes_ortho = pd.DataFrame(columns=['main', 'orthologs'])
for index, row in df_prot_orth.iterrows() : 
    query_prot = df_prot_orth.loc[index,'protein_id']
    for ortho_dict_name, prot in dictionaries.items(): 
            tax_id = str(df_species.loc[df_species['Assembly Accession'].str[:13]==ortho_dict_name, 'Organism Taxonomic ID'].iloc[0])
            query_gene = prot[query_prot] + '_' + tax_id
            df_genes_ortho.loc[index,'main'] = query_gene

    ortho_list = []
    prot_ortho_list=list(str(df_prot_orth.loc[index,'orthologs_protein_id']).split(sep=','))
    for i in range(len(prot_ortho_list)): 
        for ortho_dict_name, prot in dictionaries.items(): 
            if prot_ortho_list[i] in prot : 
                tax_id = str(df_species.loc[df_species['Assembly Accession'].str[:13]==ortho_dict_name, 'Organism Taxonomic ID'].iloc[0])
                gene_ortho = prot[prot_ortho_list[i]] + '_' + tax_id
                ortho_list.append(gene_ortho)
    df_genes_ortho.loc[index, 'orthologs']=ortho_list

pd.DataFrame.to_csv(df_genes_ortho,'/home/collot/stage/collot/out_ortho/genes_orthologs_homo_sapiens.csv' )