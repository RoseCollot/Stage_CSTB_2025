import gffutils
import pandas as pd
import os
import argparse

def create_database(path_to_gff,path_to_database):
    """create a database from a gff file"""
    if not os.path.exists(path_to_database):    #create a new database or read the existing one
        db = gffutils.create_db(path_to_gff, path_to_database, 
                                merge_strategy="create_unique", keep_order=True, checklines=10,force=False, force_gff=True )
    else : 
        db = gffutils.FeatureDB(path_to_database)
    return db
 
def dataframe_gff(database):
        """Create a dataframe drom a database and delete duplicates and exons from pseudogenes and mitochrondial genome """
        exon_id = []
        transcript_id = []
        gene_id = []
        exon_start = []
        exon_end = []
        exon_length = []
        exon_strand = []
        chromosome_id = []
        for exon in database.features_of_type('exon'):  #select only exons
                if len(list(database.parents(exon, featuretype='mRNA'))) == 0:  
                        continue
                else : 
                        exon_id.append(exon.id)
                        exon_start.append(exon.start)
                        exon_end.append(exon.end) 
                        exon_length.append(exon.end - exon.start)
                        exon_strand.append(exon.strand)
                        chromosome_id.append(exon.seqid)
                        for t in database.parents(exon, featuretype='mRNA'):   
                                transcript_id.append(t.id)      #get the name of the transcrip
                                for g in database.parents(t, featuretype='gene') : 
                                        gene_id.append(g.id)    #get the name of the gene
        df = pd.DataFrame(list(zip(chromosome_id,exon_id,transcript_id,gene_id,exon_start,exon_end,exon_length,exon_strand)), 
                          columns = ['chromosome','exon','transcript','gene','start', 'end', 'length', 'strand'])
 
        df = df.drop_duplicates(subset=['start','end','strand','chromosome','gene'])
        df = df.sort_values(by=['chromosome','start'])

        pseudogenes_id = []
        for p in database.features_of_type('pseudogene'): 
                for id in p.attributes['Dbxref']:  
                        pseudogenes_id.append(id)       #make a list of all pseudogenes id
        exons_to_drop = []
        for e in database.features_of_type('exon'): 
                if 'Dbxref' in e.attributes: 
                        for id in e.attributes['Dbxref']: 
                                if id in pseudogenes_id:
                                        exons_to_drop.append(e.id)      #make a llist of all exons id if exons are from pseudogenes
        df = df.drop(df[df['exon'].isin(exons_to_drop)].index)  #drop these exons

        MT_id= None
        for r in database.features_of_type('region') :
                if r.id.startswith('NC') :
                        for name in r.attributes['Name']:
                                if name == 'MT':
                                        MT_id = r.seqid         #get the mitochrondrial genome id 
        if MT_id is not None :
                df = df.drop(df[df['chromosome']==MT_id].index)
        return df

def merge_overlap_exons(df):
    """Merge overlapping exons into one"""
    exon_start = []
    exon_end = []
    exon_chromosome = []
    length = []
    exon_strand = []
    exons= []
    nb_overlap = []
    gene_list=[]
    for chromosome in df['chromosome'].unique():  
        df_chromosome = df[df['chromosome'] == chromosome]
        for strand in df_chromosome['strand'].unique(): 
            df_strand = df_chromosome[df_chromosome['strand'] == strand]
            start = None
            end = None
            for g in df_strand['gene'].unique():
                exons_list=[]
                df_gene = df_strand[df_strand['gene'] == g]
                for i in range(len(df_gene)):
                    if start is None :  
                        start = df_gene.iloc[i]['start']
                        end = df_gene.iloc[i]['end']
                    if df_gene.iloc[i]['start'] <= end  and df_gene.iloc[i]['end'] <= end: 
                            exons_list.append(str(df_gene.iloc[i]['exon']))
                            continue
                    if df_gene.iloc[i]['start'] <= end and df_gene.iloc[i]['end'] > end : 
                            exons_list.append(str(df_gene.iloc[i]['exon']))
                            end = df_gene.iloc[i]['end']
                    else:
                        exon_start.append(start)
                        exon_end.append(end)
                        exon_chromosome.append(chromosome)
                        length.append(end-start)
                        exon_strand.append(strand)
                        exons.append(exons_list)
                        nb_overlap.append(len(exons_list))
                        exons_list = []
                        gene_list.append(g)
                        start = df_gene.iloc[i]['start'] 
                        end = df_gene.iloc[i]['end']
                        exons_list.append(str(df_gene.iloc[i]['exon']))
                    if i == len(df_gene)-1:               
                        exon_start.append(start)
                        exon_end.append(end)
                        exon_chromosome.append(chromosome)
                        length.append(end - start)
                        exon_strand.append(strand)
                        exons.append(exons_list)
                        gene_list.append(g)
                        nb_overlap.append(len(exons_list))
                        exons_list = []
                            

    df_exons = pd.DataFrame(list(zip(exon_chromosome, exon_start, exon_end, length, exon_strand, exons, nb_overlap, gene_list)),
                            columns=('chromosome', 'start', 'end', 'length', 'strand','exons', 'nb_overlapping_exons', 'gene'))
    df_exons = df_exons.drop_duplicates(subset=['chromosome','start','end','strand'])
    df_exons = df_exons.drop(df_exons[df_exons['length']<1].index)
    return df_exons

def introns_dataframe(df_exons):
    """Create introns from the exons borders"""
    intron_start= []
    intron_end = []
    intron_length = []
    intron_gene = []
    intron_chr = []
    for gene in df_exons['gene'].unique() :
        df_gene = df_exons[df_exons['gene']==gene]
        df_gene = df_gene.reset_index(drop=True)
        if len(df_gene)>1 :
            for i in range(len(df_gene)-1):
                    j=i+1
                    start = df_gene.iloc[i]['end']      #start of the intron is the end of the last exons
                    end = df_gene.iloc[j]['start']      #end of the intron is the start of the next exons
                    intron_start.append(start)
                    intron_end.append(end)
                    intron_length.append(end-start)
                    intron_gene.append(str(gene))
                    chromosome = str(df_gene['chromosome'].unique())
                    intron_chr.append(str(chromosome))
                    

    df_introns = pd.DataFrame(list(zip(intron_start,intron_end,intron_length,intron_gene,intron_chr)), 
                              columns=['start', 'end', 'length', 'gene', 'chromosome'])
    df_introns = df_introns.drop_duplicates(subset=['start','end','gene'])
    return df_introns

def stats_gene(database,df, df_introns):
    """Create a dataframe containing statistics for each gene"""
    genes = []                                            
    for g in df['gene'].unique() :                    
        df_gene = df[df['gene']==g]                    
        gene_row = []                                 
        gene_row.append(g)   
        gene_length = (int(database[g].end) - int(database[g].start))
        gene_row.append(gene_length)                                      
        nb_exons = len(df_gene['exon'])
        gene_row.append(nb_exons)
        exons_length = sum(df_gene['length'])
        gene_row.append(exons_length)
        genes.append(gene_row)
    df_exons_genes = pd.DataFrame(genes, columns = ['gene','length','nb_exons','exons_tot_length'])
    genes = []                                            
    for g in df_introns['gene'].unique() :                    
        df_gene = df_introns[df_introns['gene']==g]                    
        gene_row = []                                  
        gene_row.append(g)                                 
        nb_introns = len(df_gene)
        gene_row.append(nb_introns)
        introns_length = sum(df_gene['length'])
        gene_row.append(introns_length)
        genes.append(gene_row)
    df_introns_genes = pd.DataFrame(genes, columns = ['gene','nb_introns', 'introns_tot_length'])
    df_stats = df_exons_genes.merge(df_introns_genes, how='left')
    df_stats['ratio_introns']=(df_stats['introns_tot_length']/df_stats['length'])
    df_stats['ratio_exons']=(df_stats['exons_tot_length']/df_stats['length'])
    return df_stats

def df_to_csv (path_to_gff,path_to_db, path_to_df): 
    """Main function, convert all dataframes into csv files"""
    if not os.path.exists(path_to_df):
        os.makedirs(path_to_df)         #create a directory with the genome's name if it doesn't exist
    db = create_database(path_to_gff, path_to_db)
    df = dataframe_gff(db)
    df_exons = merge_overlap_exons(df)
    df_introns = introns_dataframe(df_exons)
    df_stats = stats_gene(db, df,df_introns)
    df.to_csv(os.path.join(path_to_df,'df'))
    df_exons.to_csv(os.path.join(path_to_df,'df_exons'))
    df_introns.to_csv(os.path.join(path_to_df,'df_introns'))
    df_stats.to_csv(os.path.join(path_to_df,'df_stats'))

if __name__=='__main__' : 
    parser = argparse.ArgumentParser()
    parser.add_argument("path_to_gff")
    parser.add_argument("path_to_db")
    parser.add_argument("path_to_df")
    args = parser.parse_args()
    df_to_csv(args.path_to_gff,args.path_to_db, args.path_to_df)
