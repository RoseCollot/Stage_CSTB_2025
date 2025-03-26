import gffutils
import pandas as pd
import numpy as np
import os
import argparse

def create_database(path_to_gff,path_to_database):
    if not os.path.exists(path_to_database):
        db = gffutils.create_db(path_to_gff, path_to_database, 
                                merge_strategy="create_unique", keep_order=True, checklines=10,force=False, force_gff=True )
    else : 
        db = gffutils.Features(path_to_database)
    return db
 
def dataframe_gff(database):
        exon_id = []
        transcript_id = []
        gene_id = []
        exon_start = []
        exon_end = []
        exon_length = []
        exon_strand = []
        chromosome_id = []
        for exon in database.features_of_type('exon'):
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
                                transcript_id.append(t.id)
                                for g in database.parents(t, featuretype='gene') :
                                        gene_id.append(g.id)
        df = pd.DataFrame(list(zip(chromosome_id,exon_id,transcript_id,gene_id,exon_start,exon_end,exon_length,exon_strand)), 
                          columns = ['chromosome','exon','transcript','gene','start', 'end', 'length', 'strand'])
  
        #delete duplicates 
        df = df.drop_duplicates(subset=['start','end','strand','chromosome','gene'])
        df = df.sort_values(by=['chromosome','start'])

        #delete exons froms pseudogenes
        pseudogenes_id = []
        for p in database.features_of_type('pseudogene'): #pour tous les pseudogènes
                for id in p.attributes['Dbxref']: #chercher id dbref 
                        pseudogenes_id.append(id)
        exons_to_drop = []
        for e in database.features_of_type('exon'): #pour tous les exons
                if 'Dbxref' in e.attributes: 
                        for id in e.attributes['Dbxref']: 
                                if id in pseudogenes_id:
                                        exons_to_drop.append(e.id)
        df = df.drop(df[df['exon'].isin(exons_to_drop)].index)

        #delete exons from MT genome 
        MT_id= None
        for r in database.features_of_type('region') :
                if r.id.startswith('NC') :
                        for name in r.attributes['Name']:
                                if name == 'MT':
                                        MT_id = r.seqid
        if MT_id is not None :
                df = df.drop(df[df['chromosome']==MT_id].index)
        return df

def merge_overlap_exons(df):
    exon_start = []
    exon_end = []
    exon_chromosome = []
    length = []
    exon_strand = []
    exons= []
    nb_overlap = []
    gene=[]
    for chromosome in df['chromosome'].unique():  
        df_chromosome = df[df['chromosome'] == chromosome]
        for strand in df_chromosome['strand'].unique():  
            df_strand = df_chromosome[df_chromosome['strand'] == strand]
            exons_list=[]
            start = 1
            end = 1
            for i in range (len(df_strand)):
                    if df_strand.iloc[i]['start'] <= end  and df_strand.iloc[i]['end'] <= end: 
                            exons_list.append(str(df_strand.iloc[i]['exon']))
                    if df_strand.iloc[i]['start'] <= end and df_strand.iloc[i]['end'] > end : 
                            end = df_strand.iloc[i]['end']
                            exons_list.append(str(df_strand.iloc[i]['exon']))
                    else:
                            if i >1 : 
                                exon_start.append(start)
                                exon_end.append(end)
                                exon_chromosome.append(chromosome)
                                length.append(end-start)
                                exon_strand.append(strand)
                                exons.append(exons_list)
                                nb_overlap.append(len(exons_list))
                                exons_list.clear()
                                gene.append(df_strand.iloc[i-1]['gene'])
                                start = df_strand.iloc[i]['start']
                                end = df_strand.iloc[i]['end']
                                exons_list.append(str(df_strand.iloc[i]['exon']))
                    if i > 1:
                        exon_start.append(start)
                        exon_end.append(end)
                        exon_chromosome.append(chromosome)
                        length.append(end-start)
                        exon_strand.append(strand)
                        exons.append(exons_list)
                        nb_overlap.append(len(exons_list))
                        exons_list.clear()
                        gene.append(df_strand.iloc[i-1]['gene'])
                            

    df_exons = pd.DataFrame(list(zip(exon_chromosome, exon_start, exon_end, length, exon_strand, exons, nb_overlap, gene)),
                            columns=('chromosome', 'start', 'end', 'length', 'strand','exons', 'nb_overlapping_exons', 'gene'))
    df_exons = df_exons.drop_duplicates(subset=['chromosome','start','end','strand'])
    df_exons = df_exons.drop(df_exons[df_exons['length']<1].index)
    return df_exons

def introns_dataframe(df_exons):
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
                    start = df_gene.iloc[i]['end']
                    end = df_gene.iloc[j]['start']
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

def stats_gene(df, df_introns):
    genes = []                                            #contient gene_uniq
    for g in df['gene'].unique() :                    #goes through each unique gene
        df_gene = df[df['gene']==g]                    
        gene_row = []                                  #row of the database (one for each gene --> 20 000)
        gene_row.append(g)                          #row get gene id        
        nb_exons = len(df_gene['exon'])
        gene_row.append(nb_exons)
        exons_length = sum(df_gene['length'])
        gene_row.append(exons_length)
        genes.append(gene_row)
    df_exons_genes = pd.DataFrame(genes, columns = ['gene', 'nb_exons', 'exons_tot_length'])
    genes = []                                            #contient gene_uniq
    for g in df_introns['gene'].unique() :                    #goes through each unique gene
        df_gene = df_introns[df_introns['gene']==g]                    
        gene_row = []                                  #row of the database (one for each gene --> 20 000)
        gene_row.append(g)                          #row get gene id        
        nb_introns = len(df_gene)
        gene_row.append(nb_introns)
        introns_length = sum(df_gene['length'])
        gene_row.append(introns_length)
        genes.append(gene_row)
    df_introns_genes = pd.DataFrame(genes, columns = ['gene', 'nb_introns', 'introns_tot_length'])
    df_stats = df_exons_genes.merge(df_introns_genes, how='left')
    df_stats['ratio_introns']=(df_stats['introns_tot_length']/df_stats['exons_tot_length'])
    return df_stats

def df_to_csv (path_to_gff,path_to_db, path_to_df): 
    if not os.path.exists(path_to_df):
        os.makedirs(path_to_df)
    db = create_database(path_to_gff, path_to_db)
    df = dataframe_gff(db)
    df_exons = merge_overlap_exons(df)
    df_introns = introns_dataframe(df_exons)
    df_stats = stats_gene(df,df_introns)
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
