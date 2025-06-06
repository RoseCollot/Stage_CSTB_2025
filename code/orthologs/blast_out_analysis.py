import pandas as pd
from Bio import SeqIO
import argparse

def blast_to_df (path_to_blast_out, path_to_fasta_in):
    """Process the .out files of BLAST alignments"""
    blast = pd.read_csv(path_to_blast_out, sep='\t', header=None)
    blast.columns = ["query", "subject", "identity", "alignment_length", "mismatches", "gap_opens","q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
    mask = blast['query'] == blast['subject']   #create a mask to skip exons aligned with themselves
    blast = blast[~mask] 
    blast['query_species'] = blast['query'].apply(lambda x: x.split('_')[-1])   #only get the last part of the query id : the taxa
    blast['subject_species'] = blast['subject'].apply(lambda x: x.split('_')[-1])
    blast['gene'] = blast['query'].apply(lambda x: x.split('_')[0])+'_'+blast['query'].apply(lambda x: x.split('_')[-1])

    exon_length = {}
    with open(path_to_fasta_in, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            exon = record.description
            exon_length[exon]=len(record.seq)   #create the tuple (exon,total length of the exon)
    query_length = []
    for e in blast['query'].values:
        if e in exon_length:
            query_length.append(exon_length[e])
    subject_length = []
    for e in blast['subject'].values:
        if e in exon_length:
            subject_length.append(exon_length[e])

    blast['query_length']=query_length
    blast['subject_length']=subject_length
    blast['query_coverage']= (blast['alignment_length']/blast['query_length'])*100

    reciprocal_hits=[]
    #create a dataframe containign only reciprocal hits
    for index, row in blast.iterrows():
        que = row['query']
        sub = row['subject']
        if que in blast['subject'].values:
            match = blast[blast['subject']==que]
            if sub in match['query'].values:
                reciprocal_hits.append(row)
    if len(reciprocal_hits)!=len(blast):
        blast = pd.DataFrame(reciprocal_hits)

    return blast

def merge_duplicates (blast):
    """merge hits when there is multiple alignment between the same two exon"""
    dup=blast.duplicated(subset=['query', 'subject'], keep=False)   #select only the hits from the same exons
    df_dup=blast[dup]
    merged_list=[]
    grouped_duplicates= df_dup.groupby(['query', 'subject'])     
    for (query, subject), group in grouped_duplicates:
        row1=group.iloc[0]  #first hit
        row2=group.iloc[1]  #second hit
        merged_row = row1   #keep the values of the first hit for the other columns 
        merged_row['q_start']=min(row1['q_start'], row2['q_start'])
        merged_row['q_end']=min(row1['q_end'], row2['q_end'])
        merged_row['s_start']=min(row1['s_start'], row2['s_start'])
        merged_row['s_end']=min(row1['s_end'], row2['s_end'])
        merged_row['alignment_length']=merged_row['q_end']-merged_row['q_start']
        merged_row['query_coverage']=(merged_row['alignment_length']/merged_row['query_length'])*100
        merged_list.append(merged_row)
    merged_hit=pd.DataFrame(merged_list)
    
    blast_no_duplicates=blast.drop_duplicates(subset=['query', 'subject'], keep=False)
    blast_final=pd.concat([merged_hit, blast_no_duplicates])
    blast_final=blast_final.sort_values(by='query', ignore_index=True)

    return blast_final

def find_exon_event (blast_final):
    def check_duplicates(blast, species_column,check_column):
        duplicates = []
        for s in blast[species_column].unique():  
            dup = blast[blast[species_column] == s]
            dup = dup.sort_values(by=[check_column])
            if dup.duplicated(subset=check_column).sum() > 0:
                duplicates.append(dup[dup.duplicated(subset=check_column, keep=False)])
        if len(duplicates) > 0:
            return pd.concat(duplicates)
        else:
            return pd.DataFrame(columns=['query', 'subject'])
    df_one_to_many = pd.DataFrame()
    df_many_to_one = pd.DataFrame()
    df_one_to_many = pd.concat([df_one_to_many, check_duplicates(blast_final, species_column='subject_species', check_column='query')])
    df_many_to_one = pd.concat([df_many_to_one, check_duplicates(blast_final, species_column='query_species', check_column='subject')])

    exons=[]
    for index, row in blast_final.iterrows():
        exon_row=[]
        exon_row.append(row['query'])
        exon_row.append(row['subject'])
        if ((row['query'], row['subject']) in df_one_to_many[['query', 'subject']].itertuples(index=False, name=None)) : #if the exon is aligned with multiple exons from the same species
            exon_row.append('fusion')
        elif ((row['query'], row['subject']) in df_many_to_one[['query', 'subject']].itertuples(index=False, name=None)) : #if mutiple exons from this species are aligned against one unique exon
            overlapping=False
            other_exons = df_many_to_one[(df_many_to_one['subject'] == row['subject'])&(df_many_to_one['query_species'] == row['query_species'])] #select all exons aligned with the same one
            for index, other_row in other_exons.iterrows():
                if other_row['query'] != row['query']:  
                    overlap = max(0, min(row['s_end'], other_row['s_end']) - max(row['s_start'], other_row['s_start']))     
                    if overlap > 9 :   #arbitrary value for overlap 
                        overlapping = True
                        break
            if overlapping==True:
                exon_row.append("duplication")  #event is called duplication if the two alignment overlap 
            else:
                exon_row.append("fission")
        else : 
            exon_row.append(None)
        exon_row.append(row['query_coverage'])
        exon_row.append(row['gene'])
        exon_row.append(row['query_species'])
        exon_row.append(row['subject_species'])
        exons.append(exon_row)
    df_changes = pd.DataFrame(exons, columns=['query','subject','event','query_coverage', 'gene', 'query_species', 'subject_species'])
    return df_changes

def save_df (path_to_tsv, df_changes):
    with open(path_to_tsv, 'w') as t:
        df_changes.to_csv(t, sep='\t', index=False)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("path_to_blast_out")
    parser.add_argument("path_to_fasta_in")
    parser.add_argument("path_to_tsv")
    args = parser.parse_args()
    blast=blast_to_df(args.path_to_blast_out, args.path_to_fasta_in)
    blast_final= merge_duplicates(blast)
    df_changes = find_exon_event(blast_final)
    save_df(args.path_to_tsv, df_changes)
