import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#complete table
df_all = pd.read_csv('/home/collot/stage/collot/collot/out_ortho/all_species_df.tsv', sep='\t')

#create a gene column
df_all['gene']=df_all['gene'].apply(lambda x: x.split('_')[-2])

#associate taxonomy id with species anmes
df_table= pd.read_csv('/home/collot/stage/collot/PrimateData/PrimateInfo.csv', sep='\t')
genid_species = dict(zip(df_table['Organism Taxonomic ID'].astype(str), df_table['Organism Name'].astype(str)))  #dictionary containing tuples (genome_id/specie_name)

#turn table into explicit table
event_columns = ['fission', 'fusion', 'full_duplication', 'others'] 
df_all[event_columns] = df_all[event_columns].fillna(0) 
df_long = df_all.melt(id_vars=['gene', 'species', 'nb_exons'], value_vars=event_columns, var_name='event_type', value_name='event_count')
df_long['species'] = df_long['species'].astype(str).map(genid_species)

#order the species in the graph
species_order = df_long.groupby('species')['event_count'].sum().sort_values(ascending=False).index
df_long['species'] = pd.Categorical(df_long['species'], categories=species_order, ordered=True)

# Plot the histogram
fig, ax = plt.subplots(figsize=(10, 10))
sns.histplot(data=df_long, y="species", hue="event_type", weights="event_count", bins=35, multiple="stack", palette="Set3")
plt.title('Number of events')
plt.tight_layout
plt.savefig('/home/collot/stage/collot/collot/out_ortho/all_events.png', bbox_inches='tight')