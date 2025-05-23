import pandas as pd
import csv
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

def hierarchical_clustering(path_to_matrix_csv, path_to_new_tree):
    with open(path_to_matrix_csv) as matrix:
        reader = csv.reader(matrix)
        rows = list(reader)
        species = rows[0][1:]  
        data_values = [row[1:] for row in rows[1:]]  

    data = pd.DataFrame(data_values, columns=species, index=species)
    data = data.apply(pd.to_numeric)

    condensed_matrix = squareform(data.values)
    linkage_data = linkage(condensed_matrix, method='ward', metric='euclidean', optimal_ordering=True)


    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    dn = dendrogram(linkage_data, labels=data.index.tolist(), orientation='left', above_threshold_color='purple')
    for label in ax.get_yticklabels():
            label.set_fontstyle('italic')
    plt.ylabel('Species')
    plt.title('Hierarchical clustering of the intron ratio')
    plt.tight_layout()  
    plt.savefig(path_to_new_tree, dpi=300)
    plt.show()
