import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import seaborn.objects as so
from matplotlib.gridspec import GridSpec
from pandas.plotting import table 

def graph_distribution(df_stats_1,df_stats_2, df_stats_3, df_stats_4, path_to_fig_distribution):
    species_1 = pd.read_csv(df_stats_1)
    species_2 = pd.read_csv(df_stats_2)
    species_3 = pd.read_csv(df_stats_3)
    species_4 = pd.read_csv(df_stats_4)

    # Create the figure and grid layout
    fig = plt.figure(figsize=(12, 12))  
    gs = GridSpec(2, 2, figure=fig)  
    title_font = {'weight': 'normal','size': 16}

    # Add each plot to the grid
    ax1 = fig.add_subplot(gs[0, 0])
    sns.histplot(species_1['ratio_introns'], ax=ax1)
    ax1.set_title(r"$\it{Pan\ troglodytes}$", fontdict=title_font)
    ax1.axvline(species_1['ratio_introns'].mean(), color='k', lw=2)
    ax1.axvline(species_1['ratio_introns'].median(), color='k', ls='--', lw=2)
    ax1.tick_params(axis='both', which='major', labelsize=18)
    ax1.set_xlabel("Ratio of Introns", fontsize=16)
    ax1.set_ylabel("Frequency", fontsize=16)

    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
    sns.histplot(species_2['ratio_introns'], ax=ax2)
    ax2.set_title(r"$\it{Otolemur\ garnettii}$", fontdict=title_font)
    ax2.axvline(species_2['ratio_introns'].mean(), color='k', lw=2)
    ax2.axvline(species_2['ratio_introns'].median(), color='k', ls='--', lw=2)
    ax2.tick_params(axis='both', which='major', labelsize=18)
    ax2.set_xlabel("Ratio of Introns", fontsize=16)
    ax2.set_ylabel("Frequency", fontsize=16)

    ax3 = fig.add_subplot(gs[1, 0],sharey=ax1)
    sns.histplot(species_3['ratio_introns'], ax=ax3)
    ax3.set_title(r"$\it{Carlito\ syrichta}$", fontdict=title_font)
    ax3.axvline(species_3['ratio_introns'].mean(), color='k', lw=2)
    ax3.axvline(species_3['ratio_introns'].median(), color='k', ls='--', lw=2)
    ax3.tick_params(axis='both', which='major', labelsize=18)
    ax3.set_xlabel("Ratio of Introns", fontsize=16)
    ax3.set_ylabel("Frequency", fontsize=16)

    ax4 = fig.add_subplot(gs[1, 1],sharey=ax1)
    sns.histplot(species_4['ratio_introns'], ax=ax4)
    ax4.set_title(r"$\it{Homo\ sapiens}$", fontdict=title_font)
    ax4.axvline(species_4['ratio_introns'].mean(), color='k', lw=2)
    ax4.axvline(species_4['ratio_introns'].median(), color='k', ls='--', lw=2)
    ax4.tick_params(axis='both', which='major', labelsize=18)
    ax4.set_xlabel("Ratio of Introns", fontsize=16)
    ax4.set_ylabel("Frequency", fontsize=16)

    # Adjust layout
    plt.tight_layout()
    plt.savefig(path_to_fig_distribution, bbox_inches='tight')

def percentils (df_stats_1,df_stats_2, df_stats_3, df_stats_4, path_to_fig_percentils):
    species_1 = pd.read_csv(df_stats_1)
    species_1 = species_1.dropna()
    species_1["species"] = r"$\it{Homo\ sapiens}$"

    species_2 = pd.read_csv(df_stats_2)
    species_2 = species_2.dropna()
    species_2["species"] = r"$\it{Pan\ troglodytes}$"

    species_3 = pd.read_csv(df_stats_3)
    species_3 = species_3.dropna()
    species_3["species"] = r"$\it{Propithecus\ coquereli}$"

    species_4 = pd.read_csv(df_stats_4)
    species_4 = species_4.dropna()
    species_4["species"] = r"$\it{Otolemur\ garnettii}$"

    selected_species = pd.concat([species_1, species_2, species_3, species_4], ignore_index=True)

    p = (so.Plot(selected_species, x="species", y="ratio_introns")
        .add(so.Dot(color="g"), so.Perc(10))
        .layout(size=(8, 5))
        .label(x="Species", y="Ratio of Introns"))

    p.savefig(path_to_fig_percentils, bbox_inches='tight')
