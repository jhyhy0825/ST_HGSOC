#!/home/jhy/anaconda2/envs/cpdb/bin/python3

## /home/jhy/spatial_transcriptomics/220905/code/CellPhoneDB/CellPhoneDB_plot.py

import anndata
import pandas as pd
import sys
import os
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

#cpdb_results = cpdb_degs_analysis_method.call(cpdb_file_path='/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/cellphonedb.zip',meta_file_path='/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/Recur_meta.txt',counts_file_path='/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/Recur_count.txt',degs_file_path='/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/Recur_markers_new.txt',counts_data='hgnc_symbol',output_path='/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/output_recur_DEG_0.2/',threshold = 0.2)

cpdb_results = cpdb_statistical_analysis_method.call(cpdb_file_path = '/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/cellphonedb.zip',meta_file_path='/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/Recur_meta.txt',counts_file_path='/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/Recur_count.txt',counts_data = 'hgnc_symbol',output_path = '/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/output_recur_0.2/',threshold=0.2)

annotation = list(cpdb_results['relevant_interactions'].columns[:13])
interaction = list(cpdb_results['relevant_interactions'].columns[13:])

relevant_interactions_long = pd.melt(cpdb_results['relevant_interactions'],id_vars = annotation,var_name = 'Interacting_cell',value_vars = interaction,value_name = 'Relevance')

relevant_interactions_long[['Cell_a', 'Cell_b']] = relevant_interactions_long['Interacting_cell'] \
        .str.split('|', expand = True)\
        .rename(columns={0 : 'Cell_a', 1: 'Cell_b'})
relevant_interactions_long.head(3)

means_long = pd.melt(cpdb_results['means'],
        id_vars = annotation,
        var_name = 'Interacting_cell',
        value_vars = interaction,
        value_name = 'Mean')
means_long.head(3)

id_cp_dict = relevant_interactions_long.groupby('id_cp_interaction')['Relevance'] \
        .sum() \
        .to_dict()

relevant_interactions_long['Recurrence'] = relevant_interactions_long['id_cp_interaction'] \
        .map(id_cp_dict)

relevant_interactions_long = relevant_interactions_long.sort_values(['Recurrence'],
        ascending = True)

relevant_interactions_long.head(3)

interactions_to_plot = ['CPI-SS052ADA696','CPI-SS0DB3F5A37','CPI-SS060B82BBC','CPI-SS0B99F22A2','CPI-SS0263EFDDF']
idx = [ id in interactions_to_plot for id in relevant_interactions_long['id_cp_interaction'] ]
relevant_interactions_plot = relevant_interactions_long[idx].copy()

from scipy.stats import zscore
relevant_interactions_plot = relevant_interactions_plot.merge(means_long[['id_cp_interaction', 'Interacting_cell', 'Mean']],
        on = ['id_cp_interaction', 'Interacting_cell'],
        how = 'inner')

import itertools

list_trophoblast = ['Ovarian cancer cell','Fibroblast','Endothelial cell','T cell','B cell','Plasma cell','Myeloid cell','Dendritic cell','Mast cell']
interact_rm = [ '|'.join(i) for i in itertools.product(list_trophoblast, repeat = 2) ]
idx_keep = [ i not in interact_rm for i in relevant_interactions_plot['Interacting_cell'] ]
relevant_interactions_plot = relevant_interactions_plot[idx_keep]
relevant_interactions_plot['Mean_scaled'] = relevant_interactions_plot.groupby('id_cp_interaction', group_keys = False)['Mean'].transform(lambda x : zscore(x, ddof = 1))
relevant_interactions_plot = relevant_interactions_plot.sort_values('Interacting_cell')

relevant_interactions_plot.head(3)


import seaborn as sns
import matplotlib.pyplot as plt

g = sns.relplot(
        data = relevant_interactions_plot,
        x = "Interacting_cell",
        y = "interacting_pair",
        hue = "Relevance",
        size = "Mean_scaled",
        palette = "vlag",
        hue_norm=(-1, 1),
        height = 6,
        aspect = 4,
        sizes = (0, 200)
        facet_kws=dict(gridspec_kws={"hspace": 0.5})
)
g.set_xticklabels(rotation = 45)

#plt.figure(figsize=(10,4))
plt.savefig('/home/jhy/spatial_transcriptomics/220905/result/seurat/CellPhoneDB/test8.png')
plt.close()

