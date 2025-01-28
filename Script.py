# importing data

import pandas as pd
path_AvsB = 'C:/Users/44794/PycharmProjects/pythonProject_Differential_gene_expression_analysis/A_vs_B.deseq2.results.tsv'
path_AvsC = 'C:/Users/44794/PycharmProjects/pythonProject_Differential_gene_expression_analysis/A_vs_C.deseq2.results.tsv'

df_AvsB = pd.read_csv('A_vs_B.deseq2.results.tsv', sep='\t')
df_AvsC = pd.read_csv('A_vs_C.deseq2.results.tsv', sep='\t')

print(df_AvsB.head())
print(df_AvsC.head())

# summary stats - Number of upregulated and downregulated genes

# First need to define the threshold value for significance
threshold_pval = 0.05
threshold_log2pval = 1

# selecting genes in the padj column which are less than the threshold val (0.05)

sig_genes_AvsB = df_AvsB[df_AvsB['pvalue'] < threshold_pval]
sig_genes_AvsC = df_AvsC[df_AvsC['pvalue'] < threshold_pval]

# upregulated genes for A vs B
upreg_AvsB = sig_genes_AvsB[sig_genes_AvsB['log2FoldChange'] > threshold_log2pval]
# downregulated genes for A vs B
downreg_AvsB = sig_genes_AvsB[sig_genes_AvsB['log2FoldChange'] < threshold_log2pval]

#printing values for A vs B
print(F"Number of upregulated genes in A vs B = {len(upreg_AvsB)}")
print(F"Number of downregulated genes in A vs B = {len(downreg_AvsB)}")


# upregulated genes for A vs C
upreg_AvsC = sig_genes_AvsC[sig_genes_AvsC['log2FoldChange'] > threshold_log2pval]
# downregulated genes for A vs B
downreg_AvsC = sig_genes_AvsC[sig_genes_AvsC['log2FoldChange'] < threshold_log2pval]

#printing values for A vs C
print(F"Number of upregulated genes in A vs C = {len(upreg_AvsC)}")
print(F"Number of downregulated genes in A vs C = {len(downreg_AvsC)}")


# summary stats 1.1 - summary of p values and log fold changes across all genes
# p values summary A vs B
pval_summary_AvsB = {
    "pval_mean_AvsB": df_AvsB['pvalue'].mean(),
    "pval_median_AvsB": df_AvsB['pvalue'].median(),
    "pval_min_AvsB": df_AvsB['pvalue'].min(),
    "pval_max_AvsB": df_AvsB['pvalue'].max()
}
print(pval_summary_AvsB)

#log2fold change summary A vs B
L2FC_summary_AvsB = {
"L2FC_mean_AvsB": df_AvsB['log2FoldChange'].mean(),
"L2FC_median_AvsB": df_AvsB['log2FoldChange'].median(),
"L2FC_min_AvsB": df_AvsB['log2FoldChange'].min(),
"L2FC_max_AvsB": df_AvsB['log2FoldChange'].max()
}
print(L2FC_summary_AvsB)

# p-values summary A vs C
pval_summary_AvsC = {
    "pval_mean_AvsC": df_AvsC['pvalue'].mean(),
    "pval_median_AvsC": df_AvsC['pvalue'].median(),
    "pval_min_AvsC": df_AvsC['pvalue'].min(),
    "pval_max_AvsC": df_AvsC['pvalue'].max()
}
print(pval_summary_AvsC)

# log2 fold change summary A vs C
L2FC_summary_AvsC = {
    "L2FC_mean_AvsC": df_AvsC['log2FoldChange'].mean(),
    "L2FC_median_AvsC": df_AvsC['log2FoldChange'].median(),
    "L2FC_min_AvsC": df_AvsC['log2FoldChange'].min(),
    "L2FC_max_AvsC": df_AvsC['log2FoldChange'].max()
}
print(L2FC_summary_AvsC)

#checking data quality A vs B
expression_values_AvsB = []
for val in df_AvsB:
    try:
        expression_values_AvsB.append(float(val))
    except ValueError:
        print(f'Skipping invalid data point: {val}')

print('Processed values:', expression_values_AvsB)

#checking data quality A vs C
expression_values_AvsC = []
for val in df_AvsC:
    try:
        expression_values_AvsC.append(float(val))
    except ValueError:
        print(f'Skipping invalid data point: {val}')

print('Processed values:', expression_values_AvsC)

#converting non-numeric values and handling missing values

df_AvsB = df_AvsB.apply(pd.to_numeric, errors='coerce') #replaces any strings in the data frame to NaN
df_AvsB = df_AvsB.fillna(0)  # convert everything with NaN to 0

df_AvsC = df_AvsC.apply(pd.to_numeric, errors='coerce')
df_AvsC = df_AvsC.fillna(0)

#########
##plots##
#########

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Histogram of p values A vs B
plt.hist(df_AvsB['pvalue'], bins=10, color='blue')
plt.xlabel('p-value')
plt.ylabel('Frequency')
plt.title('Histogram p values for A vs B')
plt.show()

# Histogram of p values A vs c
plt.hist(df_AvsC['pvalue'], bins=10, color='blue')
plt.xlabel('p-value')
plt.ylabel('Frequency')
plt.title('Histogram p values for A vs C')
plt.show()

#volvano plot for A vs B
df_AvsB['neg_log_pvalue'] = -np.log10(df_AvsB['pvalue']) #convert the p value to neg log value
plt.scatter(df_AvsB['log2FoldChange'], df_AvsB['neg_log_pvalue'], c='red', alpha=0.3)
plt.title('Volcano Plot')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10(p-value)')
plt.show()
#volvano plot for A vs C
df_AvsC['neg_log_pvalue'] = -np.log10(df_AvsC['pvalue']) #convert the p value to neg log value
plt.scatter(df_AvsC['log2FoldChange'], df_AvsC['neg_log_pvalue'], c='red', alpha=0.3)
plt.title('Volcano Plot')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10(p-value)')
plt.show()

#Heatmap
#merging both conditions and adding suffix to each gene ID
merged_ABC_conditions = pd.merge(upreg_AvsB, upreg_AvsC, on="gene_id", suffixes=("_AvsB", "_AvsC"))
#filtering heatmap data for variables of interest
Heatmap_data = merged_ABC_conditions[["gene_id", "log2FoldChange_AvsB", "log2FoldChange_AvsC"]]
#setting the index to 'gene ID'
Heatmap_data.set_index("gene_id", inplace=True)
#Generate heat map
sns.heatmap(Heatmap_data, cmap="viridis", annot=True, cbar=True)
# adding axis labels
plt.title("Heatmap of Log2Fold Changes for both AvsB and AvsC conditions")
plt.xlabel("Conditions")
plt.ylabel("Gene ID")
plt.show()

#MA plot to for the relationship between log fold change and mean expression for AvsB
df_AvsB['log_baseMean_AvsB'] = np.log10(df_AvsB['baseMean']) #Making new column of log transformed base mean values
sns.scatterplot(x='log_baseMean_AvsB', y='log2FoldChange', data=df_AvsB, color='Blue', alpha=0.9)
plt.title('MA Plot for A vs B')
plt.xlabel('log Mean Expression')
plt.ylabel('Log2FoldChange')
plt.show()
#MA plot to for the relationship between log fold change and mean expression for AvsC
df_AvsC['log_baseMean_AvsC'] = np.log10(df_AvsC['baseMean']) #Making new column of log transformed base mean values
sns.scatterplot(x='log_baseMean_AvsC', y='log2FoldChange', data=df_AvsC, color='Blue', alpha=0.9)
plt.title('MA Plot A vs C')
plt.xlabel('log Mean Expression')
plt.ylabel('Log2FoldChange')
plt.show()

# Table of sig upregulated and downregulated genes

new_upreg_AvsB = upreg_AvsB[['log2FoldChange', 'pvalue', 'padj']]
new_upreg_AvsB.to_csv('upreg_AvsB.csv', index=False)
new_downreg_AvsB = downreg_AvsB[['log2FoldChange', 'pvalue', 'padj']]
new_downreg_AvsB.to_csv('downreg_AvsB.csv', index=False)
new_upreg_AvsC = upreg_AvsC[['log2FoldChange', 'pvalue', 'padj']]
new_upreg_AvsC.to_csv('upreg_AvsC.csv', index=False)
new_downreg_AvsC = downreg_AvsC[['log2FoldChange', 'pvalue', 'padj']]
new_downreg_AvsC.to_csv('downreg_AvsC.csv', index=False)


