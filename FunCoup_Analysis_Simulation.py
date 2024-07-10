#!/usr/bin/env python
# coding: utf-8

# In[78]:


# load the zebrafish functional couplings from here: https://funcoup.org/downloads/
# citations here (most recent is ok): https://funcoup.org/help/#Citation
# also check work done here for support: https://www.nature.com/articles/ncomms5071

import pandas as pd

funcoup_all = pd.read_csv("FC5.0_D.rerio_full", sep = "\t")

# filter functional couplings with probability functional coupling (pfc) >= 90%

funcoup_90 = funcoup_all[(funcoup_all["#0:PFC"] >= 0.9)]

# keep only the gene columns

funcoup_90_genes = df_filter.iloc[:, 2:4]

funcoup_90_genes.head(2)


# In[79]:


funcoup_90_genes.shape


# In[80]:


# retrieve from BioMart the mappings between Gene stable ID <-> Gene name
biomart_mappings = pd.read_csv("mart_export.txt", sep = "\t")

# load the 98 candidate genes by Gene name
candidate_genes = ["hsd17b7", "serbp1a", "il12rb2", "LOC110438375", "LOC103909477", "cyp2p10", "cyp2p7", "cyp2n13", "cyp2p6", "hook1", 
                   "si:dkey-183n20.15", "pfas", "ipo13b", "muc5.1", "LOC562098", "mmachc", "lrrc4ca", "abhd2b", "rlbp1b", "isg20", 
                   "pdia3", "ckmt1", "ticrr", "si:ch1073-281m9.1", "alkbh3", "hsd17b12a", "LOC559196", "CBX1", "nfe2l1b", "D5F01_LYC18948", 
                   "copz2", "prr15la", "LOC110437935", "hipk3b", "imp3", "slc35a3a", "fam78ba", "cmpk", "bcl6ab", "zmym4.2", 
                   "cldn1", "cpn2", "fkbp1b", "her6", "atp13a3", "ncl", "gk5", "atp1b3a", "grk7a", "fbxo36b", 
                   "agfg1b", "si:dkeyp-13a3.3", "mffa", "stk25b", "vps11", "hyou1", "si:ch73-261i21.5", "hist2h2l", "h2ax1", "hspa8b",
                   "LOC110438375", "hspa8b", "jhy", "bsx", "lim2.1", "smarcc1a", "si:ch211-215k15.4", "dpy19l1l", "map7b", "hbs1l",
                   "armc1l", "mtfr2", "slc39a9", "pde7a", "myb", "aldh8a1", "txnrd2.2", "gnb1l", "tbx1", "si:ch211-244o22.2",
                   "kat6a", "ap3m2", "tet3", "bicdl1", "rab35a", "gcn1", "vamp5", "vamp8", "ncaph", "creb1b",
                   "mettl21a", "ccnyl1", "etnk1", "pyroxd1", "iapp", "slco1e1", "slco1d1", "itgbl1"]  # is "pign" one of the significant genes ?

# filter biomart mappings to contain only the candidate genes - 95 mapped genes
biomart_mappings_candidate_genes = biomart_mappings[biomart_mappings["Gene name"].isin(candidate_genes)]

# use the ensembl identifiers as a list
biomart_mappings_candidate_genes_ensembl = list(biomart_mappings_candidate_genes["Gene stable ID"])


# In[81]:


# filter the >0.9 interactions to involve only the candidate genes - 7357 interactions

funcoup_90_candidate_genes = funcoup_90_genes[funcoup_90_genes['2:Gene1'].isin(biomart_mappings_candidate_genes_ensembl) 
                                            | funcoup_90_genes['3:Gene2'].isin(biomart_mappings_candidate_genes_ensembl)]


# In[82]:


# make a list of all the intermediate genes

intermediate_genes = []

for line in range(len(funcoup_90_candidate_genes)):
    if line % 500 == 0:
        print(line, end = " ")
    if funcoup_90_candidate_genes.iloc[line, 0] not in biomart_mappings_candidate_genes_ensembl:
        intermediate_genes.append(funcoup_90_candidate_genes.iloc[line, 0])
    if funcoup_90_candidate_genes.iloc[line, 1] not in biomart_mappings_candidate_genes_ensembl:
        intermediate_genes.append(funcoup_90_candidate_genes.iloc[line, 1])

# count the number of occurences of intermediate genes

intermediate_gene_counts = pd.Series(intermediate_genes).value_counts()
max_threshold = intermediate_gene_counts.max()


# In[142]:


# for every level of occurence of intermediate genes (for example 14, 13, 12 ... occurences)
# caclulate the (a) number of intermediate genes, and (b) number of candidate genes, and
# estimate the enrichment level of candidate genes vs. total genes (intermediate & candidate)

p_values = []
threshold_ = []

for threshold in range(max_threshold):
    
    intermediate_genes_to_keep = intermediate_gene_counts[intermediate_gene_counts >= threshold+1].index.tolist()

    funcoup_90_intermediate_genes_to_keep = funcoup_90_genes[funcoup_90_genes['2:Gene1'].isin(intermediate_genes_to_keep) 
                                                           | funcoup_90_genes['3:Gene2'].isin(intermediate_genes_to_keep)]

    candidate_genes_subnetwork = []

    for line in range(len(funcoup_90_intermediate_genes_to_keep)):
        if funcoup_90_intermediate_genes_to_keep.iloc[line, 0] in biomart_mappings_candidate_genes_ensembl:
            candidate_genes_subnetwork.append(funcoup_90_intermediate_genes_to_keep.iloc[line, 0])
        if funcoup_90_intermediate_genes_to_keep.iloc[line, 1] in biomart_mappings_candidate_genes_ensembl:
            candidate_genes_subnetwork.append(funcoup_90_intermediate_genes_to_keep.iloc[line, 1])

    from scipy.stats import fisher_exact

    # Create a 2x2 contingency table
    cand_genes = len(set(candidate_genes_subnetwork))
    tot_genes = len(set(intermediate_genes_to_keep)) + cand_genes
    cand_genes_ref = len(set(biomart_mappings_candidate_genes_ensembl))
    tot_genes_ref = len(set(intermediate_genes)) + cand_genes_ref

    contingency_table = [[cand_genes, tot_genes],
                         [cand_genes_ref, tot_genes_ref]]

    # Perform Fisher's exact test
    odds_ratio, p_value = fisher_exact(contingency_table)

    
    print(f"Threshold: {threshold+1}; Interactions: {len(funcoup_90_intermediate_genes_to_keep)}; Cand_Genes: {cand_genes}; Int_Genes: {len(set(intermediate_genes_to_keep))}; P-val: {p_value}")
#     print("Odds ratio:", odds_ratio)

    # store p-values as a numpy array
    p_values.append(p_value)
    threshold_.append(threshold)
p_values = np.array(p_values)
threshold_ = np.array(threshold_)


# In[162]:


# ==========
# SIMULATION 
# ==========

# Step#1: I need to draw a random selection of 95 choices from the "funcoup_90_genes"

import numpy as np

SIMS = 1000

p_values_sim = []
threshold_sim = []

enumerate_ = 0

while enumerate_ < SIMS:
    
    print(enumerate_, end = " ")

    random_elements_col1 = funcoup_90_genes['2:Gene1'].unique()
    random_elements_col2 = funcoup_90_genes['3:Gene2'].unique()

    random_elements = list(set(random_elements_col1).union(set(random_elements_col2)))

    random_sample = np.random.choice(random_elements, 95, replace=False)

    funcoup_90_random_genes = funcoup_90_genes[funcoup_90_genes['2:Gene1'].isin(random_sample) 
                                             | funcoup_90_genes['3:Gene2'].isin(random_sample)]

    # make a list of all the intermediate genes

    intermediate_genes_random = []

    for line in range(len(funcoup_90_random_genes)):
#         if line % 500 == 0:
#             print(line, end = " ")
        if funcoup_90_random_genes.iloc[line, 0] not in random_sample:
            intermediate_genes_random.append(funcoup_90_random_genes.iloc[line, 0])
        if funcoup_90_random_genes.iloc[line, 1] not in random_sample:
            intermediate_genes_random.append(funcoup_90_random_genes.iloc[line, 1])

    # count the number of occurences of intermediate genes

    intermediate_random_gene_counts = pd.Series(intermediate_genes_random).value_counts()
    max_threshold_permutation = intermediate_random_gene_counts.max()

    # for every level of occurence of intermediate genes (for example 14, 13, 12 ... occurences)
    # caclulate the (a) number of intermediate genes, and (b) number of candidate genes, and
    # estimate the enrichment level of candidate genes vs. total genes (intermediate & candidate)



    for threshold in range(max_threshold_permutation):

        intermediate_random_genes_to_keep = intermediate_random_gene_counts[intermediate_random_gene_counts >= threshold+1].index.tolist()

        funcoup_90_intermediate_random_genes_to_keep = funcoup_90_random_genes[funcoup_90_random_genes['2:Gene1'].isin(intermediate_random_genes_to_keep) 
                                                                      | funcoup_90_random_genes['3:Gene2'].isin(intermediate_random_genes_to_keep)]

        candidate_random_genes_subnetwork = []

        for line in range(len(funcoup_90_intermediate_random_genes_to_keep)):
            if funcoup_90_intermediate_random_genes_to_keep.iloc[line, 0] in biomart_mappings_candidate_genes_ensembl:
                candidate_random_genes_subnetwork.append(funcoup_90_intermediate_random_genes_to_keep.iloc[line, 0])
            if funcoup_90_intermediate_random_genes_to_keep.iloc[line, 1] in biomart_mappings_candidate_genes_ensembl:
                candidate_random_genes_subnetwork.append(funcoup_90_intermediate_random_genes_to_keep.iloc[line, 1])

        from scipy.stats import fisher_exact

        # Create a 2x2 contingency table
        cand_random_genes = len(set(candidate_random_genes_subnetwork))
        tot_random_genes = len(set(intermediate_random_genes_to_keep)) + cand_random_genes
        cand_genes_ref = len(set(biomart_mappings_candidate_genes_ensembl))
        tot_genes_ref = len(set(intermediate_genes)) + cand_genes_ref

        contingency_table_perm = [[cand_random_genes, tot_random_genes],
                                 [cand_genes_ref, tot_genes_ref]]

#       Perform Fisher's exact test
        odds_ratio, p_value = fisher_exact(contingency_table_perm)


#         print(f"Threshold: {threshold+1}; Interactions: {len(funcoup_90_intermediate_random_genes_to_keep)}; Cand_Genes: {cand_random_genes}; Int_Genes: {len(set(intermediate_random_genes_to_keep))}; P-val: {p_value}")
#         print("Odds ratio:", odds_ratio)

        p_values_sim.append(p_value)
        threshold_sim.append(threshold)
    p_values_sim.append(np.nan)
    threshold_sim.append(np.nan)
    enumerate_ += 1
        
p_values_sim = np.array(p_values_sim)
threshold_sim = np.array(threshold_sim)


# In[167]:


# Create the plot
plt.figure(figsize=(12, 6))

# Plot simulated p-values
plt.plot(threshold_sim, p_values_sim, label='Simulated P-values', marker='o', linestyle='--', color='green')

# Plot the given p-values
plt.plot(threshold_, p_values, label='Observed P-values', marker='o', linestyle='-', color='blue')

# Add labels and title
plt.xlabel('Degree')
plt.ylabel('P-value')
plt.title('Comparison of Observed and Simulated P-values')
plt.yscale('log')  # Log scale for better visualization of p-values
plt.legend()
plt.grid(True)

# Show the plot
plt.show()


# In[208]:


# Create the plot
fig, ax = plt.figure(figsize=(12, 8), facecolor='aliceblue'), plt.gca()
ax.set_facecolor('aliceblue')

# Plot simulated p-values
plt.plot(threshold_sim, np.log2(p_values_sim), label='Simulated P-values', marker='o', linestyle='--', color='green')

# Plot the given p-values
plt.plot(threshold_, np.log2(p_values), label='Observed P-values', marker='o', linestyle='-', color='blue', 
         linewidth=4, markersize=10, markerfacecolor='white')

# Add labels and title with increased font size
plt.xlabel('Degree', fontsize=18, fontweight='bold', labelpad=20)
plt.ylabel(r'Log$_2$(P-value)', fontsize=16, fontweight='bold', labelpad=20)
plt.title('Comparison of\nObserved and Simulated Enrichment in Target Genes', 
          fontsize=24, fontweight='bold', pad=16)
# plt.yscale('log')  # Log scale for better visualization of p-values

# Set x-axis limit to show only up to degree 14
plt.xlim(0, 14)

# Customize the legend
plt.legend(fontsize=12, title_fontsize='13')

# Add grid lines for better readability
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Customize ticks
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Remove top and right borders
sns.despine()

# Save the plot
plt.savefig('Simulations_Interactions.png', dpi=300, bbox_inches='tight')

# Show the plot
plt.tight_layout()  # Adjust the plot to ensure everything fits without overlapping
plt.show()


