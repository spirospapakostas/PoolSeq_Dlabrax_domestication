Scripts used in the analysis of the article entitled "Genomic Signatures for Domestication in European Seabass (Dicentrarchus labrax L.) 
Underscore a Key Role for Epigenetic Regulation in Adaptation to Captivity" by Moulistanos et al.

<=
Annotations.py
This script parses BioMart annotations for the European seabass (Dicentrarchus labrax) 
using the gff3_parser Python package. It filters the annotations 
to find features within a specified range around a given SNP position.
=>

<=
depth_simulation_analysis.py
This script simulates the selection of representatives from a sample pool 
and calculates the average number of unique representatives 
as well as the lower 95% confidence interval.
=>

<=
find_biallelic_SNP.py
This script identifies biallelic loci in a dataset 
by setting to zero the counts of nucleotides 
that fall below specified thresholds. 
The thresholds are based on the total depth at each position, 
allowing for more stringent criteria at higher depths.
This script then proceeds with identifying and categorizing Single Nucleotide Polymorphisms (SNPs) into biallelic, monomorphic, and triallelic/multiallelic categories. 
=>

<=
FunCoup_Analysis_Simulation.py
This script analyzes functional couplings in zebrafish genes, 
filtering couplings with a high probability of functional coupling (PFC) 
and identifying interactions involving candidate genes. 
It also performs simulations to compare observed and simulated enrichment levels of candidate genes.
=>

<=
FunCoup_Analysis_Simulation.ipynb
The Jupyter notebook of FunCoup_Analysis_Simulation.py
=>
