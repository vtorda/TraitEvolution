# ASR_analysis
Some R functions for analyze the output of ancestral state reconstruction (ASR) analyses

The count_parsimony_ASR.R contains the count_parsimony function.

Usage:

tree - a phy object

asr - a data.frame or matrix containing ASR data for every node and tip in the tree object. The order of the rows are important. It has to start with the tip states than the node states, ordered as the phy objectum. It has to be as many columns as many states used in ASR analysis.

occurrence - logical, whether to calculate the number of occurence of a state through the phylogeny

transition_n - logical, whether to calculate the number of transitions happenned between different states through the phylogeny

type - charachter, can be "all", "tips" or "inner_nodes". Choosing to wich portion of the phylogeny should be used during the calculation.


