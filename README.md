# SARS-COV-2-Transmission-Networks
These files analayse SARS-COV-2 sequence data and estimate their corresponding transmission networks. As a consequence, the mutational steps that occur during transmission can also be followed.
An agent-based model replicating the transmission network structure is developed.

To implement the analysis on GISAID data, follow these steps:

Within GISAID, download required sequences. For the analysis associated with this fileset, the sequences were chosen through:
A.	Search:
a.	Location: Oceania/Australia
b.	Variant: VOC Delta GK
c.	Submission: to 31 October 2021
d.	Complete
B.	Download sequences into separate files containing at most 5,000 sequences with “Input for the Augur pipeline” selected.
C.	Extract from each .tar file the separate sequence .fasta and metadata .tsv (.txt) files.
D.	Download reference sequence. In this analysis this is EPI_ISL_412026

File sequence
1.	github_cov_seq_delta_input.m: inputs the chosen Delta sequences and saves them as a structure seqs_delta
2.	github_seq_delta_domains.m: make a first estimate of the 12 ORF domains for each delta sequence through pairwise alignment with the ref sequence;  this will be improved  in the next few steps. Also correct some of the sequences that have missed the first ‘A’ in the start codon through the sequencing procedure.
3.	github_seq_delta_domain_orfs_calcs: estimate improved domain ends for each ORF. Also determine the problem sequences that have mostly been caused by N at one of the domain endpoints. These are listed in seq_group_probs. 
a.	This calls github_seq_delta_domain_orfs 
4.	github_seqs_delta_lineage_domain_grps: Group the sequences by lineage, and then within each lineage, by those with the same length for the ORFS.  
a.	It calls github_domain_length_grps.
5.	github_seq_delta_domain_grps_pair_clone_distance.m: Determine the clones within each group. 
a.	Calls github_ACTG_delta_clones_recur
i.	Calls github_clone_sep
ii.	And github_clone_recur
6.	github_clone_delta_group_day_prob_edges1.m: Calculate the optimal graphs for each group and for their separate clusters. 
a.	Calls github_prob_edge to determine the probability calculation between edges based on date separation and a Gamma distribution.
b.	Calls github_nt_delta_edge1 to determine the nt distance between sequences
c.	Calls github_CLE_adjacency_clone to determine the directed graph over all nodes and edges and to extract weakly connected components.
d.	Calls github_CLE_optimal to generate the optimal arborescence
i.	Calls github_CLE_calc1 that performs the iterations of the arborescence calculation
7.	github_clone_delta_group_day_prob_edges_single_plot.m. Plot individual optimal graphs over a lineage specified by the index j0
a.	Delta_lengths.xlsx. Read in lengths of each ORF for the groups from this spreadsheet, and determines if the group is missing the start of 1a. If so the amount missing will be added when the mutational graphs are drawn.
b.	github_plot_delta_mutpath1: Determines the mutations along the longest mutational pathway and tries to cope when some of the domains cannot be included for all sequences or even a few. 
i.	Calls github_mut_seqs_dom to determine the mutations.
c.	github_cov_state_plot: to show the time course of infection and which states
8.	github_delta_within_clone_Gopt_permute.m: calculate and plot the within optimal networks either without variation (flag_in_permute=0) or with nperm random draws of the dates so that it will produce random networks of this type. 
a.	Calls github_CLE_adjacency_within_clone
10.	github_clone_delta_group_day_prob_edges_single_analyse1.m. Determine the prob dist (Poisson or NBin) over all subgroups within a lineage. Gives the lambda or (R,p) values.
11.	github_outer_clone_network_delta_all: The clonal networks have been determined for the different lineage and length groups. For each of these groups extract their root nodes and multiply align them with the reference sequence. Can use seqalignviewer to show where the mutations and indels occur.


The Clone Agent Model
1.	github_clone_agent_model
