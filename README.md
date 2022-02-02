# epatools
R functions post-processing outputs of phylogenetic placement analysis

epatools contains a set of R functions (R Core Team 2021) post-processing an output of phylogenetic placement inference, performed by Evolutionary Placement Algorithm (EPA: Berger et al. 2011) or pplacer (Matsen et al. 2010).

The functions perform five tasks:
(1) read_jplace imports the output file in .jplace format (Matsen et al. 2012)
(2) root_jplace roots the tree (and modifies branch labels accordingly), which requires specifying outgroup sequences
(3) classify_jplace classifies branches as belonging to clades defined in a suitable data frame
(4) plot.jplace shows the reference tree, where the most probable placements of specified query sequences are highlighted
(5) classify_sequences calculates probabilities of query sequences belonging to pre-specified clades

The functions depend on R package 'ape' (Paradis & Schliep 2019). plot.jplace can save the plot as a .pdf file, the output of classify_sequences can be saved as a tab-delimited file.

The example pipeline ('epa_pipeline_example.R') includes not only the post-processing itself (= five above-mentioned steps), but also inference of ML tree in RAxML (Stamatakis 2014) and running EPA, so it is possible to start with just two basic inputs, namely the sequence file and the classification file. The sequence file is an alignment of both reference and query sequences in .fasta or interleaved .phylip format. The classification file is a tab-delimited file which assigns sequences to the clades of interest. The pipeline is modular, however, so the user can step in at any point or skip parts that are not relevant to her/him. For instance, you can supply your own pre-prepared set of reference sequences or pre-estimated phylogenetic tree. In such case you have to supply your own 'refseqfile' or 'treefile' names to the script and make it sure the names of sequences and tree tips match each other and sequence names in the sequence file. For other options, see comments in 'epa_pipeline_example.R' and 'func/epatools.R'. In order to run the whole pipeline, you have to place RAxML binary file into the 'bin'.

The provided data are cytochrome b sequences of silvery mole-rats (Heliophobius) presented by Uhrová et al. (2022). Here, they are classified into just three clades: northern (N), southeastern (SE) & southwestern (SW). 

References:
Berger SA, Krompass D, Stamatakis A (2011) Performance, accuracy, and web server for evolutionary placement of short sequence reads under maximum likelihood. Systematic Biology, 60: 291–302. https://doi.org/10.1093/sysbio/syr010
Matsen FA, Hoffman NG, Gallagher A, Stamatakis A (2012) A format for phylogenetic placements. PLoS ONE 7: e31009. https://doi.org/10.1371/journal.pone.0031009
Matsen FA, Kodner R, Armbrust E (2010) pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC Bioinformatics 11: 538. https://doi.org/10.1186/1471-2105-11-538
Paradis E, Schliep K (2019) ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35, 526–528. https://doi.org/10.1093/ bioinformatics/bty633
Stamatakis A (2014) RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30: 1312– 1313. https://doi.org/10.1093/bioinformatics/btu033
R Core Team (2021) R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org
Uhrová M, Mikula O, Bennett NC, Van Daele P, Piálek L, Bryja J, Visser JH, van Vuuren BJ, Šumbera R (2022) Species limits and phylogeographic structure in two genera of solitary African mole-rats Georychus and Heliophobius. Molecular Phylogenetics and Evolution, 167: 107337. https://doi.org/10.1016/j.ympev.2021.107337
