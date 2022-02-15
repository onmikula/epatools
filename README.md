# epatools
R functions post-processing outputs of phylogenetic placement analysis

**epatools** contains a set of R functions (R Core Team 2021) post-processing an output of phylogenetic placement inference, performed by *Evolutionary Placement Algorithm* (EPA, Berger et al. 2011) or *pplacer* (Matsen et al. 2010).

The functions perform five tasks:
1. read_jplace imports the output file in .jplace format (Matsen et al. 2012)
2. root_jplace roots the tree (and modifies branch labels accordingly), using specified outgroup sequences
3. classify_jplace classifies branches as belonging to clades, which are defined in a suitable data frame
4. plot.jplace shows the reference tree and highlights the most probable placements of query sequences
5. classify_sequences calculates probabilities of query sequences belonging to pre-specified clades

The functions depend on R package ape (Paradis & Schliep 2019). 'plot.jplace' can save the plot as a .pdf file, the output of 'classify_sequences' can be saved as a tab-delimited file.

Two example pipelines are provided.

The **epa_pipeline_fivesteps.R** demonstrates just the five aforementioned steps, starting with pre-estimated .jplace file. 

The **epa_pipeline_complete.R** takes a sequence file and their classification into lineages as its inputs and includes also selection of the reference sequences, inference of the reference tree and EPA in RAxML (Stamatakis 2014). The pipeline is modular, however, so you can skip parts that are not relevant for your purpose. Namely, you can supply your own:
- pre-prepared set of reference sequences
- pre-estimated reference tree
- .jplace file from externally run EPA or pplacer analysis (as in 'epa_pipeline_fivesteps.R')

Note also that the focus here is on classification of query sequences to taxa, but the very purpose of EPA / pplacer is the phylogenetic placement of query sequences to the tree. If just this placement is of interest, the file defining the taxa is not necessary as classify_jplace and classify_sequences steps are not performed. In such case, however, reference sequences or reference tree have to be supplied. 

The provided data are cytochrome b sequences of silvery mole-rats (Heliophobius) presented by Uhrová et al. (2022). Here, they are classified into three clades: northern (N), southeastern (SE) & southwestern (SW). 

**References:**
- Berger SA, Krompass D, Stamatakis A (2011) Performance, accuracy, and web server for evolutionary placement of short sequence reads under maximum likelihood. Systematic Biology, 60: 291–302. https://doi.org/10.1093/sysbio/syr010
- Matsen FA, Hoffman NG, Gallagher A, Stamatakis A (2012) A format for phylogenetic placements. PLoS ONE 7: e31009. https://doi.org/10.1371/journal.pone.0031009
- Matsen FA, Kodner R, Armbrust E (2010) pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC Bioinformatics 11: 538. https://doi.org/10.1186/1471-2105-11-538
- Paradis E, Schliep K (2019) ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35, 526–528. https://doi.org/10.1093/ bioinformatics/bty633
- Stamatakis A (2014) RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30: 1312– 1313. https://doi.org/10.1093/bioinformatics/btu033
- R Core Team (2021) R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org
- Uhrová M, Mikula O, Bennett NC, Van Daele P, Piálek L, Bryja J, Visser JH, van Vuuren BJ, Šumbera R (2022) Species limits and phylogeographic structure in two genera of solitary African mole-rats Georychus and Heliophobius. Molecular Phylogenetics and Evolution, 167: 107337. https://doi.org/10.1016/j.ympev.2021.107337
