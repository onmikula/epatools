### TOOLS
library(ape)
source("func/epatools.R")
write.delim <- function(x, file) write.table(x, file, quote=FALSE, row.names=FALSE, sep="\t")


### INPUT
# .jplace file (Matsen et al. 2012) containing results of EPA on cytochrome b sequences of silvery mole-rat (Heliophobius; see Uhrová et al. 2022)
jplacefile <- "data/Heliophobius_MPE2022.jplace"

### PIPLINE
# 1 # import of the .jplace file
jplace <- read_jplace(jplacefile)

# 2 # rooting of the tree (and modification of branch labels), which requires outgroup sequences to be specified (here, by the identifier "outgroup")
jplace <- root_jplace(jplace, outgroup="outgroup")

# 3 # classification of the reference tree branches into the lineages defined in the data frame 'lineages'
lineages <- read.delim("data/heliophobius_lineages.txt")
jplace <- classify_jplace(jplace, lineages, ancestral=TRUE)

# 4 # plotting the reference tree and highlighting branches with the most probable placement of query sequences (or, possibly, a subset of them) 
plot(jplace, cex=0.7, no.margin=TRUE)

# 5 # classification of query sequences into the lineages
classification <- classify_sequences(jplace)
write.delim(classification, ""results/heliophobius_lineages_plus.txt"")


### REFERENCES
# Matsen FA, Hoffman NG, Gallagher A, Stamatakis A (2012) A format for phylogenetic placements. PLoS ONE 7: e31009. https://doi.org/10.1371/journal.pone.0031009
# Uhrová M, Mikula O, Bennett NC, Van Daele P, Piálek L, Bryja J, Visser JH, van Vuuren BJ, Šumbera R (2022) Species limits and phylogeographic structure in two genera of solitary African mole-rats Georychus and Heliophobius. Molecular Phylogenetics and Evolution, 167: 107337. https://doi.org/10.1016/j.ympev.2021.107337
