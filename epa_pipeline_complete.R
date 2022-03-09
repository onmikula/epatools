### TOOLS
library(ape)
source("func/epatools.R")
source("func/divergentseq.R")
write.fasta <- function(x, file) ape::write.dna(x, file, format="fasta")
write.delim <- function(x, file) write.table(x, file, quote=FALSE, row.names=FALSE, sep="\t")


### FILE NAMES
# input file names
# sequence file contains sequences for ML inference of the reference phylogenetic tree as well as query sequences to be classified
# lineage file is tab-delimited file with classification of sequences (and hence tree tips) into taxa (lineages / species / clades)
# the sequences classified in the lineage file are considered to be the reference ones, 
# but another reference sequence file can be provided and it may be even necessary if calssification to taxa is not of interest
# and the goal is just to place query sequences to the reference tree
seqfile <- "data/Heliophobius_MPE2022.fasta"
lineagefile <- "data/heliophobius_lineages.txt"

# the core of output file names
# the next line extracts the string automatically from the name of sequence file, but you can supply your own, e.g., name <- "helio"
name <- tolower(substr(sub("^.*\\/", "", seqfile), 1, 5))
epaname <- paste(name, "epa", sep="_")

# output file names
treefile <- paste("raxml/RAxML_bestTree", name, sep=".")
jplacefile <- paste("raxml/RAxML_portableTree", epaname, "jplace", sep=".")
classificationfile <- sub("^data", "results", sub("\\.txt", "_plus.txt", lineagefile))


### CLADES & OUTGROUPS
# loading 'lineagefile' with definitions of clades
lineages <- read.delim(lineagefile)
# outgroup sequences defined by either a vector of sequence names or their unique identifier(s), e.g., a character string "outgroup" 
outgroup <- "outgroup"


### REFERENCE SEQUENCES
# preparing reference sequence file (containing sequences for the reference tree)
# if no specific refseqfile is provided, the reference sequence are those that are present in the sequence file and classifed in the lineage file
# if you skip this part, you have to supply your own alignmemnt of reference sequences and specify its file name,
# e.g., refseqfile <- "data/Heliophobius_MPE2022_ref.fasta"
refseqfile <- paste0(sub("\\..*$", "", seqfile), "_ref.", sub("^.*\\.", "", seqfile))
format <- ifelse(substr(readLines(seqfile, n=1), 1, 1) == ">", "fasta", "interleaved")
refseq <- ape::read.dna(seqfile, format=format)
outgroupnames <- rownames(refseq)[grep(paste(outgroup, collapse="|"), rownames(refseq))]
seqnames <- lineages$tip
if (is.null(seqnames)) seqnames <- lineages[,1]
refseq <- refseq[c(intersect(rownames(refseq), seqnames), outgroupnames),]
refseq <- refseq[!duplicated(rownames(refseq)),]
write.fasta(refseq, refseqfile)


### SUBSET OF THE MOST DIVERGENT SEQUENCES
# an optional modfication of the reference sequences
# if the number of sequences is exceedingly high which makes inference of the tree time-consuming and unreliable, its possible to retain a subset of most divergent sequences (here n=15) and move the rest among query sequences
# refseqfile_subs is a name for the file with the subset of sequences
refseqfile_subs <- sub("_ref", "_subs", refseqfile)
divseq <- divergentseq(refseqfile, n=15, outgroup=outgroup, model="raw")
write.fasta(divseq, refseqfile_subs)
refseqfile <- refseqfile_subs


### RAxML -> EPATOOLS PIPELINE
# ML tree inference
cmd <- ifelse(.Platform$OS.type == "unix", "./bin/raxmlHPC-SSE3", "./bin/raxmlHPC.exe")
system2(cmd, args=c("-f", "a", "-m", "GTRGAMMA", "--HKY85", "-N", "100", "-x", "12345", "-p", "12345", "-s", refseqfile, "-n", name))
outputs <- list.files(pattern=paste0("RAxML.*\\.", name))
invisible(file.rename(outputs, paste("raxml", outputs, sep="/")))

# EPA inference of ML phylogenetic placements of query sequences, i.e., all those in 'seqfile', but not 'refseqfile'
# if you skip ML inference and use another reference tree instead, its file name has to be specified
# e.g., treefile <- "data/MrBayes_helio.txt"
# note, you can modify RAxML parameters (e.g., --epa-prob-threshold) - for this, consult documentation of RAxML (https://cme.h-its.org/exelixis/web/software/raxml)
cmd <- ifelse(.Platform$OS.type == "unix", "./bin/raxmlHPC-SSE3", "./bin/raxmlHPC.exe")
system2(cmd, args=c("-f", "v", "-m", "GTRGAMMA", "--HKY85", "--epa-keep-placements=100", "--epa-prob-threshold=0.01", "-s", seqfile, "-t", treefile, "-n", epaname))
outputs <- list.files(pattern=paste0("RAxML.*\\.", epaname))
invisible(file.rename(outputs, paste("raxml", outputs, sep="/")))

# loading of EPA output, rooting of the tree and classification of tree branches to clades
# if an output of externally run phylogenetic placement analysis is post-processed, the name of its .jplace file has to be specified
# e.g., jplacefile <- "raxml/RAxML_portableTree.helio_epa.jplace"

# 1 # import of the .jplace file
jplace <- read_jplace(jplacefile)
# 2 # rooting of the tree (and modification of branch labels) by the specified outgroup(s)
jplace <- root_jplace(jplace, outgroup=outgroup)
# 3 # classification of the reference tree branches into the lineages defined in the data frame 'lineages' (not required, if just placement to the tree is of intrest)
jplace <- classify_jplace(jplace, lineages, ancestral=TRUE)
# $ # plotting the reference tree and highlighting branches with the most probable placement of query sequences (or, possibly, a subset of them)
plot(jplace, cex=0.7, no.margin=TRUE)
# 5 # classification of query sequences into the lineages (if the classification and not only phylogenetic placement itself is of interest)
classification <- classify_sequences(jplace)
write.delim(classification, classificationfile)

