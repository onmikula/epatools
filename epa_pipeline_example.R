### TOOLS
library(ape)
source("func/epatools.R")
source("func/divergentseq.R")
write.fasta <- function(x, file) ape::write.dna(x, file, format="fasta")
write.delim <- function(x, file) write.table(x, file, quote=FALSE, row.names=FALSE, sep="\t")


### FILE NAMES
# input file names
# seqfile contains sequences for ML inference of the reference phylogenetic tree as well as query sequences to be classified
# lineagefile is tab-delimited file with classification of tips into lineages
seqfile <- "data/Heliophobius_MPE2022.fasta"
lineagefile <- "data/heliophobius_lineages.txt"

# the core of output file names
# the next line extracts the string automatically from the name of sequence file, but you can supply your own,
# e.g., name <- "helio"
name <- tolower(substr(sub("^.*\\/", "", seqfile), 1, 5))

# output file names
epaname <- paste(name, "epa", sep="_")
treefile <- paste("raxml/RAxML_bestTree", name, sep=".")
jplacefile <- paste("raxml/RAxML_portableTree", paste0(name,"_epa"), "jplace", sep=".")
finaloutput <- sub("^data", "results", sub("\\.txt", "_plus.txt", lineagefile))

### CLADES & OUTGROUPS
# loading 'lineagefile' with definitions of clades
lineages <- read.delim(lineagefile)
# outgroup sequences defined by either a vector of sequence names or their unique identifier(s), e.g., a character string "outgroup" 
outgroup <- "outgroup"


### REFERENCE SEQUENCES
# preparing refseqfile with sequences for the reference tree
# if you skip this part, you have to supply your own alignmemnt of reference sequences and specify its file name,
# e.g., refseqfile <- "data/Heliophobius_MPE2022_ref.fasta"
refseqfile <- paste0(sub("\\..*$", "", seqfile), "_ref.", sub("^.*\\.", "", seqfile))
format <- ifelse(substr(readLines(seqfile, n=1), 1, 1) == ">", "fasta", "interleaved")
refseq <- ape::read.dna(seqfile, format=format)
outnam <- rownames(refseq)[grep(paste(outgroup, collapse="|"), rownames(refseq))]
seqnam <- lineages$tip
if (is.null(seqnam)) seqnam <- lineages[,1]
refseq <- refseq[c(intersect(rownames(refseq), seqnam), outnam),]
write.fasta(refseq, refseqfile)


### SUBSET OF THE MOST DIVERGENT SEQUENCES
# if the number of sequences is exceedingly high which makes inference of the tree time-consuming and unreliable, its possible to retain a subset of most divergent sequences and move the rest among query sequences
# refseqfile_subs is a name for the file with the subset of sequences
refseqfile_subs <- sub("_ref", "_subs", refseqfile)
divseq <- divergentseq(refseqfile, n=15, outgroup=outgroup, model="raw")
write.fasta(divseq, refseqfile_subs)
refseqfile <- refseqfile_subs


### PIPELINE
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
jplace <- read_jplace(jplacefile)
jplace <- root_jplace(jplace, outgroup=outgroup)
jplace <- classify_jplace(jplace, lineages, ancestral=TRUE)

# plotting the reference tree with branches highlighted where the most probable placement of query sequences (or a subset of them) 
plot(jplace, cex=0.7, no.margin=TRUE)

# classification of query sequences according to EPA output
classification <- classify_sequences(jplace)
write.delim(classification, finaloutput)

