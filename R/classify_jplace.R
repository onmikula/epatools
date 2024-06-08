#' Classification of `jplace` object.
#' 
#' @description
#' Classifies tree edges with any phylogenetic placements to clades (or bipartitions).
#'
#' @param jp An object of class `jplace`.
#' @param clade A data frame defining classification of tips (1st column) to clades / bipartitions
#'   (2nd column). The clades can be mutually nested.
#' @param ancestral Logical, whether to classify ancestral branches; The default option is `TRUE`,
#'   but it is relevant only if the tree is rooted and defined clades are comprehensive
#'   and non-overlapping (otherwise switched to `FALSE` and issues a warning).
#' @returns An amended `jplace` object with additional components `clade` and `position`
#'   in the `placements` data frame, which indicate classification of the branches.
#' @export
classify_jplace <- function(jp, clade, ancestral=TRUE) {
	getpaths <- function(phy, from, to) {
		lapply(to, function(x) ape::nodepath(phy, from=from, to=x))
	}
	find_common_branch <- function(tip, phy) {
		ifelse(length(tip) == 1, match(tip, phy$tip.label), ape::getMRCA(phy, tip))
	}
	tree <- jp$tree
	place <- jp$placements
	clade <- setNames(as.data.frame(clade[,1:2]), c("tip", "clade"))
	clade <- clade[clade$tip %in% jp$tree$tip.label,]
	clade <- split(clade$tip, clade$clade)
	anc <- lapply(clade, find_common_branch, phy=tree)
	off <- lapply(clade, match, table=tree$tip.label)
	paths <- Map(getpaths, list(tree), anc, off)
	paths <- lapply(lapply(lapply(paths, unlist), unique), sort)
	paths <- Map(setdiff, paths, anc)
	###
	crowns <- setNames(lapply(paths, function(x) tree$edge.label[match(x, tree$edge[,2])]), names(clade))
	crown_lengths <- sapply(crowns, length)
	crowns <- data.frame(edge=unlist(crowns), clade=rep(names(crowns), crown_lengths), position=rep("crown", sum(crown_lengths)), row.names=NULL)	
	stems <- setNames(tree$edge.label[match(unlist(anc), tree$edge[,2])], names(clade))
	stems <- data.frame(edge=stems, clade=names(stems), position="stem", row.names=NULL)
	outgroup_tips <- as.numeric(sapply(tree$outgroup, grep, x=tree$tip.label))
	overlapping <- any(duplicated(unlist(paths)))
	partial <- any(!setdiff(seq(ape::Ntip(tree)), outgroup_tips) %in% unlist(off))
	unrooted <- !ape::is.rooted(tree)
	if (ancestral & (overlapping | unrooted | partial)) {
		ancestral <- FALSE
		message <- c("defined clades are overlapping", "the classification is not comprehensive", "the tree is not rooted")
		message <- message[c(overlapping, partial, unrooted)]
		if (length(message) > 1) {
			message <- paste(paste(message[1:(length(message) - 1)], collapse=", "), message[length(message)], sep=" and ")	
		}
		message <- paste0(toupper(substr(message, 1, 1)), substr(message, 2, nchar(message)), ",")
		warning(paste(message, "'ancestral' argument was set to FALSE."))
	}
	if (ancestral) {
		outgroup_paths <- getpaths(tree, from=ape::Ntip(tree)+1, to=outgroup_tips)
		root_neighbors <- tree$edge[tree$edge[,1] == ape::Ntip(tree) + 1,2]
		excluded <- c(unlist(anc), unlist(paths), unlist(outgroup_paths), root_neighbors)
		mrca <- setdiff(seq(ape::Nnode(tree))[-1] + ape::Ntip(tree), excluded)
		if (length(mrca) > 0) {
			ancestral_tips <- lapply(mrca, function(x) ape::extract.clade(tree, node=x)$tip.label)
			ancestral_tips <- lapply(ancestral_tips, match, table=tree$tip.label)
			ancestral_paths <- Map(getpaths, list(tree), as.list(mrca), ancestral_tips)
			ancestral_paths <- lapply(lapply(lapply(ancestral_paths, unlist), unique), sort)
			superspecies <- lapply(ancestral_paths, function(x) sapply(anc, "%in%", x))
			superspecies <- lapply(superspecies, function(x) names(x[x]))
			superspecies <- sapply(superspecies, paste, collapse="_")
			ancestors <- setNames(tree$edge.label[match(mrca, tree$edge[,2])], superspecies)
			ancestors <- data.frame(edge=ancestors, clade=names(ancestors), position="ancestral", row.names=NULL)
		} else {
			ancestors <- data.frame(edge=character(0), clade=character(0), position=character(0))
		}
	} else {
		ancestors <- data.frame(edge=character(0), clade=character(0), position=character(0))
	}
	edges <- do.call(rbind, list(crowns, stems, ancestors))
	ii <- match(place$edge_num, edges$edge)
	place$clade <- edges$clade[ii]
	place$position <- edges$position[ii]
	jp$placements <- place
	attr(jp, "classified") <- TRUE
	class(jp) <- "jplace"
	return(jp)
}
