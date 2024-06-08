#' Rooting `jplace` object.
#' 
#' @description
#' Roots the tree & adjusts edge labelling accordingly.
#'
#' @param jp An object of class `jplace`.
#' @param outgroup A character vector listing outgroup sequences or their unique identifier(s),
#'   e.g., a character string "outgroup".
#' @returns A modified `jplace` object, i.e., its `tree` is rooted and `placements$edge_num` modified.
#' @export
root_jplace <- function(jp, outgroup) {
	if (attr(jp, "classified") == TRUE) {
		stop("classified 'jplace' objects are not assumed to be rooted")
	}
	tree <- jp$tree
	match_tips <- match(seq(ape::Ntip(tree)), tree$edge[,2])
	match_nodes <- match(seq(ape::Nnode(tree)) + ape::Ntip(tree), tree$edge[,2])
	tree$tip.label <- paste0(tree$tip.label, "{", tree$edge.label[match_tips], "}")
	tree$node.label <- paste0("{", tree$edge.label[match_nodes], "}")
	outgroup_tips <- tree$tip.label[sapply(outgroup, grep, x=tree$tip.label)]
	tree <- ape::root(tree, outgroup=outgroup_tips, resolve.root=TRUE)	
	tree$edge.label <- character(ape::Nedge(tree))
	match_tips <- match(seq(ape::Ntip(tree)), tree$edge[,2])
	match_nodes <- match(seq(ape::Nnode(tree)) + ape::Ntip(tree), tree$edge[,2])
	tree$edge.label[match_tips] <- tree$tip.label
	tree$edge.label[na.omit(match_nodes)] <- tree$node.label[!is.na(match_nodes)]
	tree$edge.label <- gsub("^.*\\{|}", "", tree$edge.label)
	tree$tip.label <- gsub("\\{[[:digit:]]+}$", "", tree$tip.label)
	tree$node.label <- NULL
	root_edge <- which(tree$edge[,1] == ape::Ntip(tree) + 1)
	root_edge_label <- na.omit(suppressWarnings(as.numeric(tree$edge.label[root_edge])))
	tree$edge.label[root_edge] <- NA
	tree$outgroup <- outgroup
	jp$tree <- tree
	jp$placements$edge_num[jp$placements$edge_num %in% root_edge_label] <- NA
	attr(jp, "classified") <- FALSE
	class(jp) <- "jplace"
	return(jp)
}
