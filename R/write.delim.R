#' Writing tab-delimited file.
#' 
#' @description
#' A wrapper for [utils::write.table] writing a tab-delimited text file. The counterpart of [utils::read.delim].
#'
#' @param x A matrix or data frame to be written.
#' @param file Name of the output file.
#' @param row.names Logical, whether to use row names (default is FALSE).
#' @param col.names Logical, whether to use column names (default is TRUE).
#' @param quote Logical or numeric, whether to enclose field contents into quotes (or in which columns).
#' @param sep the field separator string.
#' @export

write.delim <- function(x, file, row.names=FALSE, col.names=TRUE, quote=FALSE) {
	utils::write.table(x, file, row.names=row.names, col.names=col.names, quote=quote, sep="\t")
}
