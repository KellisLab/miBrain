
#' Pivot a dataframe using Matrix::sparseMatrix
#'
#' @param df Dataframe
#' @param row Row of dataframe
#' @param col Col of dataframe
#' @param val Value to be added
#' @param dense logical to determine whether to densify
#' @return pivoted matrix
#' @export
pivot <- function(df, row, col, val, dense=TRUE, use.last.ij=FALSE) {
    row = as.factor(df[[row]])
    col = as.factor(df[[col]])
    M = Matrix::sparseMatrix(i=as.integer(row),
                             j=as.integer(col),
                             x=df[[val]],
                             use.last.ij=use.last.ij,
                             dimnames=list(levels(row), levels(col)),
                             dims=c(length(levels(row)), length(levels(col))))
    if (dense) {
        return(as.matrix(M))
    } else {
        return(M)
    }
}
