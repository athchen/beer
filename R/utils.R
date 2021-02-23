.checkCounts <- function(object){
    is_int <- vapply(PhIPData::counts(object), is.integer, logical(1))
    if(!all(is_int)){
        stop("`counts` entries must be integers.")
    }
}

get_phat <- function(object){

    n <- PhIPData::librarySize(object, withDimnames = FALSE)
    n_matrix <- matrix(rep(n, nrow(object)), nrow = nrow(object), byrow = TRUE)
    counts_matrix <- as.matrix(PhIPData::counts(object))

    counts_matrix/n_matrix
}

