.checkCounts <- function(object){
    is_int <- vapply(PhIPData::counts(object), is.integer, logical(1))
    if(!all(is_int)){
        stop("`counts` entries must be integers.")
    }
}
