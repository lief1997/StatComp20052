#' @title findols
#' @description identify the most deviant observation for each row of x using R
#' @param x a matrix, each row of which is observations of corresponding objects
#' @return a vector contains the indices of the most deviant observation in each row
#' @examples
#' \dontrun{
#' rs <- matrix(round(rnorm(60)*10,digits = 2), nrow = 4, byrow = TRUE)
#' findols(rs)
#' }
#' @importFrom stats median
#' @export
findols <- function(x) {
  findol <- function(xrow) {
    mdn <- median(xrow)
    devs <- abs(xrow - mdn)
    return(which.max(devs))
    }
  return(apply(x,1,findol)) 
}
