#' @title exactlyone
#' @description Compute the probability of exactly one of these events occurring using R
#' @param p a vector contains probabilities pi
#' @return the probability of it
#' @examples
#' \dontrun{
#' p <- c(0.1,0.2,0.3,0.4)
#' exactlyone(p=p)
#' }
#' @export
exactlyone <- function(p) {
  notp <- 1 - p
  tot <- 0.0
  for (i in 1:length(p))
    tot <- tot + p[i] * prod(notp[-i])
  return(tot)
}