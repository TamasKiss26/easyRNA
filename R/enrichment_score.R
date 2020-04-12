#' Calculate Running Enrichment Score to Any Gene Set
#'
#' \code{getRES}
#'
#' @param target_set vactor of ensembl gene ids
#' @param ranked_universe DEseq2 output file with gene rank column
#'
#' @return Running Enrichment Score (RES) plot
#'
#' @examples
#' getRES(c(1,2), c(1,2,3,4,5,6,7))
#' getRES(c('a','b'), c('a','b','c','d','e','f','g'))

getRES <- function(target_set, ranked_universe){

  if(length(target_set) >= length(ranked_universe)){
    stop('The target set has to be shorter than the univerese!')
  }

  if(sum(target_set %in% ranked_universe) < length(target_set)){
    stop('All element of the target set has to be in the universe!')
  }

  #get target location code
  target_code <- as.numeric(ranked_universe %in% target_set)
  target_code.inv <- 1- target_code

  #calculate running enrichment score
  RES <- cumsum(target_code - target_code.inv)

  #scale running enrichment score
  scale <- sum(target_code) / sum(target_code.inv)
  RES.simple <- cumsum(target_code - target_code.inv * scale)

  #combine score with the abailable data
  univ_RES <- data.frame(univ = ranked_universe, RES = RES.simple)

  return(univ_RES)
}


#' Plot Result of RES
#'
#' \code{plotRES}
#'
#' @param
#'
#' @return
#'
#' @examples
#'

plotRES <- function(){

}
