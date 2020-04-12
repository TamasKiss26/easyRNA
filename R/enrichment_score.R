#' Plot Running Enrichment Score to Any Gene Set
#'
#' \code{getRES}
#'
#' @param target_set vactor of ensembl gene ids
#' @param universe DEseq2 output file with gene rank column
#'
#' @return Running Enrichment Score (RES) plot
#'
#' @examples
#'

getRES <- function(target_set, universe){

  if(!is.data.frame(universe)){
    stop('')
  }

  #if(!is.character(target_set)){
  #  stop('')
  #}

  #get target location code
  target_code <- as.numeric(universe$Ensembl_ID %in% target_set)
  target_code.inv <- 1- target_code

  #calculate running enrichment score
  RES <- cumsum(target_code - target_code.inv)

  #scale running enrichment score
  scale <- sum(target_code) / sum(target_code.inv)
  RES.simple <- cumsum(target_code - target_code.inv * scale)

  #combine score with the abailable data
  univ_RES <- data.frame(univ = ranked_universe, RES = RES.simple)

  #plot
  gg <- univ_RES%>%
    filter(ranked_universe %in% target_set)%>%
    ggplot2::ggplot(aes(x = GRank, y = 0, color = log2FoldChange))+
    ggplot2::geom_point()+
    ggplot2::scale_color_gradient2(
      name = expression(log[2]~FC),
      midpoint=0,
      low="darkgreen",
      mid="white",
      high="darkred",
      space ="Lab" )+
    ggplot2::ggtitle(ttl)+
    ggplot2::xlab('Genes ranked')+
    ggplot2::ylab('Enrichment score (ES)')+
    ggplot2::theme_bw()
  gg <- gg+geom_line(aes(y = RES), color = 'black')

  return(gg)
}
