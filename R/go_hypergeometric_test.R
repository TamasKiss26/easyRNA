#' Run Hypergeometr Test On a Set of Target Genes
#'
#' @param target_set vector of target ensembl ids
#' @param gene_universe vector of universe ensembl ids
#' @param output_name file name of the ouptut
#' @param output_path output folder location
#'
#' @return nothing, only write results to the disk
#'
#' @examples
#'

hgTest <- function(target_set, gene_universe, output_name, output_path = './annotation_data/GO/'){

  # get ENTREZ gene ids to the universe
  universe <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = universe,
    columns = c("ENSEMBL","ENTREZID", "SYMBOL"),
    keytype = "ENSEMBL"
  )

  # get ENTREZ gene ids to the target set
  target_set <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = target_set,
    columns = c("ENSEMBL","ENTREZID"),
    keytype = "ENSEMBL"
  )

  # run hypergeometric test with the three go domains
  res <- purrr::map(
    c('CC', 'MF', 'BP'),
    function(ontologyT){

      # GOHyperGParams object
      paramsGO <- new(
        "GOHyperGParams",
        geneIds = target_set$ENTREZID,
        universeGeneIds = universe$ENTREZID,
        annotation = 'org.Mm.eg.db',
        ontology = ontologyT,
        pvalueCutoff = 0.05,
        conditional = F,
        testDirection = 'over'
      )

      # run hypergeometric test
      Over.GO <- GOstats::hyperGTest(paramsGO)

      # summary test result
      Over.GO_summary <- Over.GO%>%
        summary()%>%
        dplyr::mutate(GOIDtype = ontologyT)

      # rename collumns of the summary
      colnames(Over.GO_summary) <- c("GOID", "Pvalue", "OddsRatio", "ExpCount", "Count", "Size", "Term", "GOIDtype")

      # get gene symbols to the enriched go terms
      genes2go <- purrr::map(
        geneIdsByCategory(Over.GO),
        function(g){
          sb <- target_set%>%
            dplyr::filter(ENTREZID %in% g)%>%
            dplyr::pull('SYMBOL')
          sb <- paste(sb, collapse = ', ')
          return(sb)
        }
      )

      # attach gene names to the summary column
      genes2go <- dplyr::tibble(
        'GOID'= names(genes2go),
        'Genes' = unlist(genes2go)
      )

      # attach geses to the summary table
      Over.GO_summary <- dplyr::left_join(Over.GO_summary, genes2go, by = 'GOID')

      return(Over.GO_summary)
    }
  )

  #output
  names(res) <- c('CC', 'MF', 'BP')
  write.xlsx(
    res,
    file = paste(output_path, output_name, '.xlsx', sep = ''),
    row.names = F
  )
}
