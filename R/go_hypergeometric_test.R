#' Run Hypergeometr Test On a Set of Target Genes
#'
#' @param target_set vector of target ensembl ids
#' @param gene_universe vector of universe ensembl ids
#' @param outputN file name of the ouptut
#'
#' @return
#'
#' @examples
#'

hgTest <- function(target_set, gene_universe, outputN){

  #get ENTREZ gene ids to the target set
  target_set <- select(
    org.Mm.eg.db,
    keys = target_set,
    columns=c("ENSEMBL","ENTREZID"),
    keytype="ENSEMBL"
  )
  target_set <- target_set$'ENTREZID'

  res <- map(
    c('CC', 'MF', 'BP'),
    function(ontologyT){

      # GOHyperGParams object
      paramsGO <- new(
        "GOHyperGParams",
        geneIds = target_set,
        universeGeneIds = gene_universe,
        annotation = 'org.Mm.eg.db',
        ontology = ontologyT,
        pvalueCutoff = 0.05,
        conditional = F,
        testDirection = 'over'
      )

      # run hypergeometric test
      Over.GO <- hyperGTest(paramsGO)

      Over.GO_summary <- Over.GO%>%
        summary()%>%
        mutate(GOIDtype = ontologyT)

      # rename collumns of the result table
      colnames(Over.GO_summary) <- c("GOID", "Pvalue", "OddsRatio", "ExpCount", "Count", "Size", "Term", "GOIDtype")

      # genes to all enriched go terms
      genes2go <- map(
        geneIdsByCategory(Over.GO),
        function(g){
          id2name <- select(
            org.Mm.eg.db,
            keys = g,
            columns=c("SYMBOL","ENTREZID"),
            keytype="ENTREZID"
          )

          id2name <- id2name%>%
            pull(SYMBOL)%>%
            paste(collapse = ', ')

          return(id2name)
        }
      )

      genes2go <- dplyr::tibble(
        'GOID'= names(genes2go),
        'Genes' = unlist(genes2go)
      )

      # attach geses to the summary table
      Over.GO_summary <- left_join(Over.GO_summary, genes2go, by = 'GOID')

      return(Over.GO_summary)
    }
  )

  #output
  names(res) <- c('CC', 'MF', 'BP')
  write.xlsx(
    res,
    file = paste('./annotation_data/GO/', outputN, '.xlsx', sep = ''),
    sheetName = 'CC',
    row.names = F
  )
}
