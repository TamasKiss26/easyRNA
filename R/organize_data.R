#' Get Differentially Expressed Genes by Custum Significance Threshold
#'
#' \code{getDE} read the DESeq output file and creat or rewrite the DE column based on the manualy set fold change and p-value thesholds.
#'
#' @param fileN file name in the data folder
#' @param fct fold change threshold
#' @param pvt p-value threshold
#'
#' @return
#'
#' @examples
#'

getDE <- function(fileN, fct, pvt){
  # read data file
  df <- read_excel(
    paste('./data/raw/', fileN, '.xlsx', sep = '')
  )

  # create new variables: DE, FCsize, GRank
  df <- df%>%
    mutate(DE = if_else(
      (log2FoldChange < -log(fct,2) | log2FoldChange > log(fct,2)) & padj < pvt,
      'DE',
      'notDE',
      missing = 'notDE'
    ))%>%
    mutate(FCsize = if_else(
      log2FoldChange < -log(fct,2) | log2FoldChange > log(fct,2),
      'FCbig',
      'FCsmall'
    ))%>%
    mutate(GRank = dense_rank(-log2FoldChange))

  # drop lines without gene name (these are 'retired' ensembl genes)
  df <- df%>%
    filter(!is.na(Gene))

  # write result
  write.xlsx(
    df,
    paste('./data/raw/', fileN, '.xlsx', sep = '')
  )

  # quality check
  return(dim(df)[1])
}
