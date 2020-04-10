#' Get Differentially Expressed Genes by Custum Significance Threshold
#'
#' \code{getDE} read the DESeq output file and creat or rewrite the DE column based on the manualy set fold change and p-value thesholds.
#'
#' @param file_name file name in the data folder
#' @param fold_change fold change threshold
#' @param p_value p-value threshold
#' @param path path to the DESeq files
#'
#' @return vector of differentially expressed genes
#'
#' @examples
#'

getDE <- function(file_name, fold_change, p_value, path = './data/raw/'){

  # check input quality
  if(!is.character(file_name)){
    stop('File name has to be character verctor!')
  }

  if(!is.numeric(fold_change)){
    stop('Fold change has to be numeric!')
  }

  if(fold_change <= 0){
    stop('Fold change has to be positive number!')
  }

  if(!is.numeric(p_value)){
    stop('p-value has to be numeric!')
  }

  if(p_value > 1| p_value < 0){
    stop('p-value has to be a number between 0 and 1')
  }

  # read data file
  df <- readxl::read_excel(
    paste(path, file_name, '.xlsx', sep = '')
  )

  # check if the data frame has the right columns
  if(
    sum(c('log2FoldChange', 'padj')%in% colnames(df)) != 2
  ){
    stop('Right columns are missing! \nThe function needs: \'log2FoldChange\' and \'padj\'')
  }

  # create new variables: DE, FCsize, GRank
  df <- df%>%
    dplyr::mutate(DE = if_else(
      (log2FoldChange < -log(fold_change,2) | log2FoldChange > log(fold_change,2)) & padj < p_value,
      'DE',
      'notDE',
      missing = 'notDE'
    ))%>%
    dplyr::mutate(FCsize = if_else(
      log2FoldChange < -log(fold_change,2) | log2FoldChange > log(fold_change,2),
      'FCbig',
      'FCsmall'
    ))%>%
    dplyr::mutate(GRank = dense_rank(-log2FoldChange))

  # drop lines without gene name (these are 'retired' ensembl genes)
  df <- df%>%
    dplyr::filter(
      !is.na(Gene)
      )

  # write result
  openxlsx::write.xlsx(
    df,
    paste(path, file_name, '.xlsx', sep = '')
  )

  # quality check
  genes <- df%>%
    dplyr::filet(DE == 'DE')%>%
    pull('Ensemb_ID')%>%
    unique()

  # return gene vector
  return(genes)
}
