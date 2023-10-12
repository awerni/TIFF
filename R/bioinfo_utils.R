# Others ----------------------------------------------------------------------

#' Filter Patient Annotation Data by Prefilter
#'
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#' @param pa patient annotation data (see \code{\link{getPatientAnnotation}})
#'
#' @return filtered patient annotation data by applying tissue panel and
#' selected entries from TiffPrefilter object.
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter("PDX", selected = "adrenocortical carcinoma")
#' pa <- getPatientAnnotation()
#' getPatientAnnotationFilter(prefilter, pa)
#' 
getPatientAnnotationFilter <- function(prefilter, pa) {
  if (prefilter$db_col == "organ") pa <- pa %>% filter(tumortype == 'normal' & is.na(tumortype_adjacent))
  
  discardCols <- setdiff(c("tumortype", "tumortype_adjacent", "organ"), prefilter$db_col)
  pa %>%
    select(-one_of(discardCols)) %>%
    filter(pa[, prefilter$db_col] %in% prefilter$tissue, tissuepanel == prefilter$tissue_panel) %>%
    select(-tissuepanel) %>%
    droplevels()
}

#' Get tissue name translation and add the tissue annotation
#' 
#' @param ti data.frame with tissuename column
#' @param ti_anno data.frame, tissue annotation
#' 
#' @return list with fields: df, missing
#' 
#' @details This is an utility function for merging a data.frame (\code{ti})
#' with tissue annotation data.frame (\code{ti_anno}). It returns a list
#' with merged data.frame and and character vector of tissue names that are
#' not available in annotation data.
#' 
#' @export
#' 
#' @examples 
#' 
#' setDbOptions(getSettings())
#' 
#' ti_anno <- getTissueAnnotation()
#' ti <- exampleClassAssigment()
#' ti$class2 <- c(ti$class2, "MISSING") # add non-existing tissuename
#' ti <- XIFF::stackClasses(ti)
#' 
#' translated <- getTissuenameTranslation(ti, ti_anno)
#' translated$missing
#' 
getTissuenameTranslation <- function(ti, ti_anno) {
  ti_col <- findItemColumn(ti, ti_anno)
  if (is.null(ti_col)) return()
  
  if (ti_col != "tissuename"){
    ti <- ti %>% 
      rename(tissuename = !!ti_col) %>%
      relocate(tissuename, .before = 1)
  }
  
  df <- ti %>% inner_join(ti_anno, by = "tissuename")
  missing <- setdiff(ti$tissuename, df$tissuename)
  
  return(list(df = df, missing = missing))
}

findItemColumn <- function(df, anno, colname = "tissuename"){
  if (!is.data.frame(df) || nrow(df) == 0) return()
  if (colname %in% names(df)) return(colname)
  
  df <- select_if(df, is.character)
  if (ncol(df) == 0) return()
  
  nCommonValues <- vapply(
    X = df, 
    FUN = n_common, 
    FUN.VALUE = integer(1),
    y = anno[[colname]],
    USE.NAMES = FALSE
  )
  
  if (any(nCommonValues > 0)) {
    names(df)[which.max(nCommonValues)]
  } # else NULL
}
