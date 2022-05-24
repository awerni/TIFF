# Annotations -----------------------------------------------------------------

#' Get Tissue Annotation
#'
#' @return data.frame with tissue annotation.
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' head(getTissueAnnotation())
#' 
getTissueAnnotation <- function() {
  panels <- c("GTEx normals", "TCGA normals", "TCGA tumors", "PDX Models")
  
  sql <- paste0(
    "SELECT t.tissuename, tumortype, tumortype_adjacent, organ, tissue_subtype, ",
    "CASE WHEN dnasequenced THEN 'YES' ELSE 'NO' END AS dnasequenced, ",
    "tumorpurity, microsatellite_stability_class, immune_environment, ",
    "gi_mol_subgroup, icluster, consmolsubtype, til_pattern, autolysis_score, ",
    "rna_integrity_number, minutes_ischemia, vendorname, tissuepanel, ",
    "histology_type, histology_subtype, lossofy ",
    "FROM tissue.tissue t ",
    "JOIN tissue.tissueassignment ta ON t.tissuename = ta.tissuename ",
    "WHERE ", getSQL_filter("tissuepanel", panels)
  )
  getPostgresql(sql)
}


#' Get Patient Annotation
#'
#' @return data.frame with patient annnotation data
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' head(getPatientAnnotation())
#' 
getPatientAnnotation <- function() {
  panels <- c("GTEx normals", "TCGA normals", "TCGA tumors", "PDX Models")
  sql <- paste0(
    "SELECT t.tissuename, tumortype, tumortype_adjacent, organ, vital_status, ",
    "round(-days_to_birth/365.25) as age, p.gender, height, weight, ",
    "weight/(height/100) as bmi, race, ethnicity, days_to_last_followup, ",
    "days_to_death, ta.tissuepanel ", 
    "FROM tissue.patient p JOIN tissue.tissue t ON t.patientname = p.patientname ",
    "JOIN tissue.tissueassignment ta ON t.tissuename = ta.tissuename WHERE ", 
    getSQL_filter("tissuepanel", panels)
  )
  pa <- getPostgresql(sql)
  if (nrow(pa) == 0) return()
  
  factorCols <- c("tumortype", "tumortype_adjacent", "organ", "ethnicity", "race", "gender")
  
  pa %>% 
    mutate(
      bmi = round(bmi, 1),
      vital_status = as.factor(ifelse(vital_status, "alive", "dead"))
    ) %>%
    mutate_at(factorCols, as.factor)
}

#' Get gene annotation
#' 
#' Get all available gene annotation for the selected species.
#' The function returns as a data.frame with fields:
#' ensg, symbol, name, location
#' 
#' @param species string containing the species name
#' 
#' @return list of data.frames
#' @export
#' 
#' @examples 
#' 
#' setDbOptions(getSettings())
#' head(getGeneAnno("human"))
#' 
getGeneAnno <- function(species) {
  sql <- paste0(
    "SELECT ensg, symbol, regexp_replace(name, ' \\[.*\\]', '') AS name, ",
    "chromosome || ':' || seqregionstart || '-' || seqregionend AS location ",
    "FROM gene WHERE species = '", species, "'"
  )
  getPostgresql(sql)
}

# Mutations -------------------------------------------------------------------
#' Get mutation information
#' 
#' Returns mutation information related to the specified gene ID
#' 
#' @param ensg character string, gene ID
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#' @return data.frame with fields: tissuename, tumortype, dnamutation, aamutation,
#' aamutated, histology_type, histology_subtype
#' @export
#' 
#' @examples 
#' 
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter("tumor", c("thymoma", "mesothelioma"))
#' ensg <- "ENSG00000133703"
#' getInfoMutation(ensg, prefilter)
#' 
getInfoMutation <- function(ensg, prefilter){
  sql <- preparePrefilterSql(
    SELECT = c(
      "ps.tissuename", "dnamutation", "aamutation", 
      "coarse(aamutation) as aamutated", "histology_type", "histology_subtype"
    ),
    FROM = "tissue.processedsequenceextended ps",
    WHERE = paste0(
      prepareConditionSql(ensg = ensg),
      "AND aamutation IS NOT NULL"
    ),
    JOIN = paste(
      "JOIN tissue.tissue USING(tissuename)",
      "JOIN tissue.tissueassignment USING(tissuename)"
    ),
    prefilter = prefilter
  )
  getPostgresql(sql)
}


#' Get Sequenced Tissues
#'
#' @param tissuenames vector of tissuenames or classAssigment object
#'
#' @return
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' head(getSequencedTissues(exampleClassAssigment()))
#' 
getSequencedTissues <- function(tissuenames){

  tissuenames <- XIFF::unlistClassAssignment(tissuenames)
  sql <- paste0(
    "SELECT tissuename FROM tissue.tissue WHERE dnasequenced AND ",
    getSQL_filter("tissuename", tissuenames)
  )
  getPostgresql(sql)$tissuename
}

#' Get Analysis Mutation Data
#'
#' @param tissuenames vector of tissues or classAssigment object
#'
#' @return data.frame with columns: ensg, tissuename
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' df <- getAnalysisDataMutation(ca)
#' head(df)
#'
getAnalysisDataMutation <- function(tissuenames){
  
  tissuenames <- XIFF::unlistClassAssignment(tissuenames)
  
  sql <- paste0(
    "SELECT ensg, tissuename from tissue.processedsequence ps ",
    "JOIN transcript t on (t.enst = ps.enst) ",
    "WHERE iscanonical AND coarse(aamutation) = 'mut' AND ",
    getSQL_filter("tissuename", tissuenames)
  )
  getPostgresql(sql)
}

#' Get tissue Mutation data by gene ID
#' 
#' Returns Mutation data for the specified tissue and gene.
#' 
#' @param ensg character string, gene ID. Use NULL to get data for all genes
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param addRownames if TRUE then adds tissuename as rownames
#' @return data.frame with fields: tissuename, ensg, aamutation, aamutated, 
#' rnazygosity, dnazygosity
#' @export
#' 
#' @examples 
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' ensg <- "ENSG00000133703"
#' 
#' getTissueDataMutationById(ensg, ca)
#' 
getTissueDataMutationById <- function(ensg, tissueClasses, addRownames = TRUE){
  sql <- paste0(
    "SELECT tissuename, ensg, aamutation, coarse(aamutation) AS aamutated ",
    "FROM tissue.processedsequenceextended"
  )
  
  getTissueDataCommonById(
    sql = sql,
    tissues = tissueClasses,
    conditionSql = prepareConditionSql(ensg = ensg),
    addRownames = addRownames
  )
}

# Mutational burden -----------------------------------------------------------
#' Get waterfall Mutational Burden data
#' 
#' SQL wrapper for Mutational Burden data retrieval. The `tissuename` column 
#' in the result table is a factor with levels sorted according to the 
#' `mutational_fraction` value (desc order).
#' 
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#' @return data.frame with columns: tissuename, tumortype, mutational_fraction
#' @export
#' 
#' @examples 
#' 
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter("tumor", c("thymoma", "mesothelioma"))
#' getWaterfallDataMutationalBurden(prefilter)
#' 
getWaterfallDataMutationalBurden <- function(prefilter){
  # tumortype column is available in tissue.mutationalburden and tissue.tissue
  prefilter[["db_col"]] <- paste0("ti.", prefilter[["db_col"]])
  prefilter[["sql"]] <- gsub("tumortype", "ti.tumortype", prefilter[["sql"]])
  
  sql <- preparePrefilterSql(
    SELECT = c(
      "mb.tissuename", 
      "mutational_fraction * 100 AS mutational_fraction_in_percent"
    ),
    FROM = "tissue.mutationalburden mb",
    WHERE = "mutational_fraction IS NOT NULL",
    JOIN = paste(
      "JOIN tissue.tissue ti USING(tissuename)",
      "JOIN tissue.tissueassignment USING(tissuename)"
    ),
    prefilter = prefilter
  )
  
  getPostgresql(sql) %>% reorderByScore(valueCol = "mutational_fraction_in_percent")
}

# Copy numbers ----------------------------------------------------------------


#' Get Waterfall Data Copy Number
#'
#' @param ensg gene ensg
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#'
#' @return data.frame with tissuename, relativecopynumber and tumortype
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter()
#' data <- getWaterfallDataCopyNumber("ENSG00000133703", prefilter)
#' head(data)
#' 
getWaterfallDataCopyNumber <- function(ensg, prefilter){
  sql <- preparePrefilterSql(
    SELECT = c(
      "pcv.tissuename", 
      "2*2^log2relativecopynumber as relativecopynumber"
    ),
    FROM = "tissue.processedcopynumber pcv",
    WHERE = paste(
      prepareConditionSql(ensg = ensg), 
      "AND log2relativecopynumber IS NOT NULL"
    ),
    JOIN = paste(
      "JOIN tissue.tissue USING(tissuename)",
      "JOIN tissue.tissueassignment USING(tissuename)"
    ),
    prefilter = prefilter
  )
  getPostgresql(sql) %>% reorderByScore(valueCol = "relativecopynumber")
}

#' Get Tissue Copy Number Data
#'
#' @param tissuenames vector of tissuenames or classAssigment
#'
#' @return data.frame with columns: tissuenames, ensg and log2relativecopynumber
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' head(getTissueDataCopyNumber(ca))
#'
getTissueDataCopyNumber <- function(tissuenames){
  ti <- getSQL_filter("tissuename", tissuenames)
  sql2 <- paste0(
    "SELECT tissuename, unnest(log2relativecopynumber) AS log2relativecopynumber ",
    "FROM tissue.processedcopynumber_array WHERE ", ti
  )
  
  expr_long <- getPostgresql(sql2)
  if (nrow(expr_long) == 0) {
    return(
      data.frame(
        tissuename = character(),
        ensg = character(),
        log2relativecopynumber = numeric()
      )
    )
  }
  
  sql1 <- "SELECT ensg FROM tissue.cnaltered_ensg WHERE species = 'human' ORDER BY ensg"
  my_ensg <- getPostgresql(sql1)$ensg
  
  expr_long[, "ensg"] <- my_ensg
  expr_long %>% select(tissuename, ensg, everything())
}

#' Get tissue Copy Number data by gene ID
#' 
#' Returns Copy Number data for the specified tissue and gene.
#' 
#' @param ensg character string, gene ID. Use NULL to get data for all genes
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param addRownames if TRUE (default), the tissue names are added as rownames
#' to the output table
#' @return data.frame with fields: tissuename, ensg, relativecopynumber
#' @export
#' 
#' @examples 
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' ensg <- "ENSG00000133703"
#' getTissueDataCopyNumberById(ensg, ca)
#' 
getTissueDataCopyNumberById <- function(ensg, tissueClasses, addRownames = TRUE){
  sql <- paste0(
    "SELECT tissuename, ensg, 2*2^log2relativecopynumber AS relativecopynumber ",
    "FROM tissue.processedcopynumber"
  )
  
  getTissueDataCommonById(
    sql = sql,
    tissues = tissueClasses,
    conditionSql = prepareConditionSql(ensg = ensg),
    addRownames = addRownames
  )
}

# Expression ------------------------------------------------------------------
#' Get waterfall Gene Expression data
#' 
#' SQL wrapper for Gene Expression data retrieval. The `tissuename` column 
#' in the result table is a factor with levels sorted according to the 
#' `log2tpm` value (desc order).
#' 
#' @param ensg character string, gene ID
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#' @return data.frame with columns: tissuename, tumortype, log2tpm
#' @export
#'
#' @examples 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' ensg <- "ENSG00000133703"
#' prefilter <- makePrefilter("tumor", c("thymoma", "mesothelioma"))
#' getWaterfallDataGeneExpression(ensg, prefilter)
#' 
getWaterfallDataGeneExpression <- function(ensg, prefilter){
  sql <- preparePrefilterSql(
    SELECT = c("tissuename", "2^log2tpm AS tpm"),
    FROM = "tissue.processedrnaseqview",
    WHERE = paste(
      prepareConditionSql(ensg = ensg), 
      "AND log2tpm IS NOT NULL"
    ),
    JOIN = paste(
      "JOIN tissue.tissue USING(tissuename)",
      "JOIN tissue.tissueassignment USING(tissuename)"
    ),
    prefilter = prefilter
  )
  getPostgresql(sql) %>% reorderByScore(valueCol = "tpm")
}

#' Get Tissue Gene Expression Data
#'
#' @param tissuenames character vector with tissuenames 
#' @param unit unit of the result, possible values are `log2tpm` (default)
#' or `counts`.
#'
#' @return data.frame with columns: tissuenames, ensg, log2tpm
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' head(getTissueDataGeneExpression(ca))
#' head(getTissueDataGeneExpression(ca, "counts"))
#'
getTissueDataGeneExpression <- function(tissuenames, unit = "log2tpm") {
  
  tissuenames <- XIFF::unlistClassAssignment(tissuenames)
  
  unit_query <- paste("unnest(", unit, ") AS ", unit, sep = "", collapse = ", ")
  ti <- getSQL_filter("tissuename", tissuenames)
  
  sql1 <- "SELECT ensg FROM tissue.expressed_ensg WHERE species = 'human' ORDER BY ensg"
  my_ensg <- getPostgresql(sql1)$ensg
  
  sql2 <- paste0(
    "SELECT tissuename, ", unit_query, 
    " FROM tissue.processedrnaseq_array WHERE ", ti
  )
  expr_long <- getPostgresql(sql2)
  expr_long[, "ensg"] <- my_ensg
  expr_long %>% 
    select(tissuename, ensg, everything())
}

#' Get tissue Gene Expression data by gene ID
#' 
#' Returns Gene Expression data for the specified tissue and gene.
#' 
#' @param ensg character string, gene ID. Use NULL to get data for all genes
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param addRownames if TRUE (default), the tissue names are added as rownames
#' to the output table
#' @param unit unit of the result, possible values are `tpm` (default)
#' or `log2tpm`.
#' 
#' @return data.frame with fields: tissuename, ensg, tpm
#' @export
#' 
#' @examples 
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' ensg <- "ENSG00000133703"
#' getTissueDataGeneExpressionById(ensg, ca)
#' 
getTissueDataGeneExpressionById <- function(ensg,
                                            tissueClasses,
                                            addRownames = TRUE,
                                            unit = c("tpm", "log2tpm", "counts")){
  
  unit <- match.arg(unit)
  
  if(unit == "tpm") {
    unit <- "2^log2tpm AS tpm"
  }
  
  sql <- paste0(
    "SELECT tissuename, ensg, ", unit, " ",
    "FROM tissue.processedrnaseqview"
  )
  
  getTissueDataCommonById(
    sql = sql,
    tissues = tissueClasses,
    conditionSql = prepareConditionSql(ensg = ensg),
    addRownames = addRownames
  )
}

# Proteins --------------------------------------------------------------------
#' Get waterfall Protein RPPA data
#' 
#' SQL wrapper for Protein RPPA data retrieval. The `tissuename` column 
#' in the result table is a factor with levels sorted according to the 
#' `x_score value` (desc order).
#' 
#' @param antibody character string, antibody ID
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#' @return data.frame with columns: tissuename, tumortype, x_score
#' @export
#' 
#' @examples 
#' 
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter("tumor", c("thymoma", "mesothelioma"))
#' getWaterfallDataRppa("CDK1", prefilter)
#' 
getWaterfallDataRppa <- function(antibody, prefilter){
  sql <- preparePrefilterSql(
    SELECT = c("tissuename", "score"),
    FROM = "tissue.processedproteinexpression",
    WHERE = paste(
      prepareConditionSql(antibody = antibody), 
      "AND score IS NOT NULL"
    ),
    JOIN = paste(
      "JOIN tissue.tissue USING(tissuename)",
      "JOIN tissue.tissueassignment USING(tissuename)"
    ),
    prefilter = prefilter
  )
  getPostgresql(sql) %>% reorderByScore()
}

#' Get tissue RPPA data
#' 
#' SQL wrapper for RPPA data retrieval.
#' 
#' @param tissuenames character vector of tissue names or classAssigment object
#' @return data.frame with columns: tissuename, antibody, score
#' @export
#' 
#' @examples 
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' getTissueDataRppa(ca)
#' 
getTissueDataRppa <- function(tissuenames){
  getTissueDataRppaById(
    antibody = NULL,
    tissueClasses = tissuenames,
    addRownames = FALSE
  )
}

#' Get tissue RPPA data by antibody ID
#' 
#' Returns RPPA data for the specified tissue and antibody
#' 
#' @param antibody character string, antibody ID. Use NULL to get data for all 
#' antibodies
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param addRownames if TRUE (default), the tissue names are added as rownames
#' to the output table
#' @return data.frame with fields: tissuename, antibody, score
#' @export
#' 
#' @examples 
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' ensg <- "ENSG00000133703"
#' getTissueDataRppaById(NULL, ca)
#' 
getTissueDataRppaById <- function(antibody, tissueClasses, addRownames = TRUE) {
  sql <- paste0(
    "SELECT tissuename, antibody, score ",
    "FROM tissue.processedproteinexpression"
  )
  
  getTissueDataCommonById(
    sql = sql,
    tissues = tissueClasses,
    conditionSql = prepareConditionSql(antibody = antibody),
    addRownames = addRownames
  )
}

# Immune cells ----------------------------------------------------------------

#' Get Waterfall Data Immune Cells
#'
#' @param cellType string with cell type
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#'
#' @return data.frame with columns 'tissuename', 'score' and tumortype
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' df <- getWaterfallDataImmuneCells("Activated dendritic cells", makePrefilter())
#' head(df)
#' 
getWaterfallDataImmuneCells <- function(cellType, prefilter){
  sql <- preparePrefilterSql(
    SELECT = c("tissuename", "score"),
    FROM = "tissue.immunecelldeconvolution",
    WHERE = paste(
      prepareConditionSql(celltype = cellType), 
      "AND score IS NOT NULL"
    ),
    JOIN = paste(
      "JOIN tissue.tissue USING(tissuename)",
      "JOIN tissue.tissueassignment USING(tissuename)"
    ),
    prefilter = prefilter
  )
  getPostgresql(sql) %>% reorderByScore()
}


#' Get Tissue Data Immune Cells
#'
#' @param tissuenames vector of tissuenames or classAssigment object
#'
#' @return data.frame with columns tissuenames, celltype and score.
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' head(getTissueDataImmuneCells(exampleClassAssigment()))
#' 
getTissueDataImmuneCells <- function(tissuenames){
  getTissueDataImmuneCellsById(
    cellType = NULL,
    tissueClasses = tissuenames,
    addRownames = FALSE
  )
}


#' Get Tissue Data Immune Cells By Id
#'
#' @param cellType \code{\link{getAvailableImmuneCellTypes()}}
#' @param tissueClasses vector of tissuenames or classAssigment object
#' @param addRownames if TRUE (default) tissuenames are added 
#' as rownames to the output table.
#'
#' @return data.frame with columns tissuename, celltype and score.
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' head(getTissueDataImmuneCellsById("Neurons", NULL))
#' head(getTissueDataImmuneCellsById("Neurons", exampleClassAssigment()))
#' 
getTissueDataImmuneCellsById <- function(cellType, tissueClasses, addRownames = TRUE){
  sql <- paste0(
    "SELECT tissuename, celltype, score ",
    "FROM tissue.immunecelldeconvolution"
  )
  
  getTissueDataCommonById(
    sql = sql,
    tissues = tissueClasses,
    conditionSql = prepareConditionSql(celltype = cellType),
    addRownames = addRownames
  )
}


#' Get Available Immune Cell Types
#'
#' @return
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' getAvailableImmuneCellTypes()
#' 
getAvailableImmuneCellTypes <- function() {
  sql <- "SELECT DISTINCT celltype FROM tissue.immunecelldeconvolution ORDER BY celltype"
  getPostgresql(sql)$celltype
}

# Metabolics ------------------------------------------------------------------

#' Get Waterfall Data for Metabolics
#'
#' @param pathway metabolic pathway (see \code{\link{getAvailableMetabolicPathways}})
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#'
#' @return data.frame with tissuename, score and tumortype.
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter()
#' head(getWaterfallDataMetabolics("gluconeogenesis", prefilter))
#' 
getWaterfallDataMetabolics <- function(pathway, prefilter){
  sql <- preparePrefilterSql(
    SELECT = c("tissuename", "score"),
    FROM = "tissue.metabolics",
    WHERE = paste(
      prepareConditionSql(metabolic_pathway = pathway), 
      "AND score IS NOT NULL"
    ),
    JOIN = paste(
      "JOIN tissue.tissue USING(tissuename)",
      "JOIN tissue.tissueassignment USING(tissuename)"
    ),
    prefilter = prefilter
  )
  getPostgresql(sql) %>% reorderByScore()
}


#' Get Tissue Data Metabolics
#'
#' @param tissuenames vector of tissuenames or classAssigment object
#'
#' @return data.frame with tissuename, metabolic_pathway and score.
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' head(getTissueDataMetabolics(exampleClassAssigment()))
#' 
getTissueDataMetabolics <- function(tissuenames){
  getTissueDataMetabolicsById(
    pathway = NULL,
    tissueClasses = tissuenames,
    addRownames = FALSE
  )
}


#' Get Tissue Data Metabolics By Id
#'
#' @param pathway metabolic pathway (see \code{\link{getAvailableMetabolicPathways}})
#' @param tissueClasses classAssignment object or vector with tissue names. If NULL
#' returns all tissues
#' @param addRownames if TRUE (default) tissuenames are added 
#' as rownames to the output table.
#'
#' @return
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter()
#' head(getTissueDataMetabolicsById("gluconeogenesis", NULL))
#' 
#' getTissueDataMetabolicsById("gluconeogenesis", c("TCGA-GU-A42R-01", "TCGA-KQ-A41Q-01"))
#' 
getTissueDataMetabolicsById <- function(pathway, tissueClasses, addRownames = TRUE){
  sql <- paste0(
    "SELECT tissuename, metabolic_pathway, score ",
    "FROM tissue.metabolics"
  )
  
  getTissueDataCommonById(
    sql = sql,
    tissues = tissueClasses,
    conditionSql = prepareConditionSql(metabolic_pathway = pathway),
    addRownames = addRownames
  )
}


#' Get Available Metabolic Pathways
#'
#' @return vector with all available metabolic pathways
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' getAvailableMetabolicPathways()
#' 
getAvailableMetabolicPathways <- function() {
  sql <- "SELECT DISTINCT metabolic_pathway FROM tissue.metabolics ORDER BY metabolic_pathway"
  getPostgresql(sql)$metabolic_pathway
}

# Signaling -------------------------------------------------------------------
#' Get Available Signaling Pathways
#'
#' @return vector with all available signaling pathways
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' getAvailableSignalingPathways()
#' 
getAvailableSignalingPathways <- function() {
  sql <- "SELECT * FROM tissue.signaling_pathway WHERE false"
  res <- getPostgresql(sql)
  
  pathways <- setdiff(names(res), "tissuename")
  labels <- gsub("_", " ", toupper(pathways))
  
  structure(
    pathways,
    names = labels
  )
}


#' Get Input Data Signaling
#'
#' @param pathway string containing pathway 
#' (see \code{\link{getAvailableSignalingPathways}} for all available options)
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#'
#' @return data.frame with tissuename, {pathway}, and tumortype.
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' pathways <- getAvailableSignalingPathways()
#' prefilter <- makePrefilter("PDX", selected = "adrenocortical carcinoma")
#' head(getInputDataSignaling(pathways[1], prefilter))
#' 
getInputDataSignaling <- function(pathway, prefilter) {
  if (is.null(pathway)) return()
  
  sql <- preparePrefilterSql(
    SELECT = c(
      "s.tissuename", 
      paste0("CASE WHEN ", pathway, " THEN 'active' ELSE 'inactive' END AS ", pathway)
    ),
    FROM = "tissue.signaling_pathway s",
    JOIN = paste(
      "JOIN tissue.tissue t ON t.tissuename = s.tissuename",
      "JOIN tissue.tissueassignment ta ON t.tissuename = ta.tissuename"
    ),
    prefilter = prefilter
  )
  getPostgresql(sql)
}


#' Get Tissue Data Signaling
#'
#' @param tissuenames vector of tissue names
#'
#' @return data.frame with columns tissuename, pathway, active.
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' getTissueDataSignaling(c("TCGA-06-2557-01", "TCGA-06-2569-01"))
#' 
getTissueDataSignaling <- function(tissuenames){
  getTissueDataCommonById(
    sql = "SELECT * FROM tissue.signaling_pathway",
    tissues = tissuenames,
    addRownames = FALSE
  ) %>% 
    pivot_longer(
      cols = -tissuename, 
      names_to = "pathway", 
      values_to = "active"
    )
}


#' Get Tissue Data Signaling By Id
#'
#' @param pathway signaling pathway (see \code{\link{getAvailableSignalingPathways}})
#' @param tissueClasses vector with selected tissue names, if NULL returns all
#' available ones
#' @param addRownames if TRUE (default), the tissue names are added as rownames
#' to the output table
#'
#' @return
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' head(getTissueDataSignalingById("cell_cycle", NULL))
#' 
getTissueDataSignalingById <- function(pathway, tissueClasses, addRownames = TRUE){
  getTissueDataCommonById(
    sql = paste("SELECT tissuename,", pathway, "as active FROM tissue.signaling_pathway"),
    tissues = tissueClasses,
    addRownames = addRownames
  ) %>%
    tibble::add_column(pathway = pathway, .before = "active")
}

# Clones ----------------------------------------------------------------------

#' Get Waterfall Data Clones
#'
#' @param score string with score name 'number_of_clones' or 'clone_tree_score'
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#'
#' @return data.frame with columns tissuename, score name and tumortype
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter()
#' df <- getWaterfallDataClones("number_of_clones", prefilter)
#' head(df)
#' 
getWaterfallDataClones <- function(score, prefilter) {
  if (!score %in% c("number_of_clones", "clone_tree_score")) return()
  
  sql <- preparePrefilterSql(
    SELECT = c("t.tissuename", score),
    FROM = "tissue.tissue t",
    WHERE = paste(score, "IS NOT NULL"),
    JOIN = "JOIN tissue.tissueassignment ta ON t.tissuename = ta.tissuename",
    prefilter = prefilter
  )
  getPostgresql(sql) %>% reorderByScore(valueCol = score)
}


#' Get Tissue Data Clones
#'
#' @param tissuenames vector of tissuenames or classAssigment object
#' @param scores vector of scores to be extracted
#'
#' @return
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' 
#' getTissueDataClones(NULL)
#' 
getTissueDataClones <- function(tissuenames, scores = "clone_tree_score"){
  scores <- paste(scores, collapse = ", ")
  
  getTissueDataCommonById(
    sql = paste0("SELECT tissuename, ", scores, " FROM tissue.tissue"),
    tissues = tissuenames,
    addRownames = FALSE
  )
}

# Gene signatures -------------------------------------------------------------
#' Get Available Gene Signatures
#'
#' @return data.frame with columns: signature, description, hyperlink
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' getAvailableGeneSignatures()
#'
getAvailableGeneSignatures <- function(){
  sql <- "SELECT signature, description, hyperlink FROM genesignature ORDER BY signature"
  getPostgresql(sql)
}


#' Get Waterfall Data Gene Signatures
#'
#' @param signature gene signatures (see \code{\link{getAvailableGeneSignatures}})
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#'
#' @return data.frame with tissuename, score and tumortype
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter()
#' head(getWaterfallDataGeneSignatures("HRD", prefilter))
#' 
getWaterfallDataGeneSignatures <- function(signature, prefilter){
  sql <- preparePrefilterSql(
    SELECT = c("tissuename", "score"),
    FROM = "tissue.tissue2genesignature",
    WHERE = paste(
      prepareConditionSql(signature = signature), 
      "AND score IS NOT NULL"
    ),
    JOIN = paste(
      "JOIN tissue.tissue USING(tissuename)",
      "JOIN tissue.tissueassignment USING(tissuename)"
    ),
    prefilter = prefilter
  )
  getPostgresql(sql) %>% reorderByScore()
}

#' Get tissue Gene Signature data
#' 
#' SQL wrapper for Gene Signature data retrieval.
#' 
#' @param tissuenames character vector of tissue names
#' @return data.frame with columns: tissuename, signature, score
#' @export
#' 
#' @examples 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' getTissueDataGeneSignature(ca)
#' 
getTissueDataGeneSignature <- function(tissuenames){
  getTissueDataGeneSignatureById(
    signature = NULL,
    tissueClasses = tissuenames,
    addRownames = FALSE
  )
}

#' Get tissue Gene Signature data by signature ID
#' 
#' Returns Gene Signature data for the specified tissue and signature ID
#' 
#' @param signature character string, signature ID. Use NULL to get data for all 
#' signatures
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param addRownames if TRUE (default), the tissue names are added as rownames
#' to the output table
#' @return data.frame with fields: tissuename , signature, score
#' @export
#' 
#' @examples 
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' geneSig <- getTissueDataGeneSignatureById(NULL, ca)
#' 
getTissueDataGeneSignatureById <- function(signature, tissueClasses, addRownames = TRUE){
  sql <- paste0(
    "SELECT tissuename, signature, score ",
    "FROM tissue.tissue2genesignature"
  )
  
  getTissueDataCommonById(
    sql = sql,
    tissues = tissueClasses,
    conditionSql = prepareConditionSql(signature = signature),
    addRownames = addRownames
  )
}

# Hallmark sets ---------------------------------------------------------------
#' Get waterfall Hallmark Set data
#' 
#' SQL wrapper for Hallmark Set data retrieval. The `tissuename` column 
#' in the result table is a factor with levels sorted according to the 
#' `score` value (desc order).
#' 
#' @param geneSet character string, gene set ID
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#' @param score character string, score to retrieve from DB. Available choices:
#' gsva (default), ssgsea
#' @return data.frame with columns: tissuename, tumortype, score
#' @export
#' 
#' @examples
#' 
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter("tumor", c("thymoma", "mesothelioma"))
#' getWaterfallDataHallmarkSet("HALLMARK_P53_PATHWAY", prefilter)
#' 
getWaterfallDataHallmarkSet <- function(geneSet, prefilter, score = c("gsva", "ssgsea")){
  
  score <- match.arg(score)
  
  sql <- preparePrefilterSql(
    SELECT = c("tissuename", paste(score, " AS score")),
    FROM = "tissue.hallmarkscore",
    WHERE = paste(
      prepareConditionSql(gene_set = geneSet), 
      glue::glue("AND {score} IS NOT NULL")
    ),
    JOIN = paste(
      "JOIN tissue.tissue USING(tissuename)",
      "JOIN tissue.tissueassignment USING(tissuename)"
    ),
    prefilter = prefilter
  )
  getPostgresql(sql) %>% reorderByScore()
}

#' Get tissue Gene Set data
#' 
#' SQL wrapper for Gene Set data retrieval.
#' 
#' @param tissuenames character vector of tissue names
#' @param score character string, score to retrieve from DB. Available choices:
#' gsva (default), ssgsea
#' @return data.frame with columns:tissuename, gene_set, score
#' @export
#' 
#' @examples 
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' getTissueDataGeneSet(ca)
#' 
getTissueDataGeneSet <- function(tissuenames, score = c("gsva", "ssgsea")){
  
  score <- match.arg(score)
  
  getTissueDataGeneSetById(
    geneSet = NULL,
    tissueClasses = tissuenames,
    score = score
  )
}

#' Get tissue Gene Set data by set ID
#' 
#' Returns Gene Set data for the specified tissue and set ID
#' 
#' @param geneSet character string, set ID. Use NULL to get data for all sets
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param score character string, score to retrieve from DB. Available choices:
#' gsva (default), ssgsea
#' @return data.frame with fields: tissuename, gene_set, !!score
#' @export
#' 
#' @examples 
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' getTissueDataGeneSetById("HALLMARK_P53_PATHWAY", ca)
#' 
getTissueDataGeneSetById <- function(geneSet, tissueClasses, score = c("gsva", "ssgsea")) {
  
  score <- match.arg(score)
  
  sql <- paste0(
    "SELECT tissuename, gene_set, ", score, " AS score ",
    "FROM tissue.hallmarkscore"
  )
  
  getTissueDataCommonById(
    sql = sql,
    tissues = tissueClasses,
    conditionSql = prepareConditionSql(gene_set = geneSet)
  )
}

#' Get Available Hallmark Gene Sets
#'
#' @return vector of available gene sets
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' getAvailableHallmarkGeneSets()
#'
getAvailableHallmarkGeneSets <- function(){
  sql <- "SELECT DISTINCT gene_set FROM tissue.hallmarkscore"
  getPostgresql(sql)$gene_set
}

#' Get Data for Hallmark Genesets Heatmap
#'
#' @param tissueClasses classAssignment object or list with two tissuename classes
#' @param score character string, score to retrieve from DB. Available choices:
#' gsva (default), ssgsea
#'
#' @return matrix with gene sets in rows and tissuename in columns
#' @export
#' 
#' @importFrom XIFF unlistClassAssignment
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' getHeatmapDataHallmark(ca)[1:5,1:3]
#'
getHeatmapDataHallmark <- function(tissueClasses, score = c("gsva", "ssgsea")){
  
  score <- match.arg(score)
  tissuenames <- unlistClassAssignment(tissueClasses)
  data <- getTissueDataGeneSet(tissuenames, score)
  
  if (nrow(data) > 0){
    data %>%
      pivot_wider(names_from = tissuename, values_from = score) %>%
      column_to_rownames("gene_set")
  }
}

# Dim red ---------------------------------------------------------------------
#' Get Data for Dimension Reduction of Gene Expression
#'
#' @param sampleClasses classAssignment object, or list with two classes.
#' @param ensg_gene_set NA (default, returns the most varying genes) or
#' vector of gene ensgs
#' @param p logical, ShinySession or FutureManager task object; used for progress
#'
#' @return list with elements:
#' - data.counts - matrix with gene expression counts
#' - data.log2tpm - matrix with gene expression log2tpm
#' - assignment - data.frame with tissuename and class assignment
#' 
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' res <- getExpressionDimRedData(ca, p = TRUE)
#' str(res)
#' 
#' resGeneSet <- getExpressionDimRedData(ca, "HALLMARK_P53_PATHWAY", p = TRUE)
#' str(resGeneSet)
#'
getExpressionDimRedData <- function(sampleClasses, ensg_gene_set = NA, p = FALSE) {
  progress <- ProcessProgress$new("Expression DimRed", p)
  
  assignment <- stackClasses(sampleClasses)
  data <- if (all(is.na(ensg_gene_set))) {
    progress$update(0.2, "retrieving expression of most varying genes...")
    getExpressionDimRedDataNGenes(assignment$tissuename)
    
  } else {
    progress$update(0.2, "retrieving expression of a gene sets...")
    getExpressionDimRedDataGeneSet(assignment$tissuename, unlist(ensg_gene_set))
  }
  
  if (nrow(data) == 0) progress$error("no data available")
  
  progress$update(0.5, "building cross tables...")
  
  assignment <- assignment %>% filter(tissuename %in% unique(data$tissuename))
  data.cross.counts <- tapply(data$counts, list(data$ensg, data$tissuename), sum)
  data.cross.log2tpm <- tapply(data$log2tpm, list(data$ensg, data$tissuename), sum)
  
  list(
    data.counts = data.cross.counts[, assignment$tissuename], 
    data.log2tpm = data.cross.log2tpm[, assignment$tissuename], 
    assignment = assignment
  )
}

#' Get Gene Expression Data in a Long Format for Most Varying Genes
#'
#' @param tissuenames classAssignment object or vector of tissue names
#' @param limit maximum number of most varying genes to be returned
#'
#' @return data.frame with columns: tissuenames, ensg, counts,
#' log2tpm
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' head(getExpressionDimRedDataNGenes(ca, 10))
#'
getExpressionDimRedDataNGenes <- function(tissuenames, limit = 1000){
  tissueFilter <- getSQL_filter("tissuename", tissuenames)
  sql <- paste0(
    "SELECT tissuename, ensg, counts, log2tpm ",
    "FROM tissue.processedrnaseqview ",
    "WHERE ", tissueFilter, " AND ensg = ANY(ARRAY",
    "(SELECT ensg FROM tissue.processedrnaseqview WHERE ", tissueFilter,
    " GROUP BY ensg ORDER BY variance(log2tpm) DESC LIMIT ", limit, "))"
  )
  getPostgresql(sql)
}

#' Get Gene Expression Data in a Long Format from Gene Set
#'
#' @param tissuenames classAssignment object or vector of tissue names
#' @param ensg vector of gene ensgs
#'
#' @return data.frame with columns: tissuename, ensg, counts,
#' log2tpm
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' head(getExpressionDimRedDataGeneSet(ca, "HALLMARK_P53_PATHWAY"))
#'
getExpressionDimRedDataGeneSet <- function(tissuenames, ensg){
  sql <- paste0(
    "SELECT tissuename, ensg, counts, log2tpm FROM tissue.processedrnaseqview",
    " WHERE ", getSQL_filter("tissuename", tissuenames),
    " AND ", getSQL_filter("ensg", ensg)
  )
  getPostgresql(sql)
}

# Other properties ------------------------------------------------------------

#' Add other proportions data to tissue annotation
#'
#' @param tissuename character vector with tissue names
#' @param tissueAnno tissue annotation data.frame
#' @param p logical, ShinySession or FutureManager task object; used for progress
#' @param mutationalFractionThresholds thresholds used to determine cutoffs for
#' \code{mutational_fraction_in_percent} column
#'
#' @return data.frame with the all the columns from \code{anno} and new:
#' 
#' - rnasequenced,
#' - copynumbered,
#' - mutational_fraction_in_percent,
#' - mutational_burden_class
#' - tumor_purity_class
#' 
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' anno <- getTissueAnnotation()
#' head(getOtherPropData(ca, anno))
#' 
getOtherPropData <- function(tissuename, tissueAnno, p = FALSE, 
                             mutationalFractionThresholds = c(0.1, 0.5)) {
  if (is.null(tissuename)) return()
  
  progress <- ProcessProgress$new("Other prop data", p)
  progress$update(0.5, "retrieving other tissue property data")
  
  tisstring <- paste(unlist(tissuename), collapse = "','")
  sql <- paste0("SELECT t.tissuename, dnasequenced, ",
                "CASE WHEN exists(SELECT tissuename IS NOT NULL FROM tissue.processedrnaseqview p ", 
                "WHERE p.tissuename = t.tissuename AND ensg = 'ENSG00000133703') THEN 'YES' ELSE 'NO' END AS rnasequenced, ",
                "CASE WHEN exists(SELECT tissuename IS NOT NULL FROM tissue.processedcopynumber p ", 
                "WHERE p.tissuename = t.tissuename AND ensg = 'ENSG00000133703') THEN 'YES' ELSE 'NO' END AS copynumbered, ",
                "mutational_fraction * 100 AS mutational_fraction_in_percent ",
                "FROM tissue.tissue t ",
                "LEFT OUTER JOIN tissue.mutationalburden mb ON t.tissuename = mb.tissuename ",
                "WHERE t.tissuename = ANY(ARRAY['", tisstring, "']) ")
  other_prop <- getPostgresql(sql)
  
  
  mutationalFractionLables <-
    sprintf("< %s%% mutations", mutationalFractionThresholds)
  mutationalFractionLables <-
    c(mutationalFractionLables, sprintf(">= %s%% mutations", tail(mutationalFractionThresholds, 1)))
  mutationalFractionThresholds <-
    c(-Inf, mutationalFractionThresholds, Inf)
  
  
  other_prop <- other_prop %>% 
    mutate(mutational_fraction_in_percent = ifelse(dnasequenced, mutational_fraction_in_percent, NA)) %>%
    select(-dnasequenced) %>%
    mutate(
      mutational_burden_class = cut(mutational_fraction_in_percent, breaks = mutationalFractionThresholds, right = FALSE)
    ) 
  
  levels(other_prop$mutational_burden_class) <- mutationalFractionLables
  
  alltissues <- unlist(tissuename)
  tissueAnno <- tissueAnno %>% filter(tissuename %in% alltissues)
  result <- tissueAnno %>% 
    left_join(other_prop, by = "tissuename") %>%
    mutate(
      tumor_purity_class = droplevels(cut(
        x = tumorpurity, 
        breaks = c(0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1), 
        labels = c("0-30%", "30-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%")
      ))
    )
  
  rownames(result) <- result$tissuename
  
  if(inherits(tissuename, "classAssignment")) {
    result <- addClassAssigmentAttribute(result, tissuename)
  }
  
  result
  
}

# Others ----------------------------------------------------------------------

#' Get Tissue Groups
#'
#' @return list with tissue groups
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' getTissueGroups()
#' 
getTissueGroups <- function() {
  sql <- paste0("SELECT DISTINCT tumortype FROM tissue.tissue t JOIN tissue.tissueassignment ta ON t.tissuename = ta.tissuename ",
                "WHERE tissuepanel = 'TCGA tumors' ORDER BY tumortype")
  tumor <- getPostgresql(sql)$tumortype
  
  sql <- paste0("SELECT DISTINCT tumortype_adjacent FROM tissue.tissue t JOIN tissue.tissueassignment ta ON t.tissuename = ta.tissuename ",
                "WHERE tissuepanel = 'TCGA normals' ORDER BY tumortype_adjacent")
  adjacent_normal <- getPostgresql(sql)$tumortype_adjacent
  
  sql <- paste0("SELECT DISTINCT organ FROM tissue.tissue t JOIN tissue.tissueassignment ta ON t.tissuename = ta.tissuename ",
                "WHERE tissuepanel = 'GTEx normals' ORDER BY organ")
  normals <- getPostgresql(sql)$organ
  
  sql <-  paste0("SELECT DISTINCT tumortype FROM tissue.tissue t JOIN tissue.tissueassignment ta ON t.tissuename = ta.tissuename ",
                 "WHERE tissuepanel = 'PDX Models' ORDER BY tumortype")
  
  pdx <-  getPostgresql(sql)$tumortype
  
  return(list(tumor = tumor, adjacent_normal = adjacent_normal, normal = normals, pdx = pdx))
}

#' Get Gene information from Symbol
#'
#' @param mySymbol gene symbol
#' @param species species name
#'
#' @return list with values: ensg, name, symbol
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' getGeneFromSymbol("KRAS", "human")
#'
getGeneFromSymbol <- function(mySymbol, species) { # DUPLICATE
  sql <- paste0("SELECT ensg, regexp_replace(name,' \\[Source.*\\]','') AS name ",
                "FROM gene WHERE upper(symbol) = ('", toupper(mySymbol), "') ",
                "AND species = '", species, "' AND ensg LIKE 'ENS%' AND length(chromosome) <= 2")
  annot <- getPostgresql(sql)
  if (nrow(annot) == 0) return()
  list(ensg = annot$ensg, name = annot$name, symbol = mySymbol)
}

#' Get Data Base Information
#'
#' @return data.frame with columns: Data Source, Version
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' getInformation()
#'
getInformation <- function() { # DUPLICATE
  getPostgresql("select description, information FROM information ORDER BY description") %>%
    rename(`Data Source` = description, Version = information)
}

# Common ----------------------------------------------------------------------
getTissueDataCommonById <- function(sql, tissues = NULL, conditionSql = NULL, 
                                    idCol = NULL, addRownames = FALSE, callback = NULL){
  
  sql <- paste(sql, "WHERE 1=1 ")
  
  ca <- NULL
  
  if (is.list(tissues)){
    if(inherits(tissues, "classAssignment")) {
      ca <- tissues
    }
    tissues <- c(tissues$class1, tissues$class2)
  }
  
  clSql <- if (length(tissues) > 0){
    paste0(" AND ", getSQL_filter("tissuename", tissues))
  }
  
  if (!is.null(conditionSql)){
    conditionSql <- paste("AND", conditionSql)
  }
  
  fullSql <- paste(sql, clSql, conditionSql)
  df <- getPostgresql(fullSql)
  
  if (is.function(callback)){
    df <- callback(df)
  }
  
  ti <- df$tissuename
  if (addRownames && anyDuplicated(ti) == 0){
    rownames(df) <- ti
  }
  
  if(!is.null(ca)) {
    # add class assigment attribute for further processing
    df <- addClassAssigmentAttribute(df, ca)
  }
  
  df
}



#' Prepare Prefilter Sql
#'
#' @param SELECT columns to be selected
#' @param FROM table to select from
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#' @param WHERE where statements to be added to the query
#' @param JOIN join statement
#'
#' @noRd
#' @return
#' @examples
#' 
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter(selected = c("x", "y"))
#' 
#' preparePrefilterSql(
#'   SELECT = "tbl.value",
#'   FROM = "tbl",
#'   WHERE = NULL,
#'   JOIN = "JOIN tissue.tissue ti USING(tissuename)",
#'   prefilter = prefilter
#' )
#' 
#' preparePrefilterSql(
#'   SELECT = "tbl.value",
#'   FROM = "tbl",
#'   WHERE = "a = 1 AND b = 2",
#'   JOIN = "JOIN tissue.tissue ti USING(tissuename)",
#'   prefilter = prefilter
#' )
#' 
preparePrefilterSql <- function(SELECT, FROM, prefilter, WHERE = NULL, JOIN = NULL){
  dbCol <- prefilter$db_col
  
  conditionArgs <- list(
    tissuepanel = prefilter$tissue_panel,
    dummy = prefilter$tissue
  )
  names(conditionArgs)[2] <- dbCol
  conditionSql <- do.call(prepareConditionSql, conditionArgs)
  
  if (!is.null(JOIN)){
    JOIN <- paste0(" ", trimws(JOIN))
  }
  
  addAnd <- function(x) {
    if(!is.null(x) && nchar(x) > 0) {
      paste(" AND", x) 
    }
    else {
      ""
    }
  }
  
  paste0(
    "SELECT ", paste(c(SELECT, dbCol), collapse = ", "),
    " FROM ", FROM, 
    JOIN,
    " WHERE 1=1", addAnd(WHERE), addAnd(conditionSql),
    prefilter$sql
  )
}

.exampleClassAssigment <- function(low = 0.15,
                                   high = 10,
                                   tt = c(
                                     "cholangiocarcinoma",
                                     "colon carcinoma",
                                     "melanoma",
                                     "ovarian carcinoma",
                                     "renal carcinoma",
                                     "rhabdomyosarcoma",
                                     "gastric carcinoma"
                                   )) {

  mutData <- getWaterfallDataMutationalBurden(list(db_col = "tumortype", sql = ""))
  
  if(length(tt) > 0) {
    mutData <- mutData %>% filter(tumortype %in% tt)
  }
  
  mutData <- mutData %>% mutate(
    class = case_when(
      mutational_fraction_in_percent > high ~ "highMutation",
      mutational_fraction_in_percent < low ~ "lowMutation",
      TRUE ~ "other"
    )
  ) %>% filter(class != "other") %>% select(tissuename, class) %>%
    distinct
  
  
  ca <-
    split(mutData$tissuename %>% as.character, mutData$class) %>%
    classAssignment(positiveClass = "highMutation")
  
  ca
}

#' Example Tissue Class Assigment Object
#'
#' @param low lower threshold for mutational_fraction_in_percent
#' @param high upper threshold for mutational_fraction_in_percent
#' @param tt tumor types
#'
#' @return ClassAssigment object
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' exampleClassAssigment()
#' 
exampleClassAssigment <- memoise::memoise(.exampleClassAssigment)

#' Makes prefilter object
#'
#' @param tissueType string denoting tissue type
#' @param selected vector of selected tissue types
#'
#' @return
#' @export
#'
#' @examples
#' 
#' makePrefilter("tumor")
#' 
makePrefilter <- function(tissueType = c("tumor", "adjacent normal", "normal", "PDX"),
                          selected = NULL) {
  
  if(!is.null(selected)) {
    selected <- XIFF::unlistClassAssignment(selected)  
  }
  
  tissueType <- match.arg(
    tissueType,
    choices = c("tumor", "adjacent normal", "normal", "PDX")
  )
  
  sql <- ""
  
  switch(
    EXPR = tissueType,
    "tumor" = {
      tissue <- selected
      db_col <- "tumortype"
      db_cols <- c(
        "tissuename",
        "tumortype",
        "Microsatellite stability",
        "Histology type",
        "Histology subtype",
        "gender",
        "race",
        "ethnicity",
        "Tumor purity",
        "Immune environment class", 
        "GI molecular subgroup",
        "iCluster", 
        "Consensus Molecular Subtype",
        "TIL pattern",
        "loss of Y chromosome",
        "DNA sequenced",
        "RNA sequenced",
        "Copynumber available"
      )
      tissue_panel <- "TCGA tumors"
    },
    
    "adjacent normal" = {
      tissue <- selected
      db_col <- "tumortype_adjacent"
      db_cols <- c("tissuename", "tumortype_adjacent", "DNA sequenced")
      tissue_panel <- "TCGA normals"
    },
    
    "normal" = {
      tissue <- selected
      db_col <- "organ"
      db_cols <- c(
        "tissuename", "organ", "Tissue subtype", "Autolysis score", 
        "RNA integrity number", "Minutes ischemia"
      )
      sql <- " AND tumortype = 'normal' AND tumortype_adjacent IS NULL"
      tissue_panel <- "GTEx normals"
    },
    
    "PDX" = {
      tissue <- selected
      db_col <- "tumortype"
      db_cols <- c(
        "tissuename", "tumortype", "Vendor name", "Histology type", "Histology subtype"
      )
      tissue_panel <- "PDX Models"
    }
  )
  
  res <- list(
    tissue_type = tissueType, 
    tissue = tissue, 
    tissue_panel = tissue_panel, 
    db_col = db_col, 
    db_cols = db_cols, 
    sql = sql
  )
  
  class(res) <- "TiffPrefilter"
  res
}


#' Print TiffPrefilter object
#'
#' @param x TiffPrefilter object
#' @param width max number of characters in a line. Used for printing. Default
#'        value should be good in most cases.
#' @param ... other arguments passed to the print method, discarded.
#'
#' @return
#' @export
#' @exportS3Method base::print 
#'
#' @examples
#' x <- makePrefilter("tumor")
#' x
#' 
print.TiffPrefilter <- function(x, width = NULL, ...) {
  
  xx <- x
  if(is.null(width)) {
    width <- max(getOption("width") * 0.8, 50)
  }
  
  
  if(length(x$tissue) == 0) {
    x$tissue <- "--empty--"
  }
  
  if(x$sql == "") {
    x$sql <- "--empty--"
  }
  
  cat("Tissue Prefilter:\n")
  cat("\n  tissue_type:", x$tissue_type)
  cat("\n  tissue_panel:", x$tissue_panel)
  cat("\n  db_col:", x$db_col)
  cat("\n  sql:", x$sql)
  cat("\n  db_cols:")
  cat(paste("\n    ", strwrap(paste(x$db_cols, collapse = ", "), width)))
  cat("\n  tissue:")
  cat(paste("\n    ", strwrap(paste(x$tissue, collapse = ", "), width)))

  invisible(xx)
}
