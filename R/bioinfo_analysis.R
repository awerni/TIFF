# Mutation --------------------------------------------------------------------
#' Get mutation association
#' 
#' Compares mutation data for 2 tissue classes.
#' 
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param anno data.frame, gene annotation
#' @param p logical, ShinySession or FutureManager task object; used for progress
#' @return tibble with columns: ensg, symbol, name, class2, class1, higher, 
#' p.value, adj.P.Value
#' @export
#' @examples
#' setDbOptions(getSettings())
#' tissueClasses <- exampleClassAssigment(high = 0.3)
#' anno <- getGeneAnno("human")
#' 
#' getMutationAssociation(tissueClasses, anno, p = TRUE)
#' 
getMutationAssociation <- function(tissueClasses, anno, p = FALSE) {
  progress <- ProcessProgress$new("Mutation association", p)
  progress$update(0.2, "fetching mutation data...")
  
  class1 <- tissueClasses$class1
  class2 <- tissueClasses$class2
  
  sequencedTissues <- getSequencedTissues(c(class1, class2))
  
  class1 <- intersect(class1, sequencedTissues)
  class2 <- intersect(class2, sequencedTissues)
  
  c1.all <- length(class1)
  c2.all <- length(class2)
  
  if (c1.all == 0 || c2.all == 0){
    progress$error("tissues from both classes need to be available")
  }
  
  data <- getAnalysisDataMutation(c(class1, class2))
  assignment <- stackClasses(tissueClasses, return_factor = TRUE)
  
  tissue_mutation <- data %>% 
    left_join(assignment, by = "tissuename") %>%
    mutate(class = paste0(class, ".mut")) %>%
    select(-tissuename) %>% 
    group_by(ensg, class) %>% 
    summarise(n = dplyr::n(), .groups = "drop") %>% 
    pivot_wider(names_from = class, values_from = n, values_fill = list(n = 0)) %>%
    mutate(
      class1.all = c1.all, 
      class2.all = c2.all,
      class1.wt = class1.all - class1.mut,
      class2.wt = class2.all - class2.mut
    ) %>%
    filter((class2.wt > 0 & class1.mut > 0) | (class2.mut > 0 & class1.wt > 0))
  
  progress$update(0.7, "calculating fisher's exact test...")
  
  tissue_mutation$p.value <- apply(
    X = tissue_mutation[, c("class1.mut", "class2.mut", "class1.wt", "class2.wt")], 
    MARGIN = 1, 
    FUN = function(x) fisher.test(matrix(x, nrow = 2))$p.value
  )
  
  tissue_mutation <- tissue_mutation %>% 
    mutate(
      adj.P.Value = p.adjust(p.value, method = "fdr"),
      class2Calc = class2.mut/class2.all, 
      class1Calc = class1.mut/class1.all,
      class2 = paste0(signif(100 * class2Calc, 3), "% (", class2.mut, "/", class2.all, ")"),
      class1 = paste0(signif(100 * class1Calc, 3), "% (", class1.mut, "/", class1.all, ")"),
      class2 = forcats::fct_reorder(class2, class2Calc), 
      class1 = forcats::fct_reorder(class1, class1Calc),
      higher = ifelse(class1Calc >= class2Calc, "class1", "class2")
    )
  
  progress$update(0.9, "attaching gene annotation...")
  
  anno %>% 
    inner_join(tissue_mutation, by = c("ensg" = "ensg")) %>% 
    select(ensg, symbol, name, class2, class1, higher, p.value, adj.P.Value) %>% 
    arrange(p.value, adj.P.Value)
}

# Copy number -----------------------------------------------------------------
#' Get copy number association
#' 
#' Compares copy number data for 2 tissue classes.
#' 
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param anno data.frame, gene annotation
#' @param p logical, ShinySession or FutureManager task object; used for progress
#' @return tibble with columns: ensg, symbol, name, location, class1.median, 
#' class2.median, log2FC, p.value, adj.p.val
#' @export
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @examples
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' anno <- getGeneAnno("human")
#' 
#' copynumberAssoc <- getCopyNumberAssociation(ca, anno, p = TRUE)
#' head(copynumberAssoc)
#' 
getCopyNumberAssociation <- function(tissueClasses, anno, p = FALSE) {
  progress <- ProcessProgress$new("Copy number association", p)
  progress$update(0.2, "fetching copy number data...")
  
  ti <- stackClasses(tissueClasses)
  
  allTissues <- ti$tissuename
  data <- getTissueDataCopyNumber(allTissues) %>%
    filter(!is.na(log2relativecopynumber))
  
  if (nrow(data) == 0) progress$error("no copy number data available")
  
  allTissues <- tibble(tissuename = data$tissuename %>% unique()) %>%
    left_join(ti, by = "tissuename")
  
  avail_tissues <- allTissues %>% 
    group_by(class) %>% 
    summarise(n = dplyr::n(), .groups = "drop")
  
  if (nrow(avail_tissues) < 2 | all(avail_tissues$n < 2)) {
    progress$error("tissues from both classes need to be available")
  }
  
  progress$update(0.4, "filtering genes with copy number alterations...")
  
  e <- data %>% 
    inner_join(allTissues, by = "tissuename") %>%
    group_by(ensg, class) %>% 
    summarise(
      ampl = sum(log2relativecopynumber > log2(3/2)), 
      del = sum(log2relativecopynumber < log2(1.5/2)), 
      .groups = "drop"
    ) %>%
    filter(ampl >= 2 | del >= 2) %>%
    select(ensg) %>%
    unique()
  
  if (nrow(e) == 0) progress$error("not enough data available")
  
  progress$update(0.5, "preparing data for statistical test...")
  
  data_filtered <- data %>%
    filter(ensg %in% e$ensg) %>%
    left_join(allTissues, by = "tissuename") 
  
  fold_change <- data_filtered %>%
    group_by(class, ensg) %>% 
    summarise(medianCN = median(log2relativecopynumber), .groups = "drop") %>%
    pivot_wider(names_from = class, values_from = medianCN) %>%
    mutate(
      log2FC = class1 - class2, 
      class1.median = 2*2^class1, 
      class2.median = 2*2^class2
    ) %>%
    select(-class1, -class2)
  
  data_wide <- purrr::map(
    data_filtered %>% split(.$class), 
    function(u) {
      u %>% 
        select(-class) %>% 
        pivot_wider(names_from = tissuename, values_from = log2relativecopynumber) %>% 
        # arrange(ensg) %>% # no need to sort, data already sorted in DB
        column_to_rownames("ensg") %>% 
        as.matrix()
    })
  
  progress$update(0.6, "calculating wilcoxon rank sum test...")
  
  data_wide <- ensureCommonRownames(
    m1 = data_wide$class1, 
    m2 = data_wide$class2, 
    outNames = c("class1", "class2")
  )
  
  suppressWarnings({
    result <- matrixTests::row_wilcoxon_twosample(
      x = data_wide$class1, 
      y = data_wide$class2, 
      exact = FALSE
    ) %>% 
      rownames_to_column("ensg") %>% 
      rename(p.value = pvalue) %>% 
      inner_join(fold_change, by = "ensg")
  })
  
  progress$update(0.85, "adjusting p-value...")
  
  result <- result %>% mutate(adj.P.Value = p.adjust(p.value, method = "fdr"))
  
  progress$update(0.9, "adding gene annotation...")
  
  anno %>% 
    inner_join(result, by = "ensg") %>% 
    select(
      ensg, symbol, name, location, class1.median, class2.median, 
      log2FC, p.value, adj.P.Value
    ) %>% 
    arrange(p.value, desc(abs(log2FC)))
}

# Expression ------------------------------------------------------------------
#' Differential gene expression
#' 
#' This function uses limma package to determine the differentially expressed
#' genes or transcripts across 2 classes.
#' 
#' If you want to report progress, please provide a value for `p` argument. See 
#' `ProcessProgress` for more details.
#' 
#' TODO: Add tpmHET. It's not supported yet.
#' 
#' As a result, this function returns a data.frame
#' contains the statistics for genes or transcripts. 
#' `correctTpms` is a column of determining if gene or transcript ID is  
#' highly expressed, i.e. tpm > `tpmHET`
#' (`tpmHET` - tpm High Expression Threshold, default 10, it can be set in the
#' `differentialGeneExpression_LimmaVoom` function or by using `cliff.tpmHET` 
#' option) for at least 1 tissue.
#' 
#' 
#' Columns available in `df`: ensg or enst, `anno` columns, P.Value, adj.p.val, 
#' logFC, auc.
#' 
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param anno data.frame, gene annotation
#' @param p logical, ShinySession or FutureManager task object; used for progress
#' 
#' @return data.frame
#' @export
#' @examples
#' \dontrun{
#' setDbOptions(getSettings())
#' tissueClasses <- classAssignment(sensitive = c("a", "b", "c"), resistant = c("d", "e", "f"))
#' gAnno <- getGeneAnno("human")
#' 
#' gRes <- differentialGeneExpression_LimmaVoom(tissueClasses, gAnno, p = TRUE)
#' }
differentialGeneExpression_LimmaVoom <- function(tissueClasses, anno, p = FALSE) {
  progress <- ProcessProgress$new("Differential gene expression", p)
  progress$update(0.2, "fetching expression count data...")
  
  class1 <- tissueClasses$class1
  class2 <- tissueClasses$class2
  
  options(stringsAsFactors = FALSE)
  exprData <- getTissueCounts(tissueClasses)
  
  if (is.character(exprData)) progress$error(exprData)
  
  expr <- exprData$expr
  design <- exprData$design
  tpm <- exprData$tpm
  
  progress$update(0.3, "normalizing counts...")
  
  dge <- edgeR::DGEList(counts = expr)
  dge <- edgeR::calcNormFactors(dge)
  
  progress$update(0.5, "voom transformation...")
  
  v <- limma::voom(dge, design, plot = FALSE)
  
  progress$update(0.7, "fitting linear model...")
  
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  tt <-  limma::topTable(fit, number = 1e9, coef = ncol(design))
  res <- data.frame(ensg = rownames(tt), tt)
  
  progress$update(0.9, "loading gene annotation...")
  
  correctTpms <- rownames(tpm)[rowSums(tpm > 30, na.rm = TRUE) > 0]
  res <- anno %>% inner_join(res, by = "ensg")
  rownames(res) <- res$ensg
  
  progress$update(1.0, "job done")
  
  df <- res %>% 
    rename(adj.p.val = adj.P.Val) %>% 
    arrange(P.Value, adj.p.val)
  
  df <- df %>% addClassAssigmentAttribute(tissueClasses)
  df[["correctTpms"]] <- df[["ensg"]] %in% correctTpms
  df
}

#' Get tissue expression data
#' 
#' This is a helper function that retrieves expression data for the specified 
#' tissue. 
#' 
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @return list with fields: expr, tpm, and design
#' @export
#' @examples
#' 
#' setDbOptions(getSettings())
#' tissueClasses <- exampleClassAssigment()
#' tissueCounts <- getTissueCounts(tissueClasses)
#' 
getTissueCounts <- function(tissueClasses) {
  class1 <- tissueClasses$class1
  class2 <- tissueClasses$class2
  
  classification <- prepareClassificationTable(class1, class2, addRownames = FALSE)
  
  tissuenames <- classification$tissuename
  counts.long <- getTissueDataGeneExpression(tissuenames, c("counts", "log2tpm"))
  
  tpm <- counts.long %>%
    select(tissuename, ensg, log2tpm) %>%
    mutate(tpm = 2^log2tpm, log2tpm = NULL)
  
  counts.long <- counts.long %>% select(-log2tpm)
  
  avail.tissues <- unique(counts.long$tissuename)
  missing.tissues <- setdiff(tissuenames, avail.tissues)
  if (length(missing.tissues) > 0) {
    classification <- classification %>% filter(tissuename %in% avail.tissues)
  }
  
  class1 <- intersect(class1, avail.tissues)
  class2 <- intersect(class2, avail.tissues)
  
  if (length(class1) == 0 || length(class2) == 0){
    return("tissues from both classes need to be available")
  }
  
  if (length(class1) == 1 && length(class2) == 1){
    return("1 vs 1 tissue comparison not possible")
  }
  
  ensg_sufficient_counts <- counts.long %>%
    group_by(ensg) %>%
    summarise(count_max_ok = max(counts, na.rm = TRUE) >= 20, .groups = "drop") %>%
    filter(count_max_ok) %>%
    .$ensg
  
  ensg_sufficient_class_counts <- counts.long %>%
    inner_join(classification, by = "tissuename") %>%
    group_by(ensg, class) %>%
    summarise(count_mean_ok = mean(counts, na.rm = TRUE) >= 5, .groups = "drop") %>%
    filter(count_mean_ok) %>%
    .$ensg %>%
    unique()
  
  ensg_ok <- intersect(ensg_sufficient_counts, ensg_sufficient_class_counts)
  
  expr <- counts.long %>%
    filter(ensg %in% ensg_ok) %>%
    pivot_wider(names_from = tissuename, values_from = counts) %>%
    column_to_rownames("ensg") %>%
    as.matrix()
  
  tpm <- tpm %>%
    filter(ensg %in% ensg_ok) %>%
    pivot_wider(names_from = tissuename, values_from = tpm) %>%
    column_to_rownames("ensg") %>%
    as.matrix()
  
  expr <- expr[, classification$tissuename]
  tpm <- tpm[, classification$tissuename]
  
  design <- model.matrix(
    object = ~ class, 
    data = classification %>% column_to_rownames("tissuename")
  )
  
  list(
    expr = expr,
    design = design,
    tpm = tpm
  )
}

# GSEA ------------------------------------------------------------------------
#' Gene Set Enrichment Analysis
#' 
#' Calculate GSEA for differential gene expression results
#' 
#' @param diffExResult data.frame, table result obtained from 
#' differentialGeneExpression_LimmaVoom() function
#' @param geneSets list of gene IDs, gene sets
#' @param rankType character string, "p.valueDir" or "logFC" 
#' @return tibble with columns: p.geomean, stat.mean, p.val, q.val, set.size, 
#' exp1, geneset
#' @export
#' @examples
#' \donttest{
#' setDbOptions(getSettings())
#' tissueClasses <- exampleClassAssigment()
#' anno <- getGeneAnno("human")
#' 
#' dge <- differentialGeneExpression_LimmaVoom(tissueClasses, anno)
#' geneSets <- getGSEAdata("human", "hallmark")
#' getGSEA(dge, geneSets, "p.valueDir", p = TRUE)
#' }
getGSEA <- function(diffExResult, geneSets, rankType, p = FALSE, minVal = 1e-300) {
  progress <- ProcessProgress$new("GSEA", p)
  progress$update(0.2, "preparing data for GSEA...")
  
  if (rankType == "logFC") {
    a <- diffExResult$logFC
  } else {
    diffExResult <- diffExResult %>% 
      mutate(
        P.Value = ifelse(P.Value < minVal, minVal, P.Value), # manage 0 p.values
        pRank = ifelse(logFC < 0, log(P.Value), -log(P.Value))
      )
    a <- diffExResult$pRank
  }
  
  names(a) <- diffExResult$ensg
  
  progress$update(0.6, "calculating GSEA...")
  
  GSEAresult <- gage::gage(
    exprs = a, 
    gsets = geneSets, 
    set.size = c(10, 500), 
    same.dir = FALSE
  )
  
  progress$update(1.0, "job done")
  
  GSEAresult$greater %>% 
    as_tibble() %>% 
    mutate(geneset = rownames(GSEAresult$greater)) %>% 
    filter(!is.nan(stat.mean))
}

# Proteins --------------------------------------------------------------------
#' Differential protein expression
#' 
#' Compares protein expression data for 2 tissue classes.
#' 
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param p logical, ShinySession or FutureManager task object; used for progress
#' @return tibble with columns: antibody, higher, log2FC, P.Value, adj.p.val
#' @export
#' @examples
#' setDbOptions(getSettings())
#' tissueClasses <- exampleClassAssigment(high = 1)
#' 
#' differentialProteinExpression(tissueClasses, p = TRUE)
#' 
differentialProteinExpression <- function(tissueClasses, p = FALSE) {
  differentialBayesCommon(
    sampleClasses = tissueClasses,
    dbDataFun = getTissueDataRppa,
    idCol = "antibody",
    scoreCol = "score",
    progressLabel = "Differential protein expression",
    itemLabel = "protein expression",
    p = p
  )
}

# Immune cells ----------------------------------------------------------------

#' Compare Immune Cells
#'
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param p logical, ShinySession or FutureManager task object; used for progress
#'
#' @return data.frame with columns `cell type`, method, log2FC, p.value, adj.P.Value
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment(high = 0.16)
#' head(compareImmuneCells(ca))
#' 
compareImmuneCells <- function(tissueClasses, p = FALSE){
  if (is.null(tissueClasses)) return()
  
  progress <- ProcessProgress$new("Immune cells comparison", p)
  progress$update(0.2, "retrieving immune cells data...")
  
  class1 <- tissueClasses$class1
  class2 <- tissueClasses$class2
  
  data <- getTissueDataImmuneCells(c(class1, class2))
  if (nrow(data) <= 2) progress$error("not enough data")
  
  progress$update(0.5, "preparing data for statistical test...")
  
  avail.tis <- unique(data$tissuename)
  class1 <- intersect(class1, avail.tis)
  class2 <- intersect(class2, avail.tis)
  
  if (length(class1) == 0 || length(class2) == 0){
    progress$error("tissues from both classes need to be available")
  }
  
  assignment <- stackClasses(tissueClasses)
  immuneCellData <- data %>%
    pivot_wider(names_from = celltype, values_from = score) %>%
    left_join(assignment, by = "tissuename")
  
  data_num <- immuneCellData %>%
    select_if(is.numeric) %>%
    select_if(function(x) isColumnValidForAnalysis(x, immuneCellData$class))
  
  if (ncol(data_num) == 0) {
    progress$error("no numeric data available")
  }
  
  progress$update(0.7, "calculating statistics...")
  
  log2fc <- data_num %>% purrr::map_dbl(~ getLog2FC(.x, immuneCellData$class))
  
  data_num %>% 
    purrr::map(~ wilcox2.test(.x, immuneCellData$class)) %>%
    purrr::map_dfr(broom::tidy, .id = "cell type") %>%
    select(-alternative) %>%
    mutate(
      adj.P.Value = p.adjust(p.value, method = "fdr"),
      log2FC = log2fc
    ) %>%
    select(`cell type`, method, log2FC, p.value, adj.P.Value) %>% 
    arrange(p.value)
}

# Metabolics ------------------------------------------------------------------
#' Compare Metabolics
#'
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param p logical, ShinySession or FutureManager task object; used for progress
#'
#' @return data.frame with columns `metabolic pathway`, method, diff, p.value, adj.P.Value
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment(high = 0.16)
#' head(compareMetabolics(ca))
#' 
compareMetabolics <- function(tissueClasses, p = FALSE){
  if (is.null(tissueClasses)) return()
  
  progress <- ProcessProgress$new("Metabolics comparison", p)
  progress$update(0.2, "retrieving metabolics data...")
  
  class1 <- tissueClasses$class1
  class2 <- tissueClasses$class2
  
  data <- getTissueDataMetabolics(c(class1, class2))
  if (nrow(data) <= 2) progress$error("not enough data")
  
  progress$update(0.5, "preparing data for statistical test...")
  
  avail.tis <- unique(data$tissuename)
  class1 <- intersect(class1, avail.tis)
  class2 <- intersect(class2, avail.tis)
  
  if (length(class1) == 0 || length(class2) == 0){
    progress$error("tissues from both classes need to be available")
  }
  
  assignment <- stackClasses(tissueClasses)
  metabolicsData <- data %>%
    pivot_wider(names_from = metabolic_pathway, values_from = score) %>%
    left_join(assignment, by = "tissuename")
  
  data_num <- metabolicsData %>%
    select_if(is.numeric) %>%
    select_if(function(x) isColumnValidForAnalysis(x, metabolicsData$class))
  
  if (ncol(data_num) == 0) {
    progress$error("no numeric data available")
  }
  
  medianDiff <- data_num %>% purrr::map_dbl(~ getMedianDiff(.x, metabolicsData$class))
  
  data_num %>%
    purrr::map(~ wilcox2.test(.x, metabolicsData$class)) %>%
    purrr::map_dfr(broom::tidy, .id = "metabolic pathway") %>%
    select(-alternative) %>%
    mutate(
      adj.P.Value = p.adjust(p.value, method = "fdr"),
      diff = medianDiff
    ) %>%
    select(`metabolic pathway`, method, diff, p.value, adj.P.Value) %>% 
    arrange(p.value)
}

# Signaling -------------------------------------------------------------------
#' Compare Signaling
#'
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param p logical, ShinySession or FutureManager task object; used for progress
#'
#' @return data.frame with columns pathway, class1.n_active, class1.n_total,
#' class1.percentage, class2.n_active, class2.n_total, class2.percentage, higher, p.value
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment(high = 0.16)
#' head(compareSignaling(ca))
#' 
compareSignaling <- function(tissueClasses, p = FALSE){
  if (is.null(tissueClasses)) return()
  
  progress <- ProcessProgress$new("Signaling comparison", p)
  progress$update(0.2, "retrieving signaling data...")
  
  class1 <- tissueClasses$class1
  class2 <- tissueClasses$class2
  
  data <- getTissueDataSignaling(c(class1, class2)) %>% filter(!is.na(active))
  if (nrow(data) <= 2) progress$error("not enough data")
  
  progress$update(0.5, "preparing data for statistical test...")
  
  avail.tis <- unique(data$tissuename)
  class1 <- intersect(class1, avail.tis)
  class2 <- intersect(class2, avail.tis)
  
  if (length(class1) == 0 || length(class2) == 0){
    progress$error("tissues from both classes need to be available")
  }
  
  signalingData <- prepareSignalingDataForAnalysis(data, tissueClasses)
  if (is.character(signalingData)) progress$error(signalingData)
  
  progress$update(0.7, "calculating statistics...")
  
  resPercentage <- calculateSignalingPercentageSummary(signalingData)
  resStats <- calculateSignalingStatistics(signalingData)
  
  resPercentage %>% 
    left_join(resStats, by = "pathway") %>%
    arrange(p.value)
}

prepareSignalingDataForAnalysis <- function(data, tissueClasses){
  assignment <- stackClasses(tissueClasses)
  df <- data %>% 
    mutate(
      active = ifelse(active, "active", "inactive"),
      active = factor(active, levels = c("active", "inactive"))
    ) %>%
    left_join(assignment, by = "tissuename") %>%
    select(-tissuename)
  
  validPathways <- getValidSignalingPathways(df)
  if (length(validPathways) == 0) return("no valid data available")
  
  df %>% filter(pathway %in% validPathways)
}

getValidSignalingPathways <- function(df){
  df %>%
    group_by(pathway) %>%
    summarise(n_classes = dplyr::n_distinct(class)) %>%
    filter(n_classes == 2) %>% # class1 and class2
    pull(pathway)
}

calculateSignalingPercentageSummary <- function(df){
  res1 <- df %>%
    group_by(class, pathway, active, .drop = FALSE) %>%
    summarise(n = dplyr::n(), .groups = "drop")
  
  res2 <- res1 %>%
    group_by(class, pathway) %>%
    summarise(n_total = sum(n), .groups = "drop")
  
  res1 %>%
    filter(active == "active") %>%
    left_join(res2, by = c("class", "pathway")) %>%
    arrange(pathway, class) %>%
    group_by(pathway) %>%
    summarise(
      class1.n_active = n[1],
      class1.n_total = n_total[1],
      class1.percentage = class1.n_active / class1.n_total,
      class2.n_active = n[2],
      class2.n_total = n_total[2],
      class2.percentage = class2.n_active / class2.n_total,
      higher = `if`(class1.percentage > class2.percentage, "class1", "class2"),
      .groups = "drop"
    )
}

calculateSignalingStatistics <- function(df, adjMethod = "fdr"){
  df %>%
    split(., .$pathway) %>%
    purrr::map(~ fisher.test(
      x = factor(.x$active),
      y = factor(.x$class), 
      simulate.p.value = TRUE, 
      conf.int = FALSE, 
      B = 10000
    )) %>%
    purrr::map_dfr(broom::tidy, .id = "pathway") %>%
    select(pathway, p.value) %>%
    mutate(adj.P.Value = p.adjust(p.value, method = adjMethod))
}

# Clones ----------------------------------------------------------------------

#' Compare Clones
#'
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param p logical, ShinySession or FutureManager task object; used for progress
#'
#' @return data.frame
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' 
#' ca <- exampleClassAssigment(high = 0.16)
#' compareClones(ca)
#' 
compareClones <- function(tissueClasses, p = FALSE){
  if (is.null(tissueClasses)) return()
  
  progress <- ProcessProgress$new("Clones comparison", p)
  progress$update(0.2, "retrieving clone data")
  
  class1 <- tissueClasses$class1
  class2 <- tissueClasses$class2
  
  data <- getTissueDataClones(c(class1, class2), c("number_of_clones", "clone_tree_score"))
  if (nrow(data) <= 2) progress$error("not enough data")
  
  progress$update(0.5, "preparing data for statistical test...")
  
  avail.tis <- unique(data$tissuename)
  class1 <- intersect(class1, avail.tis)
  class2 <- intersect(class2, avail.tis)
  
  if (length(class1) == 0 || length(class2) == 0){
    progress$error("tissues from both classes need to be available")
  }
  
  assignment <- stackClasses(tissueClasses)
  data <- data %>% left_join(assignment, by = "tissuename")
  data_num <- data %>% 
    select_if(is.numeric) %>%
    select_if(~ isColumnValidForAnalysis(.x, data$class))
  
  if (ncol(data_num) == 0) {
    progress$error("no numeric data available")
  }
  
  progress$update(0.7, "calculating clone statistics...")
  
  medianDiff <- data_num %>% purrr::map_dbl(~ getMedianDiff(.x, data$class))
  
  data_num %>%
    purrr::map(~ wilcox2.test(.x, data$class)) %>%
    purrr::map_dfr(broom::tidy, .id = "clone score") %>%
    select(-alternative) %>%
    mutate(
      adj.P.Value = p.adjust(p.value, method = "bonferroni"),
      diff = medianDiff
    ) %>%
    select(`clone score`, method, diff, p.value, adj.P.Value) %>% 
    arrange(p.value)
}

# Survival --------------------------------------------------------------------

#' Perform Survival Analysis
#'
#' @param data data.frame with days_to_death, days_to_last_followup and class 
#' columns.
#'
#' @return list with \code{df} object used to fit the model and \code{fit}
#' containing survfit object.
#' 
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' 
#' data <- getPatientAnnotation()
#' ca <- exampleClassAssigment()
#' ca  <- classAssignment2df(ca)
#' dataSurv <- inner_join(data, ca)
#' survivalAnalysis(dataSurv)
#' 
survivalAnalysis <- function(data){
  if (is.null(data)) return()
  
  df <- as_tibble(data) %>%
    mutate(
      vital_status = if_else(is.na(days_to_death), 0, 1), 
      times = coalesce(days_to_death, days_to_last_followup)
    ) %>%
    select(tissuename, class, times, vital_status)
  
  fit <- survival::survfit(survival::Surv(times, vital_status) ~ class, data = df)
  
  list(
    df = df,
    fit = fit
  )
}

# Gene signatures -------------------------------------------------------------
#' Differential gene signatures
#' 
#' Compares gene signature data for 2 tissue classes.
#' 
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param p logical, ShinySession or FutureManager task object; used for progress
#' @return tibble with columns: signature, higher, logFC, P.Value, adj.p.val
#' @export
#' @examples
#' \dontrun{
#' setDbOptions(getSettings())
#' tissueClasses <- classAssignment(sensitive = c("a", "b", "c"), resistant = c("d", "e", "f"))
#' 
#' differentialGeneSignature(tissueClasses, p = TRUE)
#' }
differentialGeneSignature <- function(tissueClasses, p = FALSE) {
  differentialBayesCommon(
    sampleClasses = tissueClasses,
    dbDataFun = getTissueDataGeneSignature,
    idCol = "signature",
    scoreCol = "score",
    progressLabel = "Differential gene signature",
    itemLabel = "gene signature",
    p = p
  )
}

# Hallmark sets ---------------------------------------------------------------
#' Differential Hallmark sets
#' 
#' Compares Hallmark sets data for 2 tissue classes.
#' 
#' @param tissueClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param p logical, ShinySession or FutureManager task object; used for progress
#' @return tibble with columns: gene_set, higher, logFC, P.Value, adj.p.val
#' @export
#' @examples
#' \dontrun{
#' setDbOptions(getSettings())
#' tissueClasses <- classAssignment(sensitive = c("a", "b", "c"), resistant = c("d", "e", "f"))
#' 
#' differentialHallmarkSets(tissueClasses, p = TRUE)
#' }
differentialHallmarkSets <- function(tissueClasses, p = FALSE) {
  differentialBayesCommon(
    sampleClasses = tissueClasses,
    dbDataFun = getTissueDataGeneSet,
    idCol = "gene_set",
    scoreCol = "score",
    progressLabel = "Differential gene sets",
    itemLabel = "gene set",
    p = p
  )
}

# Multiple features -----------------------------------------------------------


#' Get Multiple Features Axis Data
#'
#' @param sourceType data source type ('expr' or 'copy')
#' @param gene ensg or list with gene ensg and symbol
#' @param axis name of the target axis ('x' or 'y')
#' @param prefilter TiffPrefilter object (see \code{\link{makePrefilter}})
#'
#' @return list with varLabel and scale name, and data.frame with
#' tissuename {x/y}Var and tumortype.
#' 
#' @export
#' @examples
#' 
#' setDbOptions(getSettings())
#' 
#' prefilter <- makePrefilter()
#' getMultipleFeaturesAxisData("expr", "ENSG00000133703", "x", prefilter)
#' 
getMultipleFeaturesAxisData <- function(sourceType, gene, axis, prefilter){
  
  if(is.character(gene)) {
    gene <- list(ensg = gene, symbol = XIFF::getGeneSymbol(gene))
  }
  
  scalePrefix <- paste0("scale_", axis)
  scale <- "continuous"
  
  if (sourceType == "expr"){
    d <- getWaterfallDataGeneExpression(gene$ensg, prefilter) 
    
    var <- "tpm"
    varLabel <- "TPM"
    scale <- "log10"
  } else if (sourceType == "copy"){
    d <- getWaterfallDataCopyNumber(gene$ensg, prefilter) 
    
    var <- "relativecopynumber"
    varLabel <- "relative copy number"
    scale <- "sqrt"
  } else {
    stop("sourceType can only be 'expr' or 'copy'.")
  }
  
  dict <- structure(
    c(var, "tissuename"),
    names = c(paste0(axis, "Var"), "tissuename")
  )
  
  if (!is.null(d)){
    d <- d %>%
      mutate(tissuename = as.character(tissuename)) %>%
      rename(!!dict)
  }
  
  list(
    d = d,
    varLabel = sprintf("%s [%s]", gene$symbol, varLabel),
    scale = paste0(scalePrefix, "_", scale)
  )
}

# Other properties ------------------------------------------------------------
#' Get other properties statistics
#' 
#' Calculates statistics for continuous variables (Wilcoxon test) and class 
#' variables (Fisher's exact test) available in annotation.
#' 
#' @param sampleClasses classAssignment object or a 
#' list with fields: class1 and class2. Each field should
#' contain a character vector of tissue names
#' @param sAnno data.frame, tissue annotation
#' @param p logical, ShinySession or FutureManager task object; used for progress
#' @return data.frame with columns: property, p.value, method, estimate, data, 
#' statistic, higher
#' @export
#' @examples
#' setDbOptions(getSettings())
#' tissueClasses <- exampleClassAssigment(low = 0.05, tt = NULL)
#' anno <- getTissueAnnotation()
#' 
#' getOtherPropStatistics(tissueClasses, anno, p = TRUE)
#' 
getOtherPropStatistics <- function(sampleClasses, sAnno, p = FALSE) {
  if (is.null(sampleClasses)) return()
  
  progress <- ProcessProgress$new("Other prop stats", p)
  progress$update(0.35, "calculating class property statistics...")
  
  assignment <- stackClasses(sampleClasses, return_factor = TRUE)
  sAnno <- sAnno %>% inner_join(assignment, by = "tissuename")
  
  data_factor <- sAnno %>% 
    select_if(is.factor) %>% 
    droplevels() %>%
    select_if(~ nlevels(.x) > 1) %>%
    select_if(~ isColumnValidForAnalysis(.x, sAnno$class)) %>% 
    select(-class)
  
  classResult <- if (ncol(data_factor) > 0) {
    data_factor %>% 
      purrr::map(~ fisher.test(
        x = .x, 
        y = sAnno$class, 
        simulate.p.value = TRUE, 
        conf.int = FALSE, 
        B = 10000
      )) %>% 
      purrr::map_dfr(broom::tidy, .id = "property") %>% 
      select(-matches("alternative")) %>%
      mutate_at(vars(matches("estimate"), "p.value"), signif, 3) %>%
      mutate(
        data = "class", 
        method = as.character(method)
      ) %>%
      filter(property != "class")
  } # else NULL
  
  progress$update(0.65, "calculating continuous property statistics...")
  
  data_num <- sAnno %>%
    select_if(is.numeric) %>%
    select_if(~ isColumnValidForAnalysis(.x, sAnno$class))
  
  wilcox2.testExt <- function(x, class) {
    y <- split(x, class)
    res <- wilcox.test(y[[1]], y[[2]])
    
    df <- broom::tidy(res, stringsAsFactors = FALSE)
    
    m1 <- median(y[[1]], na.rm = TRUE)
    m2 <- median(y[[2]], na.rm = TRUE)
    df$higher <- if (m1 > m2) "class1" else "class2"
    df
  }
  
  numResult <- if (ncol(data_num) > 0) {  
    data_num %>%
      purrr::map_dfr(~ wilcox2.testExt(.x, sAnno$class), .id = "property") %>%
      select(-alternative) %>%
      mutate_at(c("p.value"), signif, 3)  %>%
      mutate(data = "numeric", method = as.character(method))
  } # else NULL
  
  if (!is.null(classResult) & !is.null(numResult)) {
    bind_rows(classResult, numResult) %>% arrange(p.value)
  } else if (!is.null(classResult)) {
    classResult %>% arrange(p.value)
  } else if (!is.null(numResult)) {
    numResult %>% arrange(p.value)
  } # else NULL
}

# Helpers ---------------------------------------------------------------------
wilcox2.test <- function(x, class) {
  y <- split(x, class)
  wilcox.test(y[[1]], y[[2]])
}

getMedianDiff <- function(x, class) {
  m <- split(x, class) %>% purrr::map_dbl(median, na.rm = TRUE)
  m[1] - m[2]
}

getLog2FC <- function(x, class) {
  m <- split(x, class) %>% purrr::map_dbl(median, na.rm = TRUE)
  if (m[2] == 0) return(+Inf)
  if (any(m < 0)) return(NaN)
  log2(m[1]/m[2])
}

isColumnValidForAnalysis <- function(x, class) {
  tmp <- tibble(
    x = x, 
    class = as.character(class)
  ) %>%
    filter(!is.na(x) & !is.nan(x)) %>%
    group_by(class) %>%
    summarise(n = dplyr::n())
  
  nrow(tmp) > 1
}
