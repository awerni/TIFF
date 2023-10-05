# Mutation plots --------------------------------------------------------------
#' Generate waterfall plot for mutational burden
#' 
#' Produces a waterfall plot with mutational burden data
#' 
#' @param data data.frame, output of getWaterfallDataMutationalBurden
#' @param fill character string, variable to use for bars coloring
#' @return ggplot object
#' @export
#' @examples
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter("tumor", c("thymoma", "mesothelioma"))
#' df <- getWaterfallDataMutationalBurden(prefilter)
#' generateMutationalBurdenWaterfallPlot(df)
#'
generateMutationalBurdenWaterfallPlot <- function(data, fill = "tumortype") {
  if (is.null(data)) return()
  
  generateWaterfallPlot(
    data = data,
    dataCol = "mutational_fraction_in_percent",
    ylabel = "mutational fraction [% mutations of all genes]",
    trans = "log10",
    fill = fill
  )
}

#' Generate mutation plot by type
#' 
#' Produces mutation plot that shows class comparison. Available types: 
#' "bar" (default), "pie", "coverage". Coverage plot may be used to show 
#' the data availability in each class.
#' 
#' @param df data.frame, output of getTissueDataMutationById()
#' @param plotType character string, plot type
#' @param ca classAssignment object if it cannot be inferred from \code{df} param.
#'        However if df is an output of getTissueDataMutationById, then this
#'        parameter does not need to be specified.
#' @param ensg gene ensg. In most cases there's no need to provide that value. 
#'        It will be interfered from \code{df}.
#' @param geneSymbol if null \code{XIFF::getGeneSymbol(ensg)} will be used to get 
#'        gene symbol.
#' @return ggplot object
#' @export
#' 
#' @importFrom XIFF getGeneSymbol
#' @examples
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' ensg <- "ENSG00000133703"
#' df <- getTissueDataMutationById(ensg, ca)
#' symbol <- XIFF::getGeneSymbol(ensg)
#' generateMutationPlot(df, "bar")
#' generateMutationPlot(df, "pie")
#' generateMutationPlot(df, "coverage")
#'
generateMutationPlot <- function(df, plotType = "bar", ca = NULL, ensg = NULL, geneSymbol = NULL){
  
  ca <- getClassAssigmentAttributeIfNull(ca, df)
  
  if(is.null(ensg)) {
    ensg <- df$ensg[[1]]
  }
  
  if(is.null(geneSymbol)) {
    geneSymbol <- XIFF::getGeneSymbol(ensg)
  }
  
  gene <- list(
    symbol = geneSymbol,
    ensg = ensg
  )
  
  sampleClasses <- ca
  classLabel    <- getClassLabel(ca)
    
  switch(
    EXPR = plotType,
    bar = generateMutationBarPlot(gene, df, sampleClasses, classLabel),
    pie = generateMutationPieChart(gene, df, sampleClasses, classLabel),
    coverage = generateDataCoveragePlot(df, "aamutation", ca),
    all = list(
      bar = generateMutationBarPlot(gene, df, sampleClasses, classLabel),
      pie = generateMutationPieChart(gene, df, sampleClasses, classLabel),
      coverage = generateDataCoveragePlot(df, "aamutation", ca)
    )
  )
}

generateMutationBarPlot <- function(gene, mutation, sampleClasses, classLabel) {
  if (is.null(gene) | is.null(mutation) | is.null(sampleClasses)) return()
  
  assignment <- stackClasses(sampleClasses, classLabel, return_factor = TRUE)
  data <- assignment %>% 
    inner_join(mutation, by = "tissuename") %>% 
    filter(!is.na(aamutation))
  
  if (length(unique(data$aamutation)) > 20) {
    new_mut <- forcats::fct_lump(data$aamutation[data$aamutated == "mut"], prop = 0.1, other_level = "other mut")
    data$aamutation[data$aamutated == "mut"] <- as.character(new_mut)
  }
  myLabel <- paste0(ifelse(is.na(gene$symbol), gene$ensg, gene$symbol), "\nstatus")
  
  ggplot(data, aes(x = class, fill = aamutation)) + 
    geom_bar() + 
    coord_flip() + 
    # scale_fill_viridis(discrete = TRUE) +
    ggtitle(paste("\n", gene$symbol, "-", gene$ensg)) + 
    xlab("") + 
    ylab("number of tissues") + 
    labs(fill = myLabel) + 
    theme(
      legend.position = "bottom", 
      text = element_text(size = 15)
    ) 
}

generateMutationPieChart <-  function(gene, mutation, sampleClasses, classLabel) {
  if (is.null(gene) | is.null(mutation) | is.null(sampleClasses)) return(NULL)
  
  assignment <- stackClasses(sampleClasses, classLabel, return_factor = TRUE)
  data <- assignment %>% 
    inner_join(mutation, by = "tissuename") %>% 
    filter(!is.na(aamutation))
  data_sum <- data %>% 
    group_by(class, aamutated) %>% 
    summarise(n = dplyr::n(), .groups = "drop_last") %>% 
    mutate(percent = n/sum(n))
  
  ggplot(data_sum, aes(x = 1, y = percent, fill = aamutated)) + 
    geom_bar(stat = "identity") + 
    facet_grid(facets = . ~ class) + 
    coord_polar(theta = "y") + 
    ggtitle(paste("\n", gene$symbol, "-", gene$ensg)) + 
    xlab("") + 
    # scale_color_viridis() + 
    theme(
      axis.text = element_blank(), 
      axis.ticks = element_blank(), 
      text = element_text(size = 15)
    )
}

# Mutational burden plots -----------------------------------------------------
#' Get Tissue Data Mutational Burden for generateMutationalFractionPlot
#'
#' @param ca tissue class assigment
#' @param anno tissue annotation
#'
#' @return data.frame with `tissuename` and `value` columns
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment(high = 0.5)
#' anno <- getTissueAnnotation()
#' df <- getTissueDataMutationalBurden(ca, anno)
#' 
getTissueDataMutationalBurden <- function(ca, anno) {
  result <- getOtherPropData(ca, anno)
  result %>% select(tissuename, value = mutational_fraction_in_percent) %>%
    addClassAssigmentAttribute(ca)
}

#' Generate mutational fraction plot by type
#' 
#' Produces mutational fraction plot that shows class comparison. Available types: 
#' "point" (default), "roc", "violin", "box", "coverage" and "all". Coverage plot may
#' be used to show the data availability in each class.
#' 
#' @param df data.frame, output of getTissueDataMutationalBurden()
#' @param plotType character string, plot type
#' @param title character string, the plot title
#' @param ca classAssignment object if it cannot be inferred from \code{df} param.
#'        However if df is an output of getTissueDataMutationalBurden, then this
#'        parameter does not need to be specified.
#' @return ggplot object
#' @export
#' @examples
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment(high = 0.5)
#' anno <- getTissueAnnotation()
#' df <- getTissueDataMutationalBurden(ca, anno)
#' label <- "Mutational fraction"
#' generateMutationalFractionPlot(df, "point", paste(label, "point plot"))
#' generateMutationalFractionPlot(df, "roc", paste(label, "roc plot"))
#' generateMutationalFractionPlot(df, "violin", paste(label, "violin plot"))
#' generateMutationalFractionPlot(df, "box", paste(label, "box plot"))
#' generateMutationalFractionPlot(df, "coverage", paste(label, "data coverage"))
#'
generateMutationalFractionPlot <- function(df, plotType, title = NULL, ca = NULL){
  generatePlotByType(
    data = df,
    ca = ca,
    plotType = plotType,
    dataCol = "value",
    title = title,
    ylabel = "value"
  )
}

# Copy number plots -----------------------------------------------------------
#' Generate waterfall plot for copy number
#' 
#' Produces a waterfall plot with copy number data
#' 
#' @param data data.frame, output of getWaterfallDataCopyNumber
#' @param fill character string, variable to use for bars coloring
#' @param cutoffAmplification numeric
#' @param cutoffDeletion numeric
#' @return ggplot object
#' @export
#' @examples
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter("tumor", c("thymoma", "mesothelioma"))
#' df <- getWaterfallDataCopyNumber("ENSG00000133703", prefilter)
#' generateCopynumberWaterfallPlot(df)
#'
generateCopynumberWaterfallPlot <- function(data, fill = "tumortype", 
                                            cutoffAmplification = 3, cutoffDeletion = 1.5) {
  if (is.null(data)) return()
  
  trans <- if (max(data$relativecopynumber, na.rm = TRUE) > 30){
    "sqrt"
  } else {
    "identity"
  }
  
  generateWaterfallPlot(
    data = data,
    dataCol = "relativecopynumber",
    ylabel = "relative copy number",
    trans = trans,
    fill = fill
  ) + 
    geom_hline(yintercept = cutoffAmplification, linetype = "dotted") + 
    geom_hline(yintercept = cutoffDeletion, linetype = "dotted") + 
    geom_hline(yintercept = 2, linetype = "solid")
}

generateCopynumberBarPlot <- function(data, cutoffAmplification = 3, cutoffDeletion = 1.5){
  if (is.null(data)) return()
  
  data <- data %>%
    mutate(
      x_score = factor(
        x = ifelse(
          test = relativecopynumber > cutoffAmplification, 
          yes = "amplification", 
          no = ifelse(
            test = relativecopynumber < cutoffDeletion, 
            yes = "deletion", 
            no = "normal"
          )
        ),
        levels = c("amplification", "normal", "deletion")
      )
    )
  
  generateScoreBarPlot(data, "copy number class")
}

#' Generate copy number plot by type
#' 
#' Produces copy number plot that shows class comparison. Available types: 
#' "point" (default), "roc", "violin", "box", "coverage" and "all". Coverage plot may
#' be used to show the data availability in each class.
#' 
#' @param df data.frame, output of getTissueDataCopyNumberById()
#' @param plotType character string, plot type
#' @param title character string, the plot title
#' @param ca classAssignment object if it cannot be inferred from \code{df} param.
#'        However if df is an output of getTissueDataCopyNumberById, then this
#'        parameter does not need to be specified.
#' @return ggplot object
#' @export
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' ensg <- "ENSG00000133703"
#' df <- getTissueDataCopyNumberById(ensg, ca)
#' symbol <- XIFF::getGeneSymbol(ensg)
#' generateCopynumberPlot(df, "point", paste(symbol, "point plot"))
#' generateCopynumberPlot(df, "roc", paste(symbol, "roc plot"))
#' generateCopynumberPlot(df, "violin", paste(symbol, "violin plot"))
#' generateCopynumberPlot(df, "box", paste(symbol, "box plot"))
#' generateCopynumberPlot(df, "coverage", paste(symbol, "data coverage"))
#'
generateCopynumberPlot <- function(df, plotType, title = NULL, ca = NULL){
  generatePlotByType(
    data = df,
    ca = ca,
    plotType = plotType,
    dataCol = "relativecopynumber",
    title = title,
    ylabel = "relative copy number",
    trans = "sqrt"
  )
}

#' Generate CN vs gene expression plot
#' 
#' Produces comparison plot with copy number and gene expression data.
#' 
#' @param dfCN data.frame, output of getTissueDataCopyNumberById()
#' @param dfExpr data.frame, output of getTissueDataGeneExpressionById()
#' @param title character string, the plot title
#' @param ca classAssignment object if it cannot be inferred from \code{dfCN} param.
#'        However if dfCN is an output of getTissueDataCopyNumberById, then this
#'        parameter does not need to be specified.
#' @return ggplot object
#' @export
#' @examples
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' ensg <- "ENSG00000133703"
#' dfCN <- getTissueDataCopyNumberById(ensg, ca)
#' dfExpr <- getTissueDataGeneExpressionById(ensg, ca)
#' generateCopynumberExpressionPlot(dfCN, dfExpr)
#'
generateCopynumberExpressionPlot <- function(dfCN, dfExpr, title = NULL, ca = NULL) {
  
  ca <- getClassAssigmentAttributeIfNull(ca, dfCN)
  
  assignment <- getAssignmentDf(ca)
  data <- dfCN %>% 
    inner_join(dfExpr, by = "tissuename") %>% 
    inner_join(assignment, by = "tissuename")
  
  mapping <- tooltipAes(x = tpm, y = relativecopynumber, color = class, plotFunc = geom_point)
  
  ggplot(
    data = data, 
    mapping = mapping
  ) + 
    geom_point() + 
    scale_color_manual(values = plotColors) + 
    scale_x_log10() +
    scale_y_sqrt() + 
    theme(
      text = element_text(size = 16), 
      legend.position = "none", 
      plot.title = element_text(hjust = 0.5)
    ) + 
    ggtitle(title) + 
    xlab("expr [TPM]") + 
    ylab("relative copy number")
}

# Expression plots ------------------------------------------------------------
#' Generate waterfall plot for gene expression
#' 
#' Produces a waterfall plot with gene expression data
#' 
#' @param data data.frame, output of getWaterfallDataGeneExpression
#' @param fill character string, variable to use for bars coloring
#' @return ggplot object
#' @export
#' @examples
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter("tumor", c("thymoma", "mesothelioma"))
#' df <- getWaterfallDataGeneExpression("ENSG00000133703", prefilter)
#' generateExpressionWaterfallPlot(df)
#'
generateExpressionWaterfallPlot <- function(data, fill = "tumortype") {
  if (is.null(data)) return()
  
  generateWaterfallPlot(
    data = data,
    dataCol = "tpm",
    ylabel = "TPM",
    trans = "log10",
    fill = fill
  )
}

#' Generate gene expression plot by type
#' 
#' Produces expression plot that shows class comparison. Available types: 
#' "point" (default), "roc", "violin", "box", "coverage" and "all". Coverage plot may
#' be used to show the data availability in each class.
#' 
#' @param df data.frame, output of getTissueDataGeneExpressionById()
#' @param plotType character string, plot type
#' @param title character string, the plot title
#' @param ca classAssignment object if it cannot be inferred from \code{df} param.
#'        However if df is an output of getTissueDataGeneExpressionById, then this
#'        parameter does not need to be specified.
#' @return ggplot object
#' @export
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' ensg <- "ENSG00000133703"
#' df <- getTissueDataGeneExpressionById(ensg, ca)
#' symbol <- XIFF::getGeneSymbol(ensg)
#' generateExpressionPlot(df, "point", paste(symbol, "point plot"))
#' generateExpressionPlot(df, "roc", paste(symbol, "roc plot"))
#' generateExpressionPlot(df, "violin", paste(symbol, "violin plot"))
#' generateExpressionPlot(df, "box", paste(symbol, "box plot"))
#' generateExpressionPlot(df, "coverage", paste(symbol, "data coverage"))
#' generateExpressionPlot(df, "all")
#'
generateExpressionPlot <- function(df, plotType = "point", title = NULL, ca = NULL){
  generatePlotByType(
    data = df,
    ca = ca,
    plotType = plotType,
    dataCol = "tpm",
    title = title,
    ylabel = "TPM",
    trans = "log10"
  )
}

#' Generate Expression Heatmap
#'
#' @param genes character vectors of genes
#' @param sampleAnno tissue annotation data.frame
#' @param heat_scale 
#' @param flipAxes logical. Should the axes be flipped.
#' @param sampleClasses classAsigment object. In most cases it should be inferred from 
#'        \code{data}. However if it is not possible, you may need to provide
#'        it on your own.
#'
#' @return pheatmap object
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' 
#' ca <- exampleClassAssigment()
#' gene_set <- getGSEAdata("human", gene_set = "HALLMARK_P53_PATHWAY")
#' cl_anno <- getTissueAnnotation()
#' 
#' annoVariables <- c("tumortype")
#' 
#' plot.new()
#' generateExpressionHeatmap(
#'  genes = gene_set, 
#'  sampleClasses = ca, 
#'  sampleAnno = cl_anno,
#'  annoVariables = annoVariables
#' )
#' 
#' plot.new()
#' generateExpressionHeatmap(
#'  genes = gene_set, 
#'  sampleClasses = ca, 
#'  sampleAnno = cl_anno,
#'  annoVariables = annoVariables,
#'  sampleLabel = "tumortype"
#' )
#' 
#' plot.new()
#' generateExpressionHeatmap(
#'  genes = gene_set, 
#'  sampleClasses = ca, 
#'  sampleAnno = cl_anno,
#'  annoVariables = c("tumortype", "organ"),
#'  sampleLabel = "tumortype"
#' )
#'
generateExpressionHeatmap <- function(genes, sampleClasses, sampleAnno, sampleLabel = "tissuename", annoVariables = NULL, heat_scale = "genes", flipAxes = FALSE, classLabel = NULL) {
  if (is.null(sampleClasses) || is.null(genes)) return()
  if(is.null(classLabel)) {
    classLabel <- getClassLabel(sampleClasses)
  }
  
  assignment <- stackClasses(sampleClasses, classLabel, return_factor = TRUE)
  anno.color <- plotColors[1:nlevels(assignment$class)]
  names(anno.color) = levels(assignment$class)
  
  sql <- paste0("SELECT tissuename, coalesce(symbol, pr.ensg) as gene, log2tpm FROM tissue.processedrnaseqview pr ",
                "JOIN gene g ON g.ensg = pr.ensg ",
                "WHERE pr.ensg IN ('", paste(genes, collapse = "','"), "') ",
                "AND tissuename = ANY(ARRAY['", paste(assignment$tissuename, collapse = "','"), "'])")
  data <- getPostgresql(sql)
  
  norange_gene <- data %>% group_by(gene) %>% summarise(r = max(log2tpm) - min(log2tpm), .groups = "drop") %>% filter(r == 0) %>% .$gene
  norange_ti <- data %>% group_by(tissuename) %>% summarise(r = max(log2tpm) - min(log2tpm), .groups = "drop") %>% filter(r == 0) %>% .$tissuename
  if (length(norange_gene) > 0) data <- data %>% filter(!gene %in% norange_gene)
  if (length(norange_ti) > 0) data <- data %>% filter(!tissuename %in% norange_ti)
  
  data.cross <- tapply(data$log2tpm, list(data$tissuename, data$gene), mean)
  
  class_annotation <- assignment[, "class", drop = FALSE]
  
  if("tissuename" %in% colnames(sampleAnno)) {
    rownames(sampleAnno) <- sampleAnno[["tissuename"]]
  }
  
  if(length(annoVariables) > 0) {
    class_annotation <- bind_cols(
      class_annotation,
      sampleAnno[rownames(class_annotation),annoVariables, drop = FALSE]
    )
    class_annotation <- class_annotation %>% mutate_if(is.logical, as.character)
    
    allNas <- vapply(class_annotation, function(x) all(is.na(x)), FUN.VALUE = TRUE)
    class_annotation <- class_annotation[,!allNas]
  }
  
  
  if (flipAxes) {
    s <- ifelse(heat_scale == "none", "none", ifelse(heat_scale == "tissues", "column", "row"))
    g <- getPheatmap(t(data.cross), scale = s, labels_col = sampleAnno[rownames(data.cross), sampleLabel], 
                     annotation_col = class_annotation, annotation_color = list(class = anno.color))
  } else {
    s <- ifelse(heat_scale == "none", "none", ifelse(heat_scale == "genes", "column", "row"))
    g <- getPheatmap(data.cross, scale = s, labels_row = sampleAnno[rownames(data.cross), sampleLabel], 
                     annotation_row = class_annotation, annotation_color = list(class = anno.color))
  }
  
  g
}

# Protein plots ---------------------------------------------------------------

#' Generate Protein Waterfall Plot
#'
#' @param data data.frame created by \code{\link{getWaterfallDataRppa}}
#' @param fill column name from data to be used for filling the bars
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#'
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter("tumor", c("thymoma", "mesothelioma"))
#' data <- getWaterfallDataRppa("CDK1", prefilter)
#' generateProteinWaterfallPlot(data)
#' 
generateProteinWaterfallPlot <- function(data, fill = "tumortype") {
  if (is.null(data)) return()
  
  generateWaterfallPlot(
    data = data,
    dataCol = "score",
    ylabel = "RPPA score",
    fill = fill
  )
}

#' Generate protein expression plot by type
#' 
#' Produces protein expression plot that shows class comparison. Available types: 
#' "point" (default), "roc", "violin", "box", "coverage", "all". Coverage plot may
#' be used to show the data availability in each class.
#' 
#' @param df data.frame, output of getTissueDataRppaById()
#' @param plotType character string, plot type
#' @param title character string, the plot title
#' @param ca classAssignment object if it cannot be inferred from \code{df} param.
#'        However if df is an output of getTissueDataRppaById, then this
#'        parameter does not need to be specified.
#' @return ggplot object
#' @export
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment(0.15,0.2)
#' antibody <- "MEK1_pS217_S221"
#' df <- getTissueDataRppaById(antibody, ca)
#' generateProteinExpressionPlot(df, "point", paste(antibody, "point plot"))
#' generateProteinExpressionPlot(df, "roc", paste(antibody, "roc plot"))
#' generateProteinExpressionPlot(df, "violin", paste(antibody, "violin plot"))
#' generateProteinExpressionPlot(df, "box",  paste(antibody, "box plot"))
#' generateProteinExpressionPlot(df, "coverage", paste(antibody, "data coverage"))
#' generateProteinExpressionPlot(df, "all")
#' 
generateProteinExpressionPlot <- function(df, plotType, title = NULL, ca = NULL) {
  generatePlotByType(
    data = df,
    ca = ca,
    plotType = plotType,
    dataCol = "score",
    title = title,
    ylabel = "RPPA score"
  )
}

# Immune cells plots ----------------------------------------------------------

#' Generate Immune Cell Waterfall Plot
#'
#' @param data data.frame create with \code{getWaterfallDataImmuneCells}
#' @param fill column used to fill the data on the plot
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' 
#' df <- getWaterfallDataImmuneCells("Activated dendritic cells", makePrefilter())
#' generateImmuneCellWaterfallPlot(df)
#' 
generateImmuneCellWaterfallPlot <- function(data, fill = "tumortype") {
  if (is.null(data)) return()
  
  generateWaterfallPlot(
    data = data,
    dataCol = "score",
    ylabel = "score",
    fill = fill
  )
}


#' Generate Immune Cells Plot
#'
#' @param df data.frame created with \code{\link{getTissueDataImmuneCellsById}}
#' @param plotType one of the follwoing "point", "roc", "violin", "box", "coverage"
#' or "all" to get a list with all plots
#' @param title string with the title to be added to the plot
#' @param ca classAssigment object (in most cases it will be inherited from df)
#' however in some cases it needs to be provided explicitly.
#'
#' @return plot or a list of plots if plotType is 'all'
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' data <- getTissueDataImmuneCellsById("Neurons", exampleClassAssigment(high = 0.16))
#' 
#' type <- c("point", "roc", "violin", "box", "coverage")
#' generateImmuneCellsPlot(data, "all")
#' 
generateImmuneCellsPlot <- function(df, plotType, title = NULL, ca = NULL){
  
  ca <- getClassAssigmentAttributeIfNull(ca, data)
  
  generatePlotByType(
    data = df,
    ca = ca,
    plotType = plotType,
    dataCol = "score",
    title = title,
    ylabel = "score"
  )
}

# Metabolics plots ------------------------------------------------------------

#' Generate Metabolics Waterfall Plot
#'
#' @param data data.frame created from \code{\link{getWaterfallDataMetabolics}}
#' @param fill column from data used to fill value
#'
#' @return ggplot
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter("tumor", c("thymoma"))
#' data <- getWaterfallDataMetabolics("gluconeogenesis", prefilter)
#' 
#' generateMetabolicsWaterfallPlot(data)
#' 
generateMetabolicsWaterfallPlot <- function(data, fill = "tumortype") {
  if (is.null(data)) return()
  
  generateWaterfallPlot(
    data = data,
    dataCol = "score",
    ylabel = "score",
    fill = fill
  )
}

#' Generate Metabolics Plot
#'
#' @param df data.frame created from \code{\link{getTissueDataMetabolicsById}}
#' @param plotType one of the follwoing "point", "roc", "violin", "box", "coverage"
#' or "all" to get a list with all plots
#' @param title character string, the plot title
#' @param ca classAssignment object if it cannot be inferred from \code{df} param.
#' However if df is an output of getTissueDataMetabolicsById, then this
#' parameter does not need to be specified.
#'
#' @return ggplot
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment(high = 0.16)
#' data <- getTissueDataMetabolicsById("gluconeogenesis", ca)
#' 
#' generateMetabolicsPlot(data, "all")
#' 
generateMetabolicsPlot <- function(df, plotType, title = NULL, ca = NULL){
  
  ca <- getClassAssigmentAttributeIfNull(ca, df)
  
  generatePlotByType(
    data = df,
    ca = ca,
    plotType = plotType,
    dataCol = "score",
    title = title,
    ylabel = "score"
  )
}

# Signaling plots -------------------------------------------------------------



#' Generate Signaling Input Barplot
#'
#' @param data data.frame created with \code{\link{getInputDataSignaling}}
#' @param fill column from data used to fill value
#'
#' @return ggplot2
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' 
#' selected <- c("cholangiocarcinoma", "lymphoid neoplasm diffuse large b-cell lymphoma", 
#'   "uterine carcinosarcoma", "kidney chromophobe", "adrenocortical carcinoma", 
#'   "uveal melanoma")
#' 
#' prefilter <- makePrefilter(selected = "uveal melanoma")
#' pathways <- getAvailableSignalingPathways()
#' data <- getInputDataSignaling(pathways[1], prefilter)
#' generateSignalingInputBarplot(data)
#' 
generateSignalingInputBarplot <- function(data, fill = "tumortype") {
  if (is.null(data)) return()
  
  pathway <- colnames(data)[2]
  pathway <- sym(pathway)
  fill <- sym(fill)
  
  ggplot(
    data = data, 
    mapping = aes(x = !!pathway, fill = !!fill)
  ) + 
    geom_bar() + 
    theme(
      legend.position = "bottom", 
      text = element_text(size = 16), 
      axis.text = element_text(size = 16)
    ) + 
    xlab(paste(toupper(pathway), "signaling")) + 
    ylab("# of tissues")
}

generateSignalingPlot <- function(df, ca, plotType, title = NULL){
  if (plotType == "coverage"){
    generateDataCoveragePlot(
      data = df, 
      col = "active", 
      ca = ca
    )
  } else {
    assignment <- getAssignmentDf(ca, useLabels = FALSE)
    data <- assignment %>%
      left_join(df, by = "tissuename") %>%
      select(state = active, class) %>% 
      filter(!is.na(state))
    
    data_sum <- data %>% 
      group_by(class) %>% 
      summarise(n = dplyr::n(), .groups = "drop")
    
    labels <- paste0(
      c(classIdToLabel("class1", ca), classIdToLabel("class2", ca)), 
      " (", data_sum$n, ")"
    )
    data <- data %>%
      inner_join(data_sum, by = "class") %>%
      mutate(label = factor(
        x = ifelse(class == "class1", labels[1], labels[2]),
        levels = labels
      ))
    
    if (plotType == "bar") {
      generateSignalingBarplot(data, title)
    } else if (plotType == "pie"){
      generateSignalingPieChart(data, title)
    }
  }
}

generateSignalingBarplot <- function(data, title = NULL){
  ggplot(
    data = data, 
    mapping = aes(x = label, fill = state)) + 
    geom_bar() + 
    theme(
      legend.position = "bottom",
      text = element_text(size = 16), 
      axis.text = element_text(size = 16)
    ) + 
    ggtitle(title) + 
    xlab("") + 
    ylab("# of tissues") + 
    coord_flip()
}

generateSignalingPieChart <- function(data, title = NULL){
  data_sum2 <- data %>% 
    group_by(label, state) %>% 
    summarise(n = dplyr::n(), .groups = "drop_last") %>% 
    mutate(percent = n/sum(n))
  
  ggplot(
    data = data_sum2, 
    mapping = aes(x = 1, y = percent, fill = state)
  ) + 
    geom_bar(stat = "identity") + 
    facet_grid(facets = . ~ label) + 
    coord_polar(theta = "y") + 
    ggtitle(title) + 
    xlab("") + 
    theme(
      axis.text = element_blank(), 
      axis.ticks = element_blank(), 
      text = element_text(size = 15)
    )
}

# Hallmark gene sets ----------------------------------------------------------
#' Generate waterfall plot for gene set
#' 
#' Produces a waterfall plot with gene set data
#' 
#' @param data data.frame, output of getWaterfallDataHallmarkSet
#' @param fill character string, variable to use for bars coloring
#' @param ylabel character string, y-axis label
#' @return ggplot object
#' @export
#' @examples
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter("tumor", c("thymoma", "mesothelioma"))
#' df <- getWaterfallDataHallmarkSet("HALLMARK_APOPTOSIS", prefilter)
#' generateHallmarkSetWaterfallPlot(df)
#'
generateHallmarkSetWaterfallPlot <- function(data, fill = "tumortype", ylabel = "gsva"){
  if (is.null(data)) return()
  
  generateWaterfallPlot(
    data = data,
    dataCol = "score",
    ylabel = ylabel,
    fill = fill
  )
}

#' Generate gene set plot by type
#' 
#' Produces gene set plot that shows class comparison. Available types: 
#' "point" (default), "roc", "violin", "box", "coverage" and "all". Coverage plot may
#' be used to show the data availability in each class.
#' 
#' @param df data.frame, output of getTissueDataGeneSetById()
#' @param plotType character string, plot type
#' @param title character string, the plot title
#' @param ca classAssignment object if it cannot be inferred from \code{df} param.
#'        However if df is an output of getTissueDataGeneSetById, then this
#'        parameter does not need to be specified.
#' @param ylabel character string, y-axis label
#' @return ggplot object
#' @export
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' geneSet <- "HALLMARK_APOPTOSIS"
#' df <- getTissueDataGeneSetById(geneSet, ca)
#' generateGeneSetPlot(df, "point", paste(geneSet, "point plot"))
#' generateGeneSetPlot(df, "roc", paste(geneSet, "roc plot"))
#' generateGeneSetPlot(df, "violin", paste(geneSet, "violin plot"))
#' generateGeneSetPlot(df, "box", paste(geneSet, "box plot"))
#' generateGeneSetPlot(df, "coverage", paste(geneSet, "data coverage"))
#' 
generateGeneSetPlot <- function(df, plotType = "point", title = NULL, ca = NULL, ylabel = "gsva") {
  generatePlotByType(
    data = df,
    ca = ca,
    plotType = plotType,
    dataCol = "score",
    title = title,
    ylabel = ylabel
  )
}

# Clones plots ----------------------------------------------------------------

#' Generate Clones Waterfall Plot
#'
#' @param data data.frame created with \code{getWaterfallDataClones}
#' @param fill column from data used to fill value
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' 
#' prefilter <- makePrefilter()
#' df <- getWaterfallDataClones("number_of_clones", prefilter)
#' generateClonesWaterfallPlot(df)
#' 
generateClonesWaterfallPlot <- function(data, fill = "tumortype") {
  if (is.null(data)) return()
  
  score <- colnames(data)[2]
  generateWaterfallPlot(
    data = data,
    dataCol = score,
    ylabel = score,
    fill = fill
  )
}


#' Generate Clones Plot
#'
#' @param df result of getTissueDataClones with score column
#' @param plotType plot type: point, roc, violin, box, coverage or all.
#' @param title plot title
#' @param ca classAssigment object (in most cases it can be
#' inherited from df)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' 
#' ca <- exampleClassAssigment(high = 0.16)
#' df <- getTissueDataClones(ca, "number_of_clones") %>%
#'  mutate(score = `number_of_clones`)
#' generateClonesPlot(df, "all")
#' 
generateClonesPlot <- function(df, plotType, title = NULL, ca = NULL){
  
  ca <- getClassAssigmentAttributeIfNull(ca, df)
  
  generatePlotByType(
    data = df,
    ca = ca,
    plotType = plotType,
    dataCol = "score",
    title = title,
    ylabel = "score"
  )
}

# Survival plot ---------------------------------------------------------------

#' Generate Kaplan-Meier Plot
#'
#' @param data tiff survival analysis created with \code{\link{survivalAnalysis}}
#' @param plot_p if TRUE (default), then the p-value is added to the plot
#' @param plot_ci if TRUE (defult FALSE), then the confidence intervals are added to the plot
#' @param plot_risk if TRUE (default), then the risk table is added to be plot
#'
#' @return ggsurvplot
#' @export
#'
#' @examples
#' 
#' setDbOptions(getSettings())
#' 
#' data <- getPatientAnnotation()
#' ca <- exampleClassAssigment(high = 0.16)
#' ca  <- classAssignment2df(ca)
#' dataSurv <- inner_join(data, ca)
#' 
#' surv <- survivalAnalysis(dataSurv)
#' generateKaplanMeierPlot(surv)
#' generateKaplanMeierPlot(surv, plot_ci = TRUE)
#' 
generateKaplanMeierPlot <- function(data, plot_p = TRUE, plot_ci = FALSE, plot_risk = TRUE) {
  fit <- data$fit
  if (is.null(fit$strata)) return()
  
  names(fit$strata) <- gsub("class=", "", names(fit$strata))
  survminer::ggsurvplot(
    fit = fit, 
    data = data$d, 
    risk.table = plot_risk, 
    pval = plot_p, 
    conf.int = plot_ci, 
    fontsize = 8, 
    font.x = 18, 
    font.y = 18, 
    font.tickslab = 18, 
    pval.size = 8, 
    font.legend = 18, 
    censor.size = 6
  )
}


#' Grid draw method for ggsurvplot
#'
#' @param x ggsurvplot object
#'
#' @return nothing, called for side effects.
#' 
#' @details required to properly draw ggsurvplot in shiny application.
#' @importFrom grid grid.draw
#' 
#' @export
#' 
grid.draw.ggsurvplot <- function(x){
  print(x, newpage = FALSE)
}

# Gene signatures plots -------------------------------------------------------
#' Generate waterfall plot for gene signatures
#' 
#' Produces a waterfall plot with gene signature data
#' 
#' @param data data.frame, output of getWaterfallDataGeneSignatures
#' @param fill character string, variable to use for bars coloring
#' @return ggplot object
#' @export
#' @examples
#' 
#' setDbOptions(getSettings())
#' prefilter <- makePrefilter("tumor", c("thymoma", "mesothelioma"))
#' df <- getWaterfallDataGeneSignatures("MERCK18", prefilter)
#' generateGeneSignatureWaterfallPlot(df)
#'
generateGeneSignatureWaterfallPlot <- function(data, fill = "tumortype"){
  if (is.null(data)) return()
  
  generateWaterfallPlot(
    data = data,
    dataCol = "score",
    ylabel = "score",
    fill = fill
  )
}

#' Generate gene signature plot by type
#' 
#' Produces gene signature plot that shows class comparison. Available types: 
#' "point" (default), "roc", "violin", "box", "coverage" and "all". Coverage plot may
#' be used to show the data availability in each class.
#' 
#' @param df data.frame, output of getTissueDataGeneSignatureById()
#' @param plotType character string, plot type
#' @param title character string, the plot title
#' @param ca classAssignment object if it cannot be inferred from \code{df} param.
#'        However if df is an output of getTissueDataGeneSignatureById, then this
#'        parameter does not need to be specified.
#' @return ggplot object
#' @export
#' @examples
#' 
#' setDbOptions(getSettings())
#' ca <- exampleClassAssigment()
#' signature <- "MERCK18"
#' df <- getTissueDataGeneSignatureById(signature, ca)
#' generateGeneSignaturePlot(df, "point", paste(signature, "point plot"))
#' generateGeneSignaturePlot(df, "roc", paste(signature, "roc plot"))
#' generateGeneSignaturePlot(df, "violin", paste(signature, "violin plot"))
#' generateGeneSignaturePlot(df, "box", paste(signature, "box plot"))
#' generateGeneSignaturePlot(df, "coverage", paste(signature, "data coverage"))
#'
generateGeneSignaturePlot <- function(df, plotType, title = NULL, ca = NULL) {
  generatePlotByType(
    data = df,
    ca = ca,
    plotType = plotType,
    dataCol = "score",
    title = title,
    ylabel = "Gene signature score"
  )
}

# Additional properties plot --------------------------------------------------
#' Generate Other Properties Plot
#'
#' @param sampleClasses classAssigment object
#' @param sAnno tissue annotation (see: \code{getTissueAnnotation})
#' @param property which property to be plotted
#' @param classLabel if sampleClasses is a list instead of classAssigment
#' then the classLabels needs to be provided.
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' setDbOptions(getSettings())
#' sampleClasses <- exampleClassAssigment()
#' sAnno <- getTissueAnnotation()
#' generateOtherPropertiesPlot(sampleClasses, sAnno, "lossofy")
#' generateOtherPropertiesPlot(sampleClasses, sAnno, "organ")
#'
generateOtherPropertiesPlot <- function(sampleClasses, sAnno, property, classLabel = NULL) {
  
  if (is.null(sampleClasses)) return()
  
  if(is.null(classLabel)) {
    classLabel <- XIFF::getClassLabel(sampleClasses)  
  }
  assignment <- stackClasses(sampleClasses, classLabel, return_factor = TRUE)
  sAnno <- sAnno %>% inner_join(assignment, by = "tissuename") %>% droplevels()
  
  if(is.character(property)) {
    property <- list(
      property =  property,
      data = ifelse(is.numeric(sAnno[[property]]), "numeric", "class")
    )
  }
  
  if (property$data == "class") {
    font_size <- if (length(levels(sAnno[, property$property])) > 20) 16 else 20
    
    ggplot(sAnno, aes_string(
      x = paste0("`", property$property, "`"), 
      fill = paste0("`", property$property, "`")
    )) + 
      geom_bar(na.rm = TRUE) + 
      theme(
        legend.position = "none", 
        text = element_text(size = font_size)
      ) + 
      xlab(property$property) + 
      ylab("# of tissues") + 
      coord_flip() + 
      facet_wrap(~class, scales = "free_x")
    
  } else if (property$data == "numeric") {
    trans <- if (property$property == "Mutational fraction") {
      "log10"
    } else {
      "identity"
    }
    
    ggplot(sAnno, aes_string(
      x = "class", 
      y = paste0("`", property$property, "`"), 
      fill = "class"
    )) + 
      geom_boxplot(na.rm = TRUE) + 
      scale_y_continuous(trans = trans) + 
      theme(
        text = element_text(size = 20), 
        legend.position = "none"
      ) + 
      coord_flip()
  }
}
