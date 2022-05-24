tissueOverviewTabUI <- function(id){
  ns <- NS(id)
  
  tabsetPanel(
    id = ns("tabset"), 
    columnTabPanel(
      title = "description",
      value = "description",
      tissueOverviewDescriptionTabUI_sidebar(ns("description")),
      tissueOverviewDescriptionTabUI_main(ns("description"))
    ),
    columnTabPanel(
      title = "gene expression: varying genes",
      value = "varying_genes",
      tissueOverviewGeneExpressionTabUI_sidebar(ns("varying_genes"), mode = "varying_genes"),
      tissueOverviewGeneExpressionTabUI_main(ns("varying_genes"))
    ),
    columnTabPanel(
      title = "gene expression: gene set",
      value = "gene_set",
      tissueOverviewGeneExpressionTabUI_sidebar(ns("gene_set"), mode = "gene_set"),
      tissueOverviewGeneExpressionTabUI_main(ns("gene_set"))
    ),
    columnTabPanel(
      title = "hallmark gene set heatmap",
      value = "hallmark",
      tissueOverviewHallmarkTabUI_sidebar(ns("hallmark")),
      tissueOverviewHallmarkTabUI_main(ns("hallmark"))
    )
  )
}

tissueOverviewTab <- function(input, output, session, fm, classSelection, classLabel, classStack, 
                              TissueAnnotation, TissueAnnotationFocus, gsea_data_hallmark){
  # Description subtab --------------------------------------------------------
  callModule(
    module = tissueOverviewDescriptionTab,
    id = "description",
    classSelection = classSelection, 
    classLabel = classLabel, 
    TissueAnnotation = TissueAnnotation, 
    TissueAnnotationFocus = TissueAnnotationFocus
  )
  
  # Gene expression varying genes subtab --------------------------------------
  callModule(
    module = tissueOverviewGeneExpressionTab,
    id = "varying_genes",
    fm = fm,
    classSelection = classSelection, 
    classLabel = classLabel, 
    classStack = classStack,
    mode = "varying_genes",
    label = "Varying genes",
    gsea_data_hallmark = gsea_data_hallmark, 
    TissueAnnotation = TissueAnnotation,
    TissueAnnotationFocus = TissueAnnotationFocus
  )
  
  # Gene expression gene set subtab -------------------------------------------
  callModule(
    module = tissueOverviewGeneExpressionTab,
    id = "gene_set",
    fm = fm,
    classSelection = classSelection, 
    classLabel = classLabel, 
    classStack = classStack,
    mode = "gene_set",
    label = "Gene set",
    gsea_data_hallmark = gsea_data_hallmark, 
    TissueAnnotation = TissueAnnotation,
    TissueAnnotationFocus = TissueAnnotationFocus
  )
  
  callModule(
    module = tissueOverviewHallmarkTab,
    id = "hallmark",
    classSelection = classSelection, 
    classLabel = classLabel,
    TissueAnnotationFocus = TissueAnnotationFocus,
    gsea_data_hallmark = gsea_data_hallmark
  )
}
