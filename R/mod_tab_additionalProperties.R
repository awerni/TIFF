additionalPropertiesTabUI <- function(id){
  ns <- NS(id)
  
  tabsetPanel(
    id = ns("tabset"),
    columnTabPanel(
      title = "class comparison",
      value = "class comparison",
      additionalPropertiesClassComparisonTabUI_sidebar(ns("comparison")),
      additionalPropertiesClassComparisonTabUI_main(ns("comparison"))
    ),
    columnTabPanel(
      title = "gene signatures",
      value = "gene signatures",
      additionalPropertiesGeneSignaturesTabUI_sidebar(ns("gene_signatures")),
      additionalPropertiesGeneSignaturesTabUI_main(ns("gene_signatures"))
    ),
    columnTabPanel(
      title = "hallmark gene set comparison",
      value = "hallmark gene set comparison",
      additionalPropertiesHallmarkComparisonTabUI_sidebar(ns("hallmark")),
      additionalPropertiesHallmarkComparisonTabUI_main(ns("hallmark"))
    )
  )
}

additionalPropertiesTab <- function(input, output, session, classSelection, classLabel, 
                                    TissueAnnotationFocus, Result_otherPropStatistic, 
                                    msigDBLink, geneSignatures){
  callModule(
    module = additionalPropertiesClassComparisonTab,
    id = "comparison",
    classSelection = classSelection, 
    classLabel = classLabel, 
    TissueAnnotationFocus = TissueAnnotationFocus, 
    Result_otherPropStatistic = Result_otherPropStatistic
  )
  
  callModule(
    module = additionalPropertiesGeneSignaturesTab,
    id = "gene_signatures",
    classSelection = classSelection, 
    classLabel = classLabel, 
    TissueAnnotationFocus = TissueAnnotationFocus,
    geneSignatures = geneSignatures
  )
  
  callModule(
    module = additionalPropertiesHallmarkComparisonTab,
    id = "hallmark",
    classSelection = classSelection, 
    classLabel = classLabel, 
    TissueAnnotationFocus = TissueAnnotationFocus,
    msigDBLink = msigDBLink
  )
}
