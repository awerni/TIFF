mutationTabUI <- function(id){
  ns <- NS(id)
  
  tabsetPanel(
    id = ns("tabset"),
    columnTabPanel(
      title = "per gene",
      value = "per_gene",
      mutationPerGeneTabUI_sidebar(ns("per_gene")),
      mutationPerGeneTabUI_main(ns("per_gene"))
    ),
    columnTabPanel(
      title = "mutational fraction",
      value = "mutational_fraction",
      mutationMutationalFractionTabUI_sidebar(ns("fraction")),
      mutationMutationalFractionTabUI_main(ns("fraction"))
    )
  )
}

mutationTab <- function(input, output, session, classSelection, classLabel, 
                        gene_anno, TissueAnnotationFocus, Result_otherPropStatistic, species){
  callModule(
    module = mutationPerGeneTab,
    id = "per_gene",
    classSelection = classSelection, 
    classLabel = classLabel, 
    gene_anno = gene_anno,
    species = species
  )
  
  callModule(
    module = mutationMutationalFractionTab,
    id = "fraction",
    classSelection = classSelection, 
    classLabel = classLabel, 
    TissueAnnotationFocus = TissueAnnotationFocus, 
    Result_otherPropStatistic = Result_otherPropStatistic
  )
}
