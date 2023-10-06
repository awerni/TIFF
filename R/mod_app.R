#' TIFF App UI
#'
#' @param settings settings object
#'
#' @importFrom bslib nav_item page_navbar
#' @importFrom shinythemes shinytheme
#' @importFrom XIFF columnTabPanel
#' @importFrom rlang `!!!`
#' @export
#' 
appUI <- function(settings) {
  
  appId <- "app"
  
  headTag <- tags$head(
    tags$title("TIFF (Tissue Differences)"),
    useXIFF(),
    shinyjs::useShinyjs(),
    tags$link(
      rel = "stylesheet",
      type = "text/css",
      href = withHash("www/main.css", "main.css")
    ),
    tags$script(src = withHash(filePath = "www/main.js", wwwPath = "main.js")),
    HTML(
      '<!-- Matomo -->
      <!-- End Matomo Code -->'
    )
  )
  
  theme <- XIFF::appTheme()
  
  div(
    headTag,
    bslib::page_navbar(
      id = "app-MainMenu",
      window_title = "Tissue Differences",
      theme = theme,
      collapsible = TRUE,
      bslib::nav_item(appUI_title(appId, "Tissue Differences", "TIFF.png")),
      appUI_main_input(appId),
      appUI_main_overview(appId),
      appUI_main_analysis(appId,),
      XIFF::appUI_main_machine_learning(appId),
      XIFF::appUI_main_about(
        appId,
        docuLink = settings[["links"]][["docuLink"]],
        aboutTabUIFunc = aboutTabUI_main
      ),
      bslib::nav_spacer(),
      !!!XIFF::appUI_title_right(
        id = appId,
        docuLink = settings[["links"]][["docuLink"]],
        dbName = settings[["db"]][["name"]],
        packageName = "TIFF"
      )
    )
  )
}

appUI_main_input <- function(id){
  
  ns <- NS(id)
  tabPanel(
    title = "Input",
    
    fluidRow(
      column_2(
        inputTabUI_selector(ns("input")),
        inputTabUI_sidebar(ns("input"))), 
      column_10(
        inputTabUI_main(ns("input"))
      )
    )
  )
  
}

appUI_main_overview <- function(id) {
  ns <- NS(id)
  tabPanel(title = "Overview", tissueOverviewTabUI(ns("overview")))
}

appUI_main_analysis <- function(id) {
  
  ns <- NS(id)
  tabset <- tabPanel(
    title = "Analysis",
    tabsetPanel(
      id = ns("AnalysisMenu"),
      tabPanel("mutation", mutationTabUI(ns("mutation"))),
      tabPanel("copy number", copyNumberTabUI(ns("copy_number"))),
      tabPanel("gene expression", expressionTabUI(ns("expression"))),
      tabPanel("protein expression", proteinTabUI(ns("protein"))),
      tabPanel("survival", survivalTabUI(ns("survival"))),
      tabPanel("derived data", derivedDataTabUI(ns("derived_data"))),
      tabPanel("additional properties", additionalPropertiesTabUI(ns("additional_props")))
    )
  )
  
  #tabset$children[[1]]$children[[2]]$attribs$style <- "margin-top: 0px;"
  tabset
}

#' TIFF App server function.
#'
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param fm future manager object
#' @param settings cliff settings
#' @param geneSignatures result of the \code{getAvailableGeneSignatures}
#' @param hallmarkGeneSets result of the \code{getAvailableHallmarkGeneSets}
#'
#' @details It should not be called directly by the user.
#' Please use \code{TIFF::run()}.
#'
#' @return shiny app server
#' @export
#'
app <- function(input, output, session, fm, settings, geneSignatures, hallmarkGeneSets){
  # Parse settings ------------------------------------------------------------
  msigDBLink <- settings[["links"]][["msigDBLink"]]
  species <- settings[["species"]]
  
  # Reactive values -----------------------------------------------------------
  classSelection <- reactiveValues(class1 = NULL, class2 = NULL)
  classLabel <- reactiveValues(class1_name = "class1", class2_name = "class2")
  classStack <- reactiveVal(tibble(tissuename = character(0)))
  
  # Common reactives ----------------------------------------------------------
  gene_anno <- reactive({
    getGeneAnno(species)
  })
  
  gsea_data_hallmark <- reactive({
    getGSEAdata(species, "hallmark")
  })
  
  TissueAnnotation <- reactive({
    factorCols <- c(
      "tumortype", "tumortype_adjacent", "organ", "tissue_subtype", 
      "microsatellite_stability_class", "immune_environment", 
      "gi_mol_subgroup", "icluster", "consmolsubtype", "dnasequenced", 
      "til_pattern", "autolysis_score", "vendorname",
      "histology_type", "histology_subtype", "lossofy"
    )
    
    ta <- getTissueAnnotation() %>%
      mutate_at(factorCols, as.factor) %>%
      rename(
        "Tissue subtype" = tissue_subtype, 
        "Microsatellite stability" = microsatellite_stability_class, 
        "Immune environment class" = immune_environment, 
        "GI molecular subgroup" = gi_mol_subgroup, 
        "DNA sequenced" =  dnasequenced, 
        "Tumor purity" = tumorpurity, 
        iCluster = icluster, 
        "Consensus Molecular Subtype" = consmolsubtype, 
        "TIL pattern" = til_pattern, 
        "Autolysis score" = autolysis_score, 
        "RNA integrity number" = rna_integrity_number, 
        "Minutes ischemia" = minutes_ischemia, 
        "Vendor name" = vendorname,
        "Histology type" = histology_type,
        "Histology subtype" = histology_subtype,
        "loss of Y chromosome" = lossofy
      ) %>%
      arrange(tumortype, organ)
    
    rownames(ta) <- ta$tissuename
    ta
  })
  
  PatientAnnotationFuller <- reactive({
    getPatientAnnotation()
  })
  
  TissuePatientAnnotation <- reactive({
    tissueAnno <- TissueAnnotation()
    req(tissueAnno)
    patientAnno <- PatientAnnotationFuller()
    req(patientAnno)
    
    patientAnno %>% 
      select(-tumortype, -tumortype_adjacent, -organ, -tissuepanel) %>%
      mutate_at(c("gender", "race", "ethnicity"), as.factor) %>%
      rename(
        "Vital status" = "vital_status", 
        "Days to last follow up" = "days_to_last_followup", 
        "Days to death" = "days_to_death"
      ) %>%
      left_join(tissueAnno, ., by = "tissuename")
  })
  
  TissuePatientAnnotationFull <- reactive({
    # Filter only selected tissues
    ta <- TissuePatientAnnotation() %>% rename(tumorpurity = `Tumor purity`)
    if (nrow(ta) == 0) return()
    # Add other props
    other_prop <- getOtherPropData(ta$tissuename, ta, p = session)
    if (nrow(other_prop) == 0) return()
    other_prop %>%
      rename(
        "RNA sequenced" = "rnasequenced", 
        "Copynumber available" = copynumbered, 
        "Mutational fraction" = mutational_fraction_in_percent,
        `Tumor purity` = tumorpurity,
        `Tumor purity class` = tumor_purity_class,
        `Mutational burden class` = mutational_burden_class
      )
     
  })
  
  TissuePatientAnnotationFocus <- reactive({
    # Filter only selected tissues
    ta <- TissuePatientAnnotationFull() %>% filter(tissuename %in% c(classSelection$class1, classSelection$class2))
    if (nrow(ta) == 0) return()
    ta
  })
  
  Result_otherPropStatistic <- reactive({
    anno <- TissuePatientAnnotationFocus()
    validate(
      need(nrow(anno) > 0, "no tissue data available..."),
      need(classSelection$class1, "goto input tab and select tissues for class1"),
      need(classSelection$class2, "goto input tab and select tissues for class2")
    )
    
    anno <- anno %>% select(-`DNA sequenced`, -`RNA sequenced`, -`Copynumber available`)
    
    cs <- reactiveValuesToList(classSelection)
    getOtherPropStatistics(cs, anno, p = session)
  })
  
  # Tabs ----------------------------------------------------------------------
  callModule(
    module = inputTab,
    id = "input",
    classSelection = classSelection, 
    classLabel = classLabel, 
    classStack = classStack, 
    fm = fm, 
    species = species,
    TissueAnnotation = TissuePatientAnnotationFull, 
    TissueAnnotationFocus = TissuePatientAnnotationFocus, 
    PatientAnnotationFuller = PatientAnnotationFuller,
    geneSignatures = geneSignatures
  )
  
  callModule(
    module = mutationTab,
    id = "mutation",
    classSelection = classSelection, 
    classLabel = classLabel, 
    gene_anno = gene_anno, 
    TissueAnnotationFocus = TissuePatientAnnotationFocus, 
    Result_otherPropStatistic = Result_otherPropStatistic,
    species = species
  )
  
  callModule(
    module = copyNumberTab,
    id = "copy_number",
    fm = fm, 
    classSelection = classSelection, 
    classLabel = classLabel, 
    gene_anno = gene_anno,
    species = species
  )
  
  callModule(
    module = expressionTab,
    id = "expression",
    fm = fm, 
    classSelection = classSelection, 
    classLabel = classLabel, 
    gene_anno = gene_anno, 
    TissueAnnotationFocus = TissuePatientAnnotationFocus, 
    gsea_data_hallmark = gsea_data_hallmark, 
    species = species,
    msigDBLink = msigDBLink
  )
  
  callModule(
    module = proteinTab,
    id = "protein",
    classSelection = classSelection, 
    classLabel = classLabel, 
    TissueAnnotationFocus = TissuePatientAnnotationFocus
  )
  
  callModule(
    module = derivedDataTab,
    id = "derived_data",
    classSelection = classSelection, 
    classLabel = classLabel
  )
  
  callModule(
    module = survivalTab,
    id = "survival",
    classSelection = classSelection, 
    classLabel = classLabel, 
    PatientAnnotationFuller = PatientAnnotationFuller
  )
  
  callModule(
    module = additionalPropertiesTab,
    id = "additional_props",
    classSelection = classSelection, 
    classLabel = classLabel, 
    TissueAnnotationFocus = TissuePatientAnnotationFocus, 
    Result_otherPropStatistic = Result_otherPropStatistic,
    msigDBLink = msigDBLink
  )
  
  callModule(
    module = tissueOverviewTab,
    id = "overview",
    fm = fm, 
    classSelection = classSelection, 
    classLabel = classLabel, 
    classStack = classStack, 
    TissueAnnotation = TissuePatientAnnotation, 
    TissueAnnotationFocus = TissuePatientAnnotationFocus, 
    gsea_data_hallmark = gsea_data_hallmark
  )
  
  callModule(
    module = machineLearningTab,
    id = "ml",
    fm = fm,
    classSelection = classSelection,
    classLabel = classLabel,
    gsea_data_hallmark = gsea_data_hallmark,
    gene_anno = gene_anno,
    AnnotationFocus = TissuePatientAnnotationFocus,
    Species = reactive(species)
  )
  
  callModule(
    module = aboutTab,
    id = "about"
  )
}
