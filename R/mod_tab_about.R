aboutTabUI_main <- function(id, docuLink){
  ns <- NS(id)
  
  text <- div(
    h3("Tissue Differences (TIFF)"),
    p("For cancer tissue data sets"),
    p("written by Andreas Wernitznig (2018-2022), Jesse Lipp (2019), Zygmunt Zawadzki (2021-2022) and Michał Jakubczak (2020-2022)"),
    p(
      "Gene Set Enrichment Analysis (GSEA) is based on the R-package", 
      a("gage", href = "https://bioconductor.org/packages/release/bioc/html/gage.html", target = "_blank"), 
      "and is using the Broad Institute's", 
      a("MSigDB", href = "http://software.broadinstitute.org/gsea/msigdb/", target = "_blank")
    ),
    p(
      "Immune environment classification originates from",
      a(
        "Thorsson et al., The Immune Landscape of Cancer", 
        href ="https://doi.org/10.1016/j.immuni.2018.03.023", target = "_blank"
      )
    ),
    p(
      "MSI-sensor score and classification is imported from", 
      a(
        "Ding et al., Perspective on Oncogenic Processes at the End of the Beginning of Cancer Genomics", 
        href = "https://doi.org/10.1016/j.cell.2018.03.033", target = "_blank"
      )
    ),
    p(
      "Gastro-intestinal molecular subtype is copied from", 
      a(
        "Liu et al., Comparative Molecular Analysis of Gastrointestinal Adenocarcinomas", 
        href = "https://doi.org/10.1016/j.ccell.2018.03.010", target = "_blank"
      )
    ),
    p(
      "iCluster is taken from", 
      a(
        "Hoadley et al. Cell-of-Origin Patterns Dominate the Molecular Classification of 10,000 Tumors from 33 Types of Cancer", 
        href = "https://doi.org/10.1016/j.cell.2018.03.022", target = "_blank"
      )
    ),
    p(
      "Consensus Molecular subtype is calulated according to", 
      a(
        "Guinney et al., The Consensus Molecular Subtypes of Colorectal Cancer", 
        href = "https://dx.doi.org/10.1038%2Fnm.3967", target = "_blank"
      )
    ),
    p(
      "Immune cell scores are based on ", 
      a(
        "Aran et al. xCell: digitally portraying the tissue cellular heterogeneity landscape", 
        href = "https://doi.org/10.1186/s13059-017-1349-1", target = "_blank"
      )
    ),
    p(
      "TIL pattern is the Tumor Infiltrating Lymphocytes (TIL) Map Structural Pattern from ",
      a(
        "Saltz et al. Spatial Organization and Molecular Correlation of Tumor-Infiltrating Lymphocytes Using Deep Learning on Pathology Images", 
        href = "https://doi.org/10.1016/j.celrep.2018.03.086", target = "_blank"
      )
    ),
    p(
      "Metabolics data are taken from ", 
      a(
        "Choi et al. Pan-cancer analysis of tumor metabolic landscape associated with genomic alterations", 
        href = "https://doi.org/10.1186/s12943-018-0895-9", target = "_blank"
      ), 
      "License: ", 
      a(
        "Creative Commons Attribution 4.0 International License", 
        href = "http://creativecommons.org/licenses/by/4.0/", target = "_blank"
      )
    ),
    p(
      "Signaling pathways are derived from ", 
      a(
        "Sanches-Vega et al. Oncogenic Signaling Pathways in The Cancer Genome Atlas", 
        href = "https://doi.org/10.1016/j.cell.2018.03.035", target = "_blank"
      )
    ),
    p(
      "Number of clones was calculated in ", 
      a(
        "Raynaud et al. Pan-cancer inference of intra-tumor heterogeneity reveals associations with different forms of genomic instability", 
        href = " https://doi.org/10.1371/journal.pgen.1007669", target = "_blank"
      )
    ),
    p(paste("Forked processes support:", `if`(supportsMulticore(), "yes", "no"))),
    p(
      R.Version()$version.string, 
      downloadLink(
        outputId = ns("session_info"),
        label = "Session info"
      )
    )
  )
  
  table <- tableOutput(ns("information"))
  
  fluidRow(
    column_3(
      div(table, class = "about-information-table")
    ),
    column_9(text)
  )
}

aboutTabUI_sidebar <- function(id){
  ns <- NS(id)
  NULL
}

aboutTab <- function(input, output, session){
  output$information <- renderTable({
    getInformation()
  })
  
  output$session_info <- downloadHandler(
    filename = "tiff_sessionInfo.txt",
    content = function(file){
      writeLines(capture.output(sessionInfo()), file)
    }
  )
}
