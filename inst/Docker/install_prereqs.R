install.packages(c(
  "BiocManager", "remotes", # install helpers
  "purrr", "dplyr", "tidyr", "forcats", "ggplot2", "readr", "readxl", "broom", # tidyverse (only required pieces)
  "shiny", "shinyjs", "shinyBS", "shinyWidgets", "shinythemes", "DT", "future", # shiny-related
  "DBI", "RPostgres", # DB handling
  "gggenes", "ggrepel", "plotROC", "pheatmap", "svglite", "gridExtra", # plotting
  "apcluster", "matrixTests", "ROCR", "Rtsne", "umap", "phateR", "fpc", # cluster, dimred
  "ModelMetrics", "caret", "randomForest", "e1071", # stats
  "NeuralNetTools", "Boruta", "neuralnet", "glmnet", # machine learning
  "logger", "config", "xfun", "xml2"
))

BiocManager::install(c(
  "limma", "DESeq2", "edgeR", "gage"
))
ignored <- c("gageData", "MASS", "viridis")
