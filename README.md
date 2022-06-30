# What is TIFF?

## TIFF shiny application
TIFF (Tissue Differences) is a web application for the analysis of
differences between two sets of tissues in terms of differential gene or
protein expression, DNA copy number or DNA mutations, mutational burden,
microsatellite instability, survival as well as serveral derived data like 
immune cell enrichment, metabolomics, activation of signaling pathways, 
number of clones or gene signatures.

TIFF is also capable to perform unsupervised analysis as descriptive
statistics and dimension reduction plots of gene expression data.

## TIFF Application Programming Interface (API)

The core functionality of TIFF is also available for R-programmers. A set of data 
retrieval, processing, and visualization functions enables access to functionalities,
which are usually only available for bioinformatics experts.

### Quick start

Start an R session in RStudio and run the following code:

```
library(TIFF)
library(XIFF)
library(tidyverse)

setDbOptions(getSettings())
tissue_anno <- getTissueAnnotation()

ggplot(tissue_anno, aes(x = tissuepanel, fill = tissuepanel)) +
  geom_bar() + theme(legend.position = "none")

```
This plot gives you an overview on the number of tumor tissues, normal
tissues and PDX models.

### Running the TIFF shiny app in RStudio
```
TIFF::run()
```
