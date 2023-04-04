# HISTA


<a href="https://zenodo.org/badge/latestdoi/271643615"><img src="https://zenodo.org/badge/271643615.svg" alt="DOI"></a>

[https://conradlab.shinyapps.io/HISTA/](https://conradlab.shinyapps.io/HISTA/)

Human Infertility Single-cell Testis Atlas
Developed by Eisa Mahyari Ph.D. @eisamahyari
PI: Don F. Conrad Ph.D.
Oregon Health & Science University (OHSU)
Div. Reproductive Genetics 

## Associated Manuscript

Eisa Mahyari, Jingtao Guo, Ana C. Lima, Daniel P. Lewinsohn, Alexandra M. Stendahl, Katinka A. Vigh-Conrad, Xichen Nie, Liina Nagirnaja, Nicole B. Rockweiler, Douglas T.Carrell, James M.Hotaling, Kenneth I.Aston, Donald F.Conrad “Comparative single-cell analysis of biopsies clarifies pathogenic mechanisms in Klinefelter syndrome.” The American Journal of Human Genetics 108, no. 10 (2021): 1924-1945

## Introduction

HISTA, is short for the Human Infertility Single-cell Testis Atlas. We have carefully combined, and processed 6 controls, with 2 infertile cases as well as 2 pre-pubescent juveniles and 2 Klinefelter single-cell sample data. 

In processing the amalgam data we utilized (SDA ref-5), a tensor-decomposition method that projects the high-dimensional transcriptomic space, into a set of prev. defined (150 in our case was optimal) latent factors which then can be assessed for batch or signal. With the batch-components removed, a new DGE (cells by genes) matrix is imputed that in theory, captures only the signal. This DGE was used for most downstream analysis such as tSNE/UMAP projects, cluster and condition based DE analysis, and psuedotime trajectory analysis. 


This Shiny App browser serves as a resource for our manuscript (under submission) as well as a tool for exploration and hypothesis generation of future experiments. 

### Refs:

1- Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2019). shiny: Web Application
  Framework for R. R package version 1.4.0. https://CRAN.R-project.org/package=shiny
  
2- Winston Chang and Barbara Borges Ribeiro (2018). shinydashboard: Create Dashboards with 'Shiny'. R
  package version 0.7.1. https://CRAN.R-project.org/package=shinydashboard
  
3- Daniel Wells and Victoria Hore (2017). SDAtools: SDAtools: A toolkit for SDA. R package version 1.0.

4- Stuart and Butler et al. Comprehensive integration of single cell data. bioRxiv (2018).

5- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5010142/


## Contact: 

Conrad Lab: conradlab.org
Packaged created by: @eisamahyari
Package maintained by: @eisamahyari

## Install : 

Tested on R 3.6.3 and 4.0.3
    

    devtools::install_github(repo = 'eisascience/HISTA', dependencies = T, upgrade = 'always')

## Launch : 

    HISTA::launchHISTA()
  
