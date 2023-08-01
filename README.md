# HISTA



<a href="https://zenodo.org/record/8206586"><img src="https://zenodo.org/badge/271643615.svg" alt="DOI"></a>


To use the latest version on the web:

[https://conradlab.shinyapps.io/HISTA/](https://conradlab.shinyapps.io/HISTA/)

The Human Infertility Single-cell Testis Atlas (HISTA): An interactive molecular scRNA-Seq reference of the human testis

Developed by Eisa Mahyari Ph.D. @eisamahyari

PI: Don F. Conrad Ph.D.

Oregon Health & Science University (OHSU)

Oregon National Primate Research Center (ONPRC) 

## Associated Manuscript

Eisa Mahyari, Jingtao Guo, Ana C. Lima, Daniel P. Lewinsohn, Alexandra M. Stendahl, Katinka A. Vigh-Conrad, Xichen Nie, Liina Nagirnaja, Nicole B. Rockweiler, Douglas T.Carrell, James M.Hotaling, Kenneth I.Aston, Donald F.Conrad “Comparative single-cell analysis of biopsies clarifies pathogenic mechanisms in Klinefelter syndrome.” The American Journal of Human Genetics 108, no. 10 (2021): 1924-1945

## Introduction

The Human Infertility Single-cell Testis Atlas (HISTA) is an interactive web tool and a reference for navigating the transcriptome of the human testis. It was developed using joint analyses of scRNA-Seq datasets derived from a dozen donors, including healthy adult controls, juveniles, and several infertility cases. HISTA provides visualization and hypothesis testing tools using 23429 genes measured across 26093 cells. A central feature of HISTA is to combine structured dimensionality reduction with genomic annotation to fine-map cellular characteristics and discover genomic modules. 

The testis is a complex reproductive tissue with unique microenvironments and cell types, that recently has been the focus of numerous single-cell transcriptomics (scRNA-Seq) studies to unravel both normal and pathological features1–7. Recently, we published our findings 1 with 26093 cells, derived from testis biopsies of 2 juveniles, 6 normal adults, 1 adult with azoospermia, 1 adult with ejaculatory dysfunction, and 2 adults with Klinefelter Syndrome (KS). To facilitate continued research and to give access and easy navigation of the testis transcriptomic data at the core of our work, we developed the Human Infertility Single-cell Testis Atlas (HISTA). 

HISTA was initially used to identify molecular and cellular signatures that distinguished infertility in Klinefelter Syndrome (KS) in comparison to healthy controls as well as other forms of infertility such as non-obstructive azoospermia (NOA) and ejaculatory dysfunction


More recently, Nagirnaja, et al., utilized HISTA to identify molecular subforms of NOA. They observed differences in gene expression in patients with NOA who had different histological diagnoses (such as maturation arrest (MA), Sertoli cell only (SCO), or "unknown"). Using HISTA, they found that the gene signatures were different for patients with MA compared to those with SCO or "unknown" histology

### For local runs

You will need to download the rds file with the data once you install this code base. 

https://zenodo.org/record/8206603

### Refs:

Mahyari, E. et al. Comparative single-cell analysis of biopsies clarifies pathogenic mechanisms in Klinefelter syndrome. Am. J. Hum. Genet. 108, 1924–1945 (2021).

Nagirnaja, L. et al. Diverse monogenic subforms of human spermatogenic failure. Nat. Commun. 13, 7953 (2022).

Hermann, B. P. et al. The Mammalian Spermatogenesis Single-Cell Transcriptome, from Spermatogonial Stem Cells to Spermatids. Cell Rep. 25, 1650–1667.e8 (2018).

Jung, M. et al. Unified single-cell analysis of testis gene regulation and pathology in five mouse strains. eLife vol. 8 Preprint at https://doi.org/10.7554/elife.43966 (2019).

Wells, D and Hore, V. (2017). SDAtools: SDAtools: A toolkit for SDA. R package version 1.0.

Hore, V. et al. Tensor decomposition for multiple-tissue gene expression experiments. Nat. Genet. 48, 1094–1100 (2016).

Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2019). shiny: Web Application Framework for R. R package version 1.4.0. https://CRAN.R-project.org/package=shiny
  
Winston Chang and Barbara Borges Ribeiro (2018). shinydashboard: Create Dashboards with 'Shiny'. R package version 0.7.1. https://CRAN.R-project.org/package=shinydashboard
  
Mahyari, E. & Conrad, D. HISTA. (2021) doi:https://zenodo.org/badge/latestdoi/271643615

## Contact: 

Conrad Lab: conradlab.org
Packaged created by: @eisamahyari
Package maintained by: @eisamahyari

## Install : 

Tested on R 3.6.3 and 4.0.3
    

    devtools::install_github(repo = 'eisascience/HISTA', dependencies = T, upgrade = 'always')

## Launch : 

    HISTA::launchHISTA()
  
