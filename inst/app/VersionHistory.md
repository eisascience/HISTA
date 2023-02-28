
# HISTA


## V2.9.4

### Top enriched componnents with specific Gene sets

Similar to the lncRNA specific analysis, however you can input any gene set of interest to identify enrichment of the set in the top 200 genes of each component; see the 'Top Loaded Components' tab.

Note: There is already an enrichment analyis tab which uses a statistical method to formally test components enriched with a smaller gene set. However, in this tab, we can input large sets of genes of interest, as an example we have input the genes with -AS (antisense genes) to highlight this feature. 

## V2.9.3

### lncRNA analysis added

Added a new tab that that include the lncRNA specific analysis; the "lncRNA - Expr of top-loaded" tab. 

## V2.9.2

### Pseudotime pdf download button fixed

The original code had a bug which prevented download of the pseudotime plot as a pdf.


## V2.9.1

### Gene Expr per cell type (boxplot) tab bug fix

The original code had a bug which the metadata was not properly merged so wrong labels were shown.


## V2.9 <current>

### Computation of Gene expression, fliter for selected cells prior to dot-product of SDA matrices

The original code in HISTA prior to this version, when selecting a subset of cells, it did not subset the SDA matrices prior to the dot-product, rather the dot-product was done on the entire matrix then subseting was done. see V2.8 GEX by cell type update.

### Cell cycle labels
Cell cycle (G1, G2M, or S phase) was inferred using Seurat methods on a Seurat-processed object of this data. This label is now available via the selection bottons as a meta feature. 

### Add version history
Prior to this entry, the version history reflects major updates, as minor changes were not captured. However, going forward, all changes will be captured herein.

## V2.8

### GEX by cell type 
A new tab is added to show Gene expression of a select gene, but also allow for selecting which cells to filter, as oppse to the entire set of cells. 

## V2.7

### fliter for removed components

The original code in HISTA prior to this version, did not remove the "batch/removed" components prior to the dot-product, thus the gene expression shown was closer to the original GEX (Y), however with this fix, the expression is batch-removed and thus this estimate (Y') is likely the more accurate signal.

### Zhao 2020 tab added
Parallel analysis with Laurentino et al. 2020 donors. 


## V2.6

### Addition of index tables 
A new tab was added to house the metadata associated to the SDA components and our observations

### Laurentino 2019 tab added
Parallel analysis with Laurentino et al. 2019 KS donor. 

## V2.5

### Germ cell pseudotime added
A new tab for the germ cell only pseudotime analysis is added. 

## V2.0 

### ShinyApps.io modifications
To host this app on ShinyApps.io the code was slightly updated; mainly the input/output


## V1.5

### Converted to Rpackage 
For easy depoyment in various situations, public and private, HISTA was made into an R package
