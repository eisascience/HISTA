
# HISTA

## User Manual (V2.9.4)

## Description of tabs

### Home Page:

The Home page of HISTA provides a graphical introduction into navigating HISTA. On the left is a tabbed menu bar to navigate the main features of HISTA, described in detail below.
 
### Main tab:

This is the main tab that HISTA loads. On the center top, two information boxes provide available background information for the selected SDA component and gene found in the interactive menu on the left panel of this tab. This menu has several parameters to change that alters what is being displayed.
 
* From top to bottom, SDA components can be searched via a numerical input.
* In the next input box, genes (human symbols) can be typed in to not only display the tSNE projected expression of the gene, but to also search in which of components the gene of interest is mostly found or in other words highly loaded.
* In the next selection, via radio button provided, it is possible to select which of the pre-processed tSNE plots are shown (details below); there are three options, the batch-removed SDA cell score matrix (default), the pre-SDA DropSim-normalized expression matrix, or lastly the batch-removed SDA-imputed expression matrix.
* The next set of radio buttons are to visualize available metadata such as donors, replicates, conditions, experiments, and cell type (default).
* At the bottom of this tab, the top loaded positive or negative genes are ranked and listed, however in the interactive menu, how many are shown can be inputted, with 20 being the default.
* The final items in the interactive menu are several buttons to download the top loaded gene lists for export as well as manual navigation of the SDA components.
 
Specifically, each of the 23429 genes available in HISTA can be projected on a precomputed 2D representation of the cells; by default, a tSNE computed on the batch-removed cell-score matrix. The gene expression visualized on the tSNE plot is computed internally, as the dot-product of the selected gene’s loading vector and the selected component’s cell score vector. Two other pre-computed tSNEs are also available via the radio buttons found in the navigation menu under "data origin"; one is derived using the normalized gene matrix (the input to SDA), and the other on the post-SDA batch-removed imputed gene matrix. Other menu options enable selecting which metadata is projected on the metadata tSNE plot (e.g., cell types, donors, etc.).
The cell-score tSNE visualizes the scores as predefined color bins, to visually maximize the signal relative to the overall distribution of the scores. The key takeaway from this plot is the direction of the score as it relates to the direction of the associated gene loadings. The "Cell Scores Across" plot, visualizes the amplitude of the cell scores across the cells where the y-axis is the cell scores and the cells are binned per donor across the x-axis.
We performed GO enrichment analysis of the top 150 positively and negatively (enriched genes) and the GO plots visualize the results, with an adjusted p-value threshold of 0.01 to highlight significance. These top-loaded genes are listed at the bottom of the main page, although the default number of genes selected on the menu is 20 genes.
To visualize any chromosome-specific enrichment relative to the genes, we organize the top loaded genes relative to the human chromosomes and the magnitude of the loading in the "chromosome location" plot.

### Index of Components:

A table of the SDA components and our observations summarized pertaining to each. To curate this table, iterative rounds of analysis were performed on each component and this table represents a summarized form of our SDA findings.
 
### Gene Expression per pathology (boxplot):

Gene expression boxplots of a searched gene, across pathology (CNT, INF1, INF2, KS, JUV) with Wilcox testing, with the ability to subset by cell type. Specifically, to test if the expression of a single gene is distributed similarly (computing p-value via Wilcox-rank-sum test) between the available experimental conditions in HISTA (CNT, KS, JUV, etc.); additionally, through the available radio buttons, it is possible to focus the hypothesis testing by cell type. For example, we can search XIST across all cell types and find it is significantly enriched in Klinefelter Syndrome (KS), (p<2.22e-16 compared to controls). Next, by selecting Sertoli cells (SC), we find that this significant enrichment is lost; as previously reported 1.
 
### Gene Expression per cell type (boxplot):

Gene expression boxplots of a searched gene, across cell types, with ability to subset by pathology (All, CNT, INF1, INF2, KS, JUV). This tab enables the user to quantify the expression of a specific gene, across all cell types. Additionally, using the available radio buttons, the cells can be subsetted by experimental conditions. For example, we can observe that indeed, the Sertoli cell marker SOX9 is exclusively expressed in Sertoli cells. This approach is powerful at rapidly visualizing expression distribution across the available cell types.
 
### Gene Expression per cell type (tSNE):

Gene expression of a searched gene, batch-corrected, mapped on the tSNE projection, with ability to subset by cell type.
 
### Fingerprinting (heatmap):

Two heatmaps, identifying signature/barcode pattern, annotating quantitatively each component, one for positive cell scores, and one for negative cell scores. Chi-squared analysis of the number of cells (positive or negative) per component identifies enrichment of these cells per component and when contrasted by pathology (NT, INF1, INF2, KS, JUV), we add an extra dimension to the enrichment analysis. The pair-wise hierarchical clustering, then identifies similar patterns of enrichment for positively or negatively scored cells, split by pathology. Briefly, the columns represent the SDA components and the rows are the experimental conditions. By defining thresholds on the cell score matrix, the number of cells that are scored positively or negatively are enumerated and passed through a Chi-Squared test; visualized are the transformed residuals which highlight enrichment or depletion. By performing pairwise hierarchical clustering, the components and the samples that are most similar, cluster together. Interestingly, one of the results of this analysis supports our other findings that INF2, the patient with secondary azoospermia, is in fact more similar to the adult controls than INF1, the idiopathic azoospermia patient.
 
### Cell scores per celltype (tSNE):

Scores of a searched component, mapped on the tSNE projection, with ability to subset by cell type. To provide a deeper scope of visualizing the cells scores, projected on the pre-computed tSNE, in this tab it is possible to subset the tSNE by cell type using the available radio buttons. For example, the cells of SDA component 1 with the highest absolute scores can only be found in the spermatid subset of the germ cells. By zooming in to these cells, we can better observe the distinct banding pattern that positively scores the last and early spermatids but negatively scores the spermatids in between them. Digging into the gene loadings of this component, we observe specific gene regulation patterns that explain the banding observed in the cell scores; for example, the top positively loaded genes contain SPRR4 and PRM1 whereas the top negatively loaded genes contain FSCN3 and PRM3 supporting the regulation of spermiogenesis as the spermatids finish their maturation.
 
### Metadata per celltype (tSNE):

Available metadata, mapped on the tSNE projection, with ability to subset by cell type. This tab enables the user to create parallel visualization to the focused tSNEs of the previous tab, but with metadata, providing a closer look at the origin and background of each cell.
 
### Pseudotime (germ only):

A pseudotime trajectory, as previously described, was inferred on the tSNE 2D projection of the germ cells, giving order to cells, driven by the transcriptomic landscape, that also parallels known spermatogenesis trajectory. Plotting the cells scores of SDA components (y-axis) pertaining to germ cells across the cells ordered by pseudotime (x-axis) we commonly see "wave" patterns which translate to transcriptional kinetics of that component relative to the affected cells. There is additional ability to select available metadata to identify correlating scoring patterns.
 
### Pseudotime Component Index:

An index of the SDA components pertaining to germ cells and a summary of our observations. A key feature of this table is the order given to these components relative to the cells (i.e., stages of spermatogenesis) they score with the most magnitude. This ranking is found by identifying the main peak/maxima of the density curve that is fitted to the scores (y-axis) by pseudotime (x-axis) scatter via a peak finding algorithm; the rank is defined by the position of this maxima. For example, the third ranked SDA component, SDA #149, splits early spermatogonial cells into 3 sections where the intermediate section is scored negatively opposite to the other sections, which we describe more in detail in the vignettes section.
 
### Enrichment analysis:

Given a set of genes, hypergeometric testing (K=150) is done to identify which components are highly enriched with those genes. In application, we find an adjusted p-value less than 0.01 to find significant enrichments, however, the fold-enrichment can also be highly informative when significance is not determined. The enrichment analysis tab helps to identify which components are enriched for a specific gene set; for this test we use the universe of detected genes with known GO annotation (N=8025).

### Top loaded components:

There are several figures in this tab. First, there is a histogram layered with a density plot comparing two distributions. Given a set of genes enriched in some biology of interest, e.g., antisense genes, we can compute the distribution of how many of these genes are found across the top positive and negative loaded components. This distribution is layered over a second distribution derived from an equal-length random sample of genes. By comparing these distributions, we are able to evaluate components enriched for or depleted of our selected genes relative to a random sample. 

In the next plot the components are ranked by how many genes of the input set are found in the top loaded components, split by positive and negative loading direction. Based on the first figure in this tab, we can use general statistics like the mean number of genes derived from a random set (~9 for the default gene set) to filter which components are more or less than average enriched or depleted respectively of the input gene set. 

The last figure is a correlation heatmap of the gene loading components, filtered by the input gene. This identifies structural similarity amongst the components in the defined gene space. Each component is enriched with a certain pathology and cell type which are annotated on rows of this heatmap. 

### LncRNAs:

As described in the lncRNA vignette of this manuscript, we have reproduced some of these figures. A Venn diagram shows the overlap of Ensembl lncRNA annotated genes and all genes found in HISTA. As described in the “top loaded components” tab we observe distribution overlap of the number of these genes found in the top loaded genes across components, compared to a random equal-length set of genes. Then the components are sorted and split by direction to identify lncRNA enriched components. See Fig. 2 for more details. 

For convenience, a numerical input box allows for searching and listing the found top loaded lncRNAs; which are really 1348 lncRNAs that overlap the Ensembl annotated lncRNAs and genes detected in HISTA.

### Component Correlations:

The top loaded genes really are at the core of many assessments, including translation of findings. We provide a way to explore and evaluate the relationship of the components narrowed by the number of top genes. First select a component of interest (numeric input). Then use the slider to select how many of the top genes to include in the correlation analysis. For convenience, the top loaded genes used are shown. 
        	
        	
 
### Soma only W. LN19:

In 2019 a new Klinefelter scRNA-Seq was made available through collaboration with Leurentino et al., but after our initial data freeze the SDA-HITA analysis. Furthermore, after integrating this new single KS donor using our customized Seurat pipelines, we found evidence for differences that we believe is correlated with sampling procedure and the age and health of the donor. Therefore, we plan to combine this data in later releases of HISTA, however for the current scope and release, we limited our analysis to only validation of our results with this new data. This tab is provided as support material for our KS manuscript 1.

#### Brief notes:

* As with the existing KS donors, there were no germ cells in this new patient so we could only focus on the somatic cells.
* There is a single large cluster of Sertoli cells (SC), made up of at least 3 sub clusters.
* The two, smaller in diameter subclusters are enriched with JUV SC but only one of them is enriched with the SC from other adults.
* The largest SC subcluster is mostly derived from LN19.
* We previously defined 3 major subtypes of Leydig cells (LC); the progenitors (PLC), the immature (ILC), and the mature (MLC) groups.
* The LN19 LC are found mostly in PLC and MLC clusters.
* LN19 contributes a distinct cluster of cells which we believed to be of myoid phenotype.
 
### LC only W. Zhao21:

In 2021 a new scRNA-Seq dataset was published by Zhao et al. We plan to combine this data in later releases of HISTA, however for the current scope and release, we limited our analysis to only validation of our results with this new data. In fact, more recently in attempts to download the raw sequencing files of this data, significant anomalies have been identified, which may limit our future integration of all of this data. Nonetheless, from their processed count data, we were able to validate our Leydig cell findings, as support material for our KS manuscript 1. As mentioned above, we have previously defined 3 major subtypes of Leydig cells (LC); the progenitors (PLC), the immature (ILC), and the mature (MLC) groups. By combining the LN19 and Zhao21 data, we validate we find these clusters in an integrated dataset, composed of cells from multiple donors, which is critical as these populations, such as the MLCs are fairly infrequent relative to other cells of the testis, especially in normal conditions.


## Case Studies

### Case example 1: 

You have one or more genes and you wish to learn more about them within HISTA.

* To look at the expression of each gene, search the gene in the main tab
* The figure header of the tSNE-expression plot in the main tab, also lists in order the SDA components and the associated loading, which enables identifying gene sets that correlate with the gene of interest
* If you have more than 4 genes, you can use the enrichment analysis tab to find which components they are enriched in and search those components for annotation and finding additional correlating genes.
 
### Case example 2: 

You have a cell type of interest and you wish to learn about the genetic signatures found within HISTA pertaining to the cell type.

* Identify the location of the cell type of interest in the main tab.
* Using the component index, identify which components score to this cell type and search them in the main tab to identify the expressed genes.
* Use the gene expression per cell type tab to evaluate cell type specific expression
 
### Case example 3: 

You have a hypothesis about a gene/transcript and you wish to see if it is differentially expressed across experimental conditions found in HISTA.
* Use the gene expression by pathology boxplot tab to identify significant differences
 

