
# HISTA

## User Manual (V2.9.4)

## Description of tabs


# HISTA

## User Manual (V2.9.4)

## Description of tabs

### Home Page:

        	The Home page of HISTA provides a graphical introduction to navigating HISTA. On the left is a tabbed menu bar to navigate the main features of HISTA, described in detail below.
 
### Main tab:

        	This is the 'Main’ tab that HISTA loads. On the center top, two information boxes provide available background information for the selected SDA component and gene found in the interactive menu on the left panel of this tab. This menu has several parameters to change that alter what is being displayed.

From top to bottom of the page:

* In the Inputs sections, the first input allows searching SDA components via a numerical input.
* In the next box, genes (human symbols e.g. PRM1) can be searched to display the gene's t-SNE projected expression. Additionally, HISTA will highlight the components that the gene of interest is mostly found or, in other words, highly loaded.
* In the next selection, via the radio button provided, it is possible to select which of the pre-processed t-SNE plots are shown (there are four options):
The batch-removed SDA cell score matrix t-SNE (default). This t-SNE was produced on non-batch components of the SDA score matrix i.e., cells by components
The batch-removed SDA cell score matrix UMAP.
The pre-SDA DropSim-normalized expression matrix t-SNE. This t-SNE was produced on the gene expression matrix normalized by DropSim, which includes a square-root-based transformation.
the batch-removed SDA-imputed expression matrix t-SNE. As explained earlier, by performing dot-product of the gene loading and cell score matrices, with only the non-batch components, we impute a batch-removed gene expression matrix on which t-SNE was run.
* The next set of radio buttons is to visualize available metadata such as donors, replicates, conditions, experiments, cell cycle, and cell type (default).
* The final items in the interactive menu are several buttons to download the top-loaded gene lists for export as well as manual navigation of the SDA components.
* The figures below the input section from top to bottom include 
the gene expression, cell score, and metadata 2D projection (t-SNE or UMAP).
The cell score across donors scatter plot, which aids in seeing the selected component’s score distribution
* The gene ontology (GO) plots proved a pre-computed analysis of the top loaded genes that aid in translating the general signature observed.
* The chromosome location highlights the loading weight of each gene relative to their position across the chromosomes. 
* At the bottom of this tab, the top loaded positive or negative genes are ranked and listed, but how many are shown can be inputted, with 20 being the default.

### Index of Components:

        	A table of the SDA components and our observations summarized pertaining to each. To curate this table, iterative rounds of analysis were performed on each component, representing a summarized form of our SDA findings.
 
### Fingerprinting (heatmap):

        	Two heatmaps, identifying signature/barcode pattern, annotating each component quantitatively, one for positive and one for negative cell scores. Chi-squared analysis of the number of cells (positive or negative) per component identifies the enrichment of these cells per component. When contrasted by pathology (CNT, INF1, INF2, KS, JUV), we add an extra dimension to the enrichment analysis. The pair-wise hierarchical clustering (which can be turned off using the radio buttons provided) identifies similar enrichment patterns for positively or negatively scored cells, split by available metadata options, selectable via radio buttons, e.g., cell types, donors, pathology, etc. The columns represent the SDA components, and the rows are the experimental conditions. By defining thresholds on the cell score matrix, the number of cells that are scored positively or negatively are enumerated and passed through a Chi-Squared test; visualized are the transformed residuals that highlight enrichment or depletion by performing pairwise hierarchical clustering, the components and the samples that are most similar, cluster together. Interestingly, one of the results of this analysis supports our other findings that INF2, the patient with secondary azoospermia, is more similar to the adult controls than INF1, the idiopathic azoospermia patient.

### Gene Expression per pathology (boxplot):

        	Gene expression boxplots of a searched gene across pathology (CNT, INF1, INF2, KS, JUV) with Wilcox testing, with the ability to subset by cell type. Specifically, to test if the expression of a single gene is distributed similarly (computing p-value via Wilcox-rank-sum test) between the available experimental conditions in HISTA (CNT, KS, JUV, etc.); additionally, through the available radio buttons, it is possible to focus the hypothesis testing by cell type. For example, we can search XIST across all cell types and find it is significantly enriched in Klinefelter Syndrome (KS) (p<2.22e-16 compared to controls). Next, by selecting Sertoli cells (SC), we find that this significant enrichment is lost, as previously reported1.

### Gene Expression per cell type (boxplot):

        	Gene expression boxplots of a searched gene, across cell types, with the ability to subset by pathology (All, CNT, INF1, INF2, KS, JUV). This tab enables the user to quantify the expression of a specific gene across all cell types. Unlike the per pathology boxplots, statistics are not provided due to minimizing the complexity of the figure. For example, we can observe that the Sertoli cell marker SOX9 is exclusively expressed in Sertoli cells. This approach is powerful at rapidly visualizing expression distribution across the available cell types.
        	
### Gene Expression per cell type (2D):

        	Gene expression of a searched gene, batch-corrected, mapped on the t-SNE projection, with the ability to subset by cell type. Additional radio buttons allow selecting which 2D plot to show, e.g., t-SNE or UMAP, on various scopes of our data.
        	
### Cell scores per cell type (2D):

        	Scores of a searched component, mapped on the 2D (t-SNE or UMAP) projection, with the ability to subset by cell type. To provide a deeper scope of visualizing the cell scores projected on the pre-computed 2D representation. In this tab, it is possible to subset the figure by cell type using the available radio buttons. For example, SDA component 1 demonstrates that the highest absolute scores are found in the spermatid population of the germ cells. By zooming in on these cells, we can better observe the distinct banding pattern that positively scores the last and early spermatids but negatively scores the spermatids between them. Digging into the gene loadings of this component, we observe specific gene regulation patterns that explain the banding observed in the cell scores; for example, the top positively loaded genes contain SPRR4 and PRM1, whereas the top negatively loaded genes contain FSCN3 and PRM3 supporting the regulation of spermiogenesis as the spermatids finish their maturation.
        	
### Metadata per cell type (2D):

        	In this tab, available metadata (selectable) are mapped on a 2D projection which can be subset by cell type (selectable). This tab enables the user to create parallel visualization to the focused t-SNE or UMAP of the cell score or gene expression 2D tabs, but with metadata, providing a closer look at the origin and background of each cell.
        	
### Gene correlations:

        	In this tab, the user can input a set of genes and select a cell type via the radio buttons, which visualize a heatmap of gene-gene correlations within the selected cells. 

### Component Correlations:

        	The top-loaded genes really are at the core of many assessments, including the translation of findings. We provide a way to explore and evaluate the relationship of the components narrowed by the number of top genes. First, select a component of interest (numeric input). Then use the slider to select how many top genes to include in the correlation analysis. Visualized is a heatmap of component-component correlation. For convenience, the top-loaded genes used are shown. 

### Pseudotime Meta (germ only):

          A pseudotime trajectory, as previously described, was inferred on the t-SNE 2D projection of the germ cells, giving order to cells driven by the transcriptomic landscape that also parallels the known spermatogenesis trajectory. Plotting the cell scores of SDA components (y-axis) pertaining to germ cells across the cells ordered by pseudotime (x-axis), we commonly see "wave" patterns that translate to transcriptional kinetics of that component relative to the affected cells. There is an additional ability to select available metadata to identify correlating scoring patterns.

### Pseudotime Gene (germ only):

        	This tab allows users to type in a gene of interest to visualize the expression wave across our defined pseudotime. Additionally, metadata selection is available to visualize differential expression patterns. 

### Pseudotime Component Index:

        	An index of the SDA components pertaining to germ cells and a summary of our observations. A key feature of this table is the order given to these components relative to the cells (i.e., stages of spermatogenesis) they score with the most magnitude. This ranking is found by identifying the main peak/maxima of the density curve that is fitted to the scores (y-axis) by pseudotime (x-axis) scatter via a peak finding algorithm; the position of this maxima defines the rank. For example, the third-ranked SDA component, SDA #149, splits early spermatogonial cells into 3 sections where the intermediate section is scored negatively opposite to the others, which we describe more in detail in the vignettes section.

### Enrichment analysis:

        	Given a set of genes, hypergeometric testing (K=150) identifies which components are highly enriched with those genes. In application, we find an adjusted p-value less than 0.01 to find significant enrichments. However, the fold enrichment can also be highly informative when significance is not determined. This statistical approach is most appropriate for less than ~30 more than 3 genes.

### Top loaded components:

        	This tab is another statistical approach to identify key components that parallel enrichment analysis, although this method is best suited for larger sets of genes. As a default, the input is loaded with antisense genes. There are several figures in this tab. First, there is a histogram layered with a density plot comparing two distributions. Given a set of genes enriched in some biology of interest, e.g., antisense genes, we can compute the distribution of how many of these genes are found across the top positive and negative loaded components. This distribution is layered over a second distribution derived from an equal-length random sample of genes (computed dynamically). By comparing these distributions, we can evaluate components enriched for or depleted of our selected genes relative to a random sample. 
        	In the next plot, the components are ranked by how many genes of the input set are found in the top loaded components, split by positive and negative loading direction. Based on the first figure in this tab, we can use general statistics like the mean number of genes derived from a random set (~9 for the default gene set) to filter which components are more or less than average enriched or depleted respectively of the input gene set. 
        	The last figure is a correlation heatmap of the gene loading components filtered by the input gene. This identifies structural similarity amongst the components in the defined gene space. Each component is enriched with a certain pathology and cell type annotated on this heatmap's rows. 

### lncRNAs:

          As described in this manuscript's long non-coding RNA (lncRNA) vignette, we have reproduced some of these figures. Additionally, this tab yields the lncRNAs are provided by entering the component number interest in the input box.  On the top, a Venn diagram shows the overlap of Ensembl lncRNA annotated genes and all genes found in HISTA. As described in the 'top loaded components’ tab, we observe distribution overlap of the number of these genes found in the top loaded genes across components, compared to a random equal-length set of genes. Then the components are sorted and split by direction to identify lncRNA-enriched components. 



### Soma only W. LN19:
          
          In 2019 a new Klinefelter scRNA-Seq was made available through collaboration with Leurentino et al., but after our initial data freeze, the SDA-HISTA analysis. Furthermore, after integrating this new single KS donor using our customized Seurat pipelines, we found evidence for differences that we believe are correlated with the sampling procedure and the age and health of the donor. Therefore, we plan to combine this data in later releases of HISTA. However, we limited our analysis to only validating our results with this new data for the current scope and release. This tab is provided as support material for our KS manuscript 1.
          
Brief notes:

* As with the existing KS donors, there were no germ cells in this new patient, so we could only focus on the somatic cells.
* There is a single large cluster of Sertoli cells (SC) comprising at least 3 sub-clusters.
*  The two smaller in diameter subclusters are enriched with JUV SC, but only one is enriched with the SC from other adults.
*   The largest SC subcluster is mostly derived from LN19.
* We previously defined 3 major subtypes of Leydig cells (LC); the progenitors (PLC), the immature (ILC), and the mature (MLC) groups.
*   The LN19 LC are found mostly in PLC and MLC clusters.
* LN19 contributes a distinct cluster of cells that we believed to be of myoid phenotype.
 
### LC only W. Zhao21:

        	In 2021, a new scRNA-Seq dataset was published by Zhao et al. We plan to combine this data in later releases of HISTA; however, for the current scope and release, we limited our analysis to only validating our results with this new data. More recently, in attempts to download the raw sequencing files of this data, significant anomalies have been identified, which may limit our future integration of all of this data. Many of the samples from this dataset failed bioinformatics QC processes. Nonetheless, using their processed count data, we validated our Leydig cell findings as support material for our KS manuscript1. As mentioned above, we have previously defined 3 major subtypes of Leydig cells (LC); the progenitors (PLC), the immature (ILC), and the mature (MLC) groups. By combining the LN19 and Zhao21 data, we validate we find these clusters in an integrated dataset composed of cells from multiple donors, which is critical as these populations, such as the MLCs, are fairly infrequent relative to other cells of the testis, especially in normal conditions.
 
## HISTA Case Examples:

Case example 1: You have one or more genes, and wish to learn more about them within HISTA.

* To look at the expression of each gene, search for the gene in the main tab

* The figure header of the 2D (t-SNE or UMAP) expression plot in the main tab also lists in order of the SDA components and the associated loading, which enables identifying gene sets that correlate with the gene of interest

* To test pathology differential expression, on specific cell type (or all cells), the 'gene expression per pathology’ tab visualizes this with Wilcox statistics.
*     To visualize expression across cell types but in specific pathologies, use the 'gene expression per cell type' tab. 
* If this gene is a germ cell expressed gene, try the 'pseudotime expr’ tab to observe the expression pattern across pseudotime.
* If you have more than ~3 genes, you can use the 'Enrichment analysis’ tab to find which components they are enriched in and search those components for annotation and find additional correlating genes. You can also use the 'top loaded components’ tab if your set of genes is fairly large. 


