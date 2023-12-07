# HISTA

## User Manual (V2.9.7)

## Tabs Overview

### Home Page
Graphical introduction to HISTA with a tabbed menu bar.

### Main Tab
Default tab with interactive menu for SDA components and gene exploration.

### Index of Components
Table summarizing SDA components and observations.

### Fingerprinting (Heatmap)
Heatmaps identifying signature patterns for positive and negative cell scores.

### Gene Expression per Pathology (Boxplot)
Gene expression boxplots across pathologies with Wilcox testing.

### Gene Expression per Cell Type (Boxplot)
Gene expression boxplots across cell types.

### Gene Expression per Cell Type (2D)
Gene expression on t-SNE projection, subset by cell type.

### Cell Scores per Cell Type (Boxplot)
Scores of a searched component, mapped on defined cell types.

### Cell Scores per Cell Type (2D)
Scores of a searched component, mapped on 2D projection.

### Metadata per Cell Type (2D)
Metadata mapped on 2D projection, subset by cell type.

### Gene Correlations
Heatmap of gene-gene correlations within selected cells.

### Component Correlations
Heatmap of component-component correlation based on top-loaded genes.

### Pseudotime Meta (Germ Only)
Pseudotime trajectory visualizing cell scores of SDA components across germ cells.

### Pseudotime Gene (Germ Only)
Visualize gene expression across pseudotime for a searched gene.

### Pseudotime Component Index
Index of SDA components pertaining to germ cells and observations.

### Enrichment Analysis
Hypergeometric testing to identify components enriched with a set of genes.

### Top Loaded Components
Statistical approach to identify key components enriched or depleted of a gene set.

### DEgenes
Explore Differentially Expressed (DE) Gene tables derived in the original analysis.

### lncRNAs
Explore long non-coding RNAs (lncRNAs) associated with a specific component.

### Soma Only W. LN19
Analysis limited to somatic cells with a new Klinefelter scRNA-Seq dataset.

### LC Only W. Zhao21
Validation of Leydig cell findings with a new scRNA-Seq dataset.


## Description of tabs

### Home Page:

        	The Home page of HISTA provides a graphical introduction to navigating HISTA. On the left is a tabbed menu bar to navigate the main features of HISTA, described in detail below.
 
### Main tab:

The 'Main' tab is the primary landing page when HISTA is launched. At the center-top, two information boxes present background details for the selected SDA component and gene, accessible through the interactive menu on the left panel of this tab. The menu offers various parameters to modify and customize the displayed content.

Starting from the top and progressing downwards on the page:

* **Inputs Section:**
  - The initial input field enables the search for SDA components using numerical input.
  - In the subsequent box, genes (human symbols, e.g., PRM1) can be queried to reveal the gene's t-SNE projected expression. HISTA also highlights the components where the gene of interest is predominantly found or highly loaded.
  - The following selection, through a provided radio button, allows the user to choose which pre-processed t-SNE plots to display. Options include:
    - The default batch-removed SDA cell score matrix t-SNE, produced on non-batch components of the SDA score matrix (cells by components).
    - The batch-removed SDA cell score matrix UMAP.
    - The pre-SDA DropSim-normalized expression matrix t-SNE, derived from gene expression matrix normalization by DropSim, featuring a square-root-based transformation.
    - The batch-removed SDA-imputed expression matrix t-SNE, obtained by performing a dot-product of the gene loading and cell score matrices, considering only non-batch components, resulting in a batch-removed gene expression matrix on which t-SNE was executed.

* **Visualization Options:**
  - The subsequent set of radio buttons facilitates the visualization of available metadata, such as donors, replicates, conditions, experiments, cell cycle, and cell type (default).
  - The final elements in the interactive menu include buttons for downloading top-loaded gene lists for export and manually navigating the SDA components.

* **Figures Below Input Section:**
  - Gene expression, cell score, and metadata 2D projection (t-SNE or UMAP).
  - Scatter plot of cell score across donors, providing insight into the distribution of scores for the selected component.
  - Gene ontology (GO) plots offering a pre-computed analysis of top-loaded genes, aiding in translating the observed general signature.
  - Chromosome location highlighting the loading weight of each gene relative to its position across the chromosomes.

* **Bottom of the Tab:**
  - The top-loaded positive or negative genes are ranked and listed. Users can specify the number of displayed genes, with 20 being the default.


### Index of Components:

A table of the SDA components and our observations summarized pertaining to each. To curate this table, iterative rounds of analysis were performed on each component, representing a summarized form of our SDA findings.
 
### Fingerprinting (heatmap):

Two heatmaps are utilized to identify signature/barcode patterns, annotating each component quantitatively. One heatmap is dedicated to positive cell scores, and the other to negative cell scores. A Chi-squared analysis is employed to assess the number of cells (positive or negative) per component, thus identifying the enrichment of these cells per component.

When contrasting by pathology (CNT, INF1, INF2, KS, JUV), an extra dimension is added to the enrichment analysis. The pair-wise hierarchical clustering, which can be toggled using the provided radio buttons, identifies similar enrichment patterns for positively or negatively scored cells. Users can choose from available metadata options via radio buttons, such as cell types, donors, pathology, etc.

The columns of the heatmap represent the SDA components, while the rows correspond to the experimental conditions. By defining thresholds on the cell score matrix, the number of cells scored positively or negatively is enumerated and subjected to a Chi-Squared test. Visualized are the transformed residuals that highlight enrichment or depletion. Pairwise hierarchical clustering is applied, revealing components and samples that are most similar, clustering together.

An interesting result of this analysis supports previous findings that INF2, the patient with secondary azoospermia, is more similar to the adult controls than INF1, the idiopathic azoospermia patient.

### Gene Expression per pathology (boxplot):

Gene expression boxplots depict the expression profile of a searched gene across different pathologies (CNT, INF1, INF2, KS, JUV). These boxplots are accompanied by Wilcox testing, assessing the distributional similarity of a single gene's expression between the available experimental conditions in HISTA (e.g., CNT, KS, JUV). The p-value is computed via the Wilcox rank-sum test.

Additionally, users have the ability to subset the analysis by cell type using available radio buttons. For example, searching for the XIST gene across all cell types may reveal its significant enrichment in Klinefelter Syndrome (KS) (p < 2.22e-16 compared to controls). Further refinement is possible by selecting specific cell types, such as Sertoli cells (SC), where the previously observed significant enrichment may be lost, as reported previously1.

### Gene Expression per cell type (boxplot):

Gene expression boxplots illustrate the expression pattern of a searched gene across various cell types, with the option to subset by pathology (All, CNT, INF1, INF2, KS, JUV). This tab allows users to quantify the expression of a specific gene across all available cell types.

Unlike the per pathology boxplots, statistics are not provided in this case, aiming to minimize the complexity of the figure. For instance, this visualization approach enables a quick observation, such as noting that the Sertoli cell marker SOX9 is exclusively expressed in Sertoli cells. This method proves powerful in rapidly visualizing the distribution of gene expression across the available cell types.
        	
### Gene Expression per Cell Type (2D):

This tab displays the gene expression of a searched gene, batch-corrected and mapped onto the t-SNE projection. Users have the flexibility to subset the analysis by cell type. Additional radio buttons offer the choice of displaying the 2D plot on various scopes of the data, such as t-SNE or UMAP.

This visualization provides insights into the spatial distribution of gene expression in a 2D projection, allowing users to observe how the expression of a specific gene varies across different cell types. The option to select different 2D plots enhances the exploration of the gene expression landscape in the context of the chosen dimensionality reduction method.


### Cell Scores per Cell Type (boxplot):

This tab visualizes scores of a searched component, mapped onto defined cell types, with the option to subset donor sets. It facilitates exploration into the scoring patterns of each component relative to the specified cell types.

Combined with the Gene Expression per Cell Type (Boxplot) feature, users can deeply investigate each cell type. The parallel visualization of cell scores and gene expression allows for a comprehensive understanding of how the searched component behaves across different cell types and donor sets.

### Cell Scores per Cell Type (2D):

This tab presents scores of a searched component mapped onto the 2D (t-SNE or UMAP) projection, offering the option to subset by cell type. The purpose is to provide a more in-depth visualization of cell scores projected onto the pre-computed 2D representation.

Users can subset the figure by cell type using the available radio buttons. For instance, examination of SDA component 1 reveals that the highest absolute scores are concentrated in the spermatid population of germ cells. Zooming in on these cells allows for a closer observation of the distinct banding pattern, where last and early spermatids show positive scores while the spermatids between them exhibit negative scores.

Digging into the gene loadings of this component unveils specific gene regulation patterns that explain the observed banding. For example, the top positively loaded genes include SPRR4 and PRM1, while the top negatively loaded genes consist of FSCN3 and PRM3, supporting the regulation of spermiogenesis as spermatids conclude their maturation.

### Metadata per cell type (2D):

This tab allows users to map available metadata (selectable) onto a 2D projection, with the option to subset by cell type. It enables users to create a parallel visualization to the focused t-SNE or UMAP of the cell score or gene expression 2D tabs. The inclusion of metadata provides a closer look at the origin and background of each cell, enhancing the understanding of the characteristics associated with different cell types.

### Gene correlations:

This tab allows users to input a set of genes and select a cell type via radio buttons. The visualization consists of a heatmap displaying gene-gene correlations within the selected cells. Users can explore the relationships and patterns of expression among the specified set of genes within the chosen cell type.


### Component Correlations:

The top-loaded genes play a crucial role in various assessments, influencing the translation of findings. This tab offers a way to explore and evaluate the relationship of the components, narrowed down by the number of top genes.

To begin, select a component of interest using numeric input. Then, use the slider to choose how many top genes to include in the correlation analysis. The result is visualized in a heatmap illustrating component-component correlations. For convenience, the top-loaded genes used in the analysis are displayed.

### Pseudotime Meta (germ only):

In this tab, a pseudotime trajectory is inferred on the t-SNE 2D projection of germ cells. This trajectory provides order to cells driven by the transcriptomic landscape, parallel to the known spermatogenesis trajectory.

The cell scores of SDA components (y-axis) related to germ cells are plotted across cells ordered by pseudotime (x-axis). Commonly observed "wave" patterns translate to transcriptional kinetics of that component relative to the affected cells.

Additionally, users have the ability to select available metadata, aiding in the identification of correlating scoring patterns within the pseudotime trajectory.

### Pseudotime Gene (germ only):

This tab enables users to type in a gene of interest, visualizing the expression wave across our defined pseudotime for germ cells. The plotted expression wave provides insights into how the gene's expression changes over the trajectory of pseudotime.

Furthermore, users can select metadata options to visualize any differential expression patterns associated with the specified gene across the pseudotime trajectory. This feature enhances the exploration of gene expression dynamics within the context of pseudotime in germ cells.

### Pseudotime Component Index:

This tab presents an index of the SDA components related to germ cells along with a summary of observations. A key feature of this table is the order assigned to these components relative to the cells, representing stages of spermatogenesis, where they score with the most magnitude.

The ranking is determined by identifying the main peak/maxima of the density curve fitted to the scores (y-axis) by pseudotime (x-axis) scatter using a peak-finding algorithm. The position of this maxima defines the rank. For example, the third-ranked SDA component, SDA #149, splits early spermatogonial cells into three sections, where the intermediate section is scored negatively opposite to the others. Further details are provided in the vignettes section.

### Enrichment analysis:

This tab performs hypergeometric testing (K=150) given a set of genes to identify components that are highly enriched with those genes. The analysis considers an adjusted p-value less than 0.01 for determining significant enrichments. However, the fold enrichment can also provide valuable insights even in the absence of statistical significance.

It's important to note that this statistical approach is most suitable for sets with less than ~30 genes, and it is particularly informative when dealing with more than 3 genes.

### Top loaded components:

This tab employs a statistical approach to identify key components, particularly suitable for larger sets of genes. The default input is loaded with antisense genes. The tab includes several figures:

1. **Histogram and Density Plot:**
   - Compares two distributions, depicting how many genes from the input set (e.g., antisense genes) are found across the top positive and negative loaded components.
   - The distribution is layered over a second distribution derived from an equal-length random sample of genes (computed dynamically).

2. **Ranking of Components:**
   - Components are ranked based on how many genes from the input set are found in the top loaded components, split by positive and negative loading direction.
   - This ranking provides insights into components that are more or less enriched than the average for the input gene set.

3. **Correlation Heatmap:**
   - Displays a heatmap of gene loading components filtered by the input gene.
   - Identifies structural similarity among the components in the defined gene space.
   - Each component is annotated with enrichment information for certain pathology and cell type on the heatmap's rows.

### DEgenes: 

The DEgenes tab enables users to explore the Differentially Expressed (DE) Gene tables derived from the original analysis. Key features include:

1. **Model Conditions Selection:**
   - The top drop-down menu allows users to select model conditions for comparison (e.g., control vs infertiles - CNTL_vs_INF1, CNTL_vs_INF2, CNTL_vs_KS, etc.).

2. **Unsupervised Clustering DE Selection:**
   - Users can choose DE genes derived from unsupervised clustering within subclusters of each cell type.
   - Selection is available through the second drop-down menu, allowing users to explore the DE landscape of transcripts within each cell type.
   
### lncRNAs:

The lncRNAs tab, described in this manuscript's long non-coding RNA (lncRNA) vignette, offers insights into lncRNAs. Key features include:

1. **Component Number Input:**
   - Users can explore lncRNAs by entering the component number of interest in the input box.

2. **Venn Diagram:**
   - At the top, a Venn diagram illustrates the overlap between Ensembl lncRNA annotated genes and all genes found in HISTA.

3. **Distribution Analysis:**
   - Similar to the 'Top Loaded Components' tab, the tab observes the distribution overlap of the number of lncRNA genes found in the top-loaded genes across components.
   - A comparison is made with a random equal-length set of genes.

4. **Component Sorting:**
   - Components are sorted and split by direction to identify lncRNA-enriched components.
   


### Soma only W. LN19:
          
The "Soma only W. LN19" tab provides support material for the KS manuscript 1, focusing on the new Klinefelter scRNA-Seq data from Leurentino et al. (2019). Key points include:

#### Brief Notes:
- **Germ Cell Absence:**
  - As with the existing KS donors, there were no germ cells in this new patient. Therefore, the analysis focuses exclusively on somatic cells.

- **Sertoli Cell Clusters:**
  - A single large cluster of Sertoli cells (SC) is identified, comprising at least 3 sub-clusters.
  - Two smaller subclusters with a smaller diameter are enriched with juvenile Sertoli cells (JUV SC), with only one enriched with SC from other adults.
  - The largest SC subcluster is mostly derived from LN19.

- **Leydig Cell Subtypes:**
  - Three major subtypes of Leydig cells (LC) are defined: progenitors (PLC), immature (ILC), and mature (MLC) groups.
  - LN19 LC are predominantly found in PLC and MLC clusters.

- **Distinct Cell Cluster:**
  - LN19 contributes a distinct cluster of cells believed to be of a myoid phenotype.

#### Data Integration:
- The new single KS donor data is integrated using customized Seurat pipelines.
- Differences correlated with the sampling procedure, age, and health of the donor are identified.
- The plan is to combine this data in later releases of HISTA.

#### Limitation:
- The current analysis is limited to validating results with the new data, and comprehensive integration is planned for future releases.

 
### LC only W. Zhao21:

The "Zhao et al. Validation" tab provides validation support for the KS manuscript1 based on the 2021 scRNA-Seq dataset by Zhao et al.

#### Data Overview:
- A new scRNA-Seq dataset by Zhao et al. was published in 2021.
- Future releases of HISTA plan to integrate this data; however, current limitations restrict full integration due to anomalies and failed QC processes.

#### Analysis Scope:
- The current release of HISTA is limited to validating results using processed count data from Zhao et al.'s dataset.
- Significant anomalies in raw sequencing files may limit future integration efforts.

#### Leydig Cell Findings Validation:
- Leydig cells (LC) are categorized into three major subtypes: progenitors (PLC), immature (ILC), and mature (MLC) groups.
- Validation focuses on confirming the presence of these clusters in an integrated dataset that combines LN19 and Zhao21 data.
- Integration involves cells from multiple donors, providing critical validation, especially for infrequent populations like MLCs under normal conditions.

#### Limitations:
- Challenges in downloading raw sequencing files and failed QC processes in many samples may impact future data integration efforts.
- The validation is based on processed count data, acknowledging the limitations imposed by these constraints.

 
## HISTA Case Examples:

### Case 1: Gene-Centric Exploration in HISTA

If you have one or more genes and want to learn more about them within HISTA, follow these steps:

1. **Main Tab:**
   - To examine the expression of each gene, search for the gene in the "Main" tab.

2. **2D Expression Plot (t-SNE or UMAP) in Main Tab:**
   - The figure header of the 2D expression plot in the "Main" tab provides information on SDA components and their loading.
   - This information helps identify gene sets that correlate with the gene of interest.

3. **Gene Expression Per Pathology Tab:**
   - Use the 'Gene Expression Per Pathology’ tab to visualize pathology-specific differential expression.
   - This tab allows you to specify a particular cell type or analyze expression across all cells.

4. **Gene Expression Per Cell Type Tab:**
   - Explore expression across cell types within specific pathologies using the 'Gene Expression Per Cell Type' tab.

5. **Pseudotime Expr Tab:**
   - If the gene is expressed in germ cells, check the 'Pseudotime Expr’ tab to observe expression patterns across pseudotime.

6. **Enrichment Analysis Tab:**
   - For more than ~3 genes, use the 'Enrichment Analysis’ tab to identify components enriched with these genes.
   - Search these components for annotations and additional correlating genes.

7. **Top Loaded Components Tab:**
   - If your gene set is large, explore the 'Top Loaded Components’ tab for a heatmap of component-component correlation.

These steps allow comprehensive exploration of gene expression patterns and relationships within HISTA.


### Case 2: Exploring Differential Expression Patterns in HISTA

Suppose you're interested in understanding the differential expression patterns across various pathologies and cell types in HISTA. Here's a step-by-step guide:

1. **Main Tab:**
   - Start in the "Main" tab to get an overview of available information.
   - Explore the 2D (t-SNE or UMAP) expression plot to identify potential genes of interest.

2. **Gene Expression Per Pathology Tab:**
   - Navigate to the 'Gene Expression Per Pathology’ tab.
   - Choose a specific cell type or analyze expression across all cells.
   - Utilize Wilcox statistics to visualize the differential expression of genes across different pathologies.

3. **Gene Expression Per Cell Type Tab:**
   - Move to the 'Gene Expression Per Cell Type' tab to focus on expression patterns within specific cell types.
   - This provides a detailed view of how genes vary across different cell types within selected pathologies.

4. **Pseudotime Expr Tab:**
   - If interested in germ cell expression patterns, check the 'Pseudotime Expr’ tab.
   - Observe how gene expression changes across pseudotime, gaining insights into developmental trajectories.

5. **Enrichment Analysis Tab:**
   - For a broader perspective, switch to the 'Enrichment Analysis’ tab.
   - Identify components enriched with specific gene sets, helping to understand their functional relevance.

6. **Top Loaded Components Tab:**
   - Explore the 'Top Loaded Components’ tab for a heatmap of component-component correlation.
   - This step is particularly useful when dealing with a larger set of genes, revealing relationships and potential co-regulation.

By following these steps, you can gain a comprehensive understanding of how genes are expressed and regulated across different conditions and cell types in HISTA.


### Case 3: Exploring Pseudotime Expression Dynamics in HISTA

Let's say you are interested in understanding how the expression of certain genes evolves across pseudotime, representing the progression of spermatogenesis. Here's how you can explore this within HISTA:
markdown
Copy code
* Start by searching for your genes of interest in the 'Main' tab to observe their expression patterns in the 2D (t-SNE or UMAP) projection.

* Head to the 'Pseudotime Expr' tab to visualize the expression pattern of these genes along the pseudotime trajectory.

* To gain a comprehensive view, consider using the 'Gene Expression Per Pathology' tab to evaluate differential expression across different pathologies.

* If you have a set of genes, the 'Enrichment Analysis' tab can help identify which SDA components are enriched with these genes, providing deeper insights into their functional relevance.

