# HISTA Version History

## V2.9.7 (Current)

### New Features
- Introduced 'SDA Score per Cell Types' tab, enabling selection of SDA components and visualization through box plots relative to cell types. This complements the 'Gene Expression per Cell Types' tab for comprehensive investigation of expression patterns across populations.

- Added 'DEgenes' tab for searching differentially expressed genes in various modeling conditions, such as case-control and unsupervised clustering.

### Bug Fixes
- Addressed missing figure download button.

## V2.9.6

- Resolved gene-expression t-SNE download button issue.

## V2.9.5

### Bug Fixes
- Fixed t-SNE batch-removed imputed DGE loading issue.

### Updates
- Updated color set to remove ambiguous colors.
- Computed UMAP in the Main tab, making it selectable from the data origin.
- Added UMAP to various tabs, updating names from (t-SNE) to (2D).
- Enhanced interactivity in the Fingerprinting tab.
- Reorganized menu/tab items for improved grouping.

## V2.9.4

- Introduced 'Top Enriched Components with Specific Gene Sets' tab.

## V2.9.3

- Added 'lncRNA Analysis' tab for specific analysis related to long non-coding RNAs.

## V2.9.2

- Fixed pseudotime PDF download button.

## V2.9.1

- Resolved bug in 'Gene Expr per Cell Type (Boxplot)' tab related to metadata merging.

## V2.9

### Major Updates
- Improved computation of gene expression, filtering selected cells prior to dot-product of SDA matrices.
- Added cell cycle labels (G1, G2M, or S phase).
- Initiated version history documentation.

## V2.8

- Added 'GEX by Cell Type' tab for gene expression of a select gene with cell type filtering.

## V2.7

- Implemented filtering for removed components.
- Added 'Zhao 2020' tab for parallel analysis with Laurentino et al. 2020 donors.

## V2.6

- Introduced index tables tab housing metadata associated with SDA components and observations.
- Added 'Laurentino 2019' tab for parallel analysis with Laurentino et al. 2019 KS donor.

## V2.5

- Added 'Germ Cell Pseudotime' tab for germ cell-only pseudotime analysis.

## V2.0

- Modified for ShinyApps.io hosting.

## V1.5

- Converted to R package for easy deployment.
