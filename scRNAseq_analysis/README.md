## For all scripts in this folder:
* Seurat v4.1.0
* R ≥ v4.0

## 1. Process and integrate snRNA-seq datasets from Ruf-Zamojski et al. 2021.
This script takes .h5 matrix files from snRNA-seq data, performs data filtering, normalization, merging between replicates (of the same sex) and integration between sexes following Seurat workflows. Cell types are labelled based on gene markers identified in Ruf-Zamojski et al. 2021.

* Script is run from: `20220517_integrated_process_ruf_snrnaseq.R`
* Raw data (.h5 files) for this dataset can be downloaded here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151961
* Only "FrozPit-XXX" samples were used (snRNA-seq from snap-frozen nuclei)
* Place downloaded .h5 files in: `scRNAseq_analysis/input/`
  * More computational resources are required to find the data integration anchors, therefore, we ran this step on a high-performance computing cluster.
  * We imported the data integration anchor file locally and used it to continue with the integration and analysis.
* Output should be in `scRNAseq_analysis/output/`
  * `20211214_ruf_2021_pit_GSE151961_integrated_seurat.rds` generated from this script is the Seurat object is required for KS test, RNA-seq deconvolution and scMappR scripts (below), therefore, **move** this object to `scRNAseq_analysis/input/`.




## 2. KS test of co-expression module genes in snRNA-seq data
This script tests for the enrichment of each co-expression module gene expression within snRNA-seq data in a given cell type compared to all other cell types using a Kolgomorov-Smirnov (KS) test. A copy of this analysis has been pre-compiled as `ruf_integrated_kstest_2022-05-30_rmd_ver.html`

* Script is run from: `ruf_integrated_kstest_2022-05-30_rmd_ver.R`
* `20211214_ruf_2021_pit_GSE151961_integrated_seurat.rds` generated from 1. is required for this script.




## 3. Comparison of RNA-seq deconvolution methods
This script runs RNA-seq deconvolution for six methods (Proportions in Admixture (WGCNA), DCQ, DeconRNASeq, CPM, NNLS, MuSiC) through R. It also outputs the input files required for the other two methods, CIBERSORT and CIBERSORTx, which were run using the browser versions of the tools. The performance of each method is then determined by correlating each deconvolution-estimated cell-type proportions at PD37 compared to snRNA-seq-determined cell-type proportions for females and males separately.

* Script is run from: `ruf_integrate_decon_final_20220524.R`
* `scRNAseq_analysis/input/20211214_ruf_2021_pit_GSE151961_integrated_seurat.rds` generated from 1. is required for this script.



## 4. Using scMappR to identify sex-biased genes in pituitary cell types
This script runs scMappR to calculates cell-weighted fold change (cwFC) for bulk sex-biased genes. We consider genes with |cwFC| > 0.5 are sex-biased genes in the given pituitary cell type. Cell-type proportions are estimated using Proportions in Admixture RNA-seq deconvolution method ("WGCNA" in scMappR). scMappR is run for all profiled ages, however our paper only focuses on results from PD37 as this time point is closest to the snRNA-seq reference samples (adult mice).

* Script is run from: `scmappr_revision_20220530_final.R`
* `scRNAseq_analysis/input/20211214_ruf_2021_pit_GSE151961_integrated_seurat.rds` generated from 1. is required for this script.
* scMappR runtimes are generally long and vary depending on the number of DE genes inputted. 
* Output from this script are in `scRNAseq_analysis/output/ruf_scmappr/scMappR_revision_pituitary_dXX_sex`.




