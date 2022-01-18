## To run scMappR scripts:
Ensure that you are using Seurat v3.2.3 and R ver ≥ 4.0.
* The script relies on a single-cell dataset reprocessed from Cheung et al. 2018.
Raw data for this dataset can be downloaded here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120410
* Place raw data in a folder: `scRNAseq_analysis/input_files/GSE120410_Cheung_2018/`
* Run `cheung_sc_reprocess_2021-02-17_rmd_ver.Rmd`
  * This should produce the following file: `scRNAseq_analysis/output_files/seurat_v3_cheung_etal_manual_cell_anno.rds` which is required for the scMappR scripts.
