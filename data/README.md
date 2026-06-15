<h1>data</h1>

The notebooks expect raw UMI count matrices under this folder. The data is not redistributed here; download it from the sources below and place it as shown. All four are UMI datasets, which is what molecular cross-validation requires.

| Key | Source | File(s) used | Expected path |
|---|---|---|---|
| Paul 2015 (MARS-seq) | GEO GSE72857 | <code>GSE72857_umitab.txt</code> | <code>data/GSE72857_umitab.txt</code> |
| PBMC 3k (10x v1) | 10x Genomics public dataset | <code>matrix.mtx</code>, <code>genes.tsv</code>, <code>barcodes.tsv</code> | <code>data/filtered_gene_bc_matrices/hg19/</code> |
| Tabula Muris marrow (10x droplet) | figshare droplet v2 / GEO GSE109774 | <code>matrix.mtx</code>, <code>genes.tsv</code>, <code>barcodes.tsv</code> | <code>data/droplet/Marrow-10X_P7_2/</code> |
| Tusi 2018 (inDrop) | GEO GSE89754 | <code>GSM2388072_basal_bone_marrow.raw_umifm_counts.csv</code> | <code>data/GSM2388072_basal_bone_marrow.raw_umifm_counts.csv</code> |

Smart-seq2 and read-count datasets (for example Nestorowa 2016, GSE81682, or the Tabula Muris FACS tarballs) are not used because the Poisson thinning step in MCV is only valid for UMI counts.
