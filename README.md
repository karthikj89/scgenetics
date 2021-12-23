# scgenetics
identifying disease critical cell types and programs from single cell RNAseq


Below we describe the steps to construct gene programs from scRNA-seq for heritaility analysis:
* cell type programs
* disease progression programs

```
bash <srcdir>/generateGenePrograms.sh <srcdir> diseaseprogression <datapath (ends with .h5ad)> <outdir> <tissue>,<samplekey>,<celltypekey>,<diseasestatuskey>,<healthylabel>,<diseaselabel>
```

* cellular processes programs

To construct gene programs from single cell data for heritability analysis you can perform
