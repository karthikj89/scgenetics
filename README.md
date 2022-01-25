# scgenetics
identifying disease critical cell types and programs from single cell RNAseq


Below we describe the steps to construct gene programs from scRNA-seq for heritaility analysis:
* cell type programs

```
bash <srcdir>/generateGenePrograms.sh <srcdir> celltype <datapath> <outdir> <tissue>,<samplekey>,<celltypekey>
```

* disease progression programs

```
bash <srcdir>/generateGenePrograms.sh <srcdir> diseaseprogression <datapath> <outdir> <tissue>,<samplekey>,<celltypekey>,<diseasestatuskey>,<healthylabel>,<diseaselabel>
```

* cellular processes programs
