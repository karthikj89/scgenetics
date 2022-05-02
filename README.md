# identifying disease critical cell types and programs from single cell RNAseq

Genome-wide association studies (GWAS) provide a powerful means to identify loci and genes contributing to disease, but in many cases the related cell types/states through which genes confer disease risk remain unknown. Deciphering such relationships is important for identifying pathogenic processes and developing therapeutics. Here, we introduce sc-linker, a framework for integrating single-cell RNA-seq (scRNA-seq), epigenomic maps and GWAS summary statistics to infer the underlying cell types and processes by which genetic variants influence disease. The inferred disease enrichments recapitulated known biology and highlighted novel cell-disease relationships, including GABAergic neurons in major depressive disorder (MDD), a disease-dependent M cell program in ulcerative colitis, and a disease-specific complement cascade process in multiple sclerosis. In autoimmune disease, both healthy and disease-dependent immune cell type programs were associated, whereas disease-dependent epithelial cell programs were most prominent, perhaps suggesting a role in disease response over initiation. Our framework provides a powerful approach for identifying the cell types and cellular processes by which genetic variants influence disease.

Here, we present the code required to run this analysis on a new single cell dataset.

Our manuscript (in press at Nature Genetics) will be available soon. The preprint can be found [here](https://www.biorxiv.org/content/10.1101/2021.03.19.436212v2). Raw and processed data can be found [online](https://alkesgroup.broadinstitute.org/LDSCORE/Jagadeesh_Dey_sclinker/scdata/).

## Data analysis
### generating gene programs from scRNA-seq for heritability analysis
* cell type programs

```
bash <srcdir>/generateGenePrograms.sh <srcdir> celltype <datapath> <outdir> <tissue>,<samplekey>,<celltypekey>
```

* disease progression programs

```
bash <srcdir>/generateGenePrograms.sh <srcdir> diseaseprogression <datapath> <outdir> <tissue>,<samplekey>,<celltypekey>,<diseasestatuskey>,<healthylabel>,<diseaselabel>
```

* cellular processes programs

### transforming gene programs to SNP annotations using epigenomic data

### measuring the heritability enrichment in human disease GWAS
