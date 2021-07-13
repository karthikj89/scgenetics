import scanpy as sc
import numpy as np
import pandas as pd
import sys
from collections import Counter

def write_celltypeprogram_matrix(adata, filedir, filename, celltypelabel):
    # set up the ordering of genes and cells
    genes = list(set(adata.var_names))
    gene2idx = {gene:i for i, gene in enumerate(genes)}
    
    pvalmtxs, logfoldmtxs, scoremtxs = [], [], []
    
    ctlabels = [celltypelabel]
    print(adata.obs.columns)
    print(ctlabels)
    for ctlabel in ctlabels:
        delabel = ctlabel + '_DE'
        cellsubsets = adata.uns[delabel]['names'].dtype.fields.keys()
        cell2idx = {cellsubset:i for i, cellsubset in enumerate(cellsubsets)}

        # create empty matrix
        pvalmtx = np.zeros((len(gene2idx), len(cell2idx)))

        logfoldmtx = np.zeros((len(gene2idx), len(cell2idx)))
        scoremtx = np.zeros((len(gene2idx), len(cell2idx)))

        # loop through and add the matrix with pvalue, logfold and score
        for gene, pval, logfold, score in zip(adata.uns[delabel]['names'], 
                                       adata.uns[delabel]['pvals_adj'], 
                                       adata.uns[delabel]['logfoldchanges'], 
                                       adata.uns[delabel]['scores']):
            for cell_subset in cellsubsets:
                if gene[cell_subset] in gene2idx:
                    pvalmtx[gene2idx[gene[cell_subset]], cell2idx[cell_subset]] = pval[cell_subset]
                    logfoldmtx[gene2idx[gene[cell_subset]], cell2idx[cell_subset]] = logfold[cell_subset]
                    scoremtx[gene2idx[gene[cell_subset]], cell2idx[cell_subset]] = score[cell_subset]

        # transform matrix to dataframe
        pvalmtxs.append(pd.DataFrame(pvalmtx, index=genes, columns=cellsubsets))
        logfoldmtxs.append(pd.DataFrame(logfoldmtx, index=genes, columns=cellsubsets))
        scoremtxs.append(pd.DataFrame(scoremtx, index=genes, columns=cellsubsets))
    pvalmtxs = pd.concat(pvalmtxs, axis=1)
    logfoldmtxs = pd.concat(logfoldmtxs, axis=1)
    scoremtxs = pd.concat(scoremtxs, axis=1)

    # write matrix to file
    pvalmtxs.to_csv("%s/%s_pval.csv"%(filedir, filename))
    logfoldmtxs.to_csv("%s/%s_logfold.csv"%(filedir, filename))
    scoremtxs.to_csv("%s/%s_score.csv"%(filedir, filename))
    
def compute_celltype_programs(filename, tissue, sampleid, celltypelabel):
    tissueadata = sc.read(filename)
    for ctlabel in [celltypelabel]:
        print(ctlabel)
        counts = Counter(tissueadata.obs[ctlabel])
        tissueadata.obs[ctlabel+'_counts'] = [counts[ct] for ct in tissueadata.obs[ctlabel]]
        adata = tissueadata[tissueadata.obs[ctlabel+'_counts'] > 10].copy()
        n_genes = adata.shape[1]
        sc.tl.rank_genes_groups(adata, ctlabel, key_added=ctlabel+'_DE', use_raw=False, method='wilcoxon', n_genes=n_genes)
        tissueadata.uns[ctlabel+'_DE'] = adata.uns[ctlabel+'_DE']
    write_celltypeprogram_matrix(tissueadata, filedir, tissue, celltypelabel)
    
    
if __name__=='__main__':
    filename = sys.argv[1]
    filedir = sys.argv[2]
    tissue = sys.argv[3]
    sampleid = sys.argv[4]
    celltypelabel = sys.argv[5]
    
    compute_celltype_programs(filename, tissue, sampleid, celltypelabel)