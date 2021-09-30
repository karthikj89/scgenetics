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
    
def write_diseaseprogression_matrix(adata, filedir, filename, celltypelabel):
    # set up the ordering of genes and cells
    #adata = adatas[filename]
    genes = list(set(adata.var_names))
    gene2idx = {gene:i for i, gene in enumerate(genes)}
    
    pvalmtxs, logfoldmtxs, scoremtxs = [], [], []
    ctlabels = [col for col in adata.uns if '_DE' in col]
    print(adata.obs.columns)
    print(ctlabels)
    for delabel in ctlabels:
        #delabel = ctlabel + '_DE'
        ct = delabel.split('_')[0]
        contamination = adata.uns['contamination_'+celltypelabel].get(ct, np.array([])).tolist()
        
        cellsubsets = adata.uns[delabel]['names'].dtype.fields.keys()
        cell2idx = {cellsubset:i for i, cellsubset in enumerate(cellsubsets)}

        # create empty matrix
        pvalmtx = np.zeros((len(gene2idx), len(cell2idx)))

        logfoldmtx = np.zeros((len(gene2idx), len(cell2idx)))
        scoremtx = np.zeros((len(gene2idx), len(cell2idx)))

        # loop through and fill up the matrix with pvalue, logfold and score
        for gene, pval, logfold, score in zip(adata.uns[delabel]['names'], 
                                       adata.uns[delabel]['pvals_adj'], 
                                       adata.uns[delabel]['logfoldchanges'], 
                                       adata.uns[delabel]['scores']):
            for cell_subset in cellsubsets:
                
                if gene[cell_subset] in contamination:
                    p = 1
                    l = 0
                    s = 0
                else:
                    p = pval[cell_subset]
                    l = logfold[cell_subset]
                    s = score[cell_subset]
                
                if gene[cell_subset] in gene2idx:
                    pvalmtx[gene2idx[gene[cell_subset]], cell2idx[cell_subset]] = p
                    logfoldmtx[gene2idx[gene[cell_subset]], cell2idx[cell_subset]] = l
                    scoremtx[gene2idx[gene[cell_subset]], cell2idx[cell_subset]] = s

        # transform matrix to dataframe
        cellsubsets = [ct for ct in cellsubsets]
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
        sc.tl.rank_genes_groups(adata, 
                                ctlabel, 
                                key_added=ctlabel+'_DE', 
                                use_raw=False, 
                                method='wilcoxon', 
                                n_genes=n_genes)
        tissueadata.uns[ctlabel+'_DE'] = adata.uns[ctlabel+'_DE']
    write_celltypeprogram_matrix(tissueadata, 
                                 filedir, 
                                 tissue, 
                                 celltypelabel)

def compute_diseaseprogression_programs(filename, tissue, patientkey, celltypelabel, diagnosislabel, healthylabel, diseaselabel):
    diseaselabel_mapping = {healthylabel:"Healthy", diseaselabel:"Disease"}
    adata = sc.read(filename)
    for ctlabel in [celltypelabel]:
        subset = adata[adata.obs[diagnosislabel]==diseaselabel]
        sc.tl.rank_genes_groups(subset, 
                                groupby=ctlabel, 
                                reference='rest', 
                                n_genes=subset.shape[1], 
                                method='wilcoxon')
        
        adata.uns['contamination_'+ctlabel] = compute_contamination(subset, celltypelabel)
        adata.obs['DEstatus'] = [diseaselabel_mapping.get(diagnosis, 'Unknown') + '_' + ct for diagnosis, ct in zip(adata.obs[diagnosislabel], adata.obs[ctlabel])]
        destatus_counts = Counter(adata.obs['DEstatus'])
        discard = False
        for ct in set(adata.obs[ctlabel]):
            discard = False
            for k in ['Healthy_'+ct, 'Disease_'+ct]:
                if destatus_counts.get(k, 0) < 5:
                    discard = True
            print(ct, destatus_counts.get('Healthy_'+ct), destatus_counts.get('Disease_'+ct), discard)
            if discard:
                continue
            sc.tl.rank_genes_groups(adata, 
                                    groupby='DEstatus', 
                                    reference='Healthy_'+ct, 
                                    groups=['Disease_'+ct], 
                                    key_added=ct+'_DE', 
                                    n_genes=adata.shape[1], 
                                    method='wilcoxon',
                                   )
        write_diseaseprogression_matrix(adata, 
                                          filedir,
                                          tissue,
                                          celltypelabel)

def compute_contamination(adata, ctlabel):
    contamination = {}
    for ct in set(adata.obs[ctlabel]):
        scores = pd.DataFrame(adata.uns['rank_genes_groups']['scores'])[ct]
        names = pd.DataFrame(adata.uns['rank_genes_groups']['names'])[ct]
        threshold = np.mean(scores) - 6*np.std(scores)
        contamination_genes = names[scores<threshold]
        contamination[ct] = contamination_genes
    return contamination
    
    
if __name__=='__main__':
    programtype = sys.argv[1]
    filename = sys.argv[2]
    filedir = sys.argv[3]
      
    if programtype=='celltype':
        tissue, sampleid, celltypelabel = sys.argv[4]
        compute_celltype_programs(filename,
                                  tissue,
                                  sampleid,
                                  celltypelabel
                                 )
    elif programtype=='diseaseprogression':
        tissue, sampleid, celltypelabel, diagnosislabel, healthylabel, diseaselabel = sys.argv[4].split(',')
        compute_diseaseprogression_programs(filename, 
                                            tissue, 
                                            sampleid, 
                                            celltypelabel, 
                                            diagnosislabel, 
                                            healthylabel, 
                                            diseaselabel
                                           )