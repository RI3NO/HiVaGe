def HiVaGePY(ds, flavour=None, HVGs_num=None):
    # Flavor-specific implementations
    # Dependency: pandas, anndata, scanpy, triku, numpy

    import pandas as pd
    import anndata as ad
    import scanpy as sc
    import triku as tk
    from numpy import inf
    
    cell_num = ds.shape[1]
    counts = ds.iloc[:, 1:cell_num].values
    counts = counts.T

    gene_names = ds['ProbeIDs']
    cell_names = ds.columns[1:cell_num]
    anndata = ad.AnnData(X=counts)

    anndata.var_names = gene_names
    anndata.obs_names = cell_names

    HVGs_num = int(HVGs_num)
    if flavour == "scanpy_seurat":
        sc.pp.log1p(anndata)
        sc.pp.highly_variable_genes(anndata, layer=None, n_top_genes=HVGs_num, min_disp=0.5, max_disp=inf, min_mean=0.0125,
                                    max_mean=3, span=0.3, n_bins=20, flavor='seurat', subset=False, inplace=True,
                                    batch_key=None, check_values=True)
        return pd.DataFrame(anndata.var[anndata.var['highly_variable']]['highly_variable'].index)
    elif flavour == "scanpy_cell_ranger":
        # Implement Cell Ranger flavor here
        sc.pp.log1p(anndata)
        sc.pp.highly_variable_genes(anndata, layer=None, n_top_genes=HVGs_num, min_disp=0.5, max_disp=inf,
                                    min_mean=0.0125,
                                    max_mean=3, span=0.3, n_bins=20, flavor='cell_ranger', subset=False, inplace=True,
                                    batch_key=None, check_values=True)
        return pd.DataFrame(anndata.var[anndata.var['highly_variable']]['highly_variable'].index)
    elif flavour == "scanpy_seurat_v3":
        # Implement Seurat v3 flavor here
        anndata_s3 = ad.AnnData(X=counts)

        anndata_s3.var_names = gene_names
        anndata_s3.obs_names = cell_names

        sc.pp.highly_variable_genes(anndata_s3, n_top_genes=HVGs_num, span=0.3, n_bins=20, flavor='seurat_v3', subset=True,
                                    inplace=True, batch_key=None, check_values=True)
        return pd.DataFrame(anndata_s3.var[anndata_s3.var['highly_variable']]['highly_variable'].index)

    elif flavour == "scanpy_Pearson":
        # Implement Pearson flavor here
        sc.experimental.pp.highly_variable_genes(anndata, theta=100, clip=None, n_top_genes=HVGs_num, batch_key=None,
                                                 chunksize=1000, flavor='pearson_residuals', check_values=True,
                                                 layer=None, subset=False, inplace=True)
        return pd.DataFrame(anndata.var[anndata.var['highly_variable']]['highly_variable'].index)


    elif flavour == "Triku":
        # Implement Triku flavor here
        sc.pp.pca(anndata)
        sc.pp.neighbors(anndata, metric='cosine', n_neighbors=int(0.5 * len(anndata) ** 0.5))
        tk.tl.triku(anndata)
        return anndata.var
        """pd.DataFrame(anndata.var[anndata.var["triku_highly_variable"]].index)"""



    else:
        raise ValueError("Invalid flavor. Available flavors are 'scanpy_seurat', 'scanpy_cell_ranger', 'scanpy_seurat_v3', 'scanpy_Pearson', and 'Triku'.")
