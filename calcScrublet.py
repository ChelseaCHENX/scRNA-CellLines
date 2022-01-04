from utils import *

def calcScrublet(adata):

    adata.var_names_make_unique()  # deduplicate gene names
    print(adata)    

    # remove doublets
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.05)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                            min_cells=3, 
                                                            min_gene_variability_pctl=85, 
                                                            n_prin_comps=30)

    scrub.plot_histogram()
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    scrub.plot_embedding('UMAP', order_points=True);
    print('%d scores for %d samples'%(len(predicted_doublets), adata.X.shape[0]))

    df = pd.DataFrame({'doublet_bln':predicted_doublets, 'doublet_score': doublet_scores})
    df.index = adata.obs.index
    nobs = adata.obs.join(df, how='inner') 
    adata.obs = nobs

    # calc QCs
    mito_genes = adata.var_names.str.startswith('MT-')
    adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    adata.obs['n_genes'] = np.sum(adata.X > 0, axis=1)

    return adata

if __name__ == "__main__":
    os.chdir('/bgfs/alee/chelsea/projects/10X/CellLine/codes')
    # wkdir = '/bgfs/alee/chelsea/projects/PublicData/BroadMCF72018/cellranger-count-v3.0.2'
    # sampleMap = np.loadtxt('../data/BroadMCF7/meta.txt', dtype=object)
    # samples = sampleMap[:,0]

    # for sp in samples:   
    #     path='{}/{}/{}/outs/filtered_feature_bc_matrix'.format(wkdir, sp, sp)
    #     adata = sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)   
    #     adata_raw = calcScrublet(adata)
    #     adata_raw.write('../data/BroadMCF7/%s.raw_tmp.h5ad'%sp)
    samples = [x.strip('.loom') for x in os.listdir('/bgfs/alee/chelsea/projects/10X/preAdapt/data/loom') if not x.startswith('SRR')]
    for n, sp in enumerate(samples):
        tmp = sc.read('/bgfs/alee/chelsea/projects/10X/preAdapt/data/loom/%s.loom'%sp, cache=True)   
        tmp.var_names_make_unique()
        tmp_raw = calcScrublet(tmp)
        tmp_raw.write('../data/Adapt/%s.raw_tmp.h5ad'%sp)
        

 
