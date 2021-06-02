import scanpy as sc
import os
from .scorect_api import *

def classify(args):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print("# Loading data from path : {}".format(args.path))
    adata=sc.read_h5ad(args.path)

    print("# Loading reference DB from : {}".format(args.db))
    ref_marker = read_markers_from_file(args.db)
    
    print("# rank genes group in {}".format(args.groupby))
    sc.tl.rank_genes_groups(adata,groupby=args.groupby,method='t-test',n_genes=1000,corr_method="bonferroni")
    marker_df = wrangle_ranks_from_anndata(adata)

    # Score cell types for each cluster
    # Let's set parameters first - K represents the number of genes included in the ranking
    # m represents the number of bins used to divide the top K genes.
    K = 300
    m = 5
    # Get the background genes - here, all the genes used to run the differential gene expression test
    background = adata.raw.var.index.tolist()
    # Now run the function
    print("# Run celltype score")
    ct_pval, ct_score = celltype_scores(nb_bins=m,
                                        ranked_genes=marker_df,
                                        K_top = K,
                                        marker_ref=ref_marker,
                                        background_genes=background)


    print("# Run assign celltype")
    # Now assign clusters to cell types
    cluster_assign = adata.obs[args.groupby]
    celltype_assign = assign_celltypes(cluster_assignment=cluster_assign, ct_pval_df=ct_pval, ct_score_df=ct_score, cutoff=0.1)
    # Add to anndata object
    adata.obs['predicted'] = celltype_assign
    # Let's compare with the true assignment now!
    adata.write(os.path.join(args.outdir,"adata.h5ad"),compression='gzip')
    print("Done!")





