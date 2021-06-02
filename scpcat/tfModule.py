import numpy as np
import scanpy as sc
import progeny
import dorothea
import os

def addMotif(args):
    adata=sc.read_h5ad(args.path)
    if "X_umap" not in adata.obsm.keys() or "X_tsne" not in adata.obsm.keys():
        raise ValueError("Please Run UMAP or tSNE first!")
    
    assert args.organism in  ["Human","Mouse"]
    
    print("Run progeny method to predict pathway activities ")
    model = progeny.load_model(
           organism=args.organism, # If working with mouse, set to Mouse
           top=args.ntop          # For sc we recommend ~1k target genes since there are dropouts
           )
    progeny.run(adata,        # Data to use
            model,        # PROGENy network
            center=True,  # Center gene expression by mean per cell
            num_perm=100, # Simulate m random activities
            norm=True,    # Normalize by number of edges to correct for large regulons
            scale=True,   # Scale values per feature so that values can be compared across cells
            use_raw=True, # Use raw adata, where we have the lognorm gene expression
            use_hvg=False, # Only use high variable genes for pathway estimation
            min_size=5    # Pathways with less than 5 targets will be ignored
           )
   
    print("Run dorothea method to predict TF activities")
    reglevel=["A","B","C"]
    regulons = dorothea.load_regulons(
              reglevel,   # Which levels of confidence to use (A most confident, E least confident)
              organism=args.organism # If working with mouse, set to Mouse
            )
    dorothea.run(adata,        # Data to use
             regulons,     # Dorothea network
             center=True,  # Center gene expression by mean per cell
             num_perm=100, # Simulate m random activities
             norm=True,    # Normalize by number of edges to correct for large regulons
             scale=True,   # Scale values per feature so that values can be compared across cells
             use_raw=True, # Use raw adata, where we have the lognorm gene expression
             use_hvg=False, # Only use high variable genes for TF estimation
             min_size=5,   # TF with less than 5 targets will be ignored
            )
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    dX=dorothea.extract(adata)
    pX=progeny.extract(adata)
    dX.write(args.outdir+"/"+"dorothea.h5ad")    
    pX.write(args.outdir+"/"+"progeny.h5ad")

