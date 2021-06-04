import scanpy as sc
import scrublet as scr
import pandas as pd
import numpy as np
import os

def detectDoublet(args):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print("# Laoding data from path : {}".format(args.path))
    adata=sc.read_h5ad(args.path)

    print("# detect doublet ..")
    counts_matrix=adata.raw.X
    scrub = scr.Scrublet(counts_matrix)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    cellDict={True:"singlecell",False:"doublet"}

    #df=pd.DataFrame({"doubletScores":doublet_scores,"Doublets":predicted_doublets},columns=["doubletScores","Doublets"],index=adata.obs_names)
    adata.obs["doubletScores"]=doublet_scores
    adata.obs["Doublets"]=np.array([cellDict[x] for x in predicted_doublets])
    adata.write(args.outdir+"/"+"adata.h5da")
