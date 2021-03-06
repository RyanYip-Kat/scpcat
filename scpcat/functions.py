import scanpy as sc
import numpy as np
import pandas as pd
import scanyuan as scy
import scrublet as scr
import os
import harmonypy as hm

from .tools import *
from .plotting import *
from .r_trasformer import Seurat2AnnData 

def _createObj(args):
  mtx_path=args.path
  if not os.path.exists(args.outdir):
     os.makedirs(args.outdir)
  
  print("### Loading data")
  adata=sc.read_10x_mtx(mtx_path,cache=True)
  cells=adata.obs_names.to_list()
  idents=[cell.split("-")[1] for cell in cells]
  adata.obs["idents"]=idents
  
  
  print("### Preprocess")
  adata.obs['n_counts'] =adata.X.sum(axis=1).A1
  adata.obs["n_genes"]=np.sum(adata.X>0,axis=1).A1
  adata.var['mt'] = adata.var_names.str.startswith('MT-')
  adata.var['rpl'] = adata.var_names.str.startswith('RPL')
  adata.var['rps'] = adata.var_names.str.startswith('RPS')
  sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rpl','rps'], percent_top=None, log1p=False, inplace=True)
  adata.raw=adata
  adata.write(args.outdir+"/"+"adata.h5ad")

def _preprocess(adata):
  sc.pp.normalize_total(adata, target_sum=1e4)
  sc.pp.log1p(adata)
  sc.pp.highly_variable_genes(adata,n_top_genes=2500)
  adata = adata[:, adata.var.highly_variable]
  sc.pp.scale(adata, max_value=10)
  return adata
  

def RunHarmony(adata,batch_key="status",max_iter=20):
    data=adata.copy()
    X=adata.X
    print("# The X matrix for harmony shape is : [{},{}]".format(X.shape[0],X.shape[1]))
    meta_data=data.obs
    assert batch_key in meta_data.columns
    vars_use=[batch_key]
    print("# Run harmony")
    ho = hm.run_harmony(X, meta_data, vars_use,max_iter_kmeans=max_iter)
    X_harmony=ho.Z_corr.transpose()

    print("# Add harmony ***")
    data.uns["harmony"]={}
    data.uns["harmony"]["params"]=ho.__dict__
    data.obsm["X_harmony"]=X_harmony
    return data


def _reduction(args):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print("")
    sc.settings.figdir=args.outdir
    sc.settings.autoshow=False
    sc.settings.autosave=True
    sc.settings.cachedir='./cache'

    sc.settings.cache_compression="gzip"
    sc.settings.n_jobs=args.n_jobs
    sc.settings.file_format_figs="pdf"

    use_rep="X_pca"
    print("# Laoding data from path : {}".format(args.path))
    adata=sc.read_h5ad(args.path)
    adata=_preprocess(adata)
    sc.tl.pca(adata,n_comps=50,svd_solver='arpack')
    
    if args.batch_correct:
      if args.batch_key is None:
         raise ValueError("Batch Key can not be NULL")
       
      adata=RunHarmony(adata,batch_key=args.batch_key)
      use_rep="X_harmony"
    
    sc.pp.neighbors(adata, n_neighbors=15, use_rep=use_rep,knn=True)
    print("# clustering...")
    sc.tl.leiden(adata,resolution=1.0)

    print("# run tsne...")
    sc.tl.tsne(adata,use_rep=use_rep)
    print("# run umap...")
    sc.tl.paga(adata)  # remove `plot=False` if you want to see the coarse-grained graph
    sc.pl.paga(adata, plot=False)
    sc.tl.umap(adata, init_pos='paga')
    adata.write(args.outdir+"/"+"adata.h5ad")


def _FindMarkers(args):
  if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

  print("# Laoding data from path : {}".format(args.path))
  adata=sc.read_h5ad(args.path)
  if args.subset1 is not None and args.column1 is not None:
    adata=subset_by_column(adata,args.subset1,args.column1)

    if args.subset2 is not None and args.column2 is not None:
        adata=subset_by_column(adata,args.subset2,args.column2)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

  if args.splitby is not None:
     labels=[x for x in adata.obs[args.splitby].unique()]
     for celltype in labels:
         c=[]
         c.append(celltype)
         print("Get from : {}".format(celltype))
         data=subset_by_column(adata,c,args.splitby)
         sc.pp.normalize_total(data, target_sum=1e4)

         sc.tl.rank_genes_groups(data,groupby=args.groupby,method=args.method,n_genes=500,corr_method="bonferroni")
         #sc.tl.filter_rank_genes_groups(data,groupby=args.groupby,key_added='rank_genes_groups_filtered',min_fold_change=0.25)
         out=os.path.join(args.outdir,celltype)

         if not os.path.exists(out):
             os.makedirs(out)
         get_rank_group_genes(data,pval=None,logfc=None,outdir=out,nTop=50)
         get_rank_group_genes(data,pval=None,logfc=None,outdir=out)

  else:
      sc.tl.rank_genes_groups(adata,groupby=args.groupby,method=args.method,n_genes=500,corr_method="bonferroni")
      get_rank_group_genes(adata,pval=None,logfc=None,outdir=args.outdir,nTop=50)
      get_rank_group_genes(adata,pval=None,logfc=None,outdir=args.outdir)
      
def _dotplot(args):
  if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
  
  sc.settings.figdir=args.outdir
  sc.settings.set_figure_params(figsize=[args.width,args.height])
  print("# Loading data from path : {}".format(args.path))
  adata=sc.read_h5ad(args.path)
  if args.subset1 is not None and args.column1 is not None:
    adata=subset_by_column(adata,args.subset1,args.column1)
    
    if args.subset2 is not None and args.column2 is not None:
        adata=subset_by_column(adata,args.subset2,args.column2)
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
  if args.level is not None:
    adata.obs[args.groupby]=adata.obs[args.groupby].cat.reorder_categories(args.level)
  
  df=pd.read_csv(args.genelist,sep="\t",header=None)
  marker_genes=df.loc[:,0].values.tolist()
  sc.pl.dotplot(adata,var_names=marker_genes,groupby=args.groupby,color_map=args.color,dendrogram=False,
        show=False,save=True,dot_min=0,dot_max=1,smallest_dot=40,standard_scale='var') # smallest_dot=40,[16,8]

  #scy.stacked_violin_t(adata, marker_genes, groupby=args.groupby,show=False,save="_stacked_violin_t",stripplot=False,jitter=False)
  #sc.pl.stacked_violin(adata,var_names=marker_genes,groupby=args.groupby,show=False,save=True,stripplot=True,jitter=False)



def _EmbPlot(args):
  sc.settings.autoshow=False
  if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
  sc.settings.set_figure_params(frameon=False,dpi_save=600,figsize=[args.width,args.height]) 
  sc.settings.figdir=args.outdir
  print("# Loading data from path : {}".format(args.path))
  adata=sc.read_h5ad(args.path)
  if args.subset1 is not None and args.column1 is not None:
    adata=subset_by_column(adata,args.subset1,args.column1)
    
    if args.subset2 is not None and args.column2 is not None:
        adata=subset_by_column(adata,args.subset2,args.column2)
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
  
  if args.groupby is not None:
      if args.level is not None:
          adata.obs[args.groupby]=adata.obs[args.groupby].cat.reorder_categories(args.level)
      sc.pl.embedding(adata=adata,basis=args.embedding,color=args.groupby,palette=args.theme,size=28,show=False,save="_"+args.groupby+".pdf")
  
  if args.genelist is not None:
    df=pd.read_csv(args.genelist,sep="\t",header=None)
    genes=df.loc[:,0].values.tolist()
    for gene in genes:
        sc.pl.embedding(adata=adata,basis=args.embedding,color=gene,gene_symbols=gene,color_map="Oranges",size=28,show=False,save="_"+gene+".pdf")

  
  
def _DoHeatmap(args):
  if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

  sc.settings.figdir=args.outdir
  #sc.settings.set_figure_params(figsize=[args.width,args.height])
  print("# Loading data from path : {}".format(args.path))
  adata=sc.read_h5ad(args.path)
  if args.subset1 is not None and args.column1 is not None:
    adata=subset_by_column(adata,args.subset1,args.column1)
    
    if args.subset2 is not None and args.column2 is not None:
        adata=subset_by_column(adata,args.subset2,args.column2)
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
  if args.level is not None:
    adata.obs[args.groupby]=adata.obs[args.groupby].cat.reorder_categories(args.level)
  
  sc.tl.rank_genes_groups(adata,groupby=args.groupby,method='t-test',n_genes=500,corr_method="bonferroni")
  sc.tl.filter_rank_genes_groups(adata,groupby=args.groupby,min_fold_change=0.1) # 0.1 is good
  
  sc.pl.rank_genes_groups_heatmap(adata,n_genes=args.n_genes,groupby=args.groupby,
        swap_axes=True,vmin=-3, vmax=3,key="rank_genes_groups",
        show=False,show_gene_labels=True,save="_rank_genes_groups",figsize=[24,16],dendrogram=False)

  sc.pl.rank_genes_groups_heatmap(adata,n_genes=args.n_genes,groupby=args.groupby,
        swap_axes=True,vmin=-3, vmax=3,key="rank_genes_groups_filtered",
        show=False,show_gene_labels=True,save="_rank_genes_groups_filtered",figsize=[24,16],dendrogram=False)

  #sc.pl.rank_genes_groups_stacked_violin(adata,n_genes=args.n_genes,stripplot=False,
  #      key="rank_genes_groups_filtered",groupby=args.groupby,swap_axes=True,
  #      figsize=[16,24],dendrogram=False,show=False,save="_rank_genes_groups_filtered")
  if args.genelist is not None:
     df=pd.read_csv(args.genelist,sep="\t",header=None)
     genes=df.loc[:,0].values.tolist()
     sc.pl.heatmap(adata,var_names=genes,groupby=args.groupby,
         swap_axes=True,show=False,show_gene_labels=True,cmap="viridis",
         vmin=-2.5, vmax=2.5,save="_genes_groups_viridis",figsize=[24,24])

   
def export_subset(args):
  if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

  print("# Laoding data from path : {}".format(args.path))
  adata=sc.read_h5ad(args.path)
  print("Before subset,data size id ({},{})".format(str(adata.n_obs),str(adata.n_vars)))
  if args.subset is not None and args.column is not None:
    adata=subset_by_column(adata,args.subset,args.column,invert=args.invert)
    print("After  subset,data size id ({},{})".format(str(adata.n_obs),str(adata.n_vars)))
  if args.bclist is not None:
    df=pd.read_csv(args.bclist,sep="\t",header=None)
    cells=df.loc[:,0].values.tolist()
    adata=subset_by_cell_feature(adata,cells=cells,invert=args.invert)
    print("After  subset,data size id ({},{})".format(str(adata.n_obs),str(adata.n_vars)))
  
  adata.write(args.outdir+"/"+"adata.h5ad")

def VlnPlot(args):
  if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

  sc.settings.figdir=args.outdir
  sc.set_figure_params(figsize=[args.width,args.height])
  print("# Laoding data from path : {}".format(args.path))
  adata=sc.read_h5ad(args.path)
  if args.subset1 is not None and args.column1 is not None:
    adata=subset_by_column(adata,args.subset1,args.column1)

    if args.subset2 is not None and args.column2 is not None:
        adata=subset_by_column(adata,args.subset2,args.column2)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

  if args.level is not None:
    adata.obs[args.groupby]=adata.obs[args.groupby].cat.reorder_categories(args.level)

  df=pd.read_csv(args.genelist,sep="\t",header=None)
  genes=df.loc[:,0].values.tolist()
  for gene in genes:
    if args.hue is not None:
        violin_hue(adata,groupby=args.groupby,
            gene_name=gene,
            hue=args.hue,
            hue_order=args.hue_order, #["YA","AA"],
            use_raw=True,
            figdir=args.outdir)
    diff,pval=gene_comparison(adata,gene=gene,by=args.groupby)
    ylabel=gene+"( "+str(pval)+" )"

    violin_sns(adata,gene_name=gene,save=True,groupby=args.groupby,figdir=args.outdir,orient=args.orient,strip=True)
    sc.pl.violin(adata, keys=gene, groupby=args.groupby, stripplot=False,ylabel=ylabel, inner='box',show=False,save="_innerBox_"+gene,rotation=45)
    sc.pl.violin(adata, keys=gene, groupby=args.groupby, stripplot=False, ylabel=ylabel,show=False,save="_"+gene,rotation=45)

  
def enrichGO(args):
  if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

  print("# Laoding data from path : {}".format(args.path))
  adata=sc.read_h5ad(args.path)
  
  if args.genelist is not None:
    df=pd.read_csv(args.genelist,sep="\t",header=None)
    genes=df.loc[:,0].values.tolist()
    
    m=sc.queries.enrich(genes,org=args.org)
    filename=os.path.join(args.outdir,"enrich_genes.csv")
    m.to_csv(filename,sep=",",index=False)
    
  else:
    sc.tl.rank_genes_groups(adata,groupby=args.groupby,method="t-test",n_genes=500,corr_method="bonferroni")
    groups=np.unique(adata.obs[args.groupby])
    tables=[]
    for g in groups:
        df=sc.queries.enrich(adata,group=str(g),org=args.org)
        df["cluster"]=g
        tables.append(df)
        m=pd.concat(tables,axis=0)
        filename=os.path.join(args.outdir,args.groupby+"_enrich.csv")
        m.to_csv(filename,sep=",",index=False)

def transformer(args):
    Seurat2AnnData(filename=args.filename,name=args.name,outdir=args.outdir)


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

def _extract(args):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    print("# Laoding data from path : {}".format(args.path))
    adata=sc.read_h5ad(args.path)

    if args.subset1 is not None and args.column1 is not None:
        adata=subset_by_column(adata,args.subset1,args.column1)
        if args.subset2 is not None and args.column2 is not None:
            adata=subset_by_column(adata,args.subset2,args.column2)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    metadata=adata.obs.copy()
    metadata["barcode"]=metadata.index
    for emb in args.emblist:
        embk="X_"+emb
        try:
            print("extract {} embedding".format(emb))
            mat=adata.obsm[embk]
            columns=[emb.upper()+"_"+str(i+1) for i in range(mat.shape[1])]
            mat=pd.DataFrame(mat,columns=columns,index=adata.obs_names)
            mat["Barcode"]=mat.index
            mat=mat.iloc[:,[2,0,1]]
            mat.to_csv(os.path.join(args.outdir,emb+"_projection.csv"),sep=",",index=False)
        except:
            print("invalid {} embedding slot,please check again!".format(emb))

    for label in args.labels:
        try:
            print("extract {} label".format(label))
            df=metadata[["barcode",label]]
            df.to_csv(os.path.join(args.outdir,label+"_cluster.csv"),sep=",",index=False)
        except:
            print("invalid {} column in metadata,please check again!".format(label))
