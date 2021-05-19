import scanpy as sc
import os
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.sparse import issparse
from matplotlib import pyplot as pl
from anndata import AnnData
from typing import Union, Optional, List, Sequence, Iterable
from matplotlib.colors import Colormap
from .tools import gene_comparison


sns.set(style="whitegrid", palette="pastel", color_codes=True)
def violin_hue(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    gene_name: Optional[Iterable[str]] = None,
    use_raw: Optional[bool] = None,
    hue : str="status",
    hue_order: Optional[str] = None,
    figdir:Optional[str] ="./figdir",
    split: bool = True,
    scale: str = "area",#'width',
    strip: bool = True,
    jitter: Union[int, float, bool] = False,
    size: int = 1,
    save: Optional[bool] =True,
):
    if adata.raw is not None and use_raw:
        X= adata.raw[:, gene_name].X
                
    else:
        X= adata[:, gene_name].X
    if issparse(X):
        X = X.toarray()
    df = pd.DataFrame(X, index=adata.obs_names, columns=["expression"])
    assert groupby in adata.obs.columns
    
    df["group"]=adata.obs[groupby]
    assert hue in adata.obs.columns
    df["hue"]=adata.obs[hue]
    df["hue"]=df["hue"].cat.reorder_categories(hue_order)

    _ax = sns.violinplot(x="group", y="expression", data=df, inner="box",
                             hue_order=None, hue='hue', split=split,palette="Set2",
                             scale=scale, orient='vertical')#palette={"AA": "y", "YA": "b"}
    if strip:
        _ax = sns.stripplot(x="group", y="expression", data=df,
                                hue='hue', dodge=True,hue_order=None,
                                jitter=jitter, color='black', size=size,ax=_ax)
    _ax.set_xlabel('')
    xlabels=_ax.get_xticklabels()
    _ax.set_xticklabels(xlabels,rotation=45)
    _ax.set_ylabel('expression')
    _ax.set_title(gene_name)
    _ax.legend_.remove()
    if save:
        if not os.path.exists(figdir):
            os.makedirs(figdir)
        filename=os.path.join(figdir,gene_name+"_violin_hue.pdf")
        pl.savefig(filename, dpi=80,figsize=(16,12),bbox_inches='tight')
        pl.close()

    else:
        return _ax


def violin_sns(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = "status",
    gene_name: Optional[Iterable[str]] = None,
    use_raw: Optional[bool] = True,
    order: Optional[str] = None,
    figdir:Optional[str] ="./figures",
    scale: str = "area",#'width',
    orient: str="h",
    strip: bool = False,
    jitter: Union[int, float, bool] = False,
    size: int = 1,
    save: Optional[bool] =True,
    add_pval: Optional[bool] =True
):
    if adata.raw is not None and use_raw:
        X= adata.raw[:, gene_name].X

    else:
        X= adata[:, gene_name].X
    if issparse(X):
        X = X.toarray()
    df = pd.DataFrame(X, index=adata.obs_names, columns=["expression"])
    assert groupby in adata.obs.columns

    df["group"]=adata.obs[groupby]

    pl.figure(figsize=(12, 8))

    if orient=="h":
        _ax = sns.violinplot(y="group", x="expression",data=df,
                             order=order,palette="Set2",
                             scale=scale)#palette={"AA": "y", "YA": "b"})
        if strip:
           _ax = sns.stripplot(y="group", x="expression", data=df,
                                 dodge=True,
                                jitter=jitter, color='black', size=size,ax=_ax)
    else :
        _ax = sns.violinplot(x="group",y="expression",data=df,
                             order=order,palette="Set2",
                             scale=scale)#palette={"AA": "y", "YA": "b"})
        if strip:
           _ax = sns.stripplot(x="group", y="expression", data=df,
                                 dodge=True,
                                jitter=jitter, color='black', size=size,ax=_ax)
    
    title=gene_name
    if add_pval:
        diff,pval=gene_comparison(adata,gene=gene_name,by=groupby)
        title=title+" : "+str(pval)

    _ax.set_title(title)
    if save:
        if not os.path.exists(figdir):
            os.makedirs(figdir)
        filename=os.path.join(figdir,gene_name+"_violin_sns.pdf")
        pl.savefig(filename, dpi=200,bbox_inches='tight')
        pl.close()

    else:
        return _ax

