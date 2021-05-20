import numpy as np
import pandas as pd
import os
import scanpy as sc
import argparse

from scpcat.functions import export_subset,_reduction,_createObj,_FindMarkers
from scpcat.functions import detectDoublet,transformer,enrichGO,_extract
from scpcat.functions import _dotplot,_EmbPlot,_DoHeatmap,VlnPlot
from scpcat.plotting import violin_hue,violin_sns

if __name__=="__main__":
  parser = argparse.ArgumentParser(description='Scanpy Cat: Single-Cell RNA-seq Simple Analysis Pipeline...,Just Very Simple,HA HA HA...')
  subparsers = parser.add_subparsers()
  ####################
  parser_create = subparsers.add_parser('create',help="create anndata object")
  parser_create.add_argument('--path',"-d",type=str, default=None,help="cellranger count or aggr pipeline out path")
  parser_create.add_argument('--outdir', '-o', type=str, default='output/', help='Output path')
  parser_create.set_defaults(func=_createObj)
  
  ####################
  parser_reduction = subparsers.add_parser('reduction',help="reduction and cluster")
  parser_reduction.add_argument('--path',"-d",type=str, default=None,help="anndata h5ad file")
  parser_reduction.add_argument('--outdir', '-o', type=str, default='output/', help='Output path')
  parser_reduction.add_argument('--batch_correct', action='store_true', help='whether do batch correct')
  parser_reduction.add_argument('--batch_key', type=str, default=None, help='if do batch correct,which column as batch key')
  parser_reduction.add_argument('--n_jobs', '-t', type=int, default=8, help='number of threads')
  parser_reduction.set_defaults(func=_reduction)
  
  ####################
  parser_marker= subparsers.add_parser('findmarker',help="find markers")
  parser_marker.add_argument('--path',"-d",type=str, default=None,help="anndata h5ad file")
  parser_marker.add_argument('--outdir', '-o', type=str, default='output/', help='Output path')

  parser_marker.add_argument("--subset1",nargs="+",type=str,default=None,help="subset1")
  parser_marker.add_argument("--subset2",nargs="+",type=str,default=None,help="subset2")
  parser_marker.add_argument("--column1",type=str,default="leiden",help="subset column1")
  parser_marker.add_argument("--column2",type=str,default=None,help="subset column2 ")
  parser_marker.add_argument("--groupby",type=str,default="status",help="groupby")
  parser_marker.add_argument("--method",type=str,default="t-test",help="t-test_overestim_var,wilcoxon")
  parser_marker.add_argument("--splitby",type=str,default=None,help="if splitby,will run in each with groupby(leiden,label_main,label_fine)")
  parser_marker.set_defaults(func=_FindMarkers)
  
  ####################
  parser_dotplot= subparsers.add_parser('dotplot',help="do dotplot")
  parser_dotplot.add_argument('--path',"-d",type=str, default=None,help="anndata h5ad file")
  parser_dotplot.add_argument('--outdir', '-o', type=str, default='output/', help='Output path')

  parser_dotplot.add_argument("--subset1",nargs="+",type=str,default=None,help="subset1")
  parser_dotplot.add_argument("--subset2",nargs="+",type=str,default=None,help="subset2")
  parser_dotplot.add_argument("--column1",type=str,default="leiden",help="subset column1")
  parser_dotplot.add_argument("--column2",type=str,default=None,help="subset column2 ")
  parser_dotplot.add_argument("--groupby",type=str,default="status",help="groupby")
  parser_dotplot.add_argument("--level",nargs="+",type=str,default=None,help="groupby factor level")
  parser_dotplot.add_argument("--color",type=str,default="Oranges",help="which color to use")
  parser_dotplot.add_argument("--genelist",type=str,required=True,default=None,help="gene list file,to show in dotplot")
  parser_dotplot.set_defaults(func=_dotplot)
  
  #####################
  parser_embplot= subparsers.add_parser('embplot',help="do embedding plot")
  parser_embplot.add_argument('--path',"-d",type=str, default=None,help="anndata h5ad file")
  parser_embplot.add_argument('--outdir', '-o', type=str, default='output/', help='Output path')

  parser_embplot.add_argument("--subset1",nargs="+",type=str,default=None,help="subset1")
  parser_embplot.add_argument("--subset2",nargs="+",type=str,default=None,help="subset2")
  parser_embplot.add_argument("--column1",type=str,default="leiden",help="subset column1")
  parser_embplot.add_argument("--column2",type=str,default=None,help="subset column2 ")
  parser_embplot.add_argument("--groupby",type=str,default="status",help="groupby")
  parser_embplot.add_argument("--level",nargs="+",type=str,default=None,help="groupby factor level")
  parser_embplot.add_argument("--theme",type=str,default=None,help="which theme use for embedding plot")
  parser_embplot.add_argument("--embedding",type=str,default="umap",help="which embedding use for embedding plot")
  parser_embplot.add_argument("--genelist",nargs="+",type=str,default=None,help="gene list file,to show in embdding plot")
  parser_embplot.set_defaults(func=_EmbPlot)
  
  ####################
  parser_heatmap= subparsers.add_parser('heatmap',help="do heatmap")
  parser_heatmap.add_argument('--path',"-d",type=str, default=None,help="anndata h5ad file")
  parser_heatmap.add_argument('--outdir', '-o', type=str, default='output/', help='Output path')

  parser_heatmap.add_argument("--subset1",nargs="+",type=str,default=None,help="subset1")
  parser_heatmap.add_argument("--subset2",nargs="+",type=str,default=None,help="subset2")
  parser_heatmap.add_argument("--column1",type=str,default="leiden",help="subset column1")
  parser_heatmap.add_argument("--column2",type=str,default=None,help="subset column2 ")
  parser_heatmap.add_argument("--groupby",type=str,default="status",help="groupby")
  parser_heatmap.add_argument("--level",nargs="+",type=str,default=None,help="groupby factor level")
  parser_heatmap.add_argument("--color",type=str,default="Oranges",help="which color to use")
  parser_heatmap.add_argument("--genelist",type=str,default=None,help="gene list file,to show in heatmap")
  parser_heatmap.add_argument('--n_genes', '-n', type=int, default=10, help='number of genes in echo group to show')
  parser_heatmap.set_defaults(func=_DoHeatmap)
  
  ####################
  parser_vln= subparsers.add_parser('vlnplot',help="do violin plot")
  parser_vln.add_argument('--path',"-d",type=str, default=None,help="anndata h5ad file")
  parser_vln.add_argument('--outdir', '-o', type=str, default='output/', help='Output path')

  parser_vln.add_argument("--subset1",nargs="+",type=str,default=None,help="subset1")
  parser_vln.add_argument("--subset2",nargs="+",type=str,default=None,help="subset2")
  parser_vln.add_argument("--column1",type=str,default="leiden",help="subset column1")
  parser_vln.add_argument("--column2",type=str,default=None,help="subset column2 ")
  parser_vln.add_argument("--groupby",type=str,default="status",help="groupby")
  parser_vln.add_argument("--level",nargs="+",type=str,default=None,help="groupby factor level")
  parser_vln.add_argument("--genelist",type=str,default=None,help="gene list file,to show in vlnplot")
  parser_vln.add_argument("--orient",type=str,default="h",help="h(hline) or v(vline")
  parser_vln.add_argument("--hue_order",nargs="+",type=str,default=None,help="hue order(HC vs VKH")
  parser_vln.add_argument("--hue",type=str,default=None,help="hue column,level length ==2,like (Status)")
  parser_vln.set_defaults(func=VlnPlot)
  
  
  ####################
  parser_s= subparsers.add_parser('subset',help="export subset")
  parser_s.add_argument('--path',"-d",type=str, default=None,help="anndata h5ad file")
  parser_s.add_argument('--outdir', '-o', type=str, default='output/', help='Output path')

  parser_s.add_argument("--subset",nargs="+",type=str,default=None,help="subset")
  parser_s.add_argument("--column",type=str,default="leiden",help="subset column")
  parser_s.add_argument("--bclist",type=str,default=None,help="cell barcode list file,to show in heatmap")
  parser_s.add_argument('--invert', action='store_true', help='whether invert')
  parser_s.set_defaults(func=export_subset)
  
  
  ####################
  parser_e= subparsers.add_parser("enrich",help="gene enrichment")
  parser_e.add_argument('--path',"-d",type=str, default=None,help="anndata h5ad file")
  parser_e.add_argument('--outdir', '-o', type=str, default='output/', help='Output path')
  parser_e.add_argument('--org', type=str, default='hsapiens', help='hsapiens or mmusculus')
  parser_e.add_argument("--genelist",type=str,default=None,help="gene list for enrichment")
  parser_e.add_argument('--groupby', type=str, default="leiden",help='groupby')
  parser_e.set_defaults(func=enrichGO)
  
  parser_t= subparsers.add_parser("transform",help="convert seurat to anadata")
  parser_t.add_argument('--filename',"-f",type=str, default=None,help="seurat rds filename")
  parser_t.add_argument("--name",type=str,default="seurat",help="prefix name")
  parser_t.add_argument("--outdir","-o",type=str,default="anndata",help="outpath to save h5ad")
  parser_t.set_defaults(func=transformer)


  parser_d= subparsers.add_parser("doublet",help="detect doublet")
  parser_d.add_argument('--path',"-d",type=str, default=None,help="anndata h5ad file")
  parser_d.add_argument('--outdir', '-o', type=str, default='output', help='Output path')
  parser_d.set_defaults(func=detectDoublet)
  
  #####################
  parser_ex= subparsers.add_parser('wraperLoupe',help="extract embedding and cluster for loupe software")
  parser_ex.add_argument('--path',"-d",type=str, default=None,help="anndata h5ad file")
  parser_ex.add_argument('--outdir', '-o', type=str, default='output/', help='Output path')

  parser_ex.add_argument("--subset1",nargs="+",type=str,default=None,help="subset1")
  parser_ex.add_argument("--subset2",nargs="+",type=str,default=None,help="subset2")
  parser_ex.add_argument("--column1",type=str,default="leiden",help="subset column1")
  parser_ex.add_argument("--column2",type=str,default=None,help="subset column2 ")
  parser_ex.add_argument("--emblist",nargs="+",type=str,default="umap",help="embedding names")
  parser_ex.add_argument("--labels",nargs="+",type=str,default="leiden",help="cluster names")
  parser_ex.set_defaults(func=_extract)
  
  args=parser.parse_args()
  args.func(args)
