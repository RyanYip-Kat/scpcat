Scanpy Cat: Single-Cell RNA-seq Simple Analysis Pipeline

### Installtation
before installtation,you must install *requirments" first(via pip install -r requirements.txt)

when use transform module(convert seurat inito anndata),you need to install packages : *Seurat,SeuratDisk*  first,and also must have R env in this envrioment

```python
git clone https://github.com/RyanYip-Kat/scpcat
cd scpcat 
python setup.py install
```
###  subcommand
there are several subcommand to reliaze functions,{create,reduction,findmarker,dotplot,embplot,heatmap,vlnplot,subset,enrich},and you can with --help to get usages
eg.
```python
scpCAT.py --help
scpCAT.py create --help
```

### usages
#### create
```python
scpCAT.py create -d /path/to/10X/outs/filtered_feature_bc_matrix -o output
```

#### reduction
```python
scpCAT.py reduction -d /path/to/adata.h5ad --batch_correct --batch_key "sample" -o output
```
...

May be you will meet a bug when use 
```python
scpCAT.py transform -f pbmc_small.rds --name seurat -o test
```

##### A simple Bug(May be)
/path/lib/python3.8/site-packages/scpcat-1.0.0-py3.8.egg/scpcat/R/Seurat2Anndata.R': No such file or directory

you can fix it by ** cp -r scpcat/R /path/lib/python3.8/site-packages/scpcat-1.0.0-py3.8.egg/scpcat/

