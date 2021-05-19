# scpcat
Scanpy Cat: Single-Cell RNA-seq Simple Analysis Pipeline

### Installtation

before installtation,you must install *requirments" first(via pip install -r requirements.txt)
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
