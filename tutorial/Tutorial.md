<h1><center>SpaDecon Tutorial</center></h1>


<center>Author: Kyle Coleman, Jian Hu, Amelia Schroeder, Mingyao Li*

### 0. Installation
- SpaDecon installation requires a python version over 3.5.  You can check your version of python by entering the following commands: 
```python
import platform
platform.python_version()
```

    '3.6.8'

We recommend creating a conda environment when installing the SpaDecon package. For instance, 
```
"""    
conda create -n SpaDecon python=3.7
conda activate SpaDeccon
```
"""        
    
There are mulitple ways to install SpaDecon:
    
- Install SpaDecon using PyPI:
```
"""    
pip3 install SpaDecon   
```
"""
    
- Download and install SpaDecon package from GitHub: 

```
"""
git clone https://github.com/kylepcoleman87/SpaDecon
cd SpaDecon/
python3 setup.py install --user
"""
```


    
### 1. Import python modules


```python
import numpy as np
import scanpy as sc
import SpaDecon as spd
from skimage import io
```

### 2. Read in data
The current version of SpaDecon requires four input data files:  
- Gene expression matrix for ST data 
- Gene expression matrix for scRNA-seq data (with cell type labels)
- Spatial coordinates of ST spots (optional)
- Histology image for ST data (optional)
<br>
SpaDecon requires the ST and scRNA-seq gene expression data to be stored as AnnData matrices.  The rows (observations) of the matrices are the samples and the columns (variables) are the genes.



```python
#Read annotated scRNA-seq data (rows = cells, columns = genes)
adata_sc = sc.read('../data/sc.h5ad')

#Read ST data (rows = spots, columns = genes)
adata_st = sc.read_10x_mtx('./spatial/anterior/filtered_feature_bc_matrix2', var_names='gene_symbols', cache = True)

#Read histology image for ST data
image=io.imread("../tutorial/data/histology.tif")

#Read file with spatial locations of spots in ST data
spatial=pd.read_csv("../tutorial/data/positions.txt",sep=",",header=None,na_filter=False,index_col=0) 
```


### 3. Run SpaDecon

```python
clf = spd.SpaDecon()
spd_proportions = clf.deconvolution(adata_sc, adata_st, image, spatial)
```

### 4. Plot Results (R\)
```R
library(Seurat)
library(SeuratData)
library(ggplot2)
anterior <- LoadData("stxBrain", type = "anterior1")
spd_proportions = read.csv('../data/spd_proportions.csv', row.names = 1, header= T, check.names = F)
anterior@meta.data = spd_proportions
SpatialFeaturePlot(anterior, features = 'L6', alpha = c(0, 1)) + ggplot2::scale_fill_gradientn(colours = heat.colors(10, rev = TRUE),limits = c(0, 1)) + ggtitle('Anterior1_L6) + theme(plot.title = element_text(size = 15, face = "bold"))
