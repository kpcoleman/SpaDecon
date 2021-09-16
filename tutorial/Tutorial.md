<h1><center>SpaDecon Tutorial</center></h1>


<center>Author: Kyle Coleman, Jian Hu, Amelia Schroeder, Edward B. Lee, Mingyao Li*  
  
   
The files used in this tutorial can be downloaded from: https://drive.google.com/drive/folders/1_eBGKMVYh4p1eiHhrJASnqOD9z_CP5e0?usp=sharing   
 
### 0. Package installation
- SpaDecon installation requires a python version over 3.5.  You can check your version of python by entering the following commands: 
```python
import platform
platform.python_version()
```

    '3.6.8'

We recommend creating and activating a new conda environment when installing the SpaDecon package. For instance, 
```bash
conda create -n SpaDecon python=3.7
conda activate SpaDeccon
```        
    
There are mulitple ways to install SpaDecon:
    
- Install SpaDecon using PyPI:

```bash
pip3 install SpaDecon   
```    
    
- Download and install SpaDecon package from GitHub: 

```bash
git clone https://github.com/kylepcoleman87/SpaDecon
cd SpaDecon/
python3 setup.py install --user
```


    
### 1. Import modules


```python
import numpy as np
import scanpy as sc
import SpaDecon as spd
from skimage import io
```

### 2. Load data
SpaDecon requires four input data files:  
- Gene expression matrix for SRT data 
- Gene expression matrix for scRNA-seq data (with cell type labels)
- Spatial coordinates of SRT spots (optional)
- Histology image for SRT data (optional)
<br>
SpaDecon requires the ST and scRNA-seq gene expression data to be stored as AnnData matrices.  The rows (observations) of the matrices are the samples and the columns (variables) are the genes.



```python
#Read annotated scRNA-seq GE data (rows = cells, columns = genes, cell types in adata_sc.obs.celltype)
adata_sc = sc.read('../data/sc.h5ad')

#Read SRT GE data (rows = spots, columns = genes)
adata_st = sc.read_10x_h5('spadecon_tutorial_data/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5)

#Read SRT histology image
image=io.imread("spadecon_tutorial_data/V1_Breast_Cancer_Block_A_Section_1_image.tif")

#Read file with SRT spatial locations
spatial=pd.read_csv("spadecon_tutorial_data/tissue_positions_list.csv",header=None,index_col=0) 
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
