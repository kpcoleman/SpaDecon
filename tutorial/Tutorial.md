<h1><center>SpaDecon Tutorial</center></h1>


<center>Authors: Kyle Coleman, Jian Hu, Amelia Schroeder, Edward B. Lee, Mingyao Li*  
  
   
 
### 0. Package installation
- SpaDecon installation requires a python version of at least 3.6. The version of python can be checked by: 
```python
import platform
platform.python_version()
```

    '3.7.11'

We recommend creating and activating a new conda environment when installing the SpaDecon package. For instance, 
```bash
conda create -n SpaDecon python=3.7
conda activate SpaDecon
```        
    
There are mulitple ways to install SpaDecon:
    
- Install SpaDecon using PyPI:

```bash
pip3 install SpaDecon   
```    
    
- Download and install SpaDecon package from GitHub: 

```bash
git clone https://github.com/kpcoleman/SpaDecon
cd SpaDecon/
python3 setup.py install --user
```


    
### 1. Import modules


```python
import SpaDecon as spd
import scanpy as sc
import pandas as pd
import numpy as np
from skimage import io
import os
```

### 2. Load data
Please download the spadecon_tutorial_data folder from: https://drive.google.com/drive/folders/1_eBGKMVYh4p1eiHhrJASnqOD9z_CP5e0?usp=sharing

SpaDecon requires four input data files:  
- Gene expression (GE) matrix for SRT data 
- GE matrix for scRNA-seq data (with cell type labels)
- Histology image for SRT data (optional)
- Spatial coordinates of SRT spots (optional)

The spatial locations file should contain 5 columns:  
  (0) Indicator (equal to 1 for all spots in SRT capture area)   
  (1) x coordinate  
  (2) y coordinate  
  (3) x pixel coordinate  
  (4) y pixel coordinate
  
<br>
SpaDecon requires the SRT and scRNA-seq GE data to be stored as AnnData matrices.  

<br>
(1) If reading in GE data in the form of h5/h5ad files:

```python
#set working directory to spadecon_tutorial_data using os.chdir()
  
#Read annotated scRNA-seq GE data (observations = cells, variables = genes, cell types in adata_sc.obs.celltype)
adata_sc = sc.read('bc_sc.h5ad')

#Read SRT GE data (rows = spots, columns = genes)
adata_srt = sc.read_10x_h5('V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5')
adata_srt.var_names_make_unique()
  
```
  
(2) If reading in GE data in the form of csv files:
```python
  #set working directory to spadecon_tutorial_data using os.chdir()
  
  #Read scRNA-seq GE data (rows = cells, columns = genes)
  sc_ge = pd.read_csv('bc_sc_ge.csv', index_col = 0)
  
  #Read scRNA-seq cell-type labels
  sc_types = pd.read_csv('bc_sc_types.csv', index_col = 0)
  
  #Convert scRNA-seq GE data to AnnData object
  adata_sc = sc.AnnData(sc_ge)
  
  #Insert cell-type labels into "celltype" column adata_sc.obs
  adata_sc.obs['celltype'] = sc_types['celltype']
  
  #Read SRT GE data (rows = cells, columns = genes)
  srt_ge = pd.read_csv('bc_srt_ge.csv', index_col = 0)
  
  #Convert SRT GE data to AnnData object
  adata_srt = sc.AnnData(srt_ge)
```  

  
```python  
#Read SRT histology image
histology = io.imread("V1_Breast_Cancer_Block_A_Section_1_image.tif")

#Read file with SRT spatial locations
locations = pd.read_csv("tissue_positions_list.csv",header=None,index_col=0) 
locations = locations.loc[adata_st.obs.index]
```

### 3. Run SpaDecon

```python
clf = spd.SpaDecon()
#The technology parameter is used to determine an upper bound on the number of cell types per spot (specify "Visium" or "ST")
clf.deconvolution(source_data=adata_sc, target_data=adata_srt, histology_image=histology, spatial_locations=locations, technology='Visium')
spadecon_proportions = clf.props
spadecon_proportions.to_csv('spadecon_proportions.csv')  
```

### 4. Visualization of results using Seurat (R\)
```R
library(Seurat)
library(SeuratData)
library(ggplot2)
st = Load10X_Spatial('.','V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5', assay = 'Spatial')
spadecon_proportions = read.csv('spadecon_proportions.csv', row.names = 1, header= T, check.names = F)
st@meta.data = spadecon_proportions
SpatialFeaturePlot(st, features = 'Tumor', alpha = c(0, 1)) + ggplot2::scale_fill_gradientn(colours = heat.colors(10, rev = TRUE),limits = c(0, 1)) + ggtitle('10X Visium Breast Cancer') + theme(plot.title = element_text(size = 15, face = "bold"))
```

![png]('tumor_heatmap.png')
