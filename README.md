# SpaDecon: cell-type deconvolution in spatial transcriptomics with semi-supervised learning

### Kyle Coleman, Jian Hu, Amelia Schroeder, Edward B. Lee, Mingyao Li*

SpaDecon is a semi-supervised learning-based method developed to perform cell-type deconvolution on spatially resolved transcriptomics (SRT) datasets. SpaDecon has been shown to provide accurate cell-type deconvolution results for both Spatial Transcriptomics (ST) and 10X Visium SRT datasets. Annotated scRNA-seq gene expression data from the same type of tissue as the SRT data are required for deconvolution.

![png](images/spadecon_workflow.png)

## SpaDecon Installation
- SpaDecon installation requires python version 3.7. The version of python can be checked by: 
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

## Tutorial
A markdown tutorial file can be found here: https://github.com/kpcoleman/SpaDecon/blob/main/tutorial/Tutorial.md

A tutorial in the form of a jupyter notebook can be found here: https://github.com/kpcoleman/SpaDecon/blob/main/tutorial/tutorial.ipynb 



## Software Requirements  
python==3.7  
keras==2.2.4  
pandas==1.2.4  
numpy==1.20.1  
scipy==1.6.2  
scanpy==1.7.0  
anndata==0.7.6  
sklearn  
tensorflow==1.14.0  

