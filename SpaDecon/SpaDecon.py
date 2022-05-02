import numpy as np
import pandas as pd
from scipy.sparse import issparse
from . ItClust import transfer_learning_clf
from . calculate_adj import distance
from . calculate_adj import calculate_adj_matrix
from . utils import find_l
import tensorflow as tf
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)


class SpaDecon(object):
    def __init__(self):
        super(SpaDecon, self).__init__()
    def deconvolution(self, source_data, target_data, histology_image=None, spatial_locations = None, p=0.5, histology = True, spatial = True, adj_matrix=None, adj = False, technology = 'Visium'):
        if technology=='Visium':
            threshold = 1/30
        elif technology=='ST':
            threshold = 1/200
        if issparse(target_data.X):
            target_data.X=target_data.X.A
        target_data.var_names=[i.upper() for i in list(target_data.var_names)]
        target_data.var["genename"]=target_data.var.index.astype("str")

        if adj:
            self.adj = adj_matrix
            l = find_l(p, [i*0.1 for i in range(1,20)], self.adj)
            adj_sub=np.exp(-1*self.adj/(2*(l**2)))
            target_data.X=np.matmul(adj_sub,target_data.X)
 
        elif spatial:
            target_data.obs["x1"] = np.array(spatial_locations[1])
            target_data.obs["x2"] = np.array(spatial_locations[2])
            target_data.obs["x3"] = np.array(spatial_locations[3])
            target_data.obs["x4"] = np.array(spatial_locations[4])
            target_data.obs["x5"] = np.array(spatial_locations[5])
            target_data=target_data[target_data.obs["x1"]==1]
            adj=calculate_adj_matrix(x=target_data.obs["x2"].tolist(),y=target_data.obs["x3"].to_list(), x_pixel=target_data.obs["x4"].to_list(), y_pixel=target_data.obs["x5"].to_list(), image=histology_image, histology = histology)
            #self.adj = adj
            l = find_l(p, [i*0.1 for i in range(1,20)], adj)
            adj_sub=np.exp(-1*adj/(2*(l**2)))
            target_data.X=np.matmul(adj_sub,target_data.X)
            del adj
        clf=transfer_learning_clf()
        clf.fit(source_data, target_data, tol = [0.01], threshold = threshold)
        type_pred = clf.predict(write=False)
        spad_props = type_pred[1]
#        spad_props.columns = [i[0] for i in type_pred[2].values()]
        spad_props.columns = clf.celltypes_final
        spad_props.index = [i[0:len(i)-7] for i in spad_props.index]
        self.props = spad_props
        return spad_props
