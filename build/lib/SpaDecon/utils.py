import numpy as np

def find_l(p, list_l, adj,sep=0.001, tol=0.01):
    print('...Finding optimal l for ','p=',p, '...',sep='')
    p_list = []
    for l in list_l:
        adj_exp=np.exp(-1*adj/(2*(l**2)))
        #print("l is ",str(l),"Percentage of total expression contributed by neighborhoods:",np.mean(np.sum(adj_exp,1))-1)
        p_list.append(np.mean(np.sum(adj_exp,1))-1)
    l1 = list_l[(np.abs(np.array(p_list)-p)).argmin()]
    if p_list[(np.abs(np.array(p_list)-0.5)).argmin()]<p:
        start = l1
        end = l1+0.1
    else:
        start = l1-0.1
        end = l1
    for i in np.arange(start, end, sep):
        adj_exp=np.exp(-1*adj/(2*(i**2)))
        q=np.mean(np.sum(adj_exp,1))-1
        #print("l=", str(i), "p=", str(round(q,5)))
        if np.abs(p-q)<tol:
            print("   l=", str(i), ", p=", str(round(q,5)), sep='')
            return i
    print("l not found, try bigger range or smaller sep!")


