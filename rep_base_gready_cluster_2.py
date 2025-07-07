import sys
from itertools import groupby

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import torch

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

text_loc="data/text_loc/"

cutoff=3

def main():
    r='\r'
    lreps=[]
    hf=open(f"{text_loc}pep.txt",'r')
    lines=hf.readlines()
    for line in lines:
        #print(line)
        l=list(line.split('_'))
        lreps.append(l[1])
    #print(len(lreps),lreps[0])
    npreps = np.array([list(i) for i in lreps[:]])
    npreps.shape
    i=1
    tempX=npreps.view(np.int32)
    x = torch.tensor(tempX.copy(), device=device)
    reps=torch.zeros([1,21],dtype=torch.int32,device=device)
    reps[0]=x[0]



    # reps[0]=x[0]

    clusterslab = torch.full([x.shape[0]], -1, device=device)
    clusterslab[0]=0

    #cutoff=3
    #print(reps)
    curClus=1
    while True:
        
        if i==x.shape[0]:
            break
        cur=torch.zeros([1,21],dtype=torch.int32,device=device)
        cur[0]=x[i]
        # print(x[i])
        
        mask = (cur[0] != reps).sum(1) > cutoff
        # print(mask)
        if torch.all(mask)==True:
            reps=torch.concat((reps,cur),0)
            clusterslab[i]=curClus
            curClus=curClus+1
        else:
            
            #print(cur[0],reps[clusterslab[torch.argmin(mask.int())]])

            clusterslab[i]=torch.argmin(mask.int())

        # print(reps)

        sys.stderr.write(f"progress: {i} / {x.shape[0]} | number of clusters found: {reps.shape[0]}{r}")    
        i=i+1
    n='\n'
    hf=open(f"{text_loc}clus3.txt",'w+')
    for pep,clu in zip(x,clusterslab): 
        pep="".join(pep.to("cpu").numpy().astype('uint32').view('U1').tolist())
        clu=int(clu.to("cpu").numpy())
        #print(pep,clu)
        #break
        hf.write(f"_{pep}__{clu}{n}")
    hf.close()



if __name__ == "__main__":
    main()
