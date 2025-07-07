import sys
from itertools import groupby

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import torch

import getData.DataProcess as dp


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

text_loc="data/text_loc/"
fasta_loc="data/fasta_loc/"



cutoff=7

def main():

    
    dp.sample_possilbe_easy_neg_data(fasta_loc=fasta_loc, text_loc=text_loc)
    
    
    knownpeps_set=set()
    r='\r'
    lreps=[]
    hf=open(f"{text_loc}easyPepRaw.txt",'r')
    lines=hf.readlines()
    for line in lines:
        #print(line)
        l=list(line.split('_'))
        lreps.append(l[1])
        knownpeps_set.add(l[1])
    #print(len(lreps),lreps[0])
    npreps = np.array([list(i) for i in lreps[:]])
    npreps.shape
    i=1
    tempX=npreps.view(np.int32)
    x = torch.tensor(tempX.copy(), device=device)
    knownpeps=[]
    hf=open(f"{text_loc}pep.txt",'r')
    lines=hf.readlines()
    for line in lines:
        #print(line)
        l=list(line.split('_'))
        knownpeps.append(l[1])
    #print(len(lreps),lreps[0])
    #reps=torch.zeros([len(knownpeps),21],dtype=torch.int32,device=device)
    knownpeps = np.array([list(i) for i in knownpeps[:]])
    tempkX=knownpeps.view(np.int32)
    reps = torch.tensor(tempkX.copy(), device=device)
    
    
    # reps[0]=x[0]
    clusterslab = torch.full([x.shape[0]], -1, device=device)
    clusterslab[0]=0
    #cutoff=3
    #print(reps)
    curClus=1
    newreps=torch.zeros([0,21],dtype=torch.int32,device=device)
    # reps[0]=x[0]
    # newreps=
    while True:
        
        if i==x.shape[0]:
            break
        cur=torch.zeros([1,21],dtype=torch.int32,device=device)
        cur[0]=x[i]
        # print(x[i])
        
        mask = (cur[0] != reps).sum(1) > cutoff
        # print(mask)
        if torch.all(mask)==True:
            mask = (cur[0] != newreps).sum(1) > cutoff
            if torch.all(mask)==True:
                newreps=torch.concat((newreps,cur),0)
                clusterslab[i]=curClus
                curClus=curClus+1
            else:
                
                #print(cur[0],reps[clusterslab[torch.argmin(mask.int())]])

                clusterslab[i]=torch.argmin(mask.int())

        # print(reps)

        sys.stderr.write(f"progress: {i} / {x.shape[0]} | number of clusters found: {newreps.shape[0]}{r}")    
        i=i+1
    n='\n'
    hf=open(f"{text_loc}clus7_easyneg.txt",'w+')
    for pep,clu in zip(x,clusterslab): 
        pep="".join(pep.to("cpu").numpy().astype('uint32').view('U1').tolist())
        clu=int(clu.to("cpu").numpy())
        hf.write(f"_{pep}__{clu}{n}")
    hf.close()



if __name__ == "__main__":
    main()