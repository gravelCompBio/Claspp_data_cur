import os
import sys
import tarfile 
import math
import random
import time

import numpy as np
from scipy.io import savemat
from unipressed import UniprotkbClient
from sklearn.cluster import SpectralClustering
from sklearn.model_selection import train_test_split



random.seed(42)
pwd = os.getcwd() 

ltu=["Phosphorylation",
    "Ubiquitination",
    "Acetylation",
    "N-linked-Glycosylation",
    "O-linked-Glycosylation",
    "Methylation",
    "Malonylation",
    "Sulfoxidation",
    "Sumoylation",
    "S-palmitoylation",
    "Glutathionylation",
    #"Succinylation",
    "Hydroxylation"
    ]

sltu=set(ltu)

lab2res={}
lab2res["Acetylation"]=set(['A','K','M'])
lab2res["Sulfoxidation"]=set(['M'])
lab2res["Methylation"]=set(['R','K'])
lab2res["Sumoylation"]=set(['K'])
lab2res["N-linked-Glycosylation"]=set(['N'])
lab2res["Phosphorylation"]=set(['S','T','Y'])
lab2res["O-linked-Glycosylation"]=set(['S','T'])
lab2res["Malonylation"]=set(['K'])
lab2res["S-palmitoylation"]=set(['C'])
lab2res["Ubiquitination"]=set(['K'])
#lab2res["Succinylation"]=set(['K'])
lab2res["Glutathionylation"]=set(['C'])
lab2res["Hydroxylation"]=set(['P','K'])


baselabs=[]
baselabs.append("ST-Phosphorylation")
baselabs.append("K-Ubiquitination")
baselabs.append("Y-Phosphorylation")
baselabs.append("K-Acetylation")
baselabs.append("N-N-linked-Glycosylation")
baselabs.append("ST-O-linked-Glycosylation")
baselabs.append("RK-Methylation")
baselabs.append("K-Sumoylation")
baselabs.append("K-Malonylation")
baselabs.append("M-Sulfoxidation")
baselabs.append("AM-Acetylation")
#baselabs.append("K-Succinylation")
baselabs.append("C-Glutathionylation")
baselabs.append("C-S-palmitoylation")
baselabs.append("PK-Hydroxylation")



bl2nc={}
bl2nc["ST-Phosphorylation"]=5
bl2nc["K-Ubiquitination"]=20
bl2nc["Y-Phosphorylation"]=1
bl2nc["K-Acetylation"]=10
bl2nc["N-N-linked-Glycosylation"]=1
bl2nc["ST-O-linked-Glycosylation"]=5
bl2nc["RK-Methylation"]=4
bl2nc["K-Sumoylation"]=1
bl2nc["K-Malonylation"]=1
bl2nc["M-Sulfoxidation"]=1
bl2nc["PK-Hydroxylation"]=1
bl2nc["AM-Acetylation"]=1
#bl2nc["K-Succinylation"]=1
bl2nc["C-Glutathionylation"]=1
bl2nc["C-S-palmitoylation"]=1

total=0
for bl in bl2nc.keys():
    total+=bl2nc[bl]

print(total)




#########################################################################################
#----------------------------------------------------------------------------------------
#  helper funtions 
#----------------------------------------------------------------------------------------
#########################################################################################


def getinfo(info):
        l=list(info.split("_"))
        uni=l[0]
        pos=l[1]
        acc=l[2]
        pep=l[3]
        lab=l[4]
        clu=l[5]
        h=list(acc.split('-'))
        spe=h[1]
        return spe,lab,pep


def getrand(ilist):
    output=random.choice(ilist)
    return output


def getaff(pep1,pep2):
    score=0
    for p in range(len(pep1)):
        if pep1[p]==pep2[p]:
            score+=1
    return score


#########################################################################################
#----------------------------------------------------------------------------------------
#  stage 1 
#----------------------------------------------------------------------------------------
#########################################################################################


def unzip_and_read_files_from_dbPTM( db_loc : str ="data/dbPTMloc/"  ):
    ############################### used unipressed to get all sequences 
    """
    unizip the raw DB PTM files and read the data
    Args:
        db_loc (str) : location of the dbptmData
    Returns:
        masterlist (list[str]) : master list of raw data 
        uni2name (dict[str:str]) : uniprot id 2 name
    """
    # ------------------------------------  Extract the info from each ptm type
    wack=[]
    uni2name={}
    masterlist=[]#   ---------------------  store raw info into masterlist (list)
    folder = os.listdir(f"{db_loc}")
    for infolder in folder:
        input_file = f"{db_loc}{infolder}"
        output_file = f"{db_loc}{infolder[:-3]}"
        if not os.path.isdir(input_file):
            with tarfile.open(input_file, 'r:gz') as tar:
                tar.extractall(path=output_file)
    time.sleep(5)
    folder = os.listdir(f"{db_loc}")
    for infolder in folder:
        input_file = f"{db_loc}{infolder}"
        output_file = f"{db_loc}{infolder[:-3]}"
        if ".gz" in input_file:
            continue
        infolderI = os.listdir(f"{db_loc}"+infolder)
        for file in infolderI:
            fullpath=f"{db_loc}"+infolder+'/'+file
            lab=file.replace(' ','-')
            time.sleep(1)
            hf=open(fullpath,'r')
            lines=hf.readlines()
            for line in lines:
                l=list(line[:-1].split('\t'))
                uni=l[1]
                ung=l[0].replace('_','-')
                pos=l[2]
                pep=l[5]
                uni2name[uni]=l[0]
                info=f"{uni}_{pos}_{ung}_{pep}_{lab}"
                masterlist.append(info)
                if '_' in ung:
                    wack.append(info)
    return uni2name, masterlist


def intial_hand_curation_update( masterlist : list[str] , csv_loc : str ="data/csv_loc/"):
    """
    fixes a very small number of problamatic uniprot accestion that are diffrult to automate 
    Args:
        masterlist (list[str]) : masterlist of rawPTM data
        csv_loc (int): directory where it places all data.
    Returns:
        masterlist (list[str]) : updated master list of raw data 
        uniset (set[str]) : non-redundant set of uniprot ids
    """
    unisetForF=set()
    newmasterlist=[]
    for info in masterlist:
        l=list(info.split('_'))
        if l[2]=='':
            unisetForF.add(l[0])
        else:
            newmasterlist.append(info)
    uni2acc={}
    hf=open(f"{csv_loc}uni2acc.tsv",'r')
    lines=hf.readlines()
    for line in lines[1:]:
        l=list(line.split('\t'))
        uni=l[0]
        acc=l[2]
        uni2acc[uni]=acc
    for info in masterlist:
        l=list(info.split('_'))
        if l[2]=='':
            uni=l[0]
            if uni not in uni2acc.keys():
                continue
            acc=uni2acc[uni]
            newmasterlist.append(f"{l[0]}_{l[1]}_{acc.replace('_','-')}_{l[3]}_{l[4]}")
    masterlist=newmasterlist
    #   --------------------  get a set of all unique uniprots
    uniset=set()#  ---------  set that stores uniprots 
    for info in masterlist:
        l=list(info.split("_"))
        uniset.add(l[0])
    return masterlist, uniset


def separate_old_and_new_data_v1(uni2seq: dict[str:str], masterlist : list[str], text_loc : str = "data/text_loc/",):
    """
    separate old and new data from dbPTM
    Args:
        uni2seq (dict[str:str]) : uniprot id 2 sequences 
        masterlist (list[str]) : masterlist of rawPTM data
        text_loc (int): directory where it places all  txt data.
        
    Returns:
        newunis (set[str]) : set of current uniprot ids 
        oldunis (set[str]) : set of out-of-date uniprot ids
    """
    nl='\n'
    newunis=set()
    oldunis=set()
    for info in masterlist:
        l=list(info.split('_'))
        uni=l[0]
        if uni in uni2seq.keys():
            newunis.add(uni)
        else:
            oldunis.add(uni)
    output=open(f"{text_loc}unilist.txt",'w+')

    for uni in oldunis:
        output.write(f"{uni}{nl}")
    output.close()

    return newunis, oldunis


def separate_old_and_new_data_v2(uni2seq: dict[str:str], olduni2seq: dict[str:str], masterlist : list[str]):
    """
    update the uni ids to sequnece mappings and separate old and new data from dbPTM
    Args:
        uni2seq (dict[str:str]) : uniprot id 2 sequences 
        olduni2seq (dict[str:str]) : first round of fixed uniprot id 2 sequences 
        masterlist (list[str]) : masterlist of rawPTM data
    Returns:
        newunis (set[str]) : set of current uniprot ids 
        oldunis (set[str]) : set of out-of-date uniprot ids
    """
    for uni in olduni2seq.keys():
        uni2seq[uni]=olduni2seq[uni]
    newunis=set()
    oldunis=set()
    for info in masterlist:
        l=list(info.split('_'))
        uni=l[0]
        if uni in uni2seq.keys():
            newunis.add(uni)
        else:
            oldunis.add(uni)

    return newunis, oldunis







def final_merge(uni2seq: dict[str:str], olduni2seq: dict[str:str]):
    """
    update the uni ids to sequnece mappings 
        uni2seq (dict[str:str]) : uniprot id 2 sequences 
        olduni2seq (dict[str:str]) : first round of fixed uniprot id 2 sequences 
        
    Returns:
        uni2seq (dict[str:str]) : updated uniprot id 2 sequences 
    """
    for uni in olduni2seq.keys():
        uni2seq[uni]=olduni2seq[uni]
    return uni2seq





def remove_broken_peptides(masterlist : list[str]):
    """
    remove peptides that that have an illegal char and peptides that don't have 21 long peptides
    Args:
        masterlist (list[str]) : masterlist of rawPTM data
    Returns:
        masterlist (list[str]) : updated masterlist of rawPTM data
    """
    newML=[]
    badLi=[]
    legRes={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-'}
    for info in masterlist:
        l=list(info.split('_'))
        pep=l[3]
        b=0
        for c in pep:
            if c not in legRes:
                b=1
        if b==1:
            badLi.append(info)
        else:
            newML.append(info)
    masterlist=newML
    k21=[]
    bad=[]
    for info in masterlist:
        l=list(info.split('_'))
        pep=l[3]
        if len(pep)==21:
            k21.append(info)
        else:
            bad.append(info)
    masterlist=k21
    return masterlist



def fix_missmatched_uniprot_postions(uni2seq: dict[str:str], masterlist : list[str]):
    """
    find all mismatching postions of PTMs and correts them. This assumes that the frist
    appearance of the peptide in the full sequence is the correct postions if previously
    broken 
    Args:
        uni2seq (dict[str:str]) : uniprot id 2 sequences 
        masterlist (list[str]) : masterlist of rawPTM data
    Returns:
        updatelist (list[str]) : updatelist of PTM data
        missmatchedseq (list[str]) : peptide not found in full seq
    """
    notmapped=[]
    missmatchedseq=[]
    wrongpos=[]
    correpos=[]
    wrong=[]
    for i,info in enumerate(masterlist):
        l=list(info.split("_"))
        uni=l[0]
        pos=int(l[1])
        pep=l[3]
        if uni not in uni2seq.keys():
            notmapped.append(uni)
            continue
        seq=uni2seq[uni]
        paddedseq="----------"+seq+"----------"
        if pep == paddedseq[pos-1:pos+20]:
            correpos.append(info)
        else:
            wrong.append(info)
    for i,info in enumerate(wrong):
        l=list(info.split("_"))
        uni=l[0]
        pos=int(l[1])
        pep=l[3]
        seq=uni2seq[uni]
        paddedseq="----------"+seq+"----------"
        if pep in paddedseq:
            wrongpos.append(info)
        else:
            missmatchedseq.append(info)
    fixedpos=set()

    for i,info in enumerate(wrongpos):
        l=list(info.split("_"))
        uni=l[0]
        pos=int(l[1])
        pep=l[3]
        seq=uni2seq[uni]
        paddedseq="----------"+seq+"----------"
        newpos=paddedseq.find(pep)+1
        newinfo=f"{uni}_{newpos}_{l[2]}_{pep}_{l[4]}"
        fixedpos.add(newinfo)
    fixedpos=set()
    for i,info in enumerate(wrongpos):
        l=list(info.split("_"))
        uni=l[0]
        pos=int(l[1])
        pep=l[3]
        seq=uni2seq[uni]
        paddedseq="----------"+seq+"----------"
        newpos=paddedseq.find(pep)+1
        newinfo=f"{uni}_{newpos}_{l[2]}_{pep}_{l[4]}"
        fixedpos.add(newinfo)
    scor=set(correpos)
    updatedlist=[]
    for info in correpos:
        updatedlist.append(info)
    for info in fixedpos:
        if info in scor:
            continue
        updatedlist.append(info)
    wrongpos=[]
    correpos=[]
    wrong=[]

    for i,info in enumerate(updatedlist):
        l=list(info.split("_"))
        uni=l[0]
        pos=int(l[1])
        pep=l[3]
        if uni not in uni2seq.keys():
            notmapped.append(uni)
            continue
        seq=uni2seq[uni]
        paddedseq="----------"+seq+"----------"
        if pep == paddedseq[pos-1:pos+20]:
            correpos.append(info)
        else:
            wrong.append(info)


        
    for i,info in enumerate(wrong):
        l=list(info.split("_"))
        uni=l[0]
        pos=int(l[1])
        pep=l[3]
        seq=uni2seq[uni]
        paddedseq="----------"+seq+"----------"
        if pep in paddedseq:
            wrongpos.append(info)
        else:
            missmatchedseq.append(info)
    return missmatchedseq, updatedlist
    




def fix_missmatched_uniprot_postions_v2( uni2seq: dict[str:str], mmUni2ver2seq: dict[str:str], missmatchedseq : list[str], updatedlist : list[str]):
    """
    final merge of fixed data
    Args:
        uni2seq (dict[str:str]) : uniprot id 2 sequences 
        mmUni2ver2seq (dict[str:str]) : missmatched uniprot 2 sequences
        masterlist (list[str]) : masterlist of rawPTM data
        updatelist (list[str]) : updatelist of PTM data
    Returns:
        masterlist (list[str]) : updated masterlist of rawPTM data
    """
    fixed={}
    notfixed=set()
    for info in missmatchedseq:
        l=list(info.split("_"))
        uni=l[0]
        pep=l[3]
        if uni not in mmUni2ver2seq.keys():
            continue
        seqs=mmUni2ver2seq[uni]
        
        for ver in seqs.keys():
            seq=seqs[ver]
            paddedseq="----------"+seq+"----------"
            if pep in paddedseq:
                fixed[info]=seq
            else:
                notfixed.add(info)
    fixedpos=set()
    ks=list(fixed.keys())
    for i,info in enumerate(ks):
        l=list(info.split("_"))
        uni=l[0]
        pos=int(l[1])
        pep=l[3]
        seq=fixed[info]
        paddedseq="----------"+seq+"----------"
        newpos=paddedseq.find(pep)+1
        newinfo=f"{uni}_{newpos}_{l[2]}_{pep}_{l[4]}"
        fixedpos.add(newinfo)
        fixed[newinfo]=seq
    updatedlistV2=[]
    udlv2s=set(updatedlistV2)
    for info in updatedlist:
        updatedlistV2.append(info)
    for info in fixedpos:
        if info in udlv2s:
            continue
        updatedlistV2.append(info)
    correpos=[]
    wrong=[]
    for i,info in enumerate(updatedlistV2):
        l=list(info.split("_"))
        uni=l[0]
        pos=int(l[1])
        pep=l[3]
        if uni not in uni2seq.keys():
            continue
        seq=uni2seq[uni]
        if info in fixed.keys():
            seq=fixed[info]
        paddedseq="----------"+seq+"----------"
        if pep == paddedseq[pos-1:pos+20]:
            correpos.append(info)
        else:
            if uni in uni2seq.keys():
                seq=uni2seq[uni]
                paddedseq="----------"+seq+"----------"
                if pep == paddedseq[pos-1:pos+20]:
                    correpos.append(info)
            wrong.append(info)
    masterlist=correpos
    return masterlist



def filter_for_human_and_major_res_per_PTM_type( masterlist : list[str],text_loc : str = "data/text_loc/"):
    """
    final merge of fixed data
    Args:
        masterlist (list[str]) : masterlist of rawPTM data
        text_loc (int): directory where it places all  txt data.
    Returns:
        writes the file that plugs into the greedy represetive based clustering 
    """
    labset=set()
    lab2count={}
    for info in masterlist:
        l=list(info.split('_'))
        lab=l[4]
        if lab not in lab2count.keys():
            lab2count[lab]=0
        lab2count[lab]=lab2count[lab]+1
        labset.add(lab)
    nonUsedPTMs=set()
    slabs=sorted(lab2count.items(), key=lambda x: x[1], reverse=True)
    for lab in slabs:
        if lab[0] not in ltu:
            nonUsedPTMs.add(lab[0])
    newML=[]
    for info in masterlist:
        l=list(info.split('_'))
        pep=l[3]
        res=pep[10]
        lab=l[4]
        acc=l[2]
        h=list(acc.split('-'))
        spe=h[1]
        if lab not in ltu:
            newML.append(info)
            continue
        elif spe=="HUMAN" and res in lab2res[lab]:
            newML.append(info)
    n='\n'
    masterlist=newML
    hf=open(f"{text_loc}cleanedinfoTranserfer.txt",'w+')
    for info in masterlist:
        hf.write(f"{info}{n}")
    hf.close()
    peps=set()
    for info in masterlist:
        l=list(info.split("_"))
        peps.add(l[3])
    hf=open(f"{text_loc}pep.txt",'w+')
    for pep in peps:
        hf.write(f"_{pep}_{n}")
    hf.close()

#########################################################################################
#----------------------------------------------------------------------------------------
#  stage 3
#----------------------------------------------------------------------------------------
#########################################################################################


def readin_files_v1(text_loc : str = "data/text_loc/"):
    """
    reads in files generated from prev code 
    Args:
        text_loc (int): directory where it places all  txt data.
    Returns:
        masterlist (list[str]) : reads in masterlist from prevouse code
        pep2clus (dict[str:str]) : reads in greedy rep based clustering 
    """
    masterlist=[]
    pep2clus={}
    newML=[]
    hf=open(f"{text_loc}cleanedinfoTranserfer.txt",'r')
    lines=hf.readlines()
    for line in lines:
        masterlist.append(line[:-1])
    hf=open(f"{text_loc}clus3.txt",'r')
    lines=hf.readlines()
    for line in lines:
        l=list(line[:-1].split('_'))
        pep=l[1]
        clu=l[3]
        pep2clus[pep]=clu
    for info in masterlist:
        l=list(info.split('_'))
        pep=l[3]
        clus=pep2clus[pep]
        ninfo=info+'_'+clus
        newML.append(ninfo)
    masterlist=newML
    return masterlist, pep2clus


def track_all_human_multi_lable_peps(masterlist : list[str], text_loc : str = "data/text_loc/"):
    """
    tracks down peptide data from clustering that and find all human and multilable data 
    Args:
        masterlist (list[str]) : list of all PTM related data
        text_loc (int): directory where it places all  txt data.
    Returns:
        pep2info (dict[str:list[str]]) : merge all info for all identical peptides
        hpep2labs (dict[str:set(str)]) : dict linking each peptide to the many possible label for human data
        bs = set(str) : data that are not asssosated a speicies 
    """
    pep2info={}
    hpep2labs={}
    brokeninfo=[]
    prank={}

    for info in masterlist:
        l=list(info.split("_"))
        pep=l[3]
        if pep not in pep2info.keys():
            pep2info[pep]=[]
        pep2info[pep].append(info)
    nonhumanpeps=set()
    humanpeps=set()
    for pep in pep2info.keys():
        for info in pep2info[pep]:
            l=list(info.split('_'))
            acc=l[2]
            a=list(acc.split('-'))
            try:
                spe=a[1]
            except:
                print(info)
                brokeninfo.append(info)
            if spe=="HUMAN":
                humanpeps.add(pep)
            else:
                nonhumanpeps.add(pep)
    for pep in humanpeps:
        infos=pep2info[pep]
        for info in infos:
            l=list(info.split('_'))
            h=list(l[2].split("-"))
            try:
                spe=h[1]
            except:
                print(info)
                continue
            lab=l[4]
            if spe=="HUMAN" and lab in sltu:
                if pep not in hpep2labs.keys():
                    hpep2labs[pep]=set()
                hpep2labs[pep].add(lab)
    mlabpeps=set()
    for pep in hpep2labs.keys():
        if 1 < len(hpep2labs[pep]):
            mlabpeps.add(pep)
    bs=set(brokeninfo)
    labs=set()
    for pep in hpep2labs.keys():
        for lab in hpep2labs[pep]:
            if lab not in sltu:
                continue
            labs.add(lab)
    mlabel2pep={}
    pep2mllabel={}
    mlabel2pep["Multi_P1_Acetylation_Ubiquitination"]=set()
    mlabel2pep["Multi_P2_O-linked-Glycosylation_Phosphorylation"]=set()
    mlabel2pep["Multi_P3_Acetylation_Malonylation_Ubiquitination"]=set()
    mlabel2pep["Multi_P4_Methylation"]=set()
    mlabel2pep["Multi_P5_Ubiquitination"]=set()
    mlabel2pep["Multi_P6_Acetylation"]=set()
    mlabel2pep["Multi_P7_other"]=set()
    #"Multi_P1_Acetylation_Ubiquitination"
    for pep in hpep2labs.keys():
        labs=hpep2labs[pep]
        if len(labs) > 1:
            # if "Malonylation" in labs and "Succinylation" in labs:
            #     mlabel2pep["Multi_P0_Malonylation_Succinylation"].add(pep)
            #     pep2mllabel[pep]="Multi_P0_Malonylation_Succinylation"

            if "Acetylation" in labs and "Ubiquitination" in labs and len(labs)==2:
                mlabel2pep["Multi_P1_Acetylation_Ubiquitination"].add(pep)
                pep2mllabel[pep]="Multi_P1_Acetylation_Ubiquitination"

            elif "O-linked-Glycosylation" in labs and "Phosphorylation" in labs:
                mlabel2pep["Multi_P2_O-linked-Glycosylation_Phosphorylation"].add(pep)
                pep2mllabel[pep]="Multi_P2_O-linked-Glycosylation_Phosphorylation"

            elif "Acetylation" in labs and "Malonylation" in labs and "Ubiquitination" in labs:
                mlabel2pep["Multi_P3_Acetylation_Malonylation_Ubiquitination"].add(pep)
                pep2mllabel[pep]="Multi_P3_Acetylation_Malonylation_Ubiquitination"

            elif "Methylation" in labs:
                mlabel2pep["Multi_P4_Methylation"].add(pep)
                pep2mllabel[pep]="Multi_P4_Methylation"

            elif "Ubiquitination" in labs:
                mlabel2pep["Multi_P5_Ubiquitination"].add(pep)
                pep2mllabel[pep]="Multi_P5_Ubiquitination"

            elif "Acetylation" in labs:
                mlabel2pep["Multi_P6_Acetylation"].add(pep)
                pep2mllabel[pep]="Multi_P6_Acetylation"
            else:
                mlabel2pep["Multi_P7_other"].add(pep)
                pep2mllabel[pep]="Multi_P7_other"

    prank={}
    #prank["Multi_P0_Malonylation_Succinylation"]=1
    prank["Multi_P1_Acetylation_Ubiquitination"]=1
    prank["Multi_P2_O-linked-Glycosylation_Phosphorylation"]=2
    prank["Multi_P5_Ubiquitination"]=3
    prank["Multi_P3_Acetylation_Malonylation_Ubiquitination"]=4
    prank["Multi_P6_Acetylation"]=5
    prank["Multi_P4_Methylation"]=6
    hf=open(f"{text_loc}mulLab.txt",'w+')
    for pep in pep2mllabel.keys():
        labs=hpep2labs[pep]
        line=pep
        for lab in labs:
            line=line+'_'+lab
        hf.write(line+'\n')
    hf.close()
    return pep2info, hpep2labs, bs, mlabpeps, pep2mllabel, prank #################


def set_sampeling_prority(masterlist : list[str], mlabpeps : set[str], pep2mllabel : dict[str:set[str]], prank : dict[str:int], pep2info : (dict[str:list[str]]), bs : set[str], text_loc : str = "data/text_loc/"):
    """
    priotirty sampeling from the representive based greedy clustering to ensure even distrubution of sequence idenity 
    (favors multi labeled data and under representied data)
    Args:
        masterlist (list[str]) : list of all PTM related data
        text_loc (int): directory where it places all  txt data.
    Returns:
        masterlist (list[str]) : list of all PTM related data (updated)
        nonIntrestClusters (list[str]) : list of cluster that contain only ptm types that fall outside the prediction task   
        clus2info (dir[str:set[str]]) : cluster id mapped to all peptide info 
        nonHumanClus (list[str]) : ptm seq ident clusters that only have non human peptides
        rank (dir[str:int]) : ptm type ranked by abundance

    """
    nonHumanClus=[]
    humanClus=[]
    humanMLPepClus=[]
    humanNonMLPepClus=[]
    nonIntrestClusters=[]
    sampledInfo=[]
    clus2info={}
    nonUsedPTMs=set()
    rank={}
    for info in masterlist:
        l=list(info.split('_'))
        clus=l[5]
        if clus not in clus2info.keys():
            clus2info[clus]=set()
        clus2info[clus].add(info)
    labset=set()
    lab2count={}
    for info in masterlist:
        l=list(info.split('_'))
        lab=l[4]
        if lab not in lab2count.keys():
            lab2count[lab]=0
        lab2count[lab]=lab2count[lab]+1
        labset.add(lab)
    slabs=sorted(lab2count.items(), key=lambda x: x[1], reverse=True)
    for lab in slabs:
        if lab[0] not in ltu:
            nonUsedPTMs.add(lab[0])
    lab2res2count={}
    newML=[]
    for info in masterlist:
        l=list(info.split('_'))
        pep=l[3]
        res=pep[10]
        lab=l[4]
        acc=l[2]
        h=list(acc.split('-'))
        spe=h[1]
        if lab not in ltu:
            continue
        elif spe=="HUMAN":
            if lab not in lab2res2count.keys():
                lab2res2count[lab]={}
            if res not in lab2res2count[lab].keys():
                lab2res2count[lab][res]=0
            lab2res2count[lab][res]=lab2res2count[lab][res]+1
    r=0
    slabs=sorted(lab2count.items(), key=lambda x: x[1], reverse=True)
    for lab in slabs:
        if lab[0] in ltu:
            r=r+1
            rank[lab[0]]=r
    for info in masterlist:
        l=list(info.split('_'))
        clus=l[5]
        if clus not in clus2info.keys():
            clus2info[clus]=set()
        clus2info[clus].add(info)
    
    for clus in clus2info.keys():
        infos=clus2info[clus]
        humanl=[]
        for info in infos:
            if info in bs:
                continue
            spe, lab, pep=getinfo(info)
            if spe=="HUMAN":
                humanl.append(info)
        if len(humanl) > 0:
            mlset=[]
            humanClus.append(clus)
            for info in humanl:
                spe, lab, pep=getinfo(info)
                if pep in mlabpeps:
                    mlset.append(pep)
            if len(mlset) > 0:
                humanMLPepClus.append(clus)
            else:
                humanNonMLPepClus.append(clus)
        else:
            nonHumanClus.append(clus)
    for clus in humanNonMLPepClus:
        infos=clus2info[clus]
        maxRankSubClus=[]
        currank=0
        b=0
        for info in infos:
            if info in bs:
                continue
            spe, lab, pep=getinfo(info)
            if spe=="HUMAN" and lab in sltu:
                r=rank[lab]
                if r >currank:
                    maxRankSubClus=[]
                    currank=r
                    maxRankSubClus.append(info)
                    b=1
                elif rank==currank:
                    maxRankSubClus.append(info)
                else:
                    continue
        if b==1:
            sampledInfo.append(getrand(maxRankSubClus))
        else:
            nonIntrestClusters.append(clus)
    mlPriority=set()
    addsample=[]
    for i,clus in enumerate(humanMLPepClus):
        infos=clus2info[clus]
        subML=[]
        b=0
        for info in infos:
            if info in bs:
                continue
            #print(info)
            l=list(info.split('_'))
            h=list(l[2].split("-"))
            spe=h[1]
            pep=l[3]
            lab=l[4]
            if pep in pep2mllabel.keys() and spe=="HUMAN" and lab in sltu:
                subML.append(info)
        maxRankSubClus=[]
        currank=0
        for info in subML:
            l=list(info.split("_"))
            pep=l[3]
            lab=pep2mllabel[pep]
            if lab=="Multi_P7_other":
                continue
            r=prank[lab]
            if r >currank:
                maxRankSubClus=[]
                currank=r
                maxRankSubClus.append(info)
                b=1
            elif rank==currank:
                maxRankSubClus.append(info)
            else:
                continue
        if b==1:
            fsamp=getrand(maxRankSubClus)
            addsample.append(fsamp)
    mlinfokey={}
    mlinfokesI={}
    for info in addsample:
        spe, lab, pep=getinfo(info)
        ins=pep2info[pep]
        labs=set()
        setofinfo=set()
        for i in ins:
            ispe, ilab, ipep=getinfo(i)
            if ispe=="HUMAN" and ilab in sltu:
                labs.add(ilab)
                setofinfo.add(i)
        mlinfokesI[info]=setofinfo
        mlinfokey[info]=labs
    masterlist=[]
    for info in addsample:
        masterlist.append(info)
    for info in sampledInfo:
        masterlist.append(info)
    return masterlist, nonIntrestClusters, clus2info, nonHumanClus, rank , mlinfokey ##########3


def set_up_medium_negitive_data_set(nonIntrestClusters : list[str], clus2info : dict[str:list[str]], nonHumanClus : list[str], rank : dict[str:int], text_loc : str = "data/text_loc/" ):
    """
    samples medium data points (data points thata are postive PTM data but are not in predicted by this model)
    Args:
        nonIntrestClusters (list[str]) : list of cluster that contain only ptm types that fall outside the prediction task   
        clus2info (dir[str:set[str]]) : cluster id mapped to all peptide info 
        nonHumanClus (list[str]) : ptm seq ident clusters that only have non human peptides
        rank (dir[str:int]) : ptm type ranked by abundance
    Returns:
        
    """
    totalInfoInMNDS=[]
    for clus in nonIntrestClusters:
        for info in clus2info[clus]:
            totalInfoInMNDS.append(info)

    for clus in nonHumanClus:
        for info in clus2info[clus]:
            totalInfoInMNDS.append(info)
    mndsRes2count={}
    for info in totalInfoInMNDS:
        spe, lab, pep=getinfo(info)
        res=pep[10]
        if res not in mndsRes2count.keys():
            mndsRes2count[res]=0
        mndsRes2count[res]=mndsRes2count[res]+1
    resos=sorted(mndsRes2count.items(), key=lambda x: x[1], reverse=True)
    rrank={}
    for i,reso in enumerate(resos):
        rrank[reso[0]]=i+1
    sampledMNDS=[]
    for clus in nonIntrestClusters:
        maxRankSubClus=[]
        currank=0
        for info in clus2info[clus]:
            spe, lab, pep=getinfo(info)
            res=pep[10]
            r=rrank[res]
            if r >currank:
                maxRankSubClus=[]
                currank=r
                maxRankSubClus.append(info)
            elif rank==currank:
                maxRankSubClus.append(info)
            else:
                continue
        sampledMNDS.append(getrand(maxRankSubClus))
    for clus in nonHumanClus:
        maxRankSubClus=[]
        currank=0
        for info in clus2info[clus]:
            spe, lab, pep=getinfo(info)
            res=pep[10]
            r=rrank[res]
            if r >currank:
                maxRankSubClus=[]
                currank=r
                maxRankSubClus.append(info)
            elif rank==currank:
                maxRankSubClus.append(info)
            else:
                continue
        sampledMNDS.append(getrand(maxRankSubClus))
    mndsRes2count={}
    for info in sampledMNDS:
        spe, lab, pep=getinfo(info)
        res=pep[10]
        if res not in mndsRes2count.keys():
            mndsRes2count[res]=0
        mndsRes2count[res]=mndsRes2count[res]+1
    resos=sorted(mndsRes2count.items(), key=lambda x: x[1], reverse=True)
    rrank={}
    for i,reso in enumerate(resos):
        rrank[reso[0]]=i+1
    hf=open(f"{text_loc}medNegDataset.txt",'w+')
    for info in sampledMNDS:
        hf.write(info+'\n')
    hf.close()


def set_up_for_spec_clustering( masterlist : list[str], mlinfokey : dict[str:str], text_loc : str = "data/text_loc/"):
    """
    sets data for each PTM type to be clustered (SK-Learn's spectral clustering)

    Args:
        masterlist (list[str]) : list of all ptm data
        mlinfokey (dict[str:str]) : maps ptm data to label
        text_loc (str) : location of the text dir
    Returns:
    """
    lab2info={}
    lab2rescount={}
    lab2count={}
    newlab2info={}
    lab2error={}
    for info in masterlist:
        if info in mlinfokey.keys():
            for lab in mlinfokey[info]:
                if lab not in lab2info.keys():
                    lab2info[lab]=set()
                lab2info[lab].add(info)
        else:
            l=list(info.split('_'))
            lab=l[4]
            if lab not in lab2info.keys():
                lab2info[lab]=set()
            lab2info[lab].add(info)
    
    for lab in lab2info.keys():
        lab2count[lab]=len(lab2info[lab])
    slabs=sorted(lab2count.items(), key=lambda x: x[1], reverse=True)
    for lab in lab2info.keys():
        lab2rescount[lab]={}
        infos=lab2info[lab]
        for info in infos:
            l=list(info.split('_'))
            pep=l[3]
            res=pep[10]
            if res not in lab2rescount[lab].keys():
                lab2rescount[lab][res]=0
            lab2rescount[lab][res]=lab2rescount[lab][res]+1
    newlab2info["RK-Methylation"]=set()
    newlab2info["K-Ubiquitination"]=set()
    newlab2info["K-Acetylation"]=set()
    newlab2info["AM-Acetylation"]=set()
    newlab2info["K-Succinylation"]=set()
    newlab2info["K-Sumoylation"]=set()
    newlab2info["K-Malonylation"]=set()
    newlab2info["M-Sulfoxidation"]=set()
    newlab2info["ST-Phosphorylation"]=set()
    newlab2info["Y-Phosphorylation"]=set()
    newlab2info["ST-O-linked-Glycosylation"]=set()
    newlab2info["C-S-palmitoylation"]=set()
    newlab2info["C-Glutathionylation"]=set()
    newlab2info["N-N-linked-Glycosylation"]=set()
    newlab2info["PK-Hydroxylation"]=set()
    for lab in lab2info.keys():
        if "Methylation"==lab:
            for info in lab2info[lab]:
                newlab2info["RK-Methylation"].add(info)
        elif "Ubiquitination"==lab:
            for info in lab2info[lab]:
                newlab2info["K-Ubiquitination"].add(info)
        elif "Acetylation"==lab:
            for info in lab2info[lab]:
                l=list(info.split('_'))
                pep=l[3]
                res=pep[10]
                if 'A'==res or 'M'==res:
                    newlab2info["AM-Acetylation"].add(info)
                if 'K'==res:
                    newlab2info["K-Acetylation"].add(info)
        elif "Succinylation"==lab:
            for info in lab2info[lab]:
                newlab2info["K-Succinylation"].add(info)
        elif "Sumoylation"==lab:
            for info in lab2info[lab]:
                newlab2info["K-Sumoylation"].add(info)
        elif "Malonylation"==lab:
            for info in lab2info[lab]:
                newlab2info["K-Malonylation"].add(info)
        elif "Sulfoxidation"==lab:
            for info in lab2info[lab]:
                newlab2info["M-Sulfoxidation"].add(info)
        elif "Phosphorylation"==lab:
            for info in lab2info[lab]:
                l=list(info.split('_'))
                pep=l[3]
                res=pep[10]
                if 'S'==res or 'T'==res:
                    newlab2info["ST-Phosphorylation"].add(info)
                if 'Y'==res:
                    newlab2info["Y-Phosphorylation"].add(info)
        elif "O-linked-Glycosylation"==lab:
            for info in lab2info[lab]:
                newlab2info["ST-O-linked-Glycosylation"].add(info)  
        elif "S-palmitoylation"==lab:
            for info in lab2info[lab]:
                newlab2info["C-S-palmitoylation"].add(info)
        elif "Glutathionylation"==lab:
            for info in lab2info[lab]:
                newlab2info["C-Glutathionylation"].add(info)
        elif "N-linked-Glycosylation"==lab:
            for info in lab2info[lab]:
                newlab2info["N-N-linked-Glycosylation"].add(info)
        elif "Hydroxylation"==lab:
            for info in lab2info[lab]:
                newlab2info["PK-Hydroxylation"].add(info)
    lab2count={}
    for lab in newlab2info.keys():
        lab2count[lab]=len(newlab2info[lab])
    slabs=sorted(lab2count.items(), key=lambda x: x[1], reverse=True)
    lab2rescount={}
    for lab in newlab2info.keys():
        lab2rescount[lab]={}
        infos=newlab2info[lab]
        for info in infos:
            l=list(info.split('_'))
            pep=l[3]
            res=pep[10]
            if res not in lab2rescount[lab].keys():
                lab2rescount[lab][res]=0
            lab2rescount[lab][res]=lab2rescount[lab][res]+1
    for slab in slabs:
        posclus=math.floor(slab[1]/2000)
        if posclus>20:
            posclus=20
        if posclus==0:
            posclus=1
    for lab in lab2count.keys():
        ar=set()
        lab2error[lab]={}
        lab2error[lab]["good"]=[]
        lab2error[lab]["mlgood"]=[]
        lab2error[lab]["bad"]=[]
        for c in lab:
            if c=='-':
                break
            else:
                ar.add(c)
        plab=lab[len(ar)+1:]
        for info in newlab2info[lab]:
            spe, tlab, pep=getinfo(info)
            if tlab == plab and pep[10] in ar:
                lab2error[lab]["good"].append(info)
            else:
                b=0
                otherlabs=mlinfokey[info]
                for olab in otherlabs:
                    if olab==plab and pep[10] in ar:
                        b=1
                if b==1:
                    lab2error[lab]["mlgood"].append(info)
                else:
                    lab2error[lab]["bad"].append(info)
    cluster2infoC={}
    for lab in newlab2info.keys():
        for info in newlab2info[lab]:
            l=list(info.split('_'))
            clus=l[5]
            if clus not in cluster2infoC.keys():
                cluster2infoC[clus]=[]
            cluster2infoC[clus].append(info)
    b=0
    for clus in cluster2infoC.keys():
        samp=cluster2infoC[clus]
        sos=len(samp)
        if sos > 1:
            b1=0
            fsamp=samp[0]
            for info in samp:
                if fsamp!=info:
                    b1=1
            if b1==1:
                print(samp)

            if sos != len(mlinfokey[fsamp]):
                print(samp)
    if not os.path.isdir(f"{text_loc}infosforaffinitymat"):
        os.makedirs(f"{text_loc}infosforaffinitymat")
    n='\n'
    for lab in newlab2info.keys():
        hf=open(f"{text_loc}infosforaffinitymat/{lab}.txt",'w+')
        for info in newlab2info[lab]:
            hf.write(f"{info}{n}")    
        hf.close()


def prep_matlab_dist_mat_for_USSC(text_loc : str = "data/text_loc/" ):
    
    """
    sets up ST-Phos to be pluged into USSC clustered (done in matlab)
    SK-learns spec clustring runs out of ram even with 2TB of ram for ST-phos data
    to lazy to implement USSC in python or get sk-learns spec clustering to be run in hard-drive space 
    The feaure of amino acid volume size was chosen becaseue the umap proj of esm2 embeddings had the best "behavoir"
    (out of all feature tested volume was the cleanest (done by eye with umap proj sorry for no imperical evideance))
    
    Args:
        text_loc (str) : location of the text dir
    Returns:
    """
    amino="G 	A 	S 	C 	D 	P 	N 	T 	E 	V 	Q 	H 	M 	I 	L 	K 	R 	F 	Y 	W"
    sizes="60.1 	88.6 	89.0 	108.5 	111.1 	112.7 	114.1 	116.1 	138.4 	140.0 	143.8 	153.2 	162.9 	166.7 	166.7 	168.6 	173.4 	189.9 	193.6 	227.8"
    size={}
    l=list(amino.split(" 	"))
    am=l
    l=list(sizes.split(" 	"))
    hy=l
    for a,h in zip(am,hy):
        size[a]=float(h)
    size['-']=float(0)
    sizeData=[]
    ff=f"{text_loc}infosforaffinitymat/"
    folder = os.listdir(ff)
    for file in folder:
        type=file[:-4]
        if type != "ST-Phosphorylation":
            continue
        hf=open(f"{ff}{file}",'r')
        lines=hf.readlines()

        for line in lines:
            temp=[]
            l=list(line.split('_'))
            pep=l[3]
            for c in pep:
                temp.append(size[c])
            sizeData.append(temp)
    newinputdata=np.array(sizeData)
    savemat(f"{text_loc}infosforaffinitymat/data_ST-Phos_VolSize_Custom.mat", {"fea":newinputdata.astype(np.float64)})


#########################################################################################
#----------------------------------------------------------------------------------------
#  stage 4
#----------------------------------------------------------------------------------------
#########################################################################################


def prep_data_for_spec_clus(text_loc : str = "data/text_loc/"):
    """
    sets up data for each PTM type to be clustered (SK-Learn's spectral clustering) 
    creates affinity matrix for each ptm type

    Args:
        text_loc (str) : location of the text dir
    Returns:
        lab2peps (dict[str:list[str]]) : this maps ptm-type to the peptides assosated 
        lab2affm (dict[str:array[int]]) : this maps ptm-type to the affinity matrix 
        lab2maxnc (dict[str:int]) : this maps ptm-type to the max allowed number of clusters 
    """
    lab2peps={}
    lab2affm={}
    lab2maxnc={}
    folder = os.listdir(f"{text_loc}infosforaffinitymat/")
    for file in folder:
        if ".mat" in file:
            continue
        lab=file[:-4]
        lab2peps[lab]=[]
        if lab == "ST-Phosphorylation":
            continue
        hf=open(f"{text_loc}infosforaffinitymat/{file}",'r')
        lines=hf.readlines()
        for line in lines:
            l=list(line.split("_"))
            pep=l[3]
            lab2peps[lab].append(pep)
    for lab in lab2peps.keys():
        npeps=len(lab2peps[lab])
        nc=math.floor(npeps / 2000)
        if nc > 20:
            nc=20
        if nc <= 2:
            continue
        lab2maxnc[lab]=nc
    for lab in lab2maxnc.keys():
        npeps=len(lab2peps[lab])
        half_affm=np.zeros((npeps,npeps))
        peps=lab2peps[lab]
        strings=peps
        npeps=len(strings)
        arr = np.array([list(s) for s in strings])
        unique_chars, encoded = np.unique(arr, return_inverse=True)
        encoded = encoded.reshape(arr.shape)
        n = len(encoded)
        distmat = np.zeros((n, n), dtype=np.int8)
        for i in range(n):
            diffs = encoded != encoded[i]
            distmat[i] = np.sum(diffs, axis=1)
            print(f"{lab} {i} / {npeps}",end='\r')
            del diffs
        max_len = arr.shape[1]
        lab2affm[lab]=max_len - distmat
    return lab2peps, lab2affm, lab2maxnc


def spec_clus_sk_learn(lab2peps : dict[str:list[str]], lab2affm : dict[str:list[int]], lab2maxnc : dict[str:int], text_loc : str = "data/text_loc/"):
    """
    sets data for each PTM type to be clustered (SK-Learn's spectral clustering)

    Args:
        text_loc (str) : location of the text dir
        lab2peps (dict[str:list[str]]) : this maps ptm-type to the peptides assosated 
        lab2affm (dict[str:np.array[int]]) : this maps ptm-type to the affinity matrix 
        lab2maxnc (dict[str:int]) : this maps ptm-type to the max allowed number of clusters 
    Returns:
        
    """
    n="\n"
    if not os.path.isdir(f"{text_loc}speccluslabs"):
        os.makedirs(f"{text_loc}speccluslabs")
    for lab in lab2maxnc.keys():
        print(f"{lab} sk-learn SpectralClustering running might take some time",end='\r')
        peps=lab2peps[lab]
        propernc=bl2nc[lab]
        for nc in range( lab2maxnc[lab]):
            nc=nc+1
            if nc==propernc:
                if nc==1:
                    continue
                clustering = SpectralClustering(n_clusters=nc,affinity='precomputed', assign_labels='discretize', random_state=42).fit(lab2affm[lab])
                filename=f"{text_loc}speccluslabs/{lab}_nc{nc}.txt"
                hf=open(filename,'w+')
                for clus in clustering.labels_:
                    hf.write(f"{clus}{n}")
                hf.close()


#########################################################################################
#----------------------------------------------------------------------------------------
#  stage 5
#----------------------------------------------------------------------------------------
#########################################################################################


def sample_possilbe_easy_neg_data(fasta_loc : str = "data/fasta_loc/", text_loc : str = "data/text_loc/"):
    """
    creates possible easy negitve data points by sampling all 
    full protein seuqnces for poteintal 21 mers that belong to the positive data points

    Args:
        fasta_loc (str) : location of the fasta dir
        text_loc (str) : location of the text dir
    Returns:
        
    """
    uni2seq={}
    hf=open(f"{fasta_loc}uniprotPostMissingPepPreConRes.fasta",'r')
    lines=hf.readlines()
    for i,line in enumerate(lines):
        if line[0]=='>':
            uni=line[1:-1]
        elif line=='\n':            
            uni2seq[uni]=seq
        else:
            seq=line[:-1]
    hf=open(f"{fasta_loc}oldunisV2.fasta",'r')
    lines=hf.readlines()
    for i,line in enumerate(lines):
        if line[0]=='>':
            uni=line[1:-1]
        elif line=='\n':
            uni2seq[uni]=seq
        else:
            seq=line[:-1]
    hf=open(f"{fasta_loc}oldunis.fasta",'r')
    lines=hf.readlines()
    for i,line in enumerate(lines):
        if line[0]=='>':
            uni=line[1:-1]
        elif line=='\n':
            uni2seq[uni]=seq
        else:
            seq=line[:-1]
    hf=open(f"{fasta_loc}mmUniFull.fasta",'r')
    lines=hf.readlines()
    for i,line in enumerate(lines):
        if line[0]=='>':
            l=list(line.split('='))
            v=l[1][:-3]
            l=list(line.split(' '))
            uni=l[0][1:-1]
            uni=uni+'v'+v
            break
        elif line=='\n':
            uni2seq[uni]=seq
        else:
            seq=line[:-1]
    unilist=list(uni2seq.keys())
    res2c={}
    aa="ACDEFGHIKLMNPQRSTVWY"
    accetpAA=set()
    for c in aa:
        res2c[c]=set()
        accetpAA.add(c)
    r='\r'
    while True:
        b=1
        for c in aa:
            if len(res2c[c])!=220000:
                b=0   
        if b==1:
            break
        rn = random.randint(0, (len(unilist)-1))
        runi=unilist[rn]
        seq=uni2seq[runi]
        rn = random.randint(0, (len(seq)-1))
        info=runi+'_'+str(rn+1)
        paddedseq='----------'+seq+'----------'
        pep=paddedseq[rn:rn+21]
        res=pep[10]
        if res not in accetpAA:
            continue
        if len(res2c[res])== 220000:
            continue
        if len(pep)!=21:
            continue
        b1=0
        for c in pep:
            if c not in accetpAA:
                b1=1
        if b1==1:
            continue
        res2c[res].add(f"{pep}_{info}")
        sys.stdout.write(f"sampling easy neg {len(res2c['W'])} / 220000 {r}")
    hf=open(f"{text_loc}easyPepRaw.txt",'w+')
    n='\n'
    for c in aa:
        for pep in res2c[c]:
            hf.write(f"_{pep}_{n}")
            
    hf.close()


#########################################################################################
#----------------------------------------------------------------------------------------
#  stage 6
#----------------------------------------------------------------------------------------
#########################################################################################


def sample_from_spec_clus(csv_loc : str = "data/csv_loc/", text_loc : str = "data/text_loc/"):
    """
    orginizes data coming form spectral clustering 
    prioritizes multi-lable data

    Args:
        csv_loc (str) : location of the csv dir
        text_loc (str) : location of the text dir
    Returns:
        pepsused (set(str)) : petides that have multi-labesl assosated 
        mLcluslab2pd (dict[str:set(str)]) : keeps track of what peptides belong to multiple lables
        cluslab2pd (dict[str:set(str)]) : all data assosated with each spec cluster 
        cluslabs (list[str]) : list of spec clus names 
        cluslab2peps (dict[str:set(str)) : all data assosated with each spec cluster
        notclustered (list[str]) : all data not clustered
        
    """
    cluslabs=[]
    cluslab2pd={}
    cluslab2peps={}
    bl2pep2c={}
    labC2c={}
    pep2labs={}
    labset=set()
    mLcluslab2pd={}
    pepsused=set()
    pepnotused=[]
    notclustered=[]
    nuratio={}
    NRcluslab2pdF={}
    posSet=set()
    res2count={}
    for lab in baselabs:
        nc=bl2nc[lab]
        for i in range(nc):
            cluslab=f"{lab}_nc{i}_tot{nc}"
            cluslabs.append(cluslab)
    
    for cluslab in cluslabs:
        cluslab2pd[cluslab]=set()
    for bl in baselabs:
        nc=bl2nc[bl]
        if bl!='ST-Phosphorylation':
            if nc > 1:
                file=f"{text_loc}speccluslabs/{bl}_nc{nc}.txt"
                print(file)
                hf=open(file,'r')
                lines1=hf.readlines()
                file=f"{text_loc}infosforaffinitymat/{bl}.txt"
                print(file)
                hf=open(file,'r')
                lines2=hf.readlines()
                print(len(lines1),len(lines2))
                bl2pep2c[bl]={}
                for line1,line2 in zip(lines1,lines2):
                    c=line1[:-1]
                    info=line2[:-1]
                    l=list(info.split('_'))
                    bl2pep2c[bl][l[3]]=int(c)
                    cluslab=f"{bl}_nc{c}_tot{nc}"
                    if cluslab not in cluslab2peps.keys():
                        cluslab2peps[cluslab]=set()
                    cluslab2peps[cluslab].add(l[3])
                    #nc,c,bl                
            else:
                file=f"{text_loc}infosforaffinitymat/{bl}.txt"
                print(file)
                hf=open(file,'r')
                lines=hf.readlines()
                #print(file)
                bl2pep2c[bl]={}
                for line in lines:
                    #print(line)
                    l=list(line[:-1].split('_'))
                    pep=l[3]
                    bl2pep2c[bl][pep]=0
                    cluslab=f"{bl}_nc{0}_tot{nc}"
                    if cluslab not in cluslab2peps.keys():
                        cluslab2peps[cluslab]=set()
                    cluslab2peps[cluslab].add(l[3])
        else:
            file=f"{text_loc}ussc_output/testdataNC-{nc}.txt"
            hf=open(file,'r')
            lines1=hf.readlines()
            print(file)
            file=f"{text_loc}infosforaffinitymat/{bl}.txt"
            hf=open(file,'r')
            lines2=hf.readlines()
            print(file)
            print(len(lines1),len(lines2))
            bl2pep2c[bl]={}
            for line1,line2 in zip(lines1,lines2):
                    c=line1[:-1]
                    info=line2[:-1]
                    l=list(info.split('_'))
                    bl2pep2c[bl][l[3]]=int(c)-1
                    cluslab=f"{bl}_nc{int(c)-1}_tot{nc}"
                    if cluslab not in cluslab2peps.keys():
                        cluslab2peps[cluslab]=set()
                    cluslab2peps[cluslab].add(l[3]) 
    
    for lab in bl2pep2c.keys():
        for pep in bl2pep2c[lab].keys():
            c=bl2pep2c[lab][pep]
            labC=lab+'_'+str(c)
            if labC not in labC2c.keys():
                labC2c[labC]=0
            labC2c[labC]=labC2c[labC]+1
    for lab in labC2c.keys():
        print(lab,labC2c[lab])
    hf=open(f"{csv_loc}pep2clusname.csv",'w+')
    for lab in bl2pep2c.keys():
        for pep in bl2pep2c[lab].keys():
            # print(lab)
            # break
            c=bl2pep2c[lab][pep]
            nc=bl2nc[lab]
            fl=f"{lab}_nc{c}_tot{nc}".replace(' ','')
            hf.write(f"{pep},{fl}"+'\n')
    hf.close()
    hf=open(f"{text_loc}mulLab.txt",'r')
    lines=hf.readlines()
    for line in lines:
        l=list(line[:-1].split("_"))
        pep=l[0]
        res=pep[10]
        for lab in l[1:]:
            lab=res+'-'+lab
            labset.add(lab)
            if lab =='T-O-linked-Glycosylation' or lab == "S-O-linked-Glycosylation":
                lab="ST-O-linked-Glycosylation"
                    
            if lab =='M-Acetylation':
                lab="AM-Acetylation"
            if lab =='S-Phosphorylation' or lab=="T-Phosphorylation":
                lab="ST-Phosphorylation"
            if lab =='K-Methylation':
                lab="RK-Methylation"
            if lab =='K-Hydroxylation':
                lab="PK-Hydroxylation"
            if pep not in pep2labs.keys():
                pep2labs[pep]=set()
            pep2labs[pep].add(lab)
    print(labset)
    for cluslab in cluslabs:
        mLcluslab2pd[cluslab]=set()
    print(mLcluslab2pd)
    i=0
    for pep in pep2labs.keys():
        labs=pep2labs[pep]
        peppass=0
        fullLabs=[]
        for lab in labs:
            if lab=='K-Succinylation':
                continue
            nc=bl2nc[lab]
            if bl2nc[lab]==1:
                c=0
            else:
                if pep not in bl2pep2c[lab].keys():
                    peppass=1
                    break
                c=bl2pep2c[lab][pep]
            fl=f"{lab}_nc{c}_tot{nc}".replace(' ','')
            fullLabs.append(fl)
        if peppass ==1:
            notclustered.append([pep,lab,labs])
        b=0
        catlab=""
        for lab in fullLabs:
            catlab=catlab+lab
            if 500 <= len(mLcluslab2pd[lab]):
                b=1
        if "Succinylation" in catlab or "Malonylation" in catlab or "Sumoylation" in catlab:
            b=0
            for lab in fullLabs:
                if 1500<= len(mLcluslab2pd[lab]):
                    b=1
        if b==1:
            for lab in fullLabs:
                if lab not in nuratio.keys():
                    nuratio[lab]=set()
                nuratio[lab].add(pep)

            pepnotused.append(pep)
            continue
        for lab in fullLabs:
            mLcluslab2pd[lab].add(pep)
            pepsused.add(pep)
    print(len(notclustered),"notmapped")
    print(len(pepsused),"used",)
    print(len(pepnotused),"not used",)
    for lab in mLcluslab2pd.keys():
        print(lab,len(mLcluslab2pd[lab]))
    for lab in nuratio.keys():
        print(lab,len(nuratio[lab]))
    for cluslab in cluslabs:
        NRcluslab2pdF[cluslab]=[]
    for bl in bl2pep2c.keys():
        for pep in bl2pep2c[bl].keys():
            if pep in pep2labs.keys():
                continue
            c=bl2pep2c[bl][pep]
            nc=bl2nc[bl]
            fl=f"{bl}_nc{c}_tot{nc}".replace(' ','')
            NRcluslab2pdF[fl].append(pep)
    for fl in NRcluslab2pdF.keys():
        print(fl,len(NRcluslab2pdF[fl]))
    cluslab2pd={}
    for cluslab in cluslabs:
        cluslab2pd[cluslab]=set()
    for fl in NRcluslab2pdF.keys():
        samp = NRcluslab2pdF[fl]
        random.shuffle(samp)
        for pep in samp:
            if 2000==len(mLcluslab2pd[fl])+len(cluslab2pd[fl]):
                break
            cluslab2pd[fl].add(pep)
    for lab in mLcluslab2pd.keys():
        print(lab,len(mLcluslab2pd[lab])+len(cluslab2pd[lab]),len(cluslab2pd[lab]))
    for lab in mLcluslab2pd.keys():
        fullset=set()
        for pep in mLcluslab2pd[lab]:
            fullset.add(pep)
        for pep in cluslab2pd[lab]:
            fullset.add(pep)
        print(lab,len(mLcluslab2pd[lab])+len(cluslab2pd[lab]),len(cluslab2pd[lab]),len(fullset))
    for lab in mLcluslab2pd.keys():
        fullset=set()
        for pep in mLcluslab2pd[lab]:
            posSet.add(pep)
            fullset.add(pep)
        for pep in cluslab2pd[lab]:
            posSet.add(pep)
            fullset.add(pep)
        print(lab,len(fullset))
    print(len(posSet))
    for pep in posSet:
        res=pep[10]
        if res not in res2count.keys():
            res2count[res]=0
        res2count[res] = res2count[res]+1
    for res in res2count.keys():
        print(res,res2count[res])
    return pepsused, mLcluslab2pd, cluslab2pd, cluslabs, cluslab2peps, notclustered


def sample_negtive_data(csv_loc : str = "data/csv_loc/", text_loc : str = "data/text_loc/"):
    """
    orginizes data coming form spectral clustering 
    prioritizes multi-lable data

    Args:
        res2count (dict[str:int]) : maps residues type to count of abundance
        csv_loc (str) : location of the csv dir
        text_loc (str) : location of the text dir
    Returns:
        negres2peps (dict[char:set(str)]) : mapping modified residue to peptides
        negSet (set(str)) : set of all negtive data 
    """
    negSet=set()
    negres2peps={}
    MNdataset=set()
    hf=open(f"{text_loc}medNegDataset.txt",'r')
    lines=hf.readlines()
    for line in lines:
        l=list(line[:-1].split('_'))
        pep=l[3]
        res=pep[10]
        if res not in negres2peps.keys():
            negres2peps[res]=set()
        negres2peps[res].add(pep)
        MNdataset.add(pep)
        negSet.add(pep)
    print(len(negSet))
    for res in negres2peps.keys():
        print(res,len(negres2peps[res]))
    negSet=set() 
    usedSet=set()
    ENdataset=set()
    res2count={}
    hf=open(f"{text_loc}clus7_easyneg.txt",'r')
    lines=hf.readlines()
    for line in lines:
        l=list(line[:-1].split('_'))
        pep=l[1]
        clus=l[3]
        if clus not in usedSet:
            res=pep[10]
            if res not in res2count.keys():
                pc=0
            else:
                pc=res2count[res]
            if 110000 <= pc + len(negres2peps[res]):
                continue
            negres2peps[res].add(pep)
            usedSet.add(clus)
            ENdataset.add(pep)
            negSet.add(pep)
    return negres2peps, negSet
    

def data_split_train_test_val(pepsused : set[str], mLcluslab2pd : dict[str:set[str]], cluslab2pd : dict[str:set[str]],
                              negres2peps : dict[str:set[str]], negSet : set[str], cluslabs : list[str],
                              cluslab2peps : dict[str:set[str]], notclustered : list[str], csv_loc : str = "data/csv_loc/"):
    """
    split data into training testing and validation
    makes sure ratios of multilabel , residue type, and postive/negtive 
    is perserved throught the split 
    (prob could do this in a few lines of code with sk-learn/ pandas 
    but want it to be done manual)

    Args:
        pepsused (set(str)) : petides that have multi-labesl assosated 
        mLcluslab2pd (dict[str:set(str)]) : keeps track of what peptides belong to multiple lables
        cluslab2pd (dict[str:set(str)]) : all data assosated with each spec cluster 
        negres2peps (dict[str:set(str)]) : mapping modified residue to peptides
        negSet (set(str)) : set of all negtive data
        cluslabs (list[str]) : list of spec clus names 
        cluslab2peps (dict[str:set(str)]) : all data assosated with each spec cluster
        notclustered (list[str]) : all data not clustered
        csv_loc (str) : location of the csv dir
    Returns:
    """
    train,pretest=train_test_split( list(pepsused) , random_state=42,test_size=0.3, shuffle=True)
    test,val=train_test_split( pretest , random_state=42,test_size=0.5, shuffle=True)
    broken=set()
    telab2c={}
    tco=0
    for lab in mLcluslab2pd.keys():
        co=0
        for pep in test:
            if pep in mLcluslab2pd[lab]:
                co=co+1
                tco=tco+1
        telab2c[lab]=co
    for lab in telab2c.keys():
        print(lab,telab2c[lab])
    valab2c={}
    for lab in mLcluslab2pd.keys():
        co=0
        for pep in val:
            if pep in mLcluslab2pd[lab]:
                co=co+1
        valab2c[lab]=co
    for lab in valab2c.keys():
        print(lab,valab2c[lab])
    trlab2c={}
    for lab in mLcluslab2pd.keys():
        co=0
        for pep in train:
            if pep in mLcluslab2pd[lab]:
                co=co+1
        trlab2c[lab]=co
    for lab in trlab2c.keys():
        print(lab,trlab2c[lab])
    for lab in trlab2c.keys():
        trl=trlab2c[lab]
        tel=telab2c[lab]
        vall=valab2c[lab]
        tot=trl+tel+vall
        if tot!=0:
            print(lab,"test ratio",tel/tot,"train ratio",trl/tot,"val ratio",vall/tot)
    r2c={}
    for pep in test:
        res=pep[10]
        if res not in r2c.keys():
            r2c[res]=0
        r2c[res]=r2c[res]+1
    for res in r2c.keys():
        print("ftest2",res,r2c[res])
    r2c={}
    for pep in train:
        res=pep[10]
        if res not in r2c.keys():
            r2c[res]=0
        r2c[res]=r2c[res]+1
    for res in r2c.keys():
        print("ftrain2",res,r2c[res])
    r2c={}
    for pep in val:
        res=pep[10]
        if res not in r2c.keys():
            r2c[res]=0
        r2c[res]=r2c[res]+1
    for res in r2c.keys():
        print("fval2",res,r2c[res])
    ftrain=set(train)
    ftest=set(test)
    fval=set(val)
    for lab in cluslab2pd.keys():
        fpd=set()
        pd=cluslab2pd[lab]
        mpd=mLcluslab2pd[lab]
        for pep in mpd:
            fpd.add(pep)
        for pep in pd:
            fpd.add(pep)
        print(len(mpd),len(pd))
        atrain,pretest=train_test_split( list(pd) , random_state=42,test_size=0.3, shuffle=True)
        atest,aval=train_test_split( pretest , random_state=42,test_size=0.5, shuffle=True)   
        tot=len(atrain)+len(atest)+len(aval)
        print(lab,len(fpd),len(atrain)/tot,len(atest)/tot,len(aval)/tot)
        for pep in atrain:
            ftrain.add(pep)
        for pep in atest:
            ftest.add(pep)
        for pep in aval:
            fval.add(pep)
        print(len(fval),len(fval),len(ftrain))
    tot=len(ftrain)+len(ftest)+len(fval)
    print(len(ftrain)/tot)
    print(len(ftest)/tot)
    print(len(fval)/tot)
    r2c={}
    for pep in ftest:
        res=pep[10]
        if res not in r2c.keys():
            r2c[res]=0
        r2c[res]=r2c[res]+1
    for res in r2c.keys():
        print("ftest2",res,r2c[res])
    r2c={}
    for pep in ftrain:
        res=pep[10]
        if res not in r2c.keys():
            r2c[res]=0
        r2c[res]=r2c[res]+1
    for res in r2c.keys():
        print("ftrain2",res,r2c[res])
    r2c={}
    for pep in fval:
        res=pep[10]
        if res not in r2c.keys():
            r2c[res]=0
        r2c[res]=r2c[res]+1
    for res in r2c.keys():
        print("fval2",res,r2c[res])
    ftrain2=ftrain
    ftest2=ftest
    fval2=fval
    for res in negres2peps.keys():
        rtrain,pretest=train_test_split( list(negres2peps[res]) , random_state=42,test_size=0.3, shuffle=True)
        rtest,rval=train_test_split( pretest , random_state=42,test_size=0.5, shuffle=True)
        tot=len(rtrain)+len(rtest)+len(rval)
        print(res,len(rtrain)/tot,len(rtest)/tot,len(rval)/tot)
        for pep in rtrain:
            ftrain2.add(pep)
            negSet.add(pep)
        for pep in rtest:
            ftest2.add(pep)
            negSet.add(pep)
        for pep in rval:
            fval2.add(pep)
            negSet.add(pep)
    tot=len(ftrain2)+len(ftest2)+len(fval2)
    print(len(ftrain2)/tot)
    print(len(ftest2)/tot)
    print(len(fval2)/tot)
    r2c={}
    newnegset=set()
    for pep in ftest2:
        res=pep[10]
        if res not in r2c.keys():
            r2c[res]=0
        r2c[res]=r2c[res]+1
    for res in r2c.keys():
        print("ftest2",res,r2c[res])
    r2c={}
    for pep in ftrain2:
        res=pep[10]
        if res not in r2c.keys():
            r2c[res]=0
        r2c[res]=r2c[res]+1
    for res in r2c.keys():
        print("ftrain2",res,r2c[res])
    r2c={}
    for pep in fval2:
        res=pep[10]
        if res not in r2c.keys():
            r2c[res]=0
        r2c[res]=r2c[res]+1
    for res in r2c.keys():
        print("fval2",res,r2c[res])
    Trainpeps=[]
    Trainlabs=[]
    pc=0
    for i,pep in enumerate(ftrain2):
        Trainpeps.append(pep)
        lab=[]
        b=0
        for ll in cluslabs:
            if pep in cluslab2peps[ll] or pep in mLcluslab2pd[ll]:
                lab.append(1)
                b=1
            else:
                lab.append(0)
        if b==1:
            pc=pc+1
        if pep in negSet:
            lab.append(1)
        else:
            lab.append(0)
        Trainlabs.append(lab)
    nolabs=set()
    for labs,pep in zip(Trainlabs,Trainpeps):
        b=1
        for lab in labs:
            if lab==0:
                continue
            if lab ==1:
                b=0
        if b==1:
            nolabs.add(pep)
    print(len(nolabs))
    posnolab=set()
    for pep in nolabs:
        b=0
        if pep in pep2labs.keys():
            posnolab.add(pep)
    print(len(posnolab))
    cnc=0
    nc=set()
    for n in notclustered:
        nc.add(n[0])
    for pep in nolabs:
        if pep in nc:
            cnc+cnc
    print(cnc)
    hf=open(f"{csv_loc}train_hd3_CustBL62SeqDistSpecClus_uniResRatio_CustNegRaio.csv",'w+')
    line="pep"
    for cl in cluslabs:
        line=line+','+cl
    line=line+',NegLab\n'
    hf.write(line)
    for pep , lab in zip(Trainpeps,Trainlabs):
        line=pep
        for l in lab:
            line=line+','+str(l)
        hf.write(line+'\n')
    hf.close()
    Trainpeps=[]
    Trainlabs=[]
    pc=0
    for i,pep in enumerate(ftest2):
        Trainpeps.append(pep)
        lab=[]
        b=0
        for ll in cluslabs:
            if pep in cluslab2peps[ll] or pep in mLcluslab2pd[ll]:
                lab.append(1)
                b=1
            else:
                lab.append(0)
        if b==1:
            pc=pc+1
        if pep in negSet:
            lab.append(1)
        else:
            lab.append(0)
        Trainlabs.append(lab)        
    nolabs=set()
    for labs,pep in zip(Trainlabs,Trainpeps):
        b=1
        for lab in labs:
            if lab==0:
                continue
            if lab ==1:
                b=0
        if b==1:
            nolabs.add(pep)
    print(len(nolabs))
    hf=open(f"{csv_loc}test_hd3_CustBL62SeqDistSpecClus_uniResRatio_CustNegRaio.csv",'w+')
    line="pep"
    for cl in cluslabs:
        line=line+','+cl
    line=line+',NegLab\n'
    hf.write(line)
    for pep , lab in zip(Trainpeps,Trainlabs):
        line=pep
        for l in lab:
            line=line+','+str(l)
        hf.write(line+'\n')
    hf.close()
    Trainpeps=[]
    Trainlabs=[]
    for i,pep in enumerate(fval2):
        Trainpeps.append(pep)
        lab=[]
        for ll in cluslabs:
            if pep in cluslab2peps[ll] or pep in mLcluslab2pd[ll]:
                lab.append(1)
            else:
                lab.append(0)
        if pep in negSet:
            lab.append(1)
        else:
            lab.append(0)
        Trainlabs.append(lab)
    nolabs=set()
    for labs,pep in zip(Trainlabs,Trainpeps):
        b=1
        for lab in labs:
            if lab==0:
                continue
            if lab ==1:
                b=0
        if b==1:
            nolabs.add(pep)
    print(len(nolabs))
    hf=open(f"{csv_loc}val_hd3_CustBL62SeqDistSpecClus_uniResRatio_CustNegRaio.csv",'w+')
    line="pep"
    for cl in cluslabs:
        line=line+','+cl
    line=line+',NegLab\n'
    hf.write(line)
    for pep , lab in zip(Trainpeps,Trainlabs):
        line=pep
        for l in lab:
            line=line+','+str(l)
        hf.write(line+'\n')
    hf.close()


def reduce_neg_ratio(csv_loc : str = "data/csv_loc/"):
    """
    reduce negtive ratio from even dist bewteen residue types across all data but 1/1000 ratio pos/ratio
    to 1/55 ratio  

    Args:
        csv_loc (str) : location of the csv dir
    Returns:
    """
    hf=open(f"{csv_loc}train_hd3_CustBL62SeqDistSpecClus_uniResRatio_CustNegRaio.csv",'r')
    lines=hf.readlines()
    sline=lines[0]
    keeplines=[]
    pos=[]
    neg=[]
    posc=[0]*54
    for line in lines[1:]:
        l=list(line[:-1].split(','))
        for i,lab in enumerate(l[1:]):
            if lab=='1':
                posc[i]=posc[i]+1
        neglab=int(l[-1])
        if neglab==1:
            neg.append(line)
        else:
            pos.append(line)
    print(len(pos),len(neg))

    aa="ACDEFGHIKLMNPQRSTVWY"
    random.seed(42)
    random.shuffle(neg)
    a2c={}
    for a in aa:
        a2c[a]=0
    b=0
    x=0
    totNegExampels=1400*1
    numPerRes=totNegExampels/20
    while(b==0):
        line=neg[x]
        l=list(line.split(','))
        pep=l[0]
        res=pep[10]
        if a2c[res]!=numPerRes:
            keeplines.append(line)
            a2c[res]=a2c[res]+1
        x=x+1
        for a in aa:
            if a2c[a]==numPerRes:
                b=1
            else:
                b=0
                break
    print(len(keeplines))
    for line in pos:
        keeplines.append(line)        
    print(len(keeplines))
    hf=open(f"{csv_loc}train_hd3_CustBL62SeqDistSpecClus_uniResRatio_1to60NegRaio.csv",'w+')
    hf.write(sline)
    for line in keeplines:
        hf.write(line)
    hf.close()
    hf=open(f"{csv_loc}val_hd3_CustBL62SeqDistSpecClus_uniResRatio_CustNegRaio.csv",'r')
    lines=hf.readlines()
    sline=lines[0]
    keeplines=[]
    pos=[]
    neg=[]
    for line in lines[1:]:
        #print(line)
        l=list(line[:-1].split(','))
        
        neglab=int(l[-1])
        if neglab==1:
            neg.append(line)
        else:
            pos.append(line)
    print(len(pos),len(neg))

    aa="ACDEFGHIKLMNPQRSTVWY"
    random.seed(42)
    random.shuffle(neg)
    a2c={}
    for a in aa:
        a2c[a]=0
    b=0
    x=0
    while(b==0):
        line=neg[x]
        l=list(line.split(','))
        pep=l[0]
        res=pep[10]
        if a2c[res]!=100:
            keeplines.append(line)
            a2c[res]=a2c[res]+1
        x=x+1
        for a in aa:
            if a2c[a]==100:
                b=1
            else:
                b=0
                break
    print(len(keeplines))
    for line in pos:
        keeplines.append(line)        
    print(len(keeplines))
    hf=open(f"{csv_loc}val_hd3_CustBL62SeqDistSpecClus_uniResRatio_1to1NegRaio.csv",'w+')
    hf.write(sline)
    for line in keeplines:
        hf.write(line)
    hf.close()
    hf=open(f"{csv_loc}test_hd3_CustBL62SeqDistSpecClus_uniResRatio_CustNegRaio.csv",'r')
    lines=hf.readlines()
    sline=lines[0]
    keeplines=[]
    pos=[]
    neg=[]
    for line in lines[1:]:
        l=list(line[:-1].split(','))
        neglab=int(l[-1])
        if neglab==1:
            neg.append(line)
        else:
            pos.append(line)
    print(len(pos),len(neg))
    aa="ACDEFGHIKLMNPQRSTVWY"
    random.seed(42)
    random.shuffle(neg)
    a2c={}
    for a in aa:
        a2c[a]=0
    b=0
    x=0
    while(b==0):
        line=neg[x]
        l=list(line.split(','))
        pep=l[0]
        res=pep[10]
        if a2c[res]!=100:
            keeplines.append(line)
            a2c[res]=a2c[res]+1
        x=x+1
        for a in aa:
            if a2c[a]==100:
                b=1
            else:
                b=0
                break
    print(len(keeplines))
    for line in pos:
        keeplines.append(line)        
    print(len(keeplines))
    hf=open(f"{csv_loc}test_hd3_CustBL62SeqDistSpecClus_uniResRatio_1to1NegRaio.csv",'w+')
    hf.write(sline)
    for line in keeplines:
        hf.write(line)
    hf.close()
