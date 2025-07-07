import os
import sys
import subprocess
import tarfile 
import requests
import threading
import queue
import time
import json
import gzip
import wget

from selenium import webdriver
from selenium.webdriver.support.ui import Select
from selenium.webdriver.common.by import By

from bs4 import BeautifulSoup
from zipfile import ZipFile 

from unipressed import UniprotkbClient
from unipressed import IdMappingClient



pwd = os.getcwd()

dbptmlist=["ADP-ribosylation",
    "AMPylation",
    "Acetylation",
    "Amidation",
    "Biotinylation",
    "Blocked amino end",
    "Butyrylation",
    "C-linked Glycosylation",
    "Carbamidation",
    "Carboxyethylation",
    "Carboxylation",
    "Cholesterol ester",
    "Citrullination",
    "Crotonylation",
    "D-glucuronoylation",
    "Deamidation",
    "Deamination",
    "Decanoylation",
    "Decarboxylation",
    "Dephosphorylation",
    "Disulfide bond	",
    "Farnesylation",
    "Formation of an isopeptide bond",
    "Formylation",
    "GPI-anchor",
    "Gamma-carboxyglutamic acid",
    "Geranylgeranylation",
    "Glutarylation",
    "Glutathionylation",
    "Hydroxyceramide ester",
    "Hydroxylation",
    "Iodination",
    "Lactoylation",
    "Lactylation",
    "Lipoylation",
    "Malonylation",
    "Methylation",
    "Myristoylation",
    "N-carbamoylation",
    "N-linked Glycosylation",
    "N-palmitoylation",
    "Neddylation",
    "Nitration",
    "O-linked Glycosylation",
    "O-palmitoleoylation",
    "O-palmitoylation",
    "Octanoylation",
    "Oxidation",
    "Phosphatidylethanolamine amidation",
    "Phosphorylation",
    "Propionylation",
    "Pyrrolidone",
    "Pyrrolylation",
    "Pyruvate",
    "S-Cyanation",
    "S-archaeol",
    "S-carbamoylation",
    "S-cysteinylation",
    "S-diacylglycerol",
    "S-linked Glycosylation",
    "S-nitrosylation",
    "S-palmitoylation",
    "Serotonylation",
    "Stearoylation",
    "Succinylation",
    "Sulfation",
    "Sulfhydration",
    "Sulfoxidation",
    "Sumoylation",
    "Thiocarboxylation",
    "UMPylation",
    "Ubiquitination",]




def download_db_PTM( dbdir : str = "data/dbPTMloc/" ):
    """
    Scrapes dbPTM and downloads all data avalable by PTM type 

    Args:
        dbdir (int): directory where it places all data.

    """
    # chrome_options = webdriver.ChromeOptions()
    # prefs = {'download.default_directory' : pwd+f"/{dbdir}"}
    # chrome_options.add_experimental_option('prefs', prefs)
    # driver = webdriver.Chrome(options=chrome_options)
    # time.sleep(1)
    print("*****Downloading_from_dbPTM*****")
    for sec in dbptmlist:
        try:
            if os.path.exists(dbdir+sec+".gz"):
                continue
            print("  "+sec,end='\r')
            link="https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/"+sec+".gz"
            #print(link)
            filename = wget.download(link, out=dbdir)
        except:
            print("  ")
            print(f"{sec} was skipped due to a 403 error")
            print("  ")
            continue
        #river.get(link)
        #time.sleep(5)
    #time.sleep(5)






def find_missing_data( masterlist : list[str], csv_loc : str ="data/csv_loc/" ):
    """
    finds missing data with the intial pass through 
    Args:
        masterlist (list[str]) : masterlist of rawPTM data
    Returns:
        missingpeps (list[str]): peps that don't seem to map well
        missingunis (list[str]): uniprots that don't seem to exist
        uni2seq (dict[str:str]) : uniprot id 2 protien sequence 
        masterlist (list[str]) : updated masterlist of rawPTM data
    """
    missingpeps=[]
    missingunis=set()
    uni2seq={}

    #   ------------------------------------ find any peps that are empty (pep info == "")
    for info in masterlist:
        #print(info)
        l=list(info.split('_'))
        #print(l)
        pep=l[3]
        uni=l[0]
        if pep == "":
            missingunis.add(uni)
            missingpeps.append(info)


    #   -----------------------------------  fill in missing peps with positions retrived from unipressed(uniprot) 
    
    x=0
    temps=[]
    temp=[]
    for i,uni in enumerate(missingunis):
        if i !=0 and i%500==0:
            temps.append(temp)
            temp=[]
        temp.append(uni)
    temps.append(temp)
    #print(len(missingunis))

    for temp in temps:
        #print(len(temp))
        reqs=UniprotkbClient.fetch_many(temp)
        for req in reqs:
            try:
                uni=req['primaryAccession']
                seq=req['sequence']['value']
                uni2seq[uni]=seq
            except:
                continue
        x=x+1
        sys.stdout.write("scraping uniprot " + str(x) + "/"+str(len(temps))+"\r")
    skipped=0
    updatedlist=[]
    for iter1,info in enumerate(masterlist):
        #print(info)
        l=list(info.split('_'))
        #print(l)
        pep=l[3]
        uni=l[0]
        pos=int(l[1])
        lab=l[4]
        ung=l[2]


        if pep == "":
            if uni not in uni2seq.keys():
                #print(uni)
                skipped+=1
                continue
            seq=uni2seq[uni]
            paddedseq="----------"+seq+"----------"
            pep=paddedseq[pos-1:+pos+20],lab
            info=f"{uni}_{str(pos)}_{ung}_{pep}_{lab}"
            updatedlist.append(info)
            
        else:
            updatedlist.append(info)
        
    #print(skipped)
    #print(len(masterlist),len(updatedlist))
    masterlist=updatedlist
    return missingpeps, missingunis, uni2seq, masterlist







def download_fasta_from_unis_v1( uniset: list[str] , fasta_loc : str ="data/fasta_loc/"  ):
    ############################### used unipressed to get all sequences 
    """
    queries Uniprot Rest API and gets proten seqs for all uniprot IDs that are up-to-date

    Args:
        uniset (list[str]): non-redundant of list of uniprots
        fasta_loc (int): directory where it places all data.
    Returns:
        Writes a fasta file 
        uni2seq (dict[str:str]) : uniprot id 2 sequences 
    """
    uni2seq={}
    temps=[]
    temp=[]
    for i,uni in enumerate(uniset):
        if i!=0 and i%500==0:
            temps.append(temp)
            temp=[]
        temp.append(uni)
            
    temps.append(temp)
    #print(len(temps))
    #print(len(temps[0]))
    x=0
    for temp in temps:
        reqs=UniprotkbClient.fetch_many(temp)
        for req in reqs:
            try:
                uni=req['primaryAccession']
                seq=req['sequence']['value']
                uni2seq[uni]=seq
            except:
                continue
        x=x+1
        sys.stdout.write("scraping uniprot " + str(x) + "/"+str(len(temps))+"\r")
    nl='\n'
    
    output=open(f"{fasta_loc}uniprotPostMissingPepPreConRes.fasta",'w+')
    for uni in uni2seq.keys():
        output.write(f'>{uni}{nl}{uni2seq[uni]}{nl}{nl}')

    output.close()
    return uni2seq












def map_out_uniprotIDs_to_uniParc_v1( oldunis: list[str] , csv_loc : str ="data/csv_loc/"  ):
    ############################### used unipressed to get all sequences 
    """
    queries Uniprot Rest API and maps uniprot id to uniparc id

    Args:
        oldunis (list[str]): non-redundant of list of uniprots
        csv_loc (int): directory where it places all data.
    Returns:
        Writes a csv mapping file 
        uni2uniparc (dict[str:str]) : uniprot id 2 uniparc id
    """
    temps=[]
    temp=set()
    for i,uni in enumerate(oldunis):
        if i!=0 and i%500==0:
            temps.append(temp)
            temp=set()
        temp.add(uni)

    uni2uniparc={}
    c=0

    for temp in temps:

        c=c+1
        request = IdMappingClient.submit(
            source="UniProtKB_AC-ID", dest="UniParc", ids=temp
        )
        x=0
        b=1

        while b==1:
            try:
                x=x+1
                r='\r'
                sys.stdout.write(f"scraping uniprot {str(x)} seconds waiting on batch {str(c)} / {str(len(temps))} uniparc can be slow to respond {r}")
                time.sleep(1)
                u2u=list(request.each_result())
                for nmid in u2u:
                    uni=nmid["from"]
                    uniprac=nmid["to"]
                    uni2uniparc[uni]=uniprac
                b=0
            except:
                b=1
    hf=open(f"{csv_loc}uni2uniparc.csv","w+")
    n='\n'
    for uni in uni2uniparc.keys():
        hf.write(f"{uni},{uni2uniparc[uni]}{n}")
    hf.close()
    return uni2uniparc






def download_fasta_from_parc_v1( uni2uniparc : dict[str:str] , fasta_loc : str ="data/fasta_loc/"  ):
    ############################### used unipressed to get all sequences 
    """
    queries Uniprot Rest API and gets proten seqs for all out-of-date uniprot IDs

    Args:
        uni2uniparc  (dict[str:str]): dictionary of uniprot id to uniparc ids
        fasta_loc (int): directory where it places all data.
    Returns:
        Writes a fasta file 
        uniparc2uni (dict[str:str]) : uniparc id 2 uniprot id  (inverse of uni2uniparc)
        olduni2seq (dict[str:str]) : out-of-date uniprot to sequnces (older versions of uniprot). 
    """
    temps=[]
    temp=[]
    uniparc2uni={}
    for i,uni in enumerate(uni2uniparc.keys()):
        if i!=0 and i%500==0:
            temps.append(temp)
            temp=[]
        temp.append(uni2uniparc[uni])
        uniparc2uni[uni2uniparc[uni]]=uni
    #print(len(temps))
    r='\r'
    olduni2seq={}
    from unipressed import UniparcClient
    for i, temp in enumerate(temps):
        reqs=UniparcClient.fetch_many(temp)
        for req in reqs:
            try:
                uniparc=req["uniParcId"]
                uni=uniparc2uni[uniparc]
                seq=req["sequence"]["value"]
                olduni2seq[uni]=seq
                #print(req["uniParcId"])
            except:
                print("boooo")
        sys.stdout.write(f"scraping uniprot {i} out of {str(len(temps))}{r}")
    n='\n'
    hf=open(f"{fasta_loc}oldunis.fasta",'w+')
    for uni in olduni2seq.keys():
        hf.write(f">{uni}{n}{olduni2seq[uni]}{n}{n}")
    hf.close()
    return uniparc2uni, olduni2seq







def map_out_uniprotIDs_to_uniParc_v2( oldunis: list[str] , csv_loc : str ="data/csv_loc/"  ):
    ############################### used unipressed to get all sequences 
    """
    round 2 of queries Uniprot Rest API and maps uniprot id to uniparc id

    Args:
        oldunis (list[str]): non-redundant of list of uniprots
        csv_loc (int): directory where it places all data.
    Returns:
        Writes a csv mapping file 
        uni2uniparc (dict[str:str]) : uniprot id 2 uniparc id
    """
    
    temps=[]
    temp=set()
    for i,uni in enumerate(oldunis):
        if i!=0 and i%500==0:
            temps.append(temp)
            temp=set()
        temp.add(uni)

    uni2uniparc={}
    c=0

    for temp in temps:

        c=c+1
        request = IdMappingClient.submit(
            source="UniProtKB_AC-ID", dest="UniParc", ids=temp
        )
        x=0
        b=1

        while b==1:
            try:
                x=x+1
                r='\r'
                sys.stdout.write(f"scraping uniprot {str(x)} seconds waiting on batch {str(c)} / {str(len(temps))} uniparc can be slow to respond {r}")
                time.sleep(1)
                u2u=list(request.each_result())
                for nmid in u2u:
                    uni=nmid["from"]
                    uniprac=nmid["to"]
                    uni2uniparc[uni]=uniprac
                b=0
            except:
                b=1
    hf=open(f"{csv_loc}uni2uniparcV2.csv","w+")
    n='\n'
    for uni in uni2uniparc.keys():
        hf.write(f"{uni},{uni2uniparc[uni]}{n}")
    hf.close()
    return uni2uniparc








def download_fasta_from_parc_v2( uni2uniparc : dict[str:str] , fasta_loc : str ="data/fasta_loc/"  ):
    ############################### used unipressed to get all sequences 
    """
    round 2 of queries Uniprot Rest API and gets proten seqs for all out-of-date uniprot IDs

    Args:
        uni2uniparc  (dict[str:str]): dictionary of uniprot id to uniparc ids
        fasta_loc (int): directory where it places all data.
    Returns:
        Writes a fasta file 
        uniparc2uni (dict[str:str]) : uniparc id 2 uniprot id  (inverse of uni2uniparc)
        olduni2seq (dict[str:str]) : out-of-date uniprot to sequnces (older versions of uniprot). 
    """
    temps=[]
    temp=[]

    uniparc2uni={}
    for i,uni in enumerate(uni2uniparc.keys()):
        if i!=0 and i%500==0:
            temps.append(temp)
            temp=[]
        temp.append(uni2uniparc[uni])
        uniparc2uni[uni2uniparc[uni]]=uni
    #print(len(temps))
    r='\r'
    olduni2seq={}
    from unipressed import UniparcClient
    for i, temp in enumerate(temps):
        reqs=UniparcClient.fetch_many(temp)
        for req in reqs:
            try:
                uniparc=req["uniParcId"]
                uni=uniparc2uni[uniparc]
                seq=req["sequence"]["value"]
                olduni2seq[uni]=seq
                #print(req["uniParcId"])
            except:
                print("boooo")
        sys.stdout.write(f"scraping uniprot {i} out of {str(len(temps))}{r}")
    n='\n'
    hf=open(f"{fasta_loc}oldunisV2.fasta",'w+')
    for uni in olduni2seq.keys():
        hf.write(f">{uni}{n}{olduni2seq[uni]}{n}{n}")
    hf.close()
    return uniparc2uni, olduni2seq

def requestUniprot(uni,q):
    url=f"https://rest.uniprot.org/unisave/{uni}?format=fasta&uniqueSequences=true"
    try:
        req = requests.get(url, verify=True, timeout=5)
        stat=req.status_code
        #print(stat)
        if 200 == stat:
            q.put(req.text)
        else:
            q.put("")
    except:
        q.put("")


def download_fasta_from_parc_v3_multiThread(missmatchedseq : list[str] , fasta_loc : str ="data/fasta_loc/"  ):
    ###################--------------REQUEST CAP IN 200 PER SECOND 
    ###################--------------DO NOT EXECEED 100 PER SECOND TO BE SAFE
    ###################--------------OR YOU WILL BE BLOCKED FROM UNIPROT REST API
    ##################---------------SEE MORE DETAILS FROM https://www.ebi.ac.uk/proteins/api/doc/
    ##################---------------CHECK BEROFE RUNING THIS CODE JUST IN CASE UNIPROT CHANGED THE CAP
    """
    round 3 of queries Uniprot Rest API and gets proten seqs for all out-of-date uniprot IDs

    Args:
        uni2uniparc  (dict[str:str]): dictionary of uniprot id to uniparc ids
        fasta_loc (int): directory where it places all data.
    Returns:
        Writes a fasta file 
        mmu_uni2seq (dict[str:str]) : uniparc id 2 uniprot id  (inverse of uni2uniparc)
    """
    mmu=set()
    for i,info in enumerate(missmatchedseq):
        l=list(info.split("_"))
        uni=l[0]
        mmu.add(uni)
    mmu_uni2seq={}
    rate=100
    lmmu=list(missmatchedseq)
    index=0
    #print(len(lmmu))
    r='\r'
    while True:
        sys.stdout.write(f"scraping uniparc {index} / {len(lmmu)}{r}")
        q = queue.Queue()
        threads=[]
        if (rate+index)>len(lmmu):
            #creat thread
            for i in range(len(lmmu)-index):
                threads.append(threading.Thread(target=requestUniprot, args=(lmmu[index+i],q)))
            #run thread
            for i in range(len(lmmu)-index):
                threads[i].start()
                time.sleep(float(1/rate))
            #join thread
            for i in range(len(lmmu)-index):
                threads[i].join()


            ###read
            for i in range(len(lmmu)-index):
                temp=q.get()
                lines=list(temp.split("\n"))
                seq=""
                for line in lines:
                    if line=="\n" or line=="":
                        continue
                    if line[0]==">":
                        if seq!="":
                            mmu_uni2seq[uni]=seq
                            seq=""
                        uni=line
                    else:
                        seq=seq+line
                

                    
                
            break
        else:
            
            #creat thread
            for i in range(rate):
                threads.append(threading.Thread(target=requestUniprot, args=(lmmu[index+i],q)))
            #run thread
            for i in range(rate):
                threads[i].start()
                time.sleep(float(1/rate))
            #join thread
            for i in range(rate):
                threads[i].join()


            ###read
            for i in range(rate):
                temp=q.get()
                lines=list(temp.split("\n"))
                seq=""
                for line in lines:
                    if line=="\n" or line=="":
                        continue
                    if line[0]==">":
                        if seq!="":
                            mmu_uni2seq[uni]=seq
                            seq=""
                        uni=line
                    else:
                        seq=seq+line
            index=index+rate
    mmUni2ver2seq={}
    hf=open(f"{fasta_loc}mmUniFull.fasta",'w+')
    n='\n'
    for uni in mmu_uni2seq.keys():
        hf.write(f"{uni}{n}{mmu_uni2seq[uni]}{n}{n}")
        mmUni2ver2seq
    hf.close()
