import os
import sys


def read_downloaded_fasta_from_unis_v1(fasta_loc : str ="data/fasta_loc/"):
    """
    reads fasta files scraped from uniprot 
    Args:
        fasta_loc (int): directory where it places all data.
    Returns:
        uni2seq (dict[str:str]) : uniprot id 2 sequences 
    """
    uni2seq={}
    hf=open(f"{fasta_loc}uniprotPostMissingPepPreConRes.fasta",'r')
    lines=hf.readlines()
    for line in lines:
        if line[0]=='>':
            uni=line[1:-1]
        elif line=='\n':
            uni2seq[uni]=seq
        else:
            seq=line[:-1]
    return uni2seq



def read_uniprotIDs_to_uniParc_v1(csv_loc : str ="data/csv_loc/" ):
    """
    reads the queries Uniprot Rest API and maps uniprot id to uniparc id 
    Args:
        csv_loc (int): directory where it places all data.
    Returns:
        uni2uniparc (dict[str:str]) : uniprot id 2 uniparc id
    """
    uni2uniparc={}
    hf=open(f"{csv_loc}uni2uniparc.csv","r")
    lines=hf.readlines()
    for line in lines:
        l=list(line[:-1].split(','))
        uni=l[0]
        uniparc=l[1]
        uni2uniparc[uni]=uniparc
    return uni2uniparc




def read_fasta_from_parc_v1(uni2uniparc : dict[str:str], fasta_loc : str ="data/fasta_loc/"):
    """
    round 2 reads fasta files scraped from uniprot (out-of-date uniprot) 
    Args:
        uni2uniparc  (dict[str:str]): dictionary of uniprot id to uniparc ids
        fasta_loc (int): directory where it places all data.
    Returns:
        uniparc2uni (dict[str:str]) : uniparc id 2 uniprot id  (inverse of uni2uniparc)
        olduni2seq (dict[str:str]) : uniprot id 2 sequences 
    """
    uniparc2uni={}
    for uni in uni2uniparc.keys():
        uniparc2uni[uni2uniparc[uni]]=uni
    olduni2seq={}
    hf=open(f"{fasta_loc}oldunis.fasta","r")
    lines=hf.readlines()
    for line in lines:
        if line[0]=='>':
            uni=line[1:-1]
        elif line[0] =='\n':
            olduni2seq[uni]=seq
        else:
            seq=line[:-1]
    return uniparc2uni, olduni2seq




def read_uniprotIDs_to_uniParc_v2(csv_loc : str ="data/csv_loc/" ):
    """
    reads the queries Uniprot Rest API and maps uniprot id to uniparc id 
    Args:
        csv_loc (int): directory where it places all data.
    Returns:
        uni2uniparc (dict[str:str]) : uniprot id 2 uniparc id
    """
    uni2uniparc={}
    hf=open(f"{csv_loc}uni2uniparcV2.csv","r")
    lines=hf.readlines()
    for line in lines:
        l=list(line[:-1].split(','))
        uni=l[0]
        uniparc=l[1]
        uni2uniparc[uni]=uniparc
    return uni2uniparc



def read_downloaded_fasta_from_unis_v2(uni2uniparc : dict[str:str],fasta_loc : str ="data/fasta_loc/"):
    """
    round 3 reads fasta files scraped from uniprot (out-of-date uniprot second iter) 
    Args:
        uni2uniparc  (dict[str:str]): dictionary of uniprot id to uniparc ids
        fasta_loc (int): directory where it places all data.
    Returns:
        uniparc2uni (dict[str:str]) : uniparc id 2 uniprot id  (inverse of uni2uniparc)
        olduni2seq (dict[str:str]) : out-of-date uniprot id 2 sequences
    """
    uniparc2uni={}
    for uni in uni2uniparc.keys():
        uniparc2uni[uni2uniparc[uni]]=uni
    olduni2seq={}
    hf=open(f"{fasta_loc}oldunisV2.fasta","r")
    lines=hf.readlines()
    for line in lines:
        if line[0]=='>':
            uni=line[1:-1]
        elif line[0] =='\n':
            olduni2seq[uni]=seq
        else:
            seq=line[:-1]
    return uniparc2uni, olduni2seq





def read_downloaded_fasta_from_unis_v3(fasta_loc : str ="data/fasta_loc/"):
    """
    round 3 reads fasta files scraped from uniprot (out-of-date uniprot second iter) 
    Args:
        fasta_loc (int): directory where it places all data.
    Returns:
        mmUni2ver2seq (dict[str:str]) : uniprot id 2 sequences 
    """
    mmUni2ver2seq={}
    hf=open(f"{fasta_loc}mmUniFull.fasta",'r')
    lines=hf.readlines() 
    for line in lines:
        if line[0]=='>':
            l=list(line.split(" "))
            #print(l)
            uni=l[0][1:-1]
            ver=l[1][3:]
            #print(uni,ver)
            
        elif line=='\n':
            if uni not in mmUni2ver2seq.keys():
                mmUni2ver2seq[uni]={}
            mmUni2ver2seq[uni][ver]=seq
            
        else:
            seq=line[:-1]
    return mmUni2ver2seq



