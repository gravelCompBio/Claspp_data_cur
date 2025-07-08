# Claspp_data_cur
Data curation for the Claspp Model 


## Introduction   


<!-- This github was Made by Nathan Gravel --> 


  


  

   


  

   

For this repo we are generating a data set in the hopes to plug in for both Contrastive learning and fine tuning.  


<p align="center">
  <img width="120%" src= "https://github.com/user-attachments/assets/e2a5f335-3710-489b-800a-df5cab462488">
</p>


This repository contains the code to replicate the CLASPP model's Data Curation. This is going above and beyond what other Deep learning based projects has done and I hope this raises the standard for what is expected in the field. With the orignal code I did NOT maintain/recover any random seeds so the data will be differnt from what I used but the process should be the same. I also rely on [Ultra-Scalable Spectral Clustering](https://github.com/huangdonghere/USPEC_USENC) and it was writen in matlab. Finally for the sk-learn version of Spectral Culstering you need a minimun of 250 GB of ram to run the K-Ubiq data set (n of 80,000). The Sequence idenity clustering is GPU accelorated (using [torch](https://pytorch.org/)) so a minumum of 8 GB of vRAM is needed. This data curation should take 12 hours to run in total. 



</br> 




For the Claspp model go to this Github. For the webtool for the CLASPP can be accessed from: 

</br> 
   

## Quick overview of the dependencies 

![Python](https://img.shields.io/badge/Python-FFD43B?style=for-the-badge&logo=python&logoColor=blue)
![Anaconda](https://img.shields.io/badge/Anaconda-%2344A833.svg?style=for-the-badge&logo=anaconda&logoColor=white)
![Jupyter](https://img.shields.io/badge/Jupyter-F37626.svg?&style=for-the-badge&logo=Jupyter&logoColor=white)
![PyTorch](https://img.shields.io/badge/PyTorch-EE4C2C?style=for-the-badge&logo=pytorch&logoColor=white)

  

![Numpy](https://img.shields.io/badge/Numpy-777BB4?style=for-the-badge&logo=numpy&logoColor=white) 
![Pandas](https://img.shields.io/badge/Pandas-2C2D72?style=for-the-badge&logo=pandas&logoColor=white) 
![scikit-learn](https://img.shields.io/badge/scikit--learn-%23F7931E.svg?style=for-the-badge&logo=scikit-learn&logoColor=white) 

  

   

</br> 


  

## Installing dependencies with version info    

  

  

### From conda:    

  

![python=3.9.23](https://img.shields.io/badge/Python-3.9.23-green)  

  

![jupyterlab=4.0.0](https://img.shields.io/badge/jupyterlab-4.0.0-blue)  

  

Python == 3.9.16  

  

   

  

### From pip:  

  

   

  

  

![numpy=2.0.2](https://img.shields.io/badge/numpy-2.0.2-blue)  


![scikit-learn=1.6.1](https://img.shields.io/badge/scikitlearn-1.6.1-blue)    


![requests=2.32.4](https://img.shields.io/badge/requests-2.32.4-blue)  


![wget=3.2](https://img.shields.io/badge/wget-3.2-blue)  


![unipressed=1.4.0](https://img.shields.io/badge/unipressed-1.4.0-blue)  


![scipy=1.13.1](https://img.shields.io/badge/scipy-1.13.1-blue) 


![torch=2.0.1](https://img.shields.io/badge/torch-2.0.1-blue)      

  

### For torch/PyTorch 

  

Make sure you go to this website https://pytorch.org/get-started/locally/ 

  

Follow along with its recommendation  

  

Installing torch can be the most complex part  

  
  
  

  

</br> 

  

</br> 

  

   

  

## Downloading this repository   

  

```   
gh repo clone gravelCompBio/Phosformer-ST 
```   

  

```   
cd Phosformer-ST 
``` 

### The following step demonstrates users how to download the training weights 


  -other repositories were used because the folder's memory size is larger than the allowed space on github 

  

  

</br> 

  

### Main option) Hugging Face  

  

Then download the link found in `multitask_MHA_esm2_t30_150M_UR50D_neg_ratio_8+8_shift_30_mask_0.2_2023-03-25_90.txt` or can be found at this link https://huggingface.co/gravelcompbio/Phosformer-ST_trainging_weights/tree/main 

  

The download link should take to a page that should look like this 

  

  

![Screenshot from 2023-07-24 13-49-54](https://github.com/gravelCompBio/Phosformer-ST/assets/75225868/bd2ebb5e-6174-4695-9cd3-730b835a8664) 

  

  

  

Click the download box highlighted in picture above 

  

  

  

</br> 

  

### Alternative option) Zenodo  

  

  

Then download the link found in `multitask_MHA_esm2_t30_150M_UR50D_neg_ratio_8+8_shift_30_mask_0.2_2023-03-25_90.txt` or can be found at this link https://zenodo.org/record/8170005 

  

The download link should take to a page that should look like this 

  

  

![Screenshot from 2023-07-20 18-14-19](https://github.com/gravelCompBio/Phos-ST-temp/assets/75225868/109c898e-49cc-4849-abb6-1dcb1f3aa5c1) 

  

Click the download box highlighted in picture above 

  

</br> 

  

### After picking one of the options above to download the training weights see below 

  

  

Once downloaded, **unizip** the folder and place in the `Phosformer-ST` along with all the other files in this github repository 

  

The final `Phosformer-ST` directory orinization should have the following files/folder  

  

- file 1 `phos-ST_Example_Code.ipynb` 

  

- file 2 `modeling_esm.py` 

   

- file 3 `configuration_esm.py` 

  

- file 4 `tokenization_esm.py` 

  

- file 5 `multitask_MHA_esm2_t30_150M_UR50D_neg_ratio_8+8_shift_30_mask_0.2_2023-03-25_90.txt` 

  

- file 6 `phosST.yml` 

   

- file 7 `Readme.md`



- file 8 `LICENSE`

  

- folder 1 `multitask_MHA_esm2_t30_150M_UR50D_neg_ratio_8+8_shift_30_mask_0.2_2023-03-25_90` (make sure it is unzipped) 

  

:tada: Once you have a folder with the files/folder above you have all the required files to run the model :tada: 

  

  

</br> 

  

</br> 

  


  

## Included in this repository are the following:   

  

   
- `scrape_and_clean_data_1.py`: intial python file that downloads the PTM data from [dbPTM](https://biomics.lab.nycu.edu.tw/dbPTM/) and [uniprot](https://www.uniprot.org/). Then cleans the data 
  
- `rep_base_gready_cluster_2.py`: preforms sequence idenity thresh-holding on ptm-sites

- `post_rbgc_prep_for_spec_clustering_3.py`: orginizes data and keep tracks of sampling, multi-label events, and negtives data points

- `spec_clus_sklearn_4.py`: Cluster all  PTM types using spectral clustering (sk-Learn)

- `get_easy_neg_data_5.py`: preforms sequence idenity thresh-holding on ptm-sites and potential easy negtive sites

- `sample_spec_cluster_6.py`: sample from spectral clusters and strifies samples across all Spectral Clusters


- `data/`: directory to hold the data/temp files   

    - `csv_loc/`: where all csv and tsv files are located 

      - `uni2acc.tsv`: hand currated exeptions to data found in dbPTM that where difficult to automate (only a handfull)
    
    - `text_loc/`: where all txt files are located
 
      - `infosforaffinitymat/` : where the order of the indexes for the affinity matrix npz files are held.

        - `ST-Phosphorylation.txt` : info for the  .mat file that is already populated
     
      - `ussc_output` : where the [Ultra-Scalable Spectral Clustering](https://github.com/huangdonghere/USPEC_USENC) clustering file goes (example format should replace on your own)
     
        - `testdataNC-5.txt` : ouput of the  [Ultra-Scalable Spectral Clustering](https://github.com/huangdonghere/USPEC_USENC) (example of what the format should be plead replace)

- `getData/`: Where all funtions exist
  
  - `DataProcess.py/`: where 90% of the code resides in the form of helper fuctions
 
  - `FileReader.py/` : where all the file reader funtions exist
 
  - `WebScraper.py` : where the code for interacting with APIs and pulling data
 
     
  
    
      
    

  

   

- `README.md`: 

  

- `LICENSE`: MIT licence (open source but don't try to make money off of it)
  

  

    

  

</br> 

  

</br> 
 

  

## ![Anaconda](https://img.shields.io/badge/Anaconda-%2344A833.svg?style=for-the-badge&logo=anaconda&logoColor=white) Installing dependencies with conda  

  

### PICK ONE of the options below  

### Main Option) Utilizing the PhosformerST.yml file 

here is a step-by-step guide to set up the environment with the yml file  

  

Just type these lines of code into the terminal after you download this repository (this assumes you have anaconda already installed) 

  

```   
conda env create -f phosST.yml -n PhosST  
```   

```   
conda deactivate 
```   

```   
conda activate phosST  
```   

  

### Alternative option) Creating this environment without yml file 

(This is if torch is not working with your version of cuda or any other problem) 

Just type these lines of code into the terminal after you download this repository (this assumes you have anaconda already installed) 

```   
conda create -n phosST python=3.9  
``` 

```   
conda deactivate 
``` 

```   
conda activate phosST  
``` 

```   
conda install -c conda-forge jupyterlab 
``` 

```   
pip3 install numpy==1.24.3 
``` 

```   
pip3 install pandas==2.0.2 
``` 

```   
pip3 install matplotlib==3.7.1 
``` 

```   
pip3 install scikit-learn==1.2.2 
``` 

```   
pip3 install tqdm==4.65.0 
``` 

```   
pip3 install fair-esm==2.0.0 
``` 

```   
pip3 install transformers==4.31.0 
``` 

### **For torch you will have to download to the torch's specification if you want gpu acceleration from this website** https://pytorch.org/get-started/locally/ 

  

```   
pip3 install torch torchvision torchaudio 
``` 

  

### the terminal line above might look different for you  

  

We provided code to test Phosformer-ST (see section below) 

  

  

</br> 

  

</br> 

  

  

  

## Utilizing the Model with our example code 

All the following code examples is done inside of the `phos-ST_Example_Code.ipynb` file using jupyter lab 

  

Once you have your environment resolved just use jupyter lab to access the example code by typing the command below in your terminal (when you're in the `Phosformer-ST` folder)  

```   

jupyter lab 

``` 

Once you open the notebook on your browser, run each cell in the notebook  

  

</br> 

  

### Testing Phosformer-ST with the example code 

There should be a positive control and a negative control example code at the bottom of the `phos-ST_Example_Code.ipynb` file which can be used to test the model. 
  

**Positive Example** 

```Python 

# P17612 KAPCA_HUMAN 

kinDomain="FERIKTLGTGSFGRVMLVKHKETGNHYAMKILDKQKVVKLKQIEHTLNEKRILQAVNFPFLVKLEFSFKDNSNLYMVMEYVPGGEMFSHLRRIGRFSEPHARFYAAQIVLTFEYLHSLDLIYRDLKPENLLIDQQGYIQVTDFGFAKRVKGRTWTLCGTPEYLAPEIILSKGYNKAVDWWALGVLIYEMAAGYPPFFADQPIQIYEKIVSGKVRFPSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWF" 

# P53602_S96_LARKRRNSRDGDPLP 

substrate="LARKRRNSRDGDPLP" 

  

phosST(kinDomain,substrate).to_csv('PostiveExample.csv') 

``` 

  

  

**Negative Example** 

```Python 

# P17612 KAPCA_HUMAN 

kinDomain="FERIKTLGTGSFGRVMLVKHKETGNHYAMKILDKQKVVKLKQIEHTLNEKRILQAVNFPFLVKLEFSFKDNSNLYMVMEYVPGGEMFSHLRRIGRFSEPHARFYAAQIVLTFEYLHSLDLIYRDLKPENLLIDQQGYIQVTDFGFAKRVKGRTWTLCGTPEYLAPEIILSKGYNKAVDWWALGVLIYEMAAGYPPFFADQPIQIYEKIVSGKVRFPSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWF" 

# Q01831_T169_PVEIEIETPEQAKTR 

substrate="PVEIEIETPEQAKTR" 

  

phosST(kinDomain,substrate).to_csv('NegitiveExample.csv') 

``` 

Both scores should show up in a csv file in the current directory

  

</br> 

  

### Inputting your own data for novel predictions 

One can simply take the code from above and modify the string variables `kinDomain` and `substrate` to make predictions on any given kinase substrate pairs 

  

**Formatting of the `kinDomain` and `substrate` for input for Phosformer-ST are as follows:** 

  

  - `kinDomain` should be a human Serine/Threonine kinase domain (not the full sequence).
     

  - `substrate` should be a 15mer with the center residue/char being the target Serine or Threonine being phosphorylated 

  

Not following these rules may result in dubious predictions  

  

  

</br> 

  

### How to interpret Phosformer-ST's output 

This model outputs a prediction score between 1 and 0.


We trained the model to uses a cutoff of 0.5 to distinguish positive and negative predictions 


A score of 0.5 or above indicates a positive prediction for peptide substrate phosphorylation by the given kinase

  
  

</br> 

  
  

## Troubleshooting 

  

If torch is not installing correctly or you do not have a GPU to run Phosformer-ST on, the CPU version of torch is perfectly fine to use 

  

Using the CPU version of torch might increase your run time so for large prediction datasets GPU acceleration is suggested 

  

If you just are here to test if it Phosformer-ST works, the example code should not take too much time to run on the CPU version of torch   

  

Also depending on your GPU the `batch_size` argument might need to be adjusted 


#### 2024-05-17
- if you get an 'EsmTokenizer' object has no attribute 'all_tokens' error when loading the tokenizer
- - Make sure you have version of  transformers==4.31.0 installed



### The model has been tested on the following computers with the following specifications for troubleshooting proposes 

  

</br> 

  

**Computer 1** 



NVIDIA Quadro RTX 5000 (16 GB vRAM)(CUDA Version: 12.1)  

  

Ubuntu 22.04.2 LTS 

  

Intel(R) Xeon(R) Bronze 3204 CPU @ 1.90GHz  (6 cores) x (1 thread per core) 

  

64 GB ram 



  

</br> 

  

**Computer 2** 



NVIDIA RTX A4000 (16 GB vRAM)(CUDA Version: 12.2)  

  

Ubuntu 20.04.6 LTS 

  

Intel(R) Xeon(R) Bronze 3204 CPU @ 1.90GHz  (6 cores) x (1 thread per core) 

  

64 GB ram 

  





</br> 


## Other accessory tools and resources
A webtool for Phosformer-ST can be accessed from: https://esbg.bmb.uga.edu/phosformer/ . A huggingface repository can be downloaded from: https://huggingface.co/gravelcompbio/Phosformer-ST_with_training_weights. A huggingface spaces app is available at: https://huggingface.co/spaces/gravelcompbio/Phosformer-ST


  

 

 
