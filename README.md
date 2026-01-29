

<!-- This github was Made by Nathan Gravel --> 
# Claspp_data_cur 


## Introduction   




  


  

   


  

   



<p align="center">
  <img width="120%" src= "https://github.com/gravelCompBio/Claspp_data_cur/blob/main/github_fig/Screenshot%20from%202025-08-05%2011-50-53.png">
</p>


This repository contains the code to replicate the CLASPP model's Data Curation. With the original code I did NOT maintain/recover any random seeds so the data will be different from what I used but the process should be the same. The Sequence identity clustering is GPU accelerated (using [torch](https://pytorch.org/)) so a minimum of 8 GB of vRAM is needed. For the sk-learn version of Spectral Clustering you need a minimum of 250 GB of ram to run the K-Ubiq data set (n of 80,000). I also rely on [Ultra-Scalable Spectral Clustering](https://github.com/huangdonghere/USPEC_USENC) and it was written in matlab for larger PTM datasets specifically ST-Phosphorylation. This data curation should take 12 hours to run in total. 


the manuscript can be found here (TBA)

</br> 






For the Claspp model go to this Github. For the webtool for the CLASPP can be accessed from: [github version of model(local)](https://github.com/gravelCompBio/Claspp_forward)  or [webtool version](https://github.com/gravelCompBio/Claspp_data_cur/tree/main)

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

# Downloading this repository   

  

```   
git clone https://github.com/gravelCompBio/Claspp_data_cur.git
```   

  

```   
cd Claspp_data_cur
``` 

  

# Installing dependencies with version info    

  

  

### From conda:    

  

![python=3.9.23](https://img.shields.io/badge/Python-3.9.23-green) 



  

   

  

### From pip:  

  

   

  

  

![numpy=2.0.2](https://img.shields.io/badge/numpy-2.0.2-blue)  ![scikit-learn=1.6.1](https://img.shields.io/badge/scikitlearn-1.6.1-blue) ![scipy=1.13.1](https://img.shields.io/badge/scipy-1.13.1-blue)  ![requests=2.32.4](https://img.shields.io/badge/requests-2.32.4-blue)  
![wget=3.2](https://img.shields.io/badge/wget-3.2-blue)  ![unipressed=1.4.0](https://img.shields.io/badge/unipressed-1.4.0-blue)  ![pandas==2.3.1](https://img.shields.io/badge/pandas-2.3.1-blue)  ![torch=2.7.1](https://img.shields.io/badge/torch-2.7.1-blue)      

  

### For torch/PyTorch 

  

Make sure you go to this website [pytorch](https://pytorch.org/get-started/locally/) 

  

Follow along with its recommendation  

  

Installing torch can be the most complex part  

  
  
  

  

</br> 

  

</br> 

  

   

  




  

## ![Anaconda](https://img.shields.io/badge/Anaconda-%2344A833.svg?style=for-the-badge&logo=anaconda&logoColor=white) Installing dependencies with conda  

    

### Creating this conda environment 
(yml file is included but torch sometimes makes it not usable depending on your nvidia driver)

Just type these lines of code into the terminal after you download this repository (this assumes you have [anaconda](https://www.anaconda.com/) already installed) 

```   
conda create -n claspp_cur python=3.9.23 
``` 

```   
conda deactivate 
``` 

```   
conda activate claspp_cur  
``` 

```   
pip3 install numpy==2.0.2
```

```   
pip3 install scikit-learn==1.6.1
```

```   
pip3 install requests==2.32.4
```

```   
pip3 install wget==3.2
``` 

```   
pip3 install unipressed==1.4.0
``` 

```   
pip3 install pandas==2.3.1
``` 

### **For torch you will have to download to the torch's specification if you want gpu acceleration from this website ** [pytorch](https://pytorch.org/get-started/locally/) ** 

  

```   
pip3 install torch torchvision torchaudio 
``` 

  

### the terminal line above might look different for you  

  


  
:tada: you are now ready to run the code :tada: 
  

</br> 

  

</br> 

  

### The following step demonstrates users how to find the different versions of the data curation 


  -other repositories were used because the folder's memory size is larger than the allowed space on github 

  

  

</br> 

  

### Other related repos



| Repo  | Link | Discription|
| ------------- | ------------- |------------------------------------------|
| GitHub  | [github version Data_cur](https://github.com/gravelCompBio/Claspp_data_cur)  | This verstion contains code but but no data. It needs you to run the code to generate all the helper-files (will take some time run this code)|
| GitHub  | [github version Forward](https://github.com/gravelCompBio/Claspp_forward)  | This verstion contains code but NOT any weights (file too big for github)|
| Huggingface | [huggingface version Forward](https://huggingface.co/esbglab/Claspp_forward)  | This verstion contains code and training weights |
| Zenodo | [zenodo version training_data]([https://zenodo.org/records/16739128?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjY4MDBkYjUwLTVlZjktNDZlNy04MjZjLTgzZjA0NjZiYmZlYyIsImRhdGEiOnt9LCJyYW5kb20iOiJlMThhZGNlMWUxN2EzNjYxNzllYjg5MWRiZjhiMWYxNSJ9.7Os5ZzQLT3TJu3Clv1Sxvh8oVtFTxxoeYLACFgKwZRjCApQdfQO2-AvctQ-eIIEojKTBGHLCcHlMTDG38AKn8A](https://zenodo.org/records/17674057))  | zenodo version of training/testing/validation data|

  

</br> 


  


### Final data organization

<p align="center">
  <img width="120%" src= "https://github.com/gravelCompBio/Claspp_data_cur/blob/main/github_fig/Screenshot%20from%202025-08-05%2011-48-48.png">
</p>

</br> 

### Workflow
<p align="center">
  <img width="120%" src= "https://github.com/gravelCompBio/Claspp_data_cur/blob/main/github_fig/Screenshot%20from%202025-08-05%2011-48-06.png">
</p>

Johnson, J. L., Yaron, T. M., Huntsman, E. M., Kerelsky, A., Song, J., Regev, A., ... & Cantley, L. C. (2023). An atlas of substrate specificities for the human serine/threonine kinome. Nature, 613(7945), 759-766.

UniProt Consortium. (2019). UniProt: a worldwide hub of protein knowledge. Nucleic acids research, 47(D1), D506-D515.
</br>

## Included in this repository are the following:   
   
- `scrape_and_clean_data_1.py`: initial python file that downloads the PTM data from [dbPTM](https://biomics.lab.nycu.edu.tw/dbPTM/) and [uniprot](https://www.uniprot.org/). Then cleans the data 
  
- `rep_base_gready_cluster_2.py`: performs sequence identity thresh-holding on ptm-sites

- `post_rbgc_prep_for_spec_clustering_3.py`: organizes data and keep tracks of sampling, multi-label events, and negatives data points

- `spec_clus_sklearn_4.py`: Cluster all  PTM types using spectral clustering (sk-Learn)

- `get_easy_neg_data_5.py`: performs sequence identity thresh-holding on ptm-sites and potential easy negative sites

- `sample_spec_cluster_6.py`: sample from spectral clusters and stratifies samples across all Spectral Clusters


- `data/`: directory to hold the data/temp files   

    - `csv_loc/`: where all csv and tsv files are located 

      - `uni2acc.tsv`: hand currated exceptions to data found in dbPTM that where difficult to automate (only a handful)
    
    - `text_loc/`: where all txt files are located
 
      - `infosforaffinitymat/` : where the order of the indexes for the affinity matrix npz files are held.

        - `ST-Phosphorylation.txt` : info for the  .mat file that is already populated
     
      - `ussc_output` : where the [Ultra-Scalable Spectral Clustering](https://github.com/huangdonghere/USPEC_USENC) clustering file goes (example format should replace on your own)
     
        - `testdataNC-5.txt` : ouput of the  [Ultra-Scalable Spectral Clustering](https://github.com/huangdonghere/USPEC_USENC) (example of what the format should look like)

- `getData/`: Where all functions exist
  
  - `DataProcess.py/`: where 90% of the code resides in the form of helper functions
 
  - `FileReader.py/` : where all the file reader functions exist
 
  - `WebScraper.py` : where the code for interacting with APIs and pulling data
 
     
  
    
      
    

  

   

- `README.md`: 

  

- `LICENSE`: MIT licence (open source)
  

  

    

  





## Uses



### Direct Use


```
python scrape_and_clean_data_1.py
```


```
python rep_base_gready_cluster_2.py
```


```
python post_rbgc_prep_for_spec_clustering_3.py
```


```
python spec_clus_sklearn_4.py
```

Run [Ultra-Scalable Spectral Clustering](https://github.com/huangdonghere/USPEC_USENC) for any ptm type you can not fit into RAM space while doing sk-learn's Spectral Clustering

```
python get_easy_neg_data_5.py
```


```
python sample_spec_cluster_6.py
```




- **Developed by:** Curation was created by Nathan Gravel.
- Codebase testing and verstion controle by Austin Downes,
- **License:** [MIT]




 

 
