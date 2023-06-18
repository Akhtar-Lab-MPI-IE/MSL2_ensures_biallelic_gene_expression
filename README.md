# Allele-specific-project

This repository includes workflows and scripts used for the paper:

MSL2 ensures biallelic gene expression in mammals

Accompanying sequencing data is available in [GEO: GSE183764](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183764).

## snakePipes

Nearly all data is processed with a custom version of snakePipes.  
The dependencies can be installed using [conda](https://docs.conda.io/en/latest/miniconda.html).  
Make sure to have [mamba](https://mamba.readthedocs.io/en/latest/) installed as well (in your base environment).  

 > conda env create -f confs/snakepipes.yml  

And the snakepipes version itself can be installed afterwards:

 > conda activate snakepipesAS  
 > pip install git+https://github.com/Akhtar-Lab-MPI-IE/snakepipes_allele_specific@allele_specific  

The source code for the snakePipes version is available [here](https://github.com/Akhtar-Lab-MPI-IE/snakepipes_allele_specific).

After creating the conda environment, the environments, specific reference and cluster config need to be set up:
 
 > conda activate snakepipesAS  

Create the snakePipes environments (optionally in a specific conda directory):

 > snakePipes createEnvs --condaDir /path/to/condaenv/directory --only CONDA_SHARED_ENV CONDA_CREATE_INDEX_ENV CONDA_RNASEQ_ENV CONDA_DNA_MAPPING_ENV CONDA_CHIPSEQ_ENV CONDA_ATAC_ENV CONDA_RMD_ENV CONDA_PREPROCESSING_ENV CONDA_SAMBAMBA_ENV CONDA_pysam_ENV  

Modify the cluster.yaml file appropriately (snakemake_cluster_cmd especially). The path to this file can be found with this command

 > snakePipes info


