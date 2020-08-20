# GWAS on BGEN files and  integrating measurements across tissues 

This file introduce the method of running GWAS on UK Biobank BGEN files, harmonization and imputation of GWAS summary result and integrating measurements across tissues.

## Prerequisites

Following tools and their prerequisteds are needed: 

* [Qctool v2.0.7](https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/download.html) (v2.0.8 will be better, but this version can not be installed on some machine)
* [BGENIE v1.2](https://jmarchini.org/bgenie/)
* [summary-gwas-imputation](https://github.com/hakyimlab/summary-gwas-imputation)
* [MetaXcan](https://github.com/hakyimlab/MetaXcan)

Python 3 are recommanded since some packages in the **summary gwas imputation** only support python 3. 

## filtering

### filtering subjects 

The original BGEN files include all of the 500k individuals and were seperated by chromosomes. An `ID.csv` file is needed to filter the subjects, by using **Qctool**

```bash
$ qctool -g example.bgen -s example.sample -og filtered.bgen incl-samples ID.csv
```

Ususally, the BGEN files released by UK Biobank exclude subject IDs. Therefore, sample files (which use as a link) are required.

This process can take a quite long time. 

### filtering SPNs

We can futher filter out the SPNs whose have very small minor allele frequency (for example less than 1%). This  process is optional but it can speed up the GWAS and save some space.

There are two steps: 

1. Generate a statistical report by running:

```bash
$ qctool -g example.bgen -snp-stats -osnp snp-stats.txt
```

This can take a very long time. The output file `snp-stats.txt` is pretty large and large memory is needed.

2. Create `rsid.txt` file which include the ID of SPNs.

```py
import pandas as pd
def readSnpStats(snpStatsPath,rsidPath):
    snpStatsPath = snpStatsPath
    ## remove first few lines
    ## read
    snpStats = pd.read_table(snpStatsPath,sep=" ")
    # filter minor allele frequency > 0.01
    snpStats = snpStats[snpStats['minor_allele_frequency']>0.01]
    print(snpStats.columns)
    #snpStats = snpStats[ ['alternate_ids', 'rsid','chromosome', 'position', 'minor_allele',
       #'major_allele']]
    snpStats = snpStats['rsid']
    #SNPID, rsid, chromosome, position,
    #snpStats.rename(columns={'alternate_ids': 'SNPID'}, index={'ONE': 'one'}, inplace=True)
    snpStats.to_csv(rsidPath,sep=" ",index = False)
    return snpStats
```

* **snpStates** is the path to `snp-stats.txt` which generated from the previous step.
* **rsidPath** is the path to output file `rsid.txt`

3. filtering SNPs

Only include those SPNs in `rsid.txt`

```bash
$ qctool -g example.bgen -og subsetted.bgen -incl-rsids rsid.txt -s example.sample
```

**-s example.sample** is not neccesary when using **Qctool v2.0.8**, since it automatically include the subject IDs in the output files.

The sample file can be created from the BGEN file `filtered.bgen` by running:

 ``` bash
qctool -g example.bgen -os example.sample
 ```

If a sample file is not specified in **Qctool v2.0.7**, the subject IDs will be excluded from the BGEN files. However, the order of the subject won't change so we can used the sample file generated from `filtered.bgen`

**This can take a very long time**

## GWAS

Before running GWAS, the genotypes files should be reordered accourding to the sample files generated from the `filtered.bgen`

This can be done by:

``` python
import pandas as pd
def getReorderPheno(phenotypes, sample_file,reorder_phenotypes):
    if os.path.exists(sample_file):
        sample = pd.read_csv(sample_file,sep = " ")
        print("reading sample file:")
        if os.path.exists(phenotypes):
            pheno = pd.read_csv(phenotypes,sep=" ")
            print(sample.head())
            reorder_pheno = pd.merge(left=sample,right=pheno,how="left",left_on="ID_1",right_on="ID")
            print("reorder pheno file: ")
            print(reorder_pheno.head())
            reorder_pheno = reorder_pheno.drop(columns=["ID_1","ID_2","missing","ID"])
            reorder_pheno = reorder_pheno.drop([0])        
            reorder_pheno.to_csv(reorder_phenotypes,index = False,sep=" ",na_rep="-999")
            return reorder_pheno
        else:
            print("phenotype file {} does not exit").format(phenotypes)
    else:
        print("sample file {} does not exit").format(sample_file) 
```

* **phenotypes** is the path to phenoyptes files
* **sample_file** is the path to samples files created from the `filtered.bgen` files
* **reorder_phenotypes** is the path to output file `reorder_phenotypes.pheno`

### running GWAS

Run **BGENIE**

``` bash
bgenie --bgen example.bgen --pheno reorder_phenotypes.pheno --out example.out
```

This process is time-comsuming

## harmonization and imputation

