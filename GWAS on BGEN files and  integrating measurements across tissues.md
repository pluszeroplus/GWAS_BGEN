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

### harmonization

The support file for harmonization and imputation `hg19ToHg38.over.chain.gz` ,`variant_metadata.txt.gz` `chrN.variants.parquet`can be found [here](https://zenodo.org/record/3657902/#.XzQWQRNKgtM)

or by runing:

```bash
wget  https://zenodo.org/record/3657902/files/sample_data.tar?download=1
```



The harmonization can only run on one phenotype each time. Shell code to run harmonization:

``` bash
chromosome=$1

for pheno in "RAESV" "RVSV" "RVEF" "RAEF" "LASV" "RVESV" "LAEDV" "RAEDV" "RVEDV" "LAESV" "LVM" "LVSV" "LVEDV" "LVESV" "RASV" "LVEF" "LAEF"
do
python gwas_parsing.py \
-gwas_file gwas_"$chromosome".gz \
-output_column_map rsid variant_id \
-output_column_map a_1  effect_allele \
-output_column_map a_0  non_effect_allele \
-output_column_map pos position \
-output_column_map chr chromosome \
--chromosome_format \
-output_column_map af frequency \
-output_column_map "$pheno"_p pvalue  \
-output_column_map "$pheno"_beta effect_size  \
-output_column_map "$pheno"_se standard_error  \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error  \
-snp_reference_metadata /reference_panel_1000G/variant_metadata.txt.gz METADATA  
-liftover /liftover/hg19ToHg38.over.chain.gz 
-output \Har_chr"$chromosome"_"$pheno".txt.gz\
echo "finish "$pheno" on "$chromosome
done

```

**-output_column_map** is used to specify the column name

```bash -output_column_map column_name_in_the_file name_used_in_code(don't change this)```

**chromosome_format** used it when chromosome is numerical: 1,2, rather than chr1, ch2



**important**

In this file:

`summary-gwas-imputation/blob/master/src/genomic_tools_lib/file_formats/Parquet.py`

In line 198, this function is used to read some files but failed.

``` pytho
def _read(file, columns=None, skip_individuals=False, to_pandas=False, specific_individuals=None):
    if columns is None:
        columns = file.schema.names
    if not skip_individuals:
        columns = ["individual"]+[x for x in columns]
    if skip_individuals and specific_individuals is not None:
            raise RuntimeError("Unsupported combination")
    v = file.read(columns=columns)
    if to_pandas:
        v = v.to_pandas()
        if specific_individuals is not None:
            indexes = set(specific_individuals)
            v = v.loc[v.individual.isin(indexes)]
    else:
        if specific_individuals:
            mask = _individual_mask(v.column(0).to_pylist(), specific_individuals)
				# this line doesn't work 
        v = {c.name:(numpy.array(c.to_pylist(), dtype=numpy.float32) if c.name != "individual" else numpy.array(c.to_pylist(), dtype=numpy.str)) for c in v}
        if specific_individuals:
            v = {k:d[mask] for k,d in v.items()}
    return v
```

We changed it or we can't perform imputation

``` pyt
def _read(file, columns=None, skip_individuals=False, to_pandas=False, specific_individuals=None):
    if columns is None:
        columns = file.schema.names
    if not skip_individuals:
        columns = ["individual"]+[x for x in columns]
    if skip_individuals and specific_individuals is not None:
            raise RuntimeError("Unsupported combination")
    v = file.read(columns=columns)
    v = v.to_pandas()
    return v
```



### Imputation

Some of the parquet files are very large, such as chromosome 1 - 12. Reading them may need a large RAM. Or the parquet pacakge can not read them and it is impossible to run the imputation. 

Unlike harmonization, imputation is not neccessary to run S-PrediXcan. In some cases, 90% of the SNPs can be found in the input file (harmonization file), the improvement of imputation may be limited.  

A shell file was created to run imputation:

``` shell
#!/bin/bash
for pheno in 'LVEF' 
do
        for chromosome in {1..21}
        do
        python gwas_summary_imputation.py \
        -by_region_file eur_ld.bed.gz \
        -gwas_file Har_chr"$chromosome"_"$pheno".txt.gz \
        -parquet_genotype reference_panel_1000G/chr"$chromosome".variants.parquet \
        -parquet_genotype_metadata /reference_panel_1000G/variant_metadata.parquet \
        -window 100000 \
        -parsimony 7 \
        -chromosome $chromosome \
        -regularization 0.1  \
        -frequency_filter 0.01 \
        --standardise_dosages \
        --cache_variants \
        -output "$pheno"/Imp_chr"$chromosome"_"$pheno".txt.gz
        
        echo "=========finish chr$chromosome $pheno==========="
        if [ $? -ne 0 ];then
        break
        fi
        done
done
```

# S-PrediXcan

The S-PrediXcan can be run by:

``` bash
#!/bin/bash
for pheno in "RAESV" "RVSV" "RVEF" "RAEF" "LASV" "RVESV" "LAEDV" "RAEDV" "RVEDV" "LAESV" "LVM" "LVSV" "LVEDV" "LVESV" "RASV" "LVEF" "LAEF"
do
python SPrediXcan.py 
--gwas_file  harmonized_gwas/All/"$pheno".txt.gz 
--snp_column panel_variant_id 
--effect_allele_column effect_allele 
--non_effect_allele_column non_effect_allele 
--zscore_column zscore 
--beta_column effect_size 
--pvalue_column pvalue 
--model_db_path  mashr_Muscle_Skeletal.db 
--covariance mashr_Muscle_Skeletal.txt.gz   
--keep_non_rsid 
--overwrite 
--model_db_snp_key varID 
--throw --output_file /Muscle_Skeletal_"$pheno".csv
echo "======finish Muscle_Skeletal "$pheno"======="
done
```

`--snp_column`, `--effect_allelee_column`, `--not_effect_allele_column`, `--zsocre_column`, `--beta_column`, `--pvalue_column` are used to specify the column name.

`mashr_Muscle_Skeletal.db` , `mashr_Muscle_Skeletal.txt.gz ` and models for other tissures can be found [here](http://predictdb.org/)



