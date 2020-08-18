---
title: 'squirrel RAD pipeline processing'
disqus: hackmd
---

# Red Squirrel RAD data analysis
> [name=Tanya Lama et al.] 

Note: Nigel's demultiplexed arctic squirrel reads are in /STACKS/03_clone_filter_out_arctic_squirrel
## **Background**: 
Red squirrels in the White Mountains National Forest may be dispersing to higher elevations where they prey on threatened bird species. This may be due to climate, demography (population expansion), biennial mast, or perhaps all three. See [here](https://www.theatlantic.com/science/archive/2019/01/climate-refugia-shelter-species-changing-world/579802/) for more details on PI Dr. Toni Lyn Morelli's work with red squirrels and climate refugia

## **Goal**: 
Explore genetic structuring among red squirrels in the White Mountains National Forest using RADseq data. 

## **Objectives**: 
Use a bioinformatic pipeline (ipyrad) to process RADseq data from raw reads to filtered SNPs in VCF format. Then use the VCF as input data for desired analyses of structure, gene flow, and landscape connectivity.

As described in my USGS contract are: 
o Bioinformatics of RADseq genomic data from red squirrels
o Analysis of gene flow and population structure
o Preparation of figures and tables

---

## Table of Contents

[TOC]

---
## Beginners Guide

If you are a total beginner to this, start here!
``` 
ls
#the squirrel environment now includes ipyrad dependencies, fastqc and chainmap
```

### 1. Sign in
Nigel decided to install a Linux system. If that doesn't work for him we will use PuTTY
```
ssh tl50a@ghpcc06.umassrc.org
ssh ng56a@ghpcc06.umassrc.org #nigel's sign in
```
### 2. Loading software with modules
Navigate to your home directory (/home/tl50a) on the cluster. Determine whether anaconda is available to load as a module using the list function. Once you’ve confirmed that a version is available, load and confirm it’s been installed.
```
module load anaconda2/4.4.0 
```
*conda* is a program that it easy to install dependencies for a given project. You can see what programs are available to install using:
```
conda search
```
Or searching the web interface at https://anaconda.org/search?q=grep

If you've never set up a conda environment before this **tutorial** is very helpful: https://github.com/thw17/BIO598_Tutorial
conda info #confirm conda is installed
module avail #list of all modules available
module list #confirm which modules are loaded

### 3. Add channels and set a path for the environment and packages
We will mostly use tools available in bioconda (fastqc, samtools etc). Note: If there is not enough room in your home directory, set the path to /nobackup/username instead.
```
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels bioconda #bioinformatics
```
#For macOS and Linux: Specify directories in which packages are located. If this key is set, the root prefix pkgs_dirs is not used unless explicitly included.

CONDA_ENVS_PATH=~/my-envs:/opt/anaconda/envs
CONDA_PKGS_DIRS=/opt/anaconda/pkgs
#You now have your own **my_root environment** directories and package cache. This means you can install packages using conda install (rather than using the read-only shared version of python available through module load).
Use *conda info* to confirm you are either in myroot or a conda environment before installing packages. If you are in the shared root, you will not have permission to install packages.

> [envs directories : /home/tl50a/.conda/envs
> /share/pkg/anaconda2/4.4.0/envs
> package cache : /share/pkg/anaconda2/4.4.0/pkgs
> /home/tl50a/.conda/pkgs]

### 4. Install packages
There are several tools we need to make available in our environment. The == syntax specifies the version you want to install. Use the conda search function above to ensure packages are available fore installing.
```
conda search ipyrad
conda install -c ipyrad ipyrad
conda list #Use the conda list function to review packages installed in *squirrel*
```
The following Python packages are installed as dependencies of ipyrad:

numpy
scipy
pandas
h5py
mpi4py ##we need this is we want to run parallelization
numba
ipyparallel
pysam
cutadapt
requests

### 5. Create an environment for your project
Our project environment is called ***squirrel***. Note: Python 3.6 is recommended because the end-of-life date for Python 2 is scheduled for 2020.
```
conda create --name squirrel python==3.6 ipyrad
```
The active environment is also displayed in front of your prompt in (parentheses) or [brackets] like this:

[(squirrel) [tl50a@ghpcc06]$…]

### 6. Save your environment
We need to save our environment and all the packages we’ve loaded into a .yaml folder in our project directory. This is now a “portable” environemnt that I could use on any machine/with new data etc. #we might skip this
#need to add code here

#### Prior to running any scripts, ALWAYS:
```
module load anaconda2/4.4.0 
source activate squirrel
```
#### The environment can be deactivated with:
```
source deactivate squirrel
```
---
## Getting Started
### cksum summary for raw reads
We downloaded the raw data from the sequencing center direct to the cluster using an ftp link. 

The first thing we'll do is run cksums on the raw data and archive it safely for long-term storage (Lisa archived the data). cksums basically helps ensure that file sizes are the exact size they should be, and none have been corrupted in the transfer process. Here is a list of cksums for each of the files we acquired from the sequencing center: 

* Copied to: /project/uma_lisa_komoroske/raw_reads/red_squirrel_Morelli
* cksums confirmed
* Archived by Lisa 

8285630061      ./Rawdata/Undetermined/Undetermined_H7HH3BBXX_L7_1.fq.gz
9136401141      ./Rawdata/Undetermined/Undetermined_H7HH3BBXX_L7_2.fq.gz
8068088202      ./Rawdata/Undetermined/Undetermined_H7HH3BBXX_L8_1.fq.gz
8907563399      ./Rawdata/Undetermined/Undetermined_H7HH3BBXX_L8_2.fq.gz
34275308803     ./Rawdata/RAD_SV1/RAD_SV1_CKDL200148269-1B_H7HH3BBXX_L7_2.fq.gz
34013595369     ./Rawdata/RAD_SV1/RAD_SV1_CKDL200148269-1B_H7HH3BBXX_L8_2.fq.gz
31099364277     ./Rawdata/RAD_SV1/RAD_SV1_CKDL200148269-1B_H7HH3BBXX_L8_1.fq.gz
31316620338     ./Rawdata/RAD_SV1/RAD_SV1_CKDL200148269-1B_H7HH3BBXX_L7_1.fq.gz -- yes moved to 
185370  ./Rawdata/Rawdata_Readme.pdf
272     ./Rawdata/Undetermined/MD5.txt

### Directory structure

In all the scripts, I used relative paths instead of absolute path. The reason for this is so that it can be easily modified for subsequent analysis of new samples.

1. In the main directory, make two directories – scripts and analyses.
```mkdir scripts analyses```
2. In the subdirectories scripts and analyses, create matching folders for each step
```
cd scripts
mkdir step1_fastqc step2_ step3_ step4_ step5_ step6_ step7_

cd analyses
mkdir step1_fastqc step2_ step3_ step4_ step5_ step6_ step7_
```
3. Then, within each subdirectory, make a directory for each sample. This can be done easily in bash:

cd scripts
for dir in */; do                                                            
for i in name name name name ; do
mkdir -- "$dir/$i";    
done;
done;

#Repeat for the analyses directory:

cd /home/tl50a/analyses
for dir in */; do                                                            
for i in ....; do
mkdir -- "$dir/$i";    
done;
done;

4. cd into the step1_fastqc folder and use ls to confirm you have created the directory correctly. 
5. In the main directory, create a third directory called download.
```mkdir download```
6. Within the download directory, create two subdirectories called fastq and refs
```
mkdir fastq
mkdir refs
```
Any fastq (fq.gz or fastq.gz) data accessible by the executable .sh scripts needs to be stored in a folder under each individual name

### Fastq Data Files and File Names

Depending how and where your sequence data were generated you may receive data as one file, or in many smaller files. The files may contain data from all of your individuals mixed together, or as separate files for each sample. If they are mixed together then the data need to be **demultiplexed** based on barcodes or indices. Step 1 of ipyrad can take data of either format, and will either demultiplex the reads or simply count/load the pre-demultiplexed data.

There is a large variety of ways to generate reduced representation genomic data sets using either restriction digestion or primer sets, and ipyrad is flexible enough to handle many of these types. 

rad – This category includes data types which use a single cutter to generate DNA fragments for sequencing based on a single cut site. e.g., RAD-seq, NextRAD.

ddrad – These methods produce fragments that were digested by two different restriction enzymes which cut the fragment on either end. During assembly this type of data is analyzed differently from the rad data type by more stringent filtering that looks for occurrences of the second (usually more common) cutter. e.g., double-digest RAD-seq.

## step1_fastqc

It is generally good practice to run the program fastqc on your raw data to get an idea of the quality of your reads, and the presence of adapter contamination. You do not need to trim your reads before starting an assembly, since ipyrad includes a built-in and recommended trimming step during step 2 of assembly (using the software tool cutadapt). If you do choose to trim your data beforehand, however, it should not cause any problems.

Working directory is scripts/step1_fastqc
1. wget .sh and wrapper.sh scripts from github
vi scripts/step1_fastqc/wrapper.sh
"i" to enter insert mode
2. uncomment the individuals array line and insert a list of sample names

If you are in Insert mode you have to leave it first with ESC.
If you *have made changes and want to save* the file, use :x
If you *haven’t made any changes*, press :q to leave the file (but you must not be in Insert mode).
If you have made changes you *don't want to save* use :q!
3. change permissions to allow execution
chmod +x step1_fastqc.sh
4. Submit a job using fastqc to check the raw reads. Script used is step1_fastqc.sh. Usage is:
./step1_fastqc.sh individual_name -o /pathto/analyses/step1/

Jobs are run from the directory/environment that the job is submitted from and unless specified otherwise, output files will be written into the directory the job was submitted from. For your script to work, you will need to call files by the path relative to the directory you submit your job from
```
bsub -q short -W 1:00 -R rusage[mem=4000] -n 1 -R span\[hosts=1\] "./step1_fastqc.sh sampleid -o /analyses/step1" 
```
Usage is: submitjob -queue short -length 1hour -memory memoryusage[1000Mb] -ncores 1core -nhosts [hosts=1] "fastqc name.fq.gz -output /path/to/analyses/step1/"
```
bsub -q short -W 1:00 -R rusage[mem=4000] -n 1 -R span\[hosts=1\] "fastqc RAD_SV1_CKDL200148269-1B_H7HH3BBXX_L7_1.fq.gz -o /path/to/analyses/step1/"
```
bjobs #check status 

Note: See the wrapper script wrapper_step1_fastqc.sh for how to automate to multiple samples.

Note: We are using paired end sequencing, so (should) run fastqc on both reads for each sample. Eventually we will align all (forward and reverse) reads to the reference.

Note: We used vi command to change step1_fastqc.sh from *.fastq.gz to .fq.gz to match our file type.

5. Download the fastqc report to your home directory
To copy from the **remote** computer to the **local** one, type, in the local computer terminal:

"scp -r tl50a@ghpcc06.umassrc.org:/pathto/analyses/step1/filename.html ./"

6. Interpret the results of the fastqc report: 
[looks good, summarize with Nigel]

### create an assembly and params file for each lane
ipyrad uses a simple text file to hold all the parameters needed for a given assembly. Start by creating a new params file using the -n flag, followed by a name for your assembly. In the example we use the name squirrel, but the name can be anything at all. Once you start analysing your own data you might call your params file something more informative, like the name of your organism. We will refer to this as the “assembly_name”.
```
ipyrad -n squirrel_L7
ipyrad -n squirrel_L8
```
> New file 'params-squirrel_L7.txt' created in /project/uma_lisa_komoroske/Morelli_RAD/scripts

There are four assembly_methods options in ipyrad: denovo, reference, denovo+reference, and denovo-reference. We chose denovo here.

The params file lists on each line one parameter followed by a ## mark, then the name of the parameter, and then a short description of its purpose. Take a look at it by using the unix command ‘cat’ (or you can use any text editor you like).

```
cat params-squirrel.txt
```
> redsquirrel ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
> ./location  ## [1] [project_dir]: Project dir
>             ## [2] [raw_fastq_path]: Location of raw **non-demultiplexed** fastq files

In general the default parameter values are sensible, and we won’t mess with them for now, but there are a few parameters we **must** change. We need to set the path to the raw data we want to analyse, and we need to set the path to the barcodes file.

In your favorite text editor (nano is a popular command line editor for linux, for Mac you can use TextEdit) open params-squirrel.txt and change these two lines to look like this, and then save it.
```
/project/uma_lisa_komoroske/Morelli_RAD/download/fastq/RAD_SV1_CKDL200148269-1B_H7HH3BBXX_L7_R1.fastq.gz ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
/project/uma_lisa_komoroske/Morelli_RAD/download/fastq/barcodes.txt ## [3] [barcodes_path]: Location of barcodes file
```
### execute step 1 of ipyrad for both assemblies
Now we will start assembling the data with ipyrad. Step 1 reads in the barcodes file and the raw data. It scans through the raw data and sorts each read based on the mapping of samples to barcodes. At the end of this step we’ll have a new directory in our project_dir called ./squirrel_fastqs. Inside this directory will be individual fastq.gz files for each sample.

**Note1**: UMA2 (sample name) is duplicated in barcodes.txt - 2 unique barcodes both match to UMA2; resultant fastq file will need to be concenated
```
bsub -q short -W 2:00 -R rusage[mem=8000] -n 2 -R span\[hosts=1\] "ipyrad -p params-squirrel_L7.txt -s 1 -r -f" 

bsub -q long -W 24:00 -R rusage[mem=8000] -n 2 -R span\[hosts=1\] "ipyrad -p params-squirrel.txt -s 3 -r -f"
```

### merge the two lanes into one assembly with a hybrid name
```
bsub -q long -W 8:00 -R rusage[mem=8000] -n 2 -R span\[hosts=1\] "ipyrad -m both params-squirrel_L7.txt params-squirrel_L8.txt"
```
view stats for the merged assembly
```
ipyrad -p params-both.txt -r
```
## step2 Filtering/Trimming data
This step filters reads based on **quality scores** and/or the occurrence of barcode+adapter combinations, and can be used to detect Illumina adapters in your reads, which is a common concern with any NGS data set, and especially so for local lab library prep protocols. Here the filter is set to the default value of 0 (zero), filters are only based on quality scores of base calls, and does not search for adapters. This is a good option if your data are already pre-filtered. The resuling filtered files from step 2 are written to a new directory called squirrel_edits/. For paired-end data ipyrad will merge overlapping reads (using vsearch for denovo assembly or simply based on mapping positions for reference-mapped assembly).

### execute step 2 of ipyrad
```
bsub -q long -W 24:00 -R rusage[mem=64000] -n 1 -R span\[hosts=1\] "ipyrad -p params-squirrel.txt -s 2 -r -f"
```
Because we sequenced our individuals on two lanes of the NovaSeq, we will need to combine our raw fastq reads before proceeding with steps 3-7. 

Usage is:
```
cat file1.fastq file2.fastq > mergedfile.fastq
```

**Note1**: We want to do the merge *after* demultiplexing (step 2) and before proceeding with step 3.

**Note2**: It's somewhat unusual that we would need to combine date from two lanes like this. The reason we sequenced data on two lanes is that the red squirrel genome is so large (6+ Gb), we would have had really scant coverage using just one lane of data. For perspective, the human genome is 2 Gb and the Canada lynx is 2.4. Most fish genomes are in the 1 Gb range, and bats and turtles including green and leatherback are also in the 2.1 range, so 6+ Gb is a really big genome for a little mammal like the red squirrel. 

### check step 2 output
```
ls iptest_edits/
```
Get current stats including # raw reads and # reads after filtering. Fortunately ipyrad tracks the state of all your steps in your current assembly, so at any time you can ask for results by invoking the -r flag.
```
ipyrad -p params-squirrel.txt -r
```
The number of raw reads and filtered reads is one of the checkpoints required by Lisa, so we've saved them on github as [here](https://github.com/ECOtlama/project_red_squirrel_RAD/blob/master/reads_raw_reads_passed_filter) 
You might also take a look at the filtered reads:
```
head -n 12 ./squirrel_edits/name_R1_.fastq
```
If you want to get even more info ipyrad tracks all kinds of wacky stats and saves them to a file inside the directories it creates for each step. For instance to see full stats for step 1:
```
 cat ./squirrel_fastqs/s1_demultiplex_stats.txt
```

---
## STACKS demultiplexing instead of step2
After running step2 for ipyrad, I was suspicious that not all of the reads were getting demultiplexed properly. the biggest clue to this was the number of no_match reads listed in s1_demultiplex_stats.txt. After some troubleshooting (separate ipyrad runs, merged, etc), I decided that we should try a different protocol that demultiplexes data specifically for the bestRAD. 

Follow the instructions [here](https://hackmd.io/YdkcvtksQNijTDcA1qvd5A?view#Stacks-Pipeline---Streamlined) for STACKS demultiplexing steps 1-3. 

## step02 prelim QC
We're going to run 02_prelim_QC.sh to compile some stats for Lisa to review as well (n reads, n cut sites). Evidently the BestRAD protocol can leave the cutsite on the forward or reverse reads, so both need to be checked and that's likely why ipyrad was not able to demux successfully.

Make two directories for L7 and L8 and follow the 01.sh script for making directories
```
mkdir STACKS STACKS_L7
cd STACKS
./01_mkdirs.sh #then do the same for STACKS_L7
```
Usage is: 
```
bsub -q long -W 10:00 -R rusage[mem=16000] -n 4 -R span\[hosts=1\] "./02_prelim_QC.sh"
```
## step03 process data
We need to run step 03_process_data.sh for L7 and L8 lanes of data separately, then cat the trimmed and sorted reads and use those as "sorted" fastq files in our ipyrad assembly step3 (below).
Running step3 on L8 8629452 
```
bsub -q long -W 24:00 -R rusage[mem=16000] -n 4 -R span\[hosts=1\] "./03_process_data.sh" 
```
Running step3 on L7 in STACKSL7 8629454
```
bsub < ./03_process_data.sh 
```
## merge L7 and L8 fastq files 
Now that our reads are individuall demultiplexed and trimmed, we need to merge the data from L7 and L8. Note: R1 reads should be joined with R1 and R2 with R2, do not mix and match. 

Step 1: Grab the unique sample ID names as a txt file
```
ls -1 *.1.fq.gz | awk -F '.' '{print $1}' | sort | uniq > ID
```
Step 2: Use a for loop to merge L7 and L8 files one record at a time. 
```
for i in `cat ./ID`; do echo cat /project/uma_lisa_komoroske/Morelli_RAD/STACKS/03_clone_filter_out/$i\.2.fq.gz /project/uma_lisa_komoroske/Morelli_RAD/STACKS_L7/03_clone_filter_out/$i\.2.fq.gz \> /project/uma_lisa_komoroske/Morelli_RAD/STACKS/03_clone_filter_out_merged/$i\.2.fq.gz; done > merge2.sh
```
Step 3: Run all of the provided cat commands
```
cat /project/uma_lisa_komoroske/Morelli_RAD/STACKS_L7/02_process_out/JF2.1.fq.gz /project/uma_lisa_komoroske/Morelli_RAD/STACKS/02_process_out/JF2.1.fq.gz > JF2_R1.fq.gz
```
Step 4: Run again on R2 reads

## remove arctic squirrel samples
Move all of the PE* acrtic squirrel individuals to /03_clone_filter_out_arctic_squirrel

## Proceed with step4 of STACKS
```
bsub < ./04_run_ustacks.sh
```
### Proceed with step5 parameter opt of STACKS
nano hs.txt and see sample selection in instructions (highest coverage based on R analyses)
```
bsub < ./05_parameter_test.sh
```


---
## (back to ipyrad) step3_clustering_mapping_reads
Step 3 first dereplicates the sequences from step 2, recording the number of times each unique read is observed. If the data are paired-end, it then uses vsearch to merge paired reads which overlap. The resulting data are then either de novo clustered (using vsearch) or mapped to a reference genome (using bwa and bedtools), depending on the selected assembly method. In either case, reads are matched together on the basis of sequence similarity and the resulting clusters are aligned using muscle.
```
bsub -q short -W 1:00 -R rusage[mem=8000] -n 2 -R span\[hosts=1\] "ipyrad -p params-squirrel.txt -s 3 -r -f"
```
bsub -q long -W 48:00 -R rusage[mem=8000] -n 2 -R span\[hosts=1\] "ipyrad -p params-squirrel_v2.txt -s 1234567 -r -f"
## STACKS parameter opt 8733738
## STACKS (all samples, real VCF) 8735431

## ipyrad step 567 on ipyrad-squirrel 8817254 #almost there!

## ipyrad 1234567 on ipyrad-squirrelv2 8733555 #still on step 3

## trial2 8735639
bsub -q long -W 72:00 -R rusage[mem=64000] -n 2 -R span\[hosts=1\] "ipyrad -p params-trial4.txt -s 34567 -r -f"
## trial3 8735640

## trial4 8736595

## trial5 8735647

## step4_het_and_error_rates 
Note: runtime was <1h and 10432 MB
Step4 jointly estimates sequencing error rate and heterozygosity based on counts of site patterns across clustered reads. These estimates are used in step5 for consensus base calling. If the max_alleles_consens is set to 1 (haploid) then heterozygosity is fixed to 0 and only error rate is estimated. For all other settings of max_alleles_consens a diploid model is used (i.e., two alleles are expected to occur equally).
## step5_concensus_base_call_filtering
Step5 estimates consensus allele sequences from clustered reads given the estimated parameters from step 4 and a binomial model. During this step we filter for maximum number of undetermined sites (Ns) per locus (max_Ns_consens). The number of alleles at each locus is recorded, but a filter for max_alleles is not applied until step7. Read depth information is also stored at this step for the VCF output in step7.
## step6_clustering_mapping_alignment
Step6 clusters consensus sequences across Samples using the same assembly method as in step 3. One allele is randomly sampled before clustering so that ambiguous characters have a lesser effect on clustering, but the resulting data retain information for heterozygotes. The clustered sequences are then aligned using muscle.
## step7_filtering_formatting_output
Step7 applies filters to the final alignments and saves the final data in a number of possible output formats. This step is most often repeated at several different settings for the parameter 21. min_samples_locus to create different assemblies with different proportions of missing data (see Assembly: Branching and Merging).

---
## Additional filtering with GATK
## Additional filtering with LEA scripts in R
## Clustering analyses sNMF and PCA
## Estimated Effective Migration Surface (EEMS) modeling
## Visualizations
## Interpretation and Summary

## Appendix and FAQ
### Documenting command-line syntax (Nigel)
$ ssh [USERNAMME][HOST]: protocol used to establish a secure connection with host system
$ pwd : print the current working directory
$ ls : list files and directories
$ cd : change the current directory (i.e., the directory in which the user is currently working)
$ vi [FILE TYPE]: open a file (i.e. vi dna.xls) and visual editor
    * $ :q : quits without saving (i.e. return)
    * $ i : insert text (used in command-line vi)
    * $ esc : turns off insert mode
    * $ :x : quit vi and save modified file
    * $ :q! : quit vi without saving
    * $ :./ : whereever I am    
$ cat: short for "concatenate", can display contents of file, multiple files in terminal 
$ less: used to view (but not hange) contents of text file one screen at a time
$ cksum : used to ensure that files transferred have not been corrupted
$ scp : allows files to be copied to, from, or between different hosts (e.g. scp #your_username@remotehost.edu:foobar.txt /some/local/directory)
$ bsub : submit a job or batch script
$ module avail : which packages are avilable
$ module load : load this package
$ bjobs : show me jobs I'm running
$ bqueues : show jobs on cluster
$ bpeek: monitors the progress of a job and identifying errors
$ mkdir: to create a directory; can create multiple directories hierarchy (e.g. mkdir folder1\folder2\folder3)
$ wget: downloads files from the internet (e.g. wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/426/925/GCF_003426925.1_ASM342692v1/GCF_003426925.1_ASM342692v1_genomic.fna.gz ./)
$ stat: displays detailed information of file
$ mv: moves files or directories from one place to another
$ source activate: activates an Anaconda environment
$ conda list: list all packages and versions installed in active environment
$ ctrl c: interrupt current job
$ nano: text editor; easier to create new documents
$ */: all files with this type

## Tanya Tracking Hours
see 1-LAMA Time Sheet 04MAY20 thru 14AUG20.xls

:::info
**Find this document incomplete?** Leave a comment!
:::

###### tags: `RAD` `red squirrel` `Genomics` 
