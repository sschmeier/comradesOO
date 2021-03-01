# comradesOO 

# COMRADES experiment


The COMRADES experimental protocol for the prediction of RNA structure in vivo was first published in 2018 (Ziv et al., 2019) where they predicted the structure of the Zika virus. The protocol has subsequently been use to predict the structure of SARS-CoV-2 (Ziv et al., 2020).

To gain a better understanding of the protocol see:


* COMRADES determines in vivo RNA structures and interactions. (2018). Omer Ziv, Marta Gabryelska, Aaron Lun, Luca Gebert. Jessica Sheu-Gruttadauria and Luke Meredith, Zhong-Yu Liu,  Chun Kit Kwok, Cheng-Feng Qin, Ian MacRae, Ian Goodfellow , John Marioni, Grzegorz Kudla, Eric Miska.  Nature Methods. Volume 15. https://doi.org/10.1038/s41592-018-0121-0   

* The Short- and Long-Range RNA-RNA Interactome of SARS-CoV-2. (2020). Omer Ziv, Jonathan Price, Lyudmila Shalamova, Tsveta Kamenova, Ian Goodfellow, Friedemann Weber, Eric A. Miska. Molecular Cell,
Volume 80
    https://doi.org/10.1016/j.molcel.2020.11.004


![Figure from Ziv et al., 2020. Virus-inoculated cells are crosslinked using clickable psoralen. Viral RNA is pulled down from the cell lysate using an array of biotinylated DNA probes, following digestion of the DNA probes and fragmentation of the RNA. Biotin is attached to crosslinked RNA duplexes via click chemistry, enabling pulling down crosslinked RNA using streptavidin beads. Half of the RNA duplexes are proximity-ligated, following reversal of the crosslinking to enable sequencing. The other half serves as a control, in which crosslink reversal proceeds the proximity ligation](https://github.com/JLP-BioInf/comradesOO/tree/master/vignettes/comradesProtocol.jpg){width=100%}


After sequencing, short reads are produced where one half of the read corresponds to one half of an RNA duplex and the other half of the reads corresponds to the other half of the RNA duplex. This package has been designed to analyse this data. The short reads need to be processed is a specific way, see the next section. 


# COMRADES data pre-processing

The quickest way to get going with the analysis of your RNA interest is to download a dataset that has previously been processed and stored in GEO:

* https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154662

Download the files that end in `.hyb` or `.hyb.gz` or similar. These are the the files produces by step 5 below and can be used as input to this packagae, if you do this, skip the data pre-processing section.

If you would prefer to analyse the raw data see the steps below:


1 Raw data ( .fastq ) can be downloded from: 

    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154662 or from
    the other experiment mentioned above. 

2 Cutadapt - Trimming adapters

    parameters:    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTC 
                   -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
                   --minimum-length 10 
    https://cutadapt.readthedocs.io/en/stable/
    
3 Pear - Joining Paired-end reads

    parameters:    -f R1.fastq
                   -r R2.fastq
    https://cme.h-its.org/exelixis/web/software/pear/
    
4 collapse.py - UMI Removal 

    parameters: --minreads 1
    https://readthedocs.org/projects/umi-tools/downloads/pdf/stable/
    
5  Hyb  -  Finding structural duplexes 
    
    https://github.com/gkudla/hyb



## Installation

Current install is via GitHub. 

```{r}
devtools::install_github("JLP-BioInf/comradesOO")
```



## Prerequisits 



* The vienna Packagae (if performing folding)


* R packages (will install automatically):
     
    + ggplot2
    + reshape2
    + GenomicRanges
    + igraph
    + heatmap3
    + MASS
    + mixtools
    + RColorBrewer
    + foreach
    + doParallel
    + R4RNA
    
    

## Inputs

* Sample Table
* Hyb output files
* RNA of interest (as annotated in the Hyb output files)
* RNA length
* RNA Sequence (if performing folding)



# comradesOO Package

* figures of what the package does  - clustering - folding 
* links to each of the classes 
* explanation of slots 
* methods


# 3 classes

The package revolves around 3 connected objects. The comradesDataSet, comradesClusteredDataSet and comradesFoldedDataSet. 

## comradesDataSet

detailed explanation of the class with slots etc

## comradesClusteredDataSet

detailed explanation of the class with slots etc

## comradesFoldedDataSet

detailed explanation of the class with slots etc

# Standard Workflow

code for a standard workflow. 
