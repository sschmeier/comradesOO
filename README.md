# comradesOO Package

# COMRADES experiment


The COMRADES experimental protocol for the prediction of RNA structure in vivo was first published in 2018 (Ziv et al., 2019) where they predicted the structure of the Zika virus. The protocol has subsequently been use to predict the structure of SARS-CoV-2 (Ziv et al., 2020).

To gain a better understanding of the protocol see:


* COMRADES determines in vivo RNA structures and interactions. (2018). Omer Ziv, Marta Gabryelska, Aaron Lun, Luca Gebert. Jessica Sheu-Gruttadauria and Luke Meredith, Zhong-Yu Liu,  Chun Kit Kwok, Cheng-Feng Qin, Ian MacRae, Ian Goodfellow , John Marioni, Grzegorz Kudla, Eric Miska.  Nature Methods. Volume 15. https://doi.org/10.1038/s41592-018-0121-0   

* The Short- and Long-Range RNA-RNA Interactome of SARS-CoV-2. (2020). Omer Ziv, Jonathan Price, Lyudmila Shalamova, Tsveta Kamenova, Ian Goodfellow, Friedemann Weber, Eric A. Miska. Molecular Cell,
Volume 80
    https://doi.org/10.1016/j.molcel.2020.11.004



After sequencing, short reads are produced where one half of the read corresponds to one half of an RNA duplex and the other half of the reads corresponds to the other half of the RNA duplex. This package has been designed to analyse this data. The short reads need to be processed is a specific way, see the next section. 

# COMRADES data pre-processing

* Cutadapt (parameters: )

    https://cutadapt.readthedocs.io/en/stable/
* 2 Pear (parameters: )

    https://cme.h-its.org/exelixis/web/software/pear/
* 3 UMI Removal (parameters: )

    https://readthedocs.org/projects/umi-tools/downloads/pdf/stable/
* 4 Hyb (parameters: )

    https://github.com/gkudla/hyb


## Prerequisits 

* A number of R libraries
* The vienna Packagae (if performing folding)

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
