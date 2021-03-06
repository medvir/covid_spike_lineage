# Identify interesting lineages of SARS-CoV-2 virus based on a part of Spike gene. 

Detect the following interesting lineages from a part of Spike sequence considering mutations in the range of (aa403 to 760) using ab1 or fasta files:


- B.1.1.7
- SA:B.1.351
- Brazil:P.1
- Brazil:P.2
- P.3 
- B.1.525
- B.1.526
- Cal B.1.427/B.1.429
- B.1.617.1
- B.1.617.2
- B.1.617.3


The sequences should be in ab1 or fasta file format. It is possible to have several overlapping sequences for the same patient. The files belonging to the same patient should start with a patient number in the following format:
`[PATINET_NUMBER]_*.ab1` or `[PATINET_NUMBER]_*.fasta`


## Requirements:

- emboss seqret
- seqtk
- biopython
- pandas > 1.0.0
- mafft
- tabulate

## Installation
The installations ``` python setup.py install ``` or ```pip install``` works. 

## Usage
Run 


```
covid_spike_lineage
```  


in the ab1 files directory.

or run the following command 


```
covid_spike_lineage [OPTIONS]

Options:
    -d --directory  Folder to input file and output files including the pdf report. [Default: The currecnt working directory]
    -e --extention  The input files extensions. options are fasta or ab1. [Default: ab1]
```


 
