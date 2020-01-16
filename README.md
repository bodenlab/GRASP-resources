<p align="center">
	<img src="/Files/Repository/reslogo.png?raw=true" alt="GRASP-resources"/>

</p>
This is the GRASP-resources page, where you will find -

[1. Supplementary material for the GRASP paper](#supplementary-material)

[2. Python notebooks designed to help with curating and analysing ancestral sequence reconstruction datasets](#notebooks)

You can also [learn more about GRASP at the GRASP-suite website](http://bodenlab.github.io/GRASP-suite) and [use GRASP now](http://grasp.scmb.uq.edu.au)

# Supplementary material

All of the supplementary material for the GRASP paper is stored in the Supplementary Material folder. Refer to the README within that folder for further information.

# Notebooks 

These Jupyter notebooks are split into two sections.

Curation - aligning, curating, and handling files before ancestral inference 

Post Inference Analysis - analysing data sets after ancestral inference

## How to use this repository ##

1. Clone this repository to your desktop

```
git clone https://github.com/bodenlab/GRASP-resources.git

```

2. Install Jupyter Notebook (python >=3.5)

[Here are the instructions to install Jupyter Notebook](https://jupyter.readthedocs.io/en/latest/install.html)

3. Install the required Python modules as specified in requirements.txt

```
pip install -r requirements.txt
```

Some notebooks require additional code that is stored in the /src folder. As long as you keep the src folder in the same relative location to the notebooks this will run correctly.

For Curation 5, the <u>standard package</u> of MAFFT is required for multiple sequence alignment.

[Here are the instructions to install MAFFT](https://mafft.cbrc.jp/alignment/software/)

4. Now you can start a Jupyter notebook from the main folder

```
jupyter notebook
```

And you will be able to navigate to the different notebooks and run the Python code within them.


## Notebooks Table of Contents ##

* **Curation 1 - Basic file handling**

This notebook shows ways to read FASTA files into Python and perform basic operations on them.

* **Curation 2 - Sequence curation**

This notebook shows how to filter sequence data sets on basis of their headers and how to summarise the species information within them.

* **Curation 3 - Checking exon counts**

This notebook shows how to query NCBI database to retrieve exon counts for a sequence data set.

* **Curation 4 - Mapping exon structure**

This notebook shows how to map the exon structure information onto a multiple sequence alignment.

* **Curation 5 - Sequence curation for ancestral sequence reconstruction**

This notebook shows how to automatically and iteratively remove sequences from a data set on the basis of length, bad characters, motifs, and internal deletions.


* **Post inference analysis 1 - Analysis of fractional distance**

This notebook allows you to analyse how the amino acid sequence at equivalent nodes changes as we increase data set size. You can specify nodes of interest in the smallest data set, which are then mapped to the equivalent nodes in the larger data sets, and then the fractional distance is calculated and plotted for all given nodes. This analysis was performed in the GRASP paper (see Figure 3).

The default notebook uses the DHAD and CYP2 data sets and recreates figures from the GRASP paper, however it can easily be adapted to your own data sets.
