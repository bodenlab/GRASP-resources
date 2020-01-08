This is the GRASP-resources page, where you will find -

1. Python notebooks designed to help with curating and analysing ancestral sequence reconstruction datasets
2. Supplementary material for the GRASP paper

You can also [learn more about GRASP at the GRASP-suite website](http://bodenlab.github.io/GRASP-suite) and [use GRASP now](http://grasp.scmb.uq.edu.au)


## How to use this repository ##

Clone this repository to your desktop

You will now have

## Notebooks ##

These Jupyter notebooks are split into two sections.

Curation - aligning, curating, and handling files before ancestral inference 

Post Inference Analysis - analysing data sets after ancestral inference

## Notebooks Table of Contents ##

* Curation 1 - Basic file handling

* Post inference analysis 1 - Analysis of fractional distance

This notebook allows you to analyse how the amino acid sequence at equivelant nodes changes as we increase data set size. You can specifiy nodes of interest in the smallest data set, which are then mapped to the equivelant nodes in the larger data sets, and then the fractional distance is calculated and plotted for all given nodes. This analysis was performed in the GRASP paper (see Figure 3).

The default notebook uses the DHAD and CYP2 data sets and recreates figures from the GRASP paper, however it can easily be adapted to your own data sets.