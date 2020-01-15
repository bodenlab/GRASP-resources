#!/usr/bin/env python3

import seqdata, sys

"""
getclosestgene

expects as input arguments: 
         [1] filepath to chip-seq peaks (.bed)
         [2] filepath to reference genome annotation file (of gene locations in genome) (.bed)
         [optional]
         --dist maximum distance to consider binding to be significant to a gene (default, 300bp)        
"""
if len(sys.argv) < 3:
    print('Usage: getclosestgene <chip-peak-bed-file> <reference-genome-annotation-file> [--dist <promoter region bp, default 300bp>]')
    print('[] represents optional parameters')
    sys.exit(2)

# parse the input arguments
chipbedfile = sys.argv[1]
refgenelocfile = sys.argv[2]
dist = 300
if len(sys.argv) > 3 and sys.argv[3] == '--dist':
    dist = sys.argv[4]

# load the chip seq peak bed file
sites = seqdata.BedFile(chipbedfile, 'Peaks')

# load the reference genome annotation file
genes = seqdata.BedFile(refgenelocfile, 'Optional')

# find the closest gene to each binding site and write to a bed output file
f = open('tss_gene.bed', 'w')
for site in sites:
    geneInfo = genes.closest(site) # finds the closest gene to the binding site, returns a tuple 
                                   # detailing the distance and the BED file entry
    # ignore if a null entry or if the distance is greater than the specified promoter region (dist)
    if geneInfo == None or int(geneInfo[0]) > int(dist):
        continue
    # record target gene if the binding site is upstream of the gene (this relies on strand information
    # from the annotation file)
    gene = geneInfo[1]
    if (gene.strand == '+' and gene.chromStart >= site.chromStart) or (gene.strand == '-' and gene.chromStart <= site.chromStart): 
        entry = (site.chrom + '\t' + str(site.chromStart) + '\t' + str(site.chromEnd) + '\t' 
                + gene.name + '\t' + str(geneInfo[0]) + '\t' + gene.strand + '\n')
        f.write(entry)
f.close()
