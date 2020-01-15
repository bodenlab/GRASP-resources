import ival

class BedEntry():

    def __init__(self, chrom, chromStart, chromEnd):
        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd
        self.usestrand = False
        self.strand = None
        self.name = ''
        self.blocks = None # interval tree with blocks

    def addOption(self,
                  name = None,
                  score = None,
                  strand = None,
                  thickStart = None,
                  thickEnd = None,
                  itemRgb = None,
                  blockCount = None,
                  blockSizes = None,
                  blockStarts = None,
                  signalValue = None,
                  pValue = None,
                  qValue = None,
                  peak = None,
                  tags = None,
                  summit = None,
                  fold = None,
                  fdr = None,
                  zscore = None,
                  bg = None):
        if name: self.name = name
        if score != None: self.score = score
        if strand:
            self.strand = strand
            self.usestrand = True # use reverse complement when sequence is requested from genome
        if thickStart != None: self.thickStart = thickStart
        if thickEnd != None: self.thickEnd = thickEnd
        if itemRgb != None: self.itemRgb = [int(color) for color in itemRgb.split(',')]
        if blockCount != None:
            if blockCount > 0:
                blockSizes = [int(sizeword) for sizeword in blockSizes.rstrip(',').split(',')]
                blockStarts = [int(startword) for startword in blockStarts.rstrip(',').split(',')]
                if len(blockSizes) != blockCount or len(blockStarts) != blockCount:
                    raise RuntimeError('Blockcount is incorrect in BED entry \"%s\"' % str(self))
                for i in range(blockCount):
                    self.addBlock(blockStarts[i], blockSizes[i])
        if signalValue != None: self.signalValue = signalValue
        if pValue != None: self.pValue = pValue
        if qValue != None: self.qValue = qValue
        if peak != None: self.peak = peak
        if tags != None: self.tags = tags
        if summit != None: self.summit = summit
        if fold != None: self.fold = fold
        if fdr != None: self.fdr = fdr
        if bg != None: self.bg = bg
        if zscore != None: self.zscore = zscore

    def __len__(self):
        if self.blocks:
            return len(self.blocks)
        else:
            return 0

    def addBlock(self, relative_start, size):
        if not self.blocks:
            self.blocks = ival.IntervalTree()
        self.blocks.put(ival.Interval(self.chromStart + relative_start, self.chromStart + relative_start + size))

    def getBlocks(self):
        return self.blocks

    def __str__(self):
        if self.strand == '+' or self.strand == '-':
           return self.chrom + ':' + str(self.chromStart) + '-' + str(self.chromEnd) + self.strand
        return self.chrom + ':' + str(self.chromStart) + '-' + str(self.chromEnd)

    def __iter__(self):
        if self.blocks:
            for b in self.blocks:
                yield (self.chrom, b.ival.min, b.ival.max)

    def isBlockOverlap(self, entry):
        if not self.blocks:
            return None
        if isinstance(entry, BedEntry):
            if (entry.chrom == self.chrom):
                return self.blocks.isect(entry.getInterval()) != None
        elif isinstance(entry, ival.Interval):
            return self.blocks.isect(entry) != None

    def loc(self, genome = None, fixedwidth = None, usesummit = False, useshift = None):
        """ Retrieve the genomic location for BED entry, or sequence if genome is provided
            genome: a dictionary with keys for sequence names, e.g. 'chr1', 'chrX', etc, and values with indexed/sliceable strings
            fixedwidth: the width of the location/sequence if the width in the BED entry is ignored, and only its centre is used
            usesummit: centre a fixedwidth window around an assigned "summit"
            useshift: centre a fixedwidth window around a shifted centre point, e.g. useshift=-125 will shiftcentre point 125bp upstream,
            to say capture a fixedwidth=350bp window with 350/2-125=50bp downstream
        """
        otherstrand = False
        if (self.usestrand):
            if (self.strand == '-'):
                otherstrand = True

        if (otherstrand == False):
            end = self.chromEnd
            start = self.chromStart
            mywidth = fixedwidth or (self.chromEnd - self.chromStart)
            mycentre = start + (self.chromEnd - self.chromStart) // 2
            if usesummit:
                mycentre = self.summit
            if useshift:
                mycentre = mycentre + useshift
            if fixedwidth: # we need to re-calculate start and end
                if genome:
                    end = min(len(genome[self.chrom]), mycentre + (mywidth // 2))
                else:
                    end = mycentre + (mywidth // 2)
                start = max(0, mycentre - (mywidth // 2))

        else: # other strand
            start = self.chromEnd
            end = self.chromStart
            mywidth = fixedwidth or (self.chromEnd - self.chromStart)
            mycentre = self.chromStart + (self.chromEnd - self.chromStart) // 2
            if usesummit:
                mycentre = self.summit
            if useshift:
                mycentre = mycentre - useshift # shift is reversed on other strand
            if fixedwidth: # we need to re-calculate start and end
                end = max(0, mycentre - (mywidth // 2))
                if genome:
                    start = min(len(genome[self.chrom]), mycentre + (mywidth // 2))
                else:
                    start = mycentre + (mywidth // 2)

        if genome: # refer to the genome sequence
            return genome[self.chrom][start : end]
        else:
            return (self.chrom, start, end)

    def setwidth(self, fixedwidth = None, usesummit = False):
        if fixedwidth:
            if usesummit:
                diff = self.summit - fixedwidth // 2
            else:
                diff = (self.chromEnd - self.chromStart) // 2 - fixedwidth // 2
            self.chromStart += diff
            self.chromStart += diff + fixedwidth
        return (self.chrom, self.chromStart, self.chromEnd)

    def getInterval(self):
        return ival.Interval(self.chromStart, self.chromEnd)

def dist(entry1, entry2, signed = False, centre2centre = False):
    """ Calculate and return the BedEntry with the closest distance (from one end of the interval of this to the end of the interval of that).
        If centre2centre is True, use the centre-to-centre distance instead.
        If signed is True, the distance is negative if this interval is after the that.
    """
    if isinstance(entry1, BedEntry) and isinstance(entry2, BedEntry):
        if (entry1.chrom == entry2.chrom):
            return ival.dist(entry1.getInterval(), entry2.getInterval(), signed, centre2centre)
    return None

class BedFile:
    """ Read BED file.

        See http://genome.ucsc.edu/FAQ/FAQformat#format1

        The first three required BED fields are (part of all supported sub-formats):

        chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
        chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
        chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.

        The 9 additional optional BED fields are (part of sub-format "Optional"):

        name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
        score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
        shade
        strand - Defines the strand - either '+' or '-'.
        thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
        thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
        itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
        blockCount - The number of blocks (exons) in the BED line.
        blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
        blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.

        ENCODE also defines broadpeaks and narrowpeaks format (part of our "Peaks" sub-format):

        name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
        score - Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were '0' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
        strand - +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
        signalValue - Measurement of overall (usually, average) enrichment for the region.
        pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
        qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
        peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.

        MACS also defines a "summit" peaks format (part of our "Summit" sub-format)
        It contains the peak summits locations for every peaks. The 5th column in this file is the .
        In addition to the required three, the following fields follow:
        length         [redundant, ignored]
        summit         summit height of fragment pileup
        tags
        pValue         [-10*log10(pvalue)]
        fold           [enrichment]
        FDR            [%; optional]

        "CCAT" BED-like file format:
        chromosome,
        peakcenter    [converted to summit],
        regionstart,
        regionend,
        tags          [tagcount],
        bg            [bgcount],
        zscore,
        fdr

    """

    def __init__(self, entries, format = 'Limited'):
        """
        Create a BedFile instance.
        :param entries: an iterable of entries or a filename
        :param format: the format of the BED file
        """
        self.format = format
        if isinstance(entries, str): # filename
            self.chroms = readBedFile(entries, format)
        else:
            self.chroms = dict()
            for entry in entries:
                # check if the chromosome has been seen before
                tree = self.chroms.get(entry.chrom)
                if not tree:
                    tree = ival.IntervalTree()
                    self.chroms[entry.chrom] = tree
                # put the entry in the interval tree for the appropriate chromosome
                iv = ival.Interval(entry.chromStart, entry.chromEnd)
                tree.put(iv, entry)

    def __len__(self):
        n = 0
        for c in self.chroms:
            n += len(self.chroms[c])
        return n

    def generate(self, chrom):
        mytree = self.chroms.get(chrom)
        if mytree != None:
            for e in mytree:
                for entry in e.values:
                    yield entry

    def __iter__(self):
        self.chromqueue = ival.Stack()
        for c in sorted(self.chroms.keys())[::-1]:
            self.chromqueue.push(self.generate(c))
        self.current = self.chromqueue.pop()
        return self

    def __next__(self):
        try:
            ret = next(self.current)
        except StopIteration:
            if not self.chromqueue.isEmpty():
                self.current = self.chromqueue.pop()
                ret = next(self.current)
            else:
                raise StopIteration
        return ret

    def __contains__(self, item):
        if isinstance(item, BedEntry):
            tree = self.chroms.get(item.chrom)
            if tree == None: return False
            else: return ival.Interval(item.chromStart, item.chromEnd) in tree
        else:
            return False

    def getOverlap(self, item):
        if isinstance(item, BedEntry):
            tree = self.chroms.get(item.chrom)
            if tree == None: return None
            else:
                iv = ival.Interval(item.chromStart, item.chromEnd)
                res = tree.isectall(iv)
                ret = []
                for r in res:
                    ret.extend(r.values)
                return ret
        else: return None

    def getClosest(self, item):
        if isinstance(item, BedEntry):
            tree = self.chroms.get(item.chrom)
            if tree == None: return None
            else:
                iv = ival.Interval(item.chromStart, item.chromEnd)
                node = tree.closest(iv)
                if node != None: return node.values
                else: return None
        else: return None


    def getOneOfClosest(self, item):
        all = self.getClosest(item)
        if all == None: return None
        else: return next(iter(all))

    def getOneOfOverlap(self, item):
        all = self.getOverlap(item)
        if all == None: return None
        elif len(all) == 0: return None
        else: return next(iter(all))

def readBedFile(filename, format = 'Limited'):
    """ Read a BED file.
        format: specifies the format of the file,
        "Limited", e.g.
            chr22 1000 5000
            chr22 2000 6000
        "Optional", e.g.
            track name=pairedReads description="Clone Paired Reads" useScore=1
            chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
            chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
            ...
            (also handles the Limited + score, and BED6 format)
        "Peaks", e.g.
            chr1    569780    569930    .    0    .    19    6.07811    -1    -1
            chr1    713300    713450    .    0    .    54    49.1167    -1    -1
        "Strand", e.g.
            chr4    185772359    185772424    -
            chr18    20513381    20513401    +
        also supports a 5th label field
            chr5    20611949        20611949        +       ENSG00000251629_20611949
            chr3    42187863        42187863        -       ENSG00000234562_42187863
        "Summit", e.g.
            # d = 130
            chr      start    end   length summit  tags -10*log10(pvalue)    fold_enrichment    FDR(%)
            chr1     8250     8671    422    286    46    145.84    11.68    0.51
            chr1    36382    36984    603    405    46    315.23    27.05    0.24
        "CCAT", e.g.
            chr8    94747805    94747070    94749250    525     3    21.519196    0.002000
            chr17   55277895    55277070    55279280    560    18    21.283333    0.002000
        "Cropped", e.g.
            chr1    851602    10
            chr1    921184    18
            chr1    931838    9
    """
    f = open(filename)
    row = 0
    acceptHeaderRows = 1
    headerRow = None
    chroms = dict()
    for line in f:
        row += 1
        words = line.strip().split()
        if len(words) == 0:
            continue # ignore empty lines
        if words[0].strip().startswith('#'):
            continue # comment
        if words[0].strip().startswith('browser'):
            continue # ignore
        if words[0].strip().startswith('track'):
            continue # ignore
        try:
            chrom = words[0]
            if format.lower().startswith('ccat'):
                chromStart = int(words[2])
                chromEnd = int(words[3])
            else: # all other standard BED formats
                chromStart = int(words[1])
                chromEnd = int(words[2])
            entry = BedEntry(chrom, chromStart, chromEnd)
            if format.lower().startswith('opt') or format.lower().startswith('bed12'):
                if len(words) >= 12:
                    entry.addOption(name = words[3], score = float(words[4]), strand = words[5], thickStart = int(words[6]), thickEnd = int(words[7]), itemRgb = words[8], blockCount = int(words[9]), blockSizes = words[10], blockStarts = words[11])
                elif len(words) >= 9:
                    entry.addOption(name = words[3], score = float(words[4]), strand = words[5], thickStart = int(words[6]), thickEnd = int(words[7]), itemRgb = words[8])
                elif len(words) >= 6:
                    entry.addOption(name = words[3], score = float(words[4]), strand = words[5])
                elif len(words) >= 5:
                    entry.addOption(name = words[3], score = float(words[4]))
                elif len(words) >= 4:
                    entry.addOption(name = words[3])
                else:
                    entry.addOption(name = '.', score = int(words[3]), strand = '.')
            elif format.lower().startswith('bed6'):
                entry.addOption(name=words[3], score=float(words[4]), strand=words[5])
            elif format.lower().startswith('strand'):
                if len(words) >= 4: # properly formatted
                    entry.addOption(strand = words[3])
                if len(words) >= 5:
                    entry.addOption(name = words[4])
            elif format.lower().startswith('peak'):
                if len(words) >= 10: # narrowpeaks
                    entry.addOption(name = words[3], score = int(words[4]), strand = words[5], signalValue = float(words[6]), pValue = float(words[7]), qValue = float(words[8]), peak = int(words[9]))
                else: # broadpeaks
                    entry.addOption(name = words[3], score = int(words[4]), strand = words[5], signalValue = float(words[6]), pValue = float(words[7]), qValue = float(words[8]))
            elif format.lower().startswith('summit'):
                if len(words) >= 9:
                    entry.addOption(summit = int(words[4]), tags = int(words[5]), pValue = float(words[6]), fold = float(words[7]), fdr = float(words[8]))
                else:
                    entry.addOption(summit = int(words[4]), tags = int(words[5]), pValue = float(words[6]), fold = float(words[7]))
            elif format.lower().startswith('ccat'):
                entry.addOption(summit = int(words[1]) - entry.chromStart, tags = int(words[4]), bg = int(words[5]), zscore = float(words[6]), fdr = float(words[7]), name = '.', score = int(words[4]), strand = '.')
            elif format.lower().startswith('crop'):
                entry.addOption(score = int(words[2]), name = '.', strand = '.')
                entry.chromEnd = entry.chromStart + 1
            # check if the chromosome has been seen before
            tree = chroms.get(chrom)
            if not tree:
                tree = ival.IntervalTree()
                chroms[chrom] = tree
            # put the entry in the interval tree for the appropriate chromosome
            iv = ival.Interval(entry.chromStart, entry.chromEnd)
            tree.put(iv, entry)
        except RuntimeError as e:
            if not acceptHeaderRows:
                raise RuntimeError('Error in BED file at row %d (%s)' % (row, e.strerror))
            else:
                headerRow = words
                acceptHeaderRows -= 1 # count down the number of header rows that can occur
    f.close()
    return chroms

def writeBedFile(entries, filename, format = 'BED6', header = None):
    """ Save the BED entries to a BED file.
        format - the format to use for WRITING, currently only BED6 ('Optional' 6-col format) is supported.
    """
    f = open(filename, 'w')
    if header:
        f.write(header + '\n')
    for row in entries:
        if row.blocks: # BED12
            f.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t" % (row.chrom, row.chromStart, row.chromEnd, row.name, row.score, row.strand, row.thickStart, row.thickEnd))
            if row.itemRgb:
                if len(row.itemRgb) == 3:
                    f.write("%d,%d,%d\t" % (row.itemRgb[0],row.itemRgb[1],row.itemRgb[2]))
                else:
                    f.write("0\t")
            else:
                f.write("0\t")
            f.write("%d\t" % (len(row)))
            blockStarts = []
            blockSizes = []
            for b in row.blocks:
                blockStarts.append(b.ival.min - row.chromStart)
                blockSizes.append(len(b.ival))
            for b in blockSizes:
                f.write("%d," % (b))
            f.write("\t")
            for b in blockStarts:
                f.write("%d," % (b))
            f.write("\n")
        else:
            if format == 'Peaks':
                f.write("%s\t%d\t%d\t%s\t%d\t%s\t%f" % (row.chrom, row.chromStart, row.chromEnd, row.name, row.score, row.strand, row.signalValue))
            elif format == 'Limited':
                f.write("%s\t%d\t%d" % (row.chrom, row.chromStart, row.chromEnd))
            elif format == 'Strand':
                f.write("%s\t%d\t%d" % (row.chrom, row.chromStart, row.chromEnd, row.strand, row.name))
            else:
                f.write("%s\t%d\t%d\t%s\t%d\t%s" % (row.chrom, row.chromStart, row.chromEnd, row.name, row.score, row.strand))
            f.write("\n")
    f.close()

if __name__ == '__main__':
#    bf = BedFile('/Users/mikael/binfpy/BIOL3014/Week7/mm10_genes.bed', 'optional')
    bf = BedFile('/Volumes/Share/ARCDP19/Analysis/Young/Young_flat.bed', 'optional')
    print(bf.chroms.keys())
    g = bf.generate('chr1')
    print(next(g))
    print(next(g))
    print(next(g))
    cnt = 0
    collect = []
    for entry in bf:
        cnt += 1
        print(str(cnt) + '\t' + str(entry))
        collect.append(entry)
        if cnt == 7:
            for b in entry:
                print('\t', b)
        if cnt == 10:
            break
    writeBedFile(collect, '/Users/mikael/Desktop/test.bed')
    bf2 = BedFile('/Users/mikael/Desktop/test.bed', 'opt')
    q = ival.Interval(3805000, 3806000)
    t2 = ival.IntervalTree()
    t2.put(q, "blah")
    for entry in bf2:
        if entry.isBlockOverlap(q):
            print('Found:', entry)
            tree = entry.getBlocks()
            t2.putAll(tree)
            for t in t2:
                print(t)



    entry1 = BedEntry('chrX', 3858266, 3858530)
    print(entry1 in bf)
    entry2 = BedEntry('chrX', 10047550, 10067694)
    for x in bf.getOverlap(entry2):
        print(x)
    entry3 = BedEntry('chr9', 102699903, 102700167)
    for x in bf.getClosest(entry3):
        print(x)
        for y in x:
            print(y)
