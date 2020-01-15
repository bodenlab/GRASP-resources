import shlex
import ival

class GtfEntry():

    '''
    GFF fields:
    seqname - The name of the sequence. Must be a chromosome or scaffold.
    source - The program that generated this feature.
    feature - The name of this type of feature. Some examples of standard feature types are "CDS" "start_codon" "stop_codon" and "exon"li>
    start - The starting position of the feature in the sequence. The first base is numbered 1.
    end - The ending position of the feature (inclusive).
    score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ":.":.
    strand - Valid entries include "+", "-", or "." (for don't know/don't care).
    frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be ".".
    group - All lines with the same group are linked together into a single item.
    '''
    def __init__(self, chrom, start, end, feature, score = ".", source = "unknown", strand = ".", frame = ".", group = None):
        self.seqname = chrom
        self.start = start
        self.end = end
        self.feature = feature
        self.score = score
        self.strand = strand
        self.source = source
        self.frame = frame
        self.group = group
        self.attr = {}
        if self.group:
            fields = self.group.split(';')
            for f in fields:
                pair = shlex.split(f.strip())
                if len(pair) == 2:
                    self.attr[pair[0]] = pair[1]

    def __getitem__(self, item):
        return self.attr[item]

    def __contains__(self, item):
        return item in self.attr

    def __str__(self):
        return str((self.seqname, self.start, self.end))

    def __len__(self):
        return self.end - self.start

    def getInterval(self):
        return ival.Interval(self.start, self.end)

def dist(entry1, entry2, signed = False, centre2centre = False):
    """ Calculate and return the BedEntry with the closest distance (from one end of the interval of this to the end of the interval of that).
        If centre2centre is True, use the centre-to-centre distance instead.
        If signed is True, the distance is negative if this interval is after the that.
    """
    if isinstance(entry1, GtfEntry) and isinstance(entry2, GtfEntry):
        if (entry1.seqname == entry2.seqname):
            return ival.dist(entry1.getInterval(), entry2.getInterval(), signed, centre2centre)
    return None

class GtfFile:
    """ Read GTF/GFF file.

        See http://genome.ucsc.edu/FAQ/FAQformat#format1
    """

    def __init__(self, entries, filter_feature = None):
        """
        Create a GtfFile instance.
        :param entries: an iterable of entries or a filename
        """
        if isinstance(entries, str): # filename
            self.chroms = readGtfFile(entries, filter_feature)
        else:
            self.chroms = dict()
            for entry in entries:
                # check if the chromosome has been seen before
                tree = self.chroms.get(entry.chrom)
                if not tree:
                    tree = ival.IntervalTree()
                    self.chroms[entry.chrom] = tree
                # put the entry in the interval tree for the appropriate chromosome
                iv = ival.Interval(entry.start, entry.end)
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
        if isinstance(item, GtfEntry):
            tree = self.chroms.get(item.chrom)
            if tree == None: return False
            else: return ival.Interval(item.start, item.end) in tree
        else:
            return False

    def getOverlap(self, item):
        if isinstance(item, GtfEntry):
            tree = self.chroms.get(item.chrom)
            if tree == None: return None
            else:
                iv = ival.Interval(item.start, item.end)
                res = tree.isectall(iv)
                ret = []
                for r in res:
                    ret.extend(r.values)
                return ret
        else: return None

    def getClosest(self, item):
        if isinstance(item, GtfEntry):
            tree = self.chroms.get(item.chrom)
            if tree == None: return None
            else:
                iv = ival.Interval(item.start, item.end)
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

    def getIndex(self, group_attribute):
        """
        Create a dictionary for a specified attribute in the group list, e.g. "gene_id", "gene_name" or "transcript_id"
        :param group_attribute:
        :return: a dictionary keyed by the values of the nominated group_attribute,
        pointing to the entry in the chromosome interval-tree
        """
        d = {}
        for entry in self:
            if group_attribute in entry.attr: # is the entry indexed with the attribute
                if not entry.attr[group_attribute] in d: # only add the first entry for this particular value
                    d[entry.attr[group_attribute]] = entry
        return d

def readGtfFile(filename, filter_feature = None):
    """ Read a GTF/GFF file.
    filename: name of file
    filter_feature: name of feature to be selected, all others ignored; None means anything
    """
    f = open(filename)
    row = 0
    acceptHeaderRows = 1
    headerRow = None
    chroms = dict()
    for line in f:
        row += 1
        words = line.strip().split('\t')
        if len(words) == 0:
            continue # ignore empty lines
        if words[0].strip().startswith('#'):
            continue # comment
        if words[0].strip().startswith('browser'):
            continue # ignore
        if words[0].strip().startswith('track'):
            continue # ignore
        try:
            seqname = words[0]
            source = words[1]
            feature = words[2]
            if filter_feature and filter_feature != feature:
                continue
            start = int(words[3])
            end = int(words[4])
            score = None
            if words[5].isnumeric():
                score = int(words[5])
            strand = '.'
            if words[6] == '+' or words[6] == '-':
                strand = words[6]
            frame = None
            if words[7].isdigit():
                frame = int(words[7])
            group = None
            if len(words) > 8:
                group = words[8]
            entry = GtfEntry(seqname, start, end, feature, score, source, strand, frame, group)
            # check if the chromosome has been seen before
            tree = chroms.get(seqname)
            if not tree:
                tree = ival.IntervalTree()
                chroms[seqname] = tree
            # put the entry in the interval tree for the appropriate chromosome
            iv = ival.Interval(entry.start, entry.end)
            tree.put(iv, entry)
        except RuntimeError as e:
            if not acceptHeaderRows:
                raise RuntimeError('Error in GTF/GFF file at row %d (%s)' % (row, e.strerror))
            else:
                headerRow = words
                acceptHeaderRows -= 1 # count down the number of header rows that can occur
    f.close()
    return chroms

def writeGtfFile(entries, filename, header = None):
    """ Save the GTF entries to a file.
    """
    f = open(filename, 'w')
    if header:
        f.write(header + '\n')
    for row in entries:
        f.write("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s" % (row.chrom, row.source, row.feature, row.start, row.end, row.score, row.strand, row.frame, row.group))
        f.write("\n")
    f.close()

if __name__ == '__main__':
    bf = GtfFile('/Users/mikael/simhome/NFIX/WT1689.gtf')
    print(bf.chroms.keys())
    g = bf.generate('chr12')
    print(next(g))
    print(next(g))
    print(next(g))
    cnt = 0
    for entry in bf:
        cnt += 1
        print(str(cnt) + '\t' + str(entry))
        if cnt == 100:
            break
