'''
Module for classes and functions that are representing and processing basic probabilities.
Also includes Markov chain and hidden Markov model.
Uses and depends on "Alphabet" that is used to define discrete random variables.
'''
import random
from sym import *
import math

#################################################################################################
# Generic utility functions
#################################################################################################

def _getMeTuple(alphas, str):
    """ Handy function that resolves what entries that are being referred to in the case
    of written wildcards etc.
    Example y = _getMeTuple([DNA_Alphabet, Protein_Alphabet], '*R') gives y = (None, 'R')
    alphas: the alphabets
    str: the string that specifies entries (may include '*' and '-' signifying any symbol) """
    assert len(str) == len(alphas), "Entry invalid"
    if not type(str) is tuple:
        list = []
        for ndx in range(len(alphas)):
            if str[ndx] == '*' or str[ndx] == '-':
                list.append(None)
            else:
                list.append(str[ndx])
        return tuple(list)
    else:
        return str

#################################################################################################
# Distrib class
#################################################################################################

class Distrib():
    """ A class for a discrete probability distribution, defined over a specified "Alphabet"
        TODO: Fix pseudo counts
              Exclude from counts, specify in constructor,
              include only when computing probabilities by standard formula (n_a + pseudo_a * N^(1/2)) / (N + N^(1/2))
              Exclude from filesaves, include with filereads (optional)
    """
    def __init__(self, alpha, pseudo = 0.0):
        """ Construct a new distribution for a specified alphabet, using an optional pseudo-count.
        alpha: alphabet
        pseudo: either a single "count" that applies to all symbols, OR a distribution/dictionary with counts.
        """
        self.pseudo = pseudo or 0.0
        self.alpha = alpha
        self.cnt = [0.0 for _ in alpha]
        try: # assume pseudo is a dictionary or a Distrib itself
            self.tot = 0
            symndx = 0
            for sym in alpha:
                cnt = float(pseudo[sym])
                self.cnt[symndx] = cnt
                self.tot = self.tot + cnt
                symndx += 1
        except TypeError: # assume pseudo is a single count for each symbol
            self.cnt = [float(self.pseudo) for _ in alpha]
            self.tot = float(self.pseudo) * len(alpha) # track total counts (for efficiency)

    def observe(self, sym, cntme = 1.0):
        """ Make an observation of a symbol
        sym: symbol that is being observed
        cntme: number/weight of observation (default is 1)
        """
        ndx = self.alpha.symbols.index(sym)
        self.cnt[ndx] = self.cnt[ndx] + cntme
        self.tot = self.tot + cntme
        return

    def reset(self):
        """ Re-set the counts of this distribution. Pseudo-counts are re-applied. """
        try:
            self.tot = 0
            symndx = 0
            for sym in self.alpha: # assume it is a Distribution
                cnt = float(self.pseudo[sym])
                self.cnt[symndx] = cnt
                self.tot = self.tot + cnt
                symndx += 1
        except TypeError: # assume pseudo is a single count for each symbol
            self.cnt = [float(self.pseudo) for _ in self.alpha]
            self.tot = float(self.pseudo) * len(self.alpha) # track total counts (for efficiency)

    def reduce(self, new_alpha):
        """ Create new distribution from self, using (smaller) alphabet new_alpha. """
        d = Distrib(new_alpha, self.pseudo)
        for sym in new_alpha:
            d.observe(sym, self.cnt[self.alpha.index(sym)])
        return d

    def count(self, sym = None):
        """ Return the absolute count(s) of the distribution
            or the count for a specified symbol. """
        if sym != None:
            ndx = self.alpha.symbols.index(sym)
            return self.cnt[ndx]
        else:
            d = {}
            index = 0
            for a in self.alpha:
                d[a] = self.cnt[index]
                index += 1
            return d

    def add(self, distrib):
        """ Add the counts for the provided distribution to the present. """
        for i in range(len(self.cnt)):
            cnt = distrib.count(self.alpha[i])
            self.cnt[i] += cnt
            self.tot += cnt

    def subtract(self, distrib):
        """ Subtract the counts for the provided distribution from the present. """
        for i in range(len(self.cnt)):
            cnt = distrib.count(self.alpha[i])
            self.cnt[i] -= cnt
            self.tot -= cnt

    def getSymbols(self):
        return self.alpha.symbols

    def __getitem__(self, sym):
        """ Retrieve the probability of a symbol (ascertained by counts incl pseudo-counts) """
        if self.tot > 0.0:
            return self.count(sym) / self.tot
        else:
            return 1.0 / len(self.alpha) # uniform

    def prob(self, sym = None):
        """ Retrieve the probability of a symbol OR the probabilities of all symbols
        (listed in order of the alphabet index). """
        if sym != None:
            return self.__getitem__(sym)
        elif self.tot > 0:
            return [ s / self.tot for s in self.cnt ]
        else:
            return [ 1.0 / len(self.alpha) for _ in self.cnt ]

    def __iter__(self):
        return self.alpha

    def __str__(self):
        """ Return a readable representation of the distribution """
        str = '< '
        for s in self.alpha:
            str += (s + ("=%4.2f " % self[s]))
        return str + ' >'

    def swap(self, sym1, sym2):
        """ Swap the entries for specified symbols. Useful for reverse complement etc.
            Note that changes are made to the current instance. Use swapxcopy if you
            want to leave this instance intact. """
        sym1ndx = self.alpha.index(sym1)
        sym2ndx = self.alpha.index(sym2)
        tmpcnt = self.cnt[sym1ndx]
        self.cnt[sym1ndx] = self.cnt[sym2ndx]
        self.cnt[sym2ndx] = tmpcnt

    def swapxcopy(self, sym1, sym2):
        """ Create a new instance with swapped entries for specified symbols.
            Useful for reverse complement etc.
            Note that changes are NOT made to the current instance.
            Use swap if you want to modify this instance. """
        newdist = Distrib(self.alpha, self.count())
        newdist.swap(sym1, sym2)
        return newdist

    def writeDistrib(self, filename = None):
        """ Write the distribution to a file or string.
            Note that the total number of counts is also saved, e.g.
            * 1000 """
        str = ''
        for s in self.alpha:
            str += (s + ("\t%f\n" % self[s]))
        str += "*\t%d\n" % self.tot
        if filename != None:
            fh = open(filename, 'w')
            fh.write(str)
            fh.close()
        return str

    def generate(self):
        """ Generate and return a symbol from the distribution using assigned probabilities. """
        alpha = self.alpha
        p = random.random() # get a random value between 0 and 1
        q = 0.0
        for sym in alpha: # pick a symbol with a frequency proportional to its probability
            q = q + self[sym]
            if p < q:
                return sym
        return alpha[len(alpha)]

    def getmax(self):
        """ Generate the symbol with the largest probability. """
        maxprob = 0.0
        maxsym = None
        for sym in self.alpha:
            if self[sym] > maxprob or maxprob == 0.0:
                maxsym = sym
                maxprob = self[sym]
        return maxsym

    def getsort(self):
        """ Return the list of symbols, in order of their probability. """
        symlist = [sym for (sym, _) in self.getProbsort()]
        return symlist

    def getProbsort(self):
        """ Return the list of symbol-probability pairs, in order of their probability. """
        s = [(sym, self.prob(sym)) for sym in self.alpha]
        ss = sorted(s, key=lambda y: y[1], reverse=True)
        return ss

    def divergence(self, distrib2):
        """ Calculate the Kullback-Leibler divergence between two discrete distributions.
            Note that when self.prob(x) is 0, the divergence for x is 0.
            When distrib2.prob(x) is 0, it is replaced by 0.0001.
        """
        assert self.alpha == distrib2.alpha
        sum = 0.0
        base = len(self.alpha)
        for sym in self.alpha:
            if self[sym] > 0:
                if distrib2[sym] > 0:
                    sum += math.log(self[sym] / distrib2[sym]) * self[sym]
                else:
                    sum += math.log(self[sym] / 0.0001) * self[sym]
        return sum

    def entropy(self):
        """ Calculate the information (Shannon) entropy of the distribution.
            Note that the base is the size of the alphabet, so maximum entropy is by definition 1.
            Also note that if the probability is exactly zero, it is replaced by a small value to
            avoid numerical issues with the logarithm. """
        sum = 0.0
        base = len(self.alpha)
        for sym in self.alpha:
            p = self.__getitem__(sym)
            if p == 0:
                p = 0.000001
            sum +=  p * math.log(p, base)
        return -sum

def writeDistribs(distribs, filename):
    """ Write a list/set of distributions to a single file. """
    str = ''
    k = 0
    for d in distribs:
        str += "[%d]\n%s" % (k, d.writeDistrib())
        k += 1
    fh = open(filename, 'w')
    fh.write(str)
    fh.close()

def _readDistrib(linelist):
    """ Extract distribution from a pre-processed list if strings. """
    symstr = ''
    d = {}
    for line in linelist:
        line = line.strip()
        if len(line) == 0 or line.startswith('#'):
            continue
        sections = line.split()
        sym, value = sections[0:2]
        if len(sym) == 1:
            if sym != '*':
                symstr += sym
        else:
            raise RuntimeError("Invalid symbol in distribution: " + sym)
        try:
            d[sym] = float(value)
        except ValueError:
            raise RuntimeError("Invalid value in distribution for symbol " + sym + ": " + value)
    if len(d) == 0:
        return None
    alpha = Alphabet(symstr)
    if '*' in list(d.keys()): # tot provided
        for sym in d:
            if sym != '*':
                d[sym] = d[sym] * d['*']
    distrib = Distrib(alpha, d)
    return distrib

def readDistribs(filename):
    """ Load a list of distributions from file.
    Note that if a row contains '* <number>' then it is assumed that each probability
    associated with the specific distribution is based on <number> counts. """
    fh = open(filename)
    string = fh.read()
    distlist = []
    linelist = []
    for line in string.splitlines():
        line = line.strip()
        if line.startswith('['):
            if len(linelist) != 0:
                distlist.append(_readDistrib(linelist))
            linelist = []
        elif len(line) == 0 or line.startswith('#'):
            pass # comment or blank line --> ignore
        else:
            linelist.append(line)
    # end for-loop, reading the file
    if len(linelist) != 0:
        distlist.append(_readDistrib(linelist))
    fh.close()
    return distlist

def readDistrib(filename):
    """ Load a distribution from file.
    Note that if a row contains '* <number>' then it is assumed that each probability
    is based on <number> counts. """
    dlist = readDistribs(filename)
    if len(dlist) > 0:  # if at least one distribution was in the file...
        return dlist[0] # return the first

import re

def _readMultiCount(linelist, format = 'JASPAR'):
    ncol = 0
    symcount = {}
    if format == 'JASPAR2010':
        for line in linelist:
            line = line.strip()
            if len(line) > 0:
                name = line.split()[0]
                counts = []
                for txt in re.findall(r'\w+', line):
                    try:
                        y = float(txt)
                        counts.append(y)
                    except ValueError:
                        pass # ignore non-numeric entries
                if len(counts) != ncol and ncol != 0:
                    raise RuntimeError('Invalid row in file: ' + line)
                ncol = len(counts)
                if len(name) == 1: # proper symbol
                    symcount[name] = counts
        alpha = Alphabet(''.join(list(symcount.keys())))
        distribs = []
        for col in range(ncol):
            d = dict([(sym, symcount[sym][col]) for sym in symcount])
            distribs.append(Distrib(alpha, d))
    elif format == 'JASPAR':
        alpha_str = 'ACGT'
        alpha = Alphabet(alpha_str)
        cnt = 0
        for sym in alpha_str:
            line = linelist[cnt].strip()
            counts = []
            for txt in re.findall(r'\w+', line):
                try:
                    y = float(txt)
                    counts.append(y)
                except ValueError:
                    pass # ignore non-numeric entries
            if len(counts) != ncol and ncol != 0:
                raise RuntimeError('Invalid row in file: ' + line)
            ncol = len(counts)
            symcount[sym] = counts
            cnt += 1
        distribs = []
        for col in range(ncol):
            d = dict([(sym, symcount[sym][col]) for sym in symcount])
            distribs.append(Distrib(alpha, d))
    else:
        raise RuntimeError('Unsupported format: ' + format)
    return distribs

def readMultiCounts(filename, format = 'JASPAR'):
    """ Read a file of raw counts for multiple distributions over the same set of symbols
        for (possibly) multiple (named) entries.
        filename: name of file
        format: format of file, default is 'JASPAR' exemplified below
        >MA0001.1 SEP4
        0    3    79    40    66    48    65    11    65    0
        94    75    4    3    1    2    5    2    3    3
        1    0    3    4    1    0    5    3    28    88
        2    19    11    50    29    47    22    81    1    6
        returns a dictionary of Distrib's, key:ed by entry name (e.g. MA001.1)
    """
    fh = open(filename)
    linelist = []
    entryname = ''
    entries = {}
    for row in fh:
        row = row.strip()
        if len(row) < 1: continue
        if row.startswith('>'):
            if len(linelist) > 0:
                entries[entryname] = _readMultiCount(linelist, format=format)
                linelist = []
            entryname = row[1:].split()[0]
        else:
            linelist.append(row)
    if len(linelist) > 0:
        entries[entryname] = _readMultiCount(linelist, format=format)
    fh.close()
    return entries

def readMultiCount(filename, format = 'JASPAR'):
    """ Read a file of raw counts for multiple distributions over the same set of symbols.
        filename: name of file
        format: format of file, default is 'JASPAR' exemplified below
        0    3    79    40    66    48    65    11    65    0
        94    75    4    3    1    2    5    2    3    3
        1    0    3    4    1    0    5    3    28    88
        2    19    11    50    29    47    22    81    1    6
        returns a list of Distrib's
    """
    d = readMultiCounts(filename, format=format)
    if len(d) > 0:
        return list(d.values())[0]

#################################################################################################
# Joint class
#################################################################################################

class Joint(object):
    """ A joint probability class.
        The JP is represented as a distribution over n-tuples where n is the number of variables.
        Variables can be for any defined alphabet. The size of each alphabet determine the
        number of entries in the table (with probs that add up to 1.0) """

    def __init__(self, alphas):
        """ A distribution of n-tuples.
        alphas: Alphabet(s) over which the distribution is defined
        """
        if type(alphas) is Alphabet:
            self.alphas = tuple( [alphas] )
        elif type(alphas) is tuple:
            self.alphas = alphas
        else:
            self.alphas = tuple( alphas )
        self.store = TupleStore(self.alphas)
        self.totalCnt = 0

    def getN(self):
        """ Retrieve the number of distributions/random variables. """
        return len(self.alphas)

    def __iter__(self):
        return self.store.__iter__()

    def reset(self):
        """ Re-set the counts of this joint distribution. Pseudo-counts are re-applied. """
        for entry in self.store:
            self.store[entry] = None
        self.totalCnt = 0

    def observe(self, key, cnt = 1):
        """ Make an observation of a tuple/key
        key: tuple that is being observed
        cnt: number/weight of observation (default is 1)
        """
        key = _getMeTuple(self.alphas, key)
        if not None in key:
            score = self.store[key]
            if (score == None):
                score = 0
            self.totalCnt += cnt
            self.store[key] = score + cnt
        else: # there are wildcards in the key
            allkeys = [mykey for mykey in self.store.getAll(key)]
            mycnt = float(cnt)/float(len(allkeys))
            self.totalCnt += cnt
            for mykey in allkeys:
                score = self.store[mykey]
                if (score == None):
                    score = 0
                self.store[mykey] = score + mycnt
        return

    def count(self, key):
        """ Return the absolute count that is used for the joint probability table. """
        key = _getMeTuple(self.alphas, key)
        score = self.store[key]
        if (score == None):
            score = 0.0
            for match in self.store.getAll(key):
                y = self.store[match]
                if y != None:
                    score += y
        return score

    def __getitem__(self, key):
        """ Determine and return the probability of a specified expression of the n-tuple
        which can involve "wildcards"
        Note that no assumptions are made regarding independence. """
        key = _getMeTuple(self.alphas, key)
        score = self.store[key]
        if (score == None):
            score = 0.0
            for match in self.store.getAll(key):
                y = self.store[match]
                if y != None:
                    score += y
        if self.totalCnt == 0:
            return 0.0
        return float(score) / float(self.totalCnt)

    def __str__(self):
        """ Return a textual representation of the JP. """
        str = '< '
        if self.totalCnt == 0.0:
            return str + 'None >'
        for s in self.store:
            if self[s] == None:
                y = 0.0
            else:
                y = self[s]
            str += (''.join(s) + ("=%4.2f " % y))
        return str + ' >'

    def items(self, sort = False):
        """ In a dictionary-like way return all entries as a list of 2-tuples (key, prob).
        If sort is True, entries are sorted in descending order of probability.
        Note that this function should NOT be used for big (>5 variables) tables."""
        if self.totalCnt == 0.0:
            return []
        ret = []
        for s in self.store:
            if self[s] != None:
                ret.append((s, self[s]))
        if sort:
            return sorted(ret, key=lambda v: v[1], reverse=True)
        return ret


class IndepJoint(Joint):

    def __init__(self, alphas, pseudo = 0.0):
        """ A distribution of n-tuples.
        All positions are assumed to be independent.
        alphas: Alphabet(s) over which the distribution is defined
        """
        self.pseudo = pseudo
        if type(alphas) is Alphabet:
            self.alphas = tuple( [alphas] )
        elif type(alphas) is tuple:
            self.alphas = alphas
        else:
            self.alphas = tuple( alphas )
        self.store = [Distrib(alpha, pseudo) for alpha in self.alphas]

    def getN(self):
        """ Retrieve the number of distributions/random variables. """
        return len(self.alphas)

    def __iter__(self):
        return TupleStore(self.alphas).__iter__()

    def reset(self):
        """ Re-set the counts of each distribution. Pseudo-counts are re-applied. """
        self.store = [Distrib(alpha, self.pseudo) for alpha in self.alphas]

    def observe(self, key, cnt = 1, countGaps = True):
        """ Make an observation of a tuple/key
        key: tuple that is being observed
        cnt: number/weight of observation (default is 1)
        """
        assert len(key) == len(self.store), "Number of symbols must agree with the number of positions"
        for i in range(len(self.store)):
            subkey = key[i]
            if subkey == '-' and countGaps == False:
                continue
            if subkey == '*' or subkey == '-':
                for sym in self.alphas[i]:
                    score = self.store[i][sym]
                    if (score == None):
                        score = 0
                    self.store[i].observe(sym, float(cnt)/float(len(self.alphas[i])))
            else:
                score = self.store[i][subkey]
                if (score == None):
                    score = 0
                self.store[i].observe(subkey, cnt)

    def __getitem__(self, key):
        """ Determine and return the probability of a specified expression of the n-tuple
        which can involve "wildcards"
        Note that variables are assumed to be independent. """
        assert len(key) == len(self.store), "Number of symbols must agree with the number of positions"
        prob = 1.0
        for i in range(len(self.store)):
            mykey = key[i]
            if mykey == '*' or mykey == '-':
                pass # same as multiplying with 1.0 (all symbols possible)
            else:
                prob *= self.store[i][mykey]
        return prob

    def get(self, sym, pos):
        """ Retrieve the probability of a specific symbol at a specified position. """
        mystore = self.store[pos]
        return mystore[sym]

    def getColumn(self, column, count = False):
        """ Retrieve all the probabilities (or counts) for a specified position.
            Returns values as a dictionary, with symbol as key."""
        d = {}
        for a in self.alphas[column]:
            if count: # absolute count
                d[a] = self.store[column].count(a)
            else: # probability
                d[a] = self.store[column][a]
        return d

    def getRow(self, sym, count = False):
        """ Retrieve the probabilities (or counts) for a specific symbol over all columns/positions.
            Returns a list of values in the order of the variables/alphabets supplied to the constructor. """
        d = []
        for store in self.store:
            if count: # absolute count
                d.append(store.count(sym))
            else: # probability
                d.append(store[sym])
        return d

    def getMatrix(self, count = False):
        """ Retrieve the full matrix of probabilities (or counts) """
        d = {}
        for a in self.alphas[0]:
            d[a] = self.getRow(a, count)
        return d

    def displayMatrix(self, count = False):
        """ Pretty-print matrix """
        print((" \t%s" % (''.join("\t%5d" % (i + 1) for i in range(len(self.alphas))))))
        for a in self.alphas[0]:
            if count:
                print(("%s\t%s" % (a, ''.join("\t%5d" % (y) for y in self.getRow(a, True)))))
            else:
                print(("%s\t%s" % (a, ''.join("\t%5.3f" % (y) for y in self.getRow(a)))))

    def __str__(self):
        """ Text representation of the table. Note that size is an issue so big tables
        will not be retrieved and displayed. """
        if self.alphas > 5:
            return '< ... too large to process ... >'
        tstore = TupleStore(self.alphas)
        str = '< '
        for key in tstore:
            p = 1.0
            for i in range(len(self.store)):
                value = self.store[i][key[i]]
                if value != None and value != 0.0:
                    p *= value
                else:
                    p = 0;
                    break;
            str += (''.join(key) + ("=%4.2f " % p))
        return str + ' >'

    def items(self, sort = False):
        """ In a dictionary-like way return all entries as a list of 2-tuples (key, prob).
        If sort is True, entries are sorted in descending order of probability.
        Note that this function should NOT be used for big (>5 variables) tables."""
        tstore = TupleStore(self.alphas)
        ret = []
        for key in tstore:
            p = 1.0
            for i in range(len(self.store)):
                value = self.store[i][key[i]]
                if value != None and value != 0.0:
                    p *= value
                else:
                    p = 0;
                    break;
            if p > 0.0:
                ret.append((key, p))
        if sort:
            return sorted(ret, key=lambda v: v[1], reverse=True)
        return ret

#################################################################################################
# Naive Bayes' classifier
#################################################################################################

class NaiveBayes():
    """ NaiveBayes implements a classifier: a model defined over a class variable
        and conditional on a list of discrete feature variables.
        Note that feature variables are assumed to be independent. """

    def __init__(self, inputs, output, pseudo_input = 0.0, pseudo_output = 0.0):
        """ Initialise a classifier.
            inputs: list of alphabets that define the values that input variables can take.
            output: alphabet that defines the possible values the output variable takes
            pseudo_input: pseudo-count used for each input variable (default is 0.0)
            pseudo_output: pseudo-count used for the output variable (default is 0.0) """
        if type(inputs) is Alphabet:
            self.inputs = tuple( [inputs] )
        elif type(inputs) is tuple:
            self.inputs = inputs
        else:
            self.inputs = tuple( inputs )
        self.condprobs = {}   # store conditional probabilities as a dictionary (class is key)
        for outsym in output: # GIVEN the class
            # for each input variable initialise a conditional probability
            self.condprobs[outsym] = [ Distrib(input, pseudo_input) for input in self.inputs ]
        self.classprob = Distrib(output, pseudo_output) # the class prior

    def observe(self, inpseq, outsym):
        """ Record an observation of an input sequence of feature values that belongs to a class.
            inpseq: sequence/list of feature values, e.g. 'ATG'
            outsym: the class assigned to these feature values. """
        condprob = self.condprobs[outsym]
        for i in range(len(inpseq)):
            condprob[i].observe(inpseq[i])
        self.classprob.observe(outsym)

    def __getitem__(self, key):
        """ Determine and return the class probability GIVEN a specified n-tuple of feature values
        The class probability is given as an instance of Distrib. """
        out = Distrib(self.classprob.alpha)
        for outsym in self.classprob.getSymbols():
            condprob = self.condprobs[outsym]
            prob = self.classprob[outsym]
            for i in range(len(key)):
                prob *= condprob[i][key[i]] or 0.0
            out.observe(outsym, prob)
        return out

#################################################################################################
# Markov chain
#################################################################################################

class MarkovChain():
    """ Markov Chain in a simple form; supports higher-orders and can determine (joint) probability of sequence.
    """

    def __init__(self, alpha, order = 1, startsym = '^', endsym = '$'):
        """ Construct a new Markov chain based on a given alphabet of characters.
        alpha: alphabet of allowed characters and states
        order: the number of states to include in memory (default is 1)
        startsym: the symbol to mark the first character in the internal sequence, and the first state
        endsym: the symbol to mark the termination of the internal sequence (and the last state)
        """
        self.order = order
        self.startsym = startsym
        self.endsym = endsym
        self.alpha = getTerminatedAlphabet(alpha, self.startsym, self.endsym)
        self.transit = TupleStore([self.alpha for _ in range(order)]) # transition probs, i.e. given key (prev state/s) what is the prob of current state

    def _getpairs(self, term_seq):
        """ Return a tuple of all (tuple) Markov pairs from a sequence. Used internally. """
        ret = []
        for i in range(len(term_seq) - self.order):
            past = tuple(term_seq[i:i + self.order])
            present = term_seq[i + self.order]
            ret.append(tuple([past, present]))
        return tuple(ret)

    def observe(self, wholeseq):
        """ Set parameters of Markov chain by counting transitions, as observed in the sequence.
         wholeseq: the sequence not including the termination symbols.
         """
        myseq = _terminate(wholeseq, self.order, self.startsym, self.endsym)
        for (past, present) in self._getpairs(myseq):
            d = self.transit[past]
            if not d: # no distrib
                d = Distrib(self.alpha)
                self.transit[past] = d
            d.observe(present)

    def __getitem__(self, wholeseq):
        """ Determine the log probability of a given sequence.
         wholeseq: the sequence not including the termination symbols.
         returns the joint probability
        """
        myseq = _terminate(wholeseq, self.order, self.startsym, self.endsym)
        logp = 0
        for (past, present) in self._getpairs(myseq):
            d = self.transit[past]
            if not d:
                return None
            p = d[present]
            if p == 0:
                return None
            logp += math.log(p)
        return logp


def _terminate(unterm_seq, order = 1, startsym = '^', endsym = '$'):
    """ Terminate sequence with start and end symbols """
    term_seq = [startsym for _ in range(order)]
    term_seq.extend(unterm_seq)
    term_seq.append(endsym)
    return term_seq

def getTerminatedAlphabet(alpha, startsym = '^', endsym = '$'):
    """ Amend the given alphabet with termination symbols """
    return Alphabet(alpha.symbols + tuple([startsym, endsym]))

#################################################################################################
# Hidden Markov model (HMM)
#################################################################################################

class HMM():
    """ Basic, first-order HMM.
     Has functionality to set up HMM, and query it with Viterbi and Forward algorithms."""

    def __init__(self, states, symbols, startstate = '^', endstate = '$'):
        """ Construct HMM with states and symbols, here given as strings of characters.
        > cpg_hmm = prob.HMM('HL','ACGT')
        """
        if isinstance(states, str):
            states = Alphabet(states)
        self.mystates = getTerminatedAlphabet(states, startstate, endstate)
        if isinstance(symbols, str):
            symbols = Alphabet(symbols)
        self.mysymbols = getTerminatedAlphabet(symbols, startstate, endstate)
        self.a = dict() # transition probabilities
        self.e = dict() # emission probabilities
        self.startsym = startstate
        self.endsym = endstate

    def transition(self, fromstate, distrib):
        """ Add a transition to the HMM, determining with the probability of transitions, e.g.
        > cpg_hmm.transition('^',{'^':0,'$':0,'H':0.5,'L':0.5})
        > cpg_hmm.transition('H',{'^':0,'$':0.001,'H':0.5,'L':0.5})
        > cpg_hmm.transition('L',{'^':0,'$':0.001,'H':0.4,'L':0.6})
        > cpg_hmm.transition('$',{'^':1,'$':0,'H':0,'L':0})
        """
        if not isinstance(distrib, Distrib):
            distrib = Distrib(self.mystates, distrib)
        self.a[fromstate] = distrib

    def emission(self, state, distrib):
        """ Add an emission probability to the HMM, e.g.
        > cpg_hmm.emission('^',{'^':1,'$':0,'A':0,'C':0,'G':0,'T':0})
        > cpg_hmm.emission('H',{'^':0,'$':0,'A':0.2,'C':0.3,'G':0.3,'T':0.2})
        > cpg_hmm.emission('L',{'^':0,'$':0,'A':0.3,'C':0.2,'G':0.2,'T':0.3})
        > cpg_hmm.emission('$',{'^':0,'$':1,'A':0,'C':0,'G':0,'T':0})
        """
        if not isinstance(distrib, Distrib):
            distrib = Distrib(self.mysymbols, distrib)
        self.e[state] = distrib

    def joint(self, symseq, stateseq):
        """
        Determine the joint probability of the sequence and the given path.
        :param symseq: sequence of characters
        :param stateseq: sequence of states
        :return: the probability
        """
        X = _terminate(symseq, 1, self.startsym, self.endsym)
        P = _terminate(stateseq, 1, self.startsym, self.endsym)
        p = 1
        for i in range(len(X) - 1):
            p = p * self.e[P[i]][X[i]] * self.a[P[i]][P[i + 1]]
        return p

    def viterbi(self, symseq, V = dict(), trace = dict()):
        """
        Determine the Viterbi path (the most probable sequence of states) given a sequence of symbols
        :param symseq: sequence of symbols
        :param V: the Viterbi dynamic programming variable as a matrix (optional; pass an empty dict if you need it)
        :param trace: the traceback (optional; pass an empty dict if you need it)
        :return: the Viterbi path as a string of characters
        > X = 'GGCACTGAA' # sequence of characters
        > states = cpg_hmm.viterbi(X)
        > print(states)
        """
        X = _terminate(symseq, 1, self.startsym, self.endsym) # put start and end symbols on sequence
        # Initialise state scores for each index in X
        for state in self.mystates:
            # Fill in emission probabilities in V for each index of X
            # (only the first position is really needed)
            V[state] = [self.e[state][x] for x in X]
            trace[state] = []
        # Next loop through the sequence
        for j in range(len(X) - 1):
            i = j + 1  # sequence index that we're processing start with 1, not 0
            for tostate in self.mystates: # check each state for i = 1, ...
                tracemax = 0
                beststate = None # the state v with max[Vv(i-1) * t(v,u)]
                for fromstate in self.mystates:
                    # determine the best score propagated forward from previous state
                    score = V[fromstate][i - 1] * self.a[fromstate][tostate]
                    if score > tracemax:
                        beststate = fromstate
                        tracemax = score
                # record the transition that will appear in the traceback
                trace[tostate].append(beststate)
                # finalise the dynamic programming score for current i in state u
                V[tostate][i] = self.e[tostate][X[i]] * tracemax
        # finally, assemble the string that describes the most probable path
        ret = ''
        traced = '$'
        for j in range(len(X)):
            i = len(X) - 2 - j
            traced = trace[traced][i]
            if j > 0 and j < len(X) - 2:
                ret = traced + ret
        return ret

    def forward(self, symseq, F = dict()):
        """
        Determine the probability of the sequence, summing over all possible state paths
        :param symseq: sequence of symbols
        :param F: the Forward dynamic programming variable as a matrix (optional; pass an empty dict if you need it)
        :return: the probability
        > X = 'GGCACTGAA' # sequence of characters
        > prob = cpg_hmm.forward(X)
        > print(prob)
        """
        X = _terminate(symseq, 1, self.startsym, self.endsym)
        # Initialise state scores for each index in X
        for state in self.mystates:
            # Fill in emission probabilities for each index in X
            F[state] = [self.e[state][x] for x in X]
        for j in range(len(X) - 1):
            i = j + 1  # sequence index that we're processing
            for tostate in self.mystates:
                mysum = 0
                for fromstate in self.mystates:
                    mysum += F[fromstate][i - 1] * self.a[fromstate][tostate]
                F[tostate][i] = self.e[tostate][X[i]] * mysum
        traced = '$'
        return F[traced][len(X) - 1]

    def writeHTML(self, X, Viterbi, Trace = None, filename = None):
        """ Generate HTML that displays a DP matrix from Viterbi (or Forward) algorithms.
        > from IPython.core.display import HTML
        > X = 'GGCACTGAA' # sequence of characters
        > V = dict()
        > T = dict()
        > cpg_hmm.viterbi(X, V, T)
        > HTML(cpg_hmm.writeHTML(X, V, T))
        """
        html = '''<html><head><meta content="text/html; charset=ISO-8859-1" http-equiv="Content-Type">\n<title>HMM dynamic programming matrix</title>\n</head><body><pre>\n'''
        html += '<table style=\"width:100\%\">\n'
        html += '<tr>\n'
        html += '<th>X</th>\n'
        for state in Viterbi:
            html += '<th>%s</th>\n' % str(state)
        html += '</tr>\n'
        # process each sequence symbol
        X = _terminate(X, 1, self.startsym, self.endsym)
        for row in range(len(X)):
            html += '<tr>\n'
            html += '<td>x<sub>%d</sub>=<pre>%s</td>\n' % (row, str(X[row]))
            for state in Viterbi:
                if Trace and row > 0:
                    html += '<td>e<sub>%s</sub>(%s)=%4.2f<pre>V<sub>%s</sub>=%3.1e<pre>&uarr;%s[%4.2f]</td>\n' % (str(state),X[row],self.e[state][X[row]],str(state),Viterbi[state][row],Trace[state][row - 1],self.a[Trace[state][row - 1]][state] if Trace[state][row - 1] != None else 0)
                else:
                    html += '<td>e<sub>%s</sub>(%s)=%4.2f<pre>V<sub>%s</sub>=%3.1e</td>\n' % (str(state),X[row],self.e[state][X[row]],str(state),Viterbi[state][row])
            html += '</tr>\n'
        html += '</table>\n'
        html += '</pre></body></html>'
        if filename:
            fh = open(filename, 'w')
            fh.write(html)
            fh.close()
        return html

#################################################################################################
# The univariate Gaussian density function.
#################################################################################################

class Gaussian():

    mu = None           # mean
    sigma = None        # standard deviation
    sigmaSquared = None # variance

    def ROOT_2PI(self):
        return (2 * math.pi) ** 2
    def LOG_ROOT_2PI(self):
        return 0.5 * (math.log(2) + math.log(math.pi))

    def __init__(self, mean, variance):
        """ Creates a univariate Gaussian distribution with the given fixed mean and variance. """
        self.mu = mean
        self.sigmaSquared = variance;
        self.sigma = math.sqrt(variance);
        self.normConst = self.sigma * self.ROOT_2PI();
        self.logNormConst = (0.5 * math.log(variance)) + self.LOG_ROOT_2PI();

    def __str__(self):
        return "<" + "%5.3f" % self.mu + u"\u00B1%5.3f" % self.sigma + ">"

    def getDensity(self, x):
        """ Returns the density of this Gaussian distribution at the given value. """
        return (math.exp(-math.pow((x - self.mu), 2) / (2 * self.sigmaSquared)) / self.normConst);

    def __getitem__(self, value):
        """ Get the probability density of a specified value for this Gaussian """
        return self.getDensity(value)

    def getMean(self):
        return self.mu

    def getVariance(self):
        return self.sigmaSquared

    def getLogDensity(self, x):
        """ Returns the natural log of the density of this Gaussian distribution at the given value. """
        return (-math.pow((x - self.mu), 2) / (2 * self.sigmaSquared)) - self.logNormConst;

    def sample(self):
        """ Returns a value sampled from this Gaussian distribution. The implementation uses the Box - Muller transformation
         [G.E.P.Box and M.E.Muller(1958) "A note on the generation of random normal deviates". Ann.Math.Stat 29: 610 - 611]. """
        U = random.random() # get a random value between 0 and 1
        V = random.random() # get a random value between 0 and 1
        return (self.mu + (self.sigma * math.sin(2 * math.pi * V) * math.sqrt((-2 * math.log(U)))))

def estimate(samples, count = None):
    """ Create a density based on the specified samples. Optionally, provide an iterable with the corresponding counts/weights. """
    mean = 0
    if count == None:
        for i in range(len(samples)):
            mean += samples[i] / len(samples)
        diff = 0
        for i in range(len(samples)):
            diff += (mean - samples[i]) * (mean - samples[i]);
        if (diff == 0):
            return None
        return Gaussian(mean, diff / len(samples))
    elif len(count) == len(samples):
        totcnt = 0
        for i in range(len(samples)):
            mean += samples[i] * count[i]
            totcnt += count[i]
        mean /= totcnt
        diff = 0
        for i in range(len(samples)):
            diff += (mean - samples[i]) * (mean - samples[i]) * count[i];
        if (diff == 0):
            return None
        return Gaussian(mean, diff / totcnt)

#################################################################################################
# The univariate Poisson distribution
#################################################################################################

class Poisson():

    def __init__(self, LAMBDA):
        """
           * Define a Poisson distribution
           * @ param lambda the average number of events per interval
        """
        self.LAMBDA = LAMBDA

    def p(self, k):
        """
         * The probability mass function
         * @param k the number of events in the interval
         * @return the probability of k events
        """
        return math.exp(k * math.log(self.LAMBDA) - self.LAMBDA - lgamma(k + 1))

    def cdf(self, k):
        """
         * The cumulative probability function.
         * The implementation calls the PMF for all values i from 0 to floor(k)
         * @param k the number of events in the interval
         * @return the cumulative probability of k events
         * https://en.wikipedia.org/wiki/Poisson_distribution
        """
        sum = 0
        for i in range(k + 1):
            sum += self.p(i)
            if (sum >= 1.0): # happens only due to poor numerical precision
                return 1.0
        return sum


def lgamma(x):
    """
     * Returns an approximation of the log of the Gamma function of x. Laczos
     * Approximation Reference: Numerical Recipes in C
     * http://www.library.cornell.edu/nr/cbookcpdf.html
    """
    cof = [ 76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 ]
    y = x
    tmp = x + 5.5
    tmp -= ((x + 0.5) * math.log(tmp))
    ser = 1.000000000190015
    for j in range(len(cof)):
        y += 1
        ser += (cof[j] / y)
    return (-tmp + math.log(2.5066282746310005 * ser / x))
