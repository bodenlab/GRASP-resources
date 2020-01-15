'''
Module with methods and classes for phylogeny.
Extended to handle n-ary trees (Jan 2019).
@author: mikael
'''
import sequence
from collections import defaultdict
import annotations

class PhyloTree:
    """ Rooted, n-ary tree for representing phylogenetic relationships.
        Functionality includes labelling and traversing nodes; reading and writing to Newick format;
        association with sequence alignment; maximum parsimony inference of ancestral sequence;
        generation of rooted tree by UPGMA.
        Known issues: Parsimony does not handle gaps in alignment.
        Programmers should note that almost all functionality is implemented through recursion. """

    def __init__(self, root):
        """ Create a tree from a node that is "root" in the tree."""
        self.root = root

    def putAlignment(self, aln):
        """ Associate the tree with a set of sequences/alignment.
            Involves assigning the sequence to the leaf nodes. """
        self.aln = aln
        self.root._assignAlignment(aln)

    def putAnnotations(self, nexus_annotations):
        self.nexus_annotations = nexus_annotations
        # Update the annotations dictionary so that it contains PhyloNode objects as keys, not text labels
        for node in self.getNodes():
            if node.label in self.nexus_annotations.leaf_annotations:
                self.nexus_annotations.leaf_annotations[node] = self.nexus_annotations.leaf_annotations[node.label]
                self.nexus_annotations.leaf_annotations.pop(node.label)

    def __str__(self):
        """ Produce a printable representation of the tree, specifically the root of the tree. """
        return str(self.root)

    def strSequences(self, start=None, end=None):
        """ Produce a sequence representation of the tree, specifically the root of the tree.
            Specify the start and end positions in the alignment for the sequence to be printed
            (if None the min and max positions will be used). """
        if self.aln != None:
            my_start = start or 0
            my_end = end or self.aln.alignlen
            return self.root._printSequences(my_start, my_end)

    def findLabel(self, label):
        """ Retrieve/return the node with the specified label.
            Returns None if not found."""
        return self.root._findLabel(label)

    def getNodes(self, strategy = 'DEPTH-FIRST'):
        """ Returns all nodes as a list """
        nodes = []
        queue = [self.root]
        while len(queue) > 0:
            node = queue.pop()
            nodes.append(node)
            # if strategy.upper().startswith('DEPTH'):
            if not node.isLeaf():
                queue.extend(node.children)
        return nodes

    def getLeaves(self):
        all = self.getNodes()
        leaves = []
        for n in all:
            if n.isLeaf():
                leaves.append(n)
        return leaves

    def getDescendantsOf(self, node, transitive=False):
        """ Retrieve and return the (list of) descendants (children) of a specified node.
            Node can be the label or the instance.
            transitive indicates if only the direct descendants (False) or if all descendants
            should be returned.
            If node does not exist, None is returned.
            If node has no descendants, an empty list will be returned."""
        if not isinstance(node, PhyloNode):
            node = self.findLabel(node)
        if node:
            return node.getDescendants(transitive)
        return None

    def getAncestorsOf(self, node, transitive=False):
        """ Retrieve and return the ancestor (transitive=False) or
            ancestors (transitive=True) of a specified node.
            Node can be the label or the instance.
            If node does not exist, None is returned.
            If node is the root of the tree, None is returned."""
        if not isinstance(node, PhyloNode):
            node = self.findLabel(node)
        if node:
            return node.getAncestors(transitive)

    def parsimony(self):
        """ Solve the "small parsimony problem",
            i.e. find the sequences on each of the internal nodes.
            See Jones and Pevzner, p. 368 and onwards, for details. """
        self.root._forwardParsimony(self.aln)  # setup and compute scores for all nodes
        self.root._backwardParsimony(self.aln)  # use scores to determine sequences
        return self.root.getSequence()  # return the sequence found at the root

    def swap_annotations(self, annotation_key):
        try:
            for node in self.getNodes():
                if node.isLeaf():
                    node.label = self.nexus_annotations.leaf_annotations[node][annotation_key]
        except:
            return

    def write_to_nexus(self, out_path, write_annotations=True, nexus_annotations=None, exclude_annotations=[], use_symbols=False):
        """
        Writes out the tree to NEXUS format, with any annotations stored in nexus_annotations added to the file
        :param out_path: The path to write the NEXUS file to
        :param nexus_annotations: The NexusAnnotations containing the annotations
        """
        if write_annotations and not nexus_annotations:
            if not self.nexus_annotations:
                raise RuntimeError("This tree file has no associated annotation file. Either associate or supply one as a parameter.")
            nexus_annotations = self.nexus_annotations
        if nexus_annotations:
            for node in self.getNodes():
                if node in self.nexus_annotations.node_annotations:
                    node.annotate_node(self.nexus_annotations.node_annotations, self.nexus_annotations.annotation_symbols, exclude_annotations, use_symbols)
            tree_annotation = str(self) + ";"
            self.swap_annotations("Original")
            for node in self.getNodes():
                if node in self.nexus_annotations.leaf_annotations:
                    node.annotate_node(self.nexus_annotations.leaf_annotations, exclude_annotations)
                leaves = []
                for node in self.getNodes():
                    if node.isLeaf():
                        leaves.append(node.label)
                leaf_annotation = ""
                for leaf in leaves:
                    leaf_annotation += "\n\t%s" % (leaf)
            with open(out_path, "w+") as file:
                file.write(
                    "#NEXUS\nbegin taxa;\n\tdimensions ntax=%d;\n\ttaxlabels%s\n;\nend;\n\nbegin trees;\n\ttree tree_1 = "
                    "[&R] %s\nend;" % (len(leaves), leaf_annotation, tree_annotation))

class PhyloNode:
    """ A class for a node in a rooted, n-ary tree.
        Contains pointers to multiple descendants/daughters,
        optional fields include data, label, sequence and dist.
        If parsimony is used scores and traceback pointers are available.
        A number of methods are named with a _ prefix. These can be, but
        are not intended to be used from outside the class. """

    def __init__(self, parent = None, label=''):
        """ Initialise a node.
            Set its parent (another PhyloNode), parent can be None.
            Set label to name it.
            Use field data for any type of information associated with node.
            Use dist to indicate the distance to its parent (if any).
            Other fields are used internally, including sequence for associated alignment,
            seqscores, back for maximum parsimony. """
        self.parent = parent
        self.children = None
        self.data = None
        self.label = label
        self.dist = None
        self.sequence = None  # The sequence after an alignment have been mapped (leaf) or the most parsimonous sequence (ancestral)
        self.seqscores = None # The scores propagated from leaves via children
        self.backptr = None   # Pointers back to children: what symbol rendered current/parent symbols

    def isLeaf(self):
        return self.nChildren() == 0

    def nChildren(self):
        if self.children == None:
            return 0
        else:
            return len(self.children)

    def __str__(self):
        """ Returns string with node (incl descendants) in a Newick style. """
        stubs = ['' for _ in range(self.nChildren())]
        label = dist = ''
        for i in range(self.nChildren()):
            stubs[i] = str(self.children[i])
        if self.dist or self.dist == 0.0:
            dist = ':' + str(self.dist)
        if self.label != None:
            label = str(self.label)
        if self.nChildren() == 0:
            return label + dist
        else:
            stubstr = '('
            for i in range(len(stubs) - 1):
                stubstr += stubs[i] + ','
            return stubstr + stubs[-1] + ')' + label + dist
        # there is no label
            '''
            if not self.left and self.right:
                return ',' + right
            elif self.left and not self.right:
                return left + ','
            elif self.left and self.right:
                return '(' + left + ',' + right + ')' + dist
                '''

    # def __le__(self, other):
    #     """ Returns indication of less than other node. """
    #     return other and self.__hash__() <= other.__hash__()
    #
    # def __eq__(self, other):
    #     """ Returns indication of equivalence to other node. """
    #     return other and self.__hash__() == other.__hash__()
    #
    # def __hash__(self):
    #     """ Returns hash of object. """
    #     return hash((self.label, self.dist, self.sequence))

    def _printSequences(self, start, end):
        """ Returns string with node (incl descendants) in a Newick style. """
        stubs = ['' for _ in range(self.nChildren())]
        label = dist = ''
        for i in range(self.nChildren()):
            stubs[i] = self._printSequences(self.children[i], start, end)
        if self.dist or self.dist == 0.0:
            dist = ':' + str(self.dist)
        if self.label != None:
            label = str(self.label)
        if self.nChildren() == 0:
            return label + dist
        else:
            stubstr = '('
            for i in range(len(stubs) - 1):
                stubstr += stubs[i] + ','
            return stubstr + stubs[-1] + ')' + label + dist

    def _findLabel(self, label):
        """ Find a node by label at this node or in any descendants (recursively). """
        if self.label == label:
            return self
        else:
            for i in range(self.nChildren()):
                found = self.children[i]._findLabel(label)
                if found:
                    return found
            return None

    def _propagateDistance(self, parent_dist):
        """ Convert absolute distances to relative.
            The only parameter is the absolute distance to the parent of this node. """
        travelled = self.dist  # absolute distance to this node
        self.dist = parent_dist - self.dist  # relative distance to this node
        for i in range(self.nChildren()):
            self.children[i]._propagateDistance(travelled)  # pass absolute distance to this node

    def _assignAlignment(self, aln):
        """ Assign an alignment to the node, which implies assigning a sequence to it if one is
            available in the alignment. """
        self.sequence = None
        for i in range(self.nChildren()):
            self.children[i]._assignAlignment(aln)
        for seq in aln.seqs:
            if seq.name == self.label:
                self.sequence = seq
                break

    """ # Not sure if this is required (putting nodes into a canonical ordering)
    def _canonise(self):
        if self.left == None and self.right == None:  # at leaf
            return self.label
        myleft = self.left._canonise()
        myright = self.right._canonise();
        if myleft > myright:
            tmpnode = self.left
            self.left = self.right
            self.right = tmpnode
            return myright
        return myleft
    """

    def _forwardParsimony(self, aln):
        """ Internal function that operates recursively to first initialise each node (forward),
            stopping only once a sequence has been assigned to the node,
            then to propagate scores from sequence assigned nodes to root (backward). """
        if self.sequence == None:  # no sequence has been assigned
            if self.nChildren() == 0:  # no children, so terminal, cannot propagate scores
                raise RuntimeError("No sequence assigned to leaf node:", self.label)
            scores = [None for _ in range(self.nChildren())]
            for i in range(self.nChildren()):
                scores[i] = self.children[i]._forwardParsimony(aln)
            # for each position in the alignment,
            # introduce (initially zero) score for each symbol in alphabet
            self.seqscores = [[0 for _ in aln.alphabet] for col in range(aln.alignlen)]
            # for each position in the alignment,
            # allocate a position to put the each child symbol from which each current node symbol score was determined
            self.backptr = [[[None for _ in aln.alphabet] for _ in range(aln.alignlen)] for _ in range(self.nChildren())]
            for col in range(aln.alignlen):
                for i in range(self.nChildren()):
                    # left child will contribute first
                    for a_parent in range(len(aln.alphabet)):
                        best_score = +9999999
                        best_symb = 0
                        for a in range(len(aln.alphabet)):
                            score = (scores[i][col][a] + (
                            1 if a != a_parent else 0))  # if we want to weight scores, this would need to change
                            if score < best_score:
                                best_symb = a
                                best_score = score
                        self.seqscores[col][a_parent] += best_score
                        self.backptr[i][col][a_parent] = best_symb
        else:
            self.seqscores = [[0 if a == sym else 999999 for a in aln.alphabet] for sym in
                              self.sequence]  # if we want to weight scores, this would need to change
        return self.seqscores

    def _backwardParsimony(self, aln, seq=None):
        """ Internal function that operates recursively to inspect scores to determine
            most parsimonious sequence, from root to leaves. """
        if self.sequence == None:  # no sequence has been assigned
            childbuf = [[] for _ in range(self.nChildren())]
            if self.nChildren() == 0:  # no children, so terminal, cannot propagate scores
                raise RuntimeError("No sequence assigned to leaf node:", self.label)
            if seq == None:  # Only root can do this, no parents to consider, so we pick the lowest scoring symbol
                currbuf = []
                for col in range(aln.alignlen):
                    min_score = 999999
                    min_symb = None
                    child_symb = [None for _ in range(self.nChildren())]
                    for a_parent in range(len(aln.alphabet)):
                        if self.seqscores[col][a_parent] < min_score:
                            min_score = self.seqscores[col][a_parent]
                            min_symb = a_parent
                            for i in range(self.nChildren()):
                                child_symb[i] = self.backptr[i][col][a_parent]
                    currbuf.append(aln.alphabet[min_symb])
                    for i in range(self.nChildren()):
                        childbuf[i].append(aln.alphabet[child_symb[i]])
                self.sequence = sequence.Sequence(currbuf, aln.alphabet, self.label, gappy=True)
            else:  # Non-root, but not leaf
                self.sequence = seq
                col = 0
                for sym_parent in self.sequence:
                    a_parent = aln.alphabet.index(sym_parent)
                    child_symb = [None for _ in range(self.nChildren())]
                    for i in range(self.nChildren()):
                        child_symb[i] = self.backptr[i][col][a_parent]
                        childbuf.append(aln.alphabet[child_symb[i]])
                    col += 1
            for i in range(self.nChildren()):
                self.children[i]._backwardParsimony(aln, sequence.Sequence(childbuf[i], aln.alphabet, self.label, gappy=True))
        return self.sequence

    def getSequence(self):
        """ Get the sequence for the node. Return None if no sequence is assigned.
            Requires that an alignment is associated with the tree, and that sequence names match node labels.
            If the explored node is not a leaf, the sequence can be determined by parsimony. """
        if self.sequence != None:  # a sequence has been assigned
            return self.sequence
        elif self.seqscores != None:  # inferred by parsimony but not yet assigned
            return None  # determine most parsimonous sequence, not yet implemented

    def isAncestorOf(self, node, transitive=True):
        """ Decide if this node is the ancestor of specified node.
            If transitive is True (default), all descendants are included.
            If transitive is False, only direct descendants are included. """
        for i in range(self.nChildren()):
            if node == self.children[i]:
                return True
            elif transitive:
                status = self.children[i].isAncestorOf(node, transitive)
                if status: return True
        else:
            return False

    def getAncestors(self, transitive=False):
        """ Retrieve and return (list of) parent nodes.
            If transitive is False (default), only the direct parent is included.
            If transitive is True, all parents (parents of parents etc) are included. """
        if self.parent == None:
            return []
        if not transitive:
            return [self.parent]
        else:
            parents = self.parent.getAncestors(transitive)
            parents.append(self.parent)
        return parents

    def getDescendants(self, transitive=False):
        """ Retrieve and return (list of) nodes descendant of this.
            If transitive is False (default), only direct descendants are included.
            If transitive is True, all descendants are (recursively) included. """
        children = []
        for i in range(self.nChildren()):
            children.append(self.children[i])
        if not transitive:
            return children
        else:
            grandchildren = []
            for c in children:
                d = c.getDescendants(transitive)
                if d:
                    grandchildren.extend(d)
            children.extend(grandchildren)
            return children

    def annotate_node(self, annotations, annotation_symbols= None, exclude_annotations = [], use_symbols=False ):
        annotation_string = "[&"
        for key, val_list in annotations[self].items():
            if type(val_list) != list:
                val_list = [val_list]
            if key not in exclude_annotations:
                # If we are using annotation symbols and the annotation has an associated symbol
                for val in val_list:
                    if use_symbols and val in annotation_symbols:
                        sorted_symbols = sorted([annotation_symbols[val] for val in val_list])
                        annotation_string += '%s="%s",' % (key, ' '.join(['%s' % (val,) for val in sorted_symbols]))
                    else:
                        annotation_string += '%s="%s",' % (key, ' '.join(['%s' % (val,) for val in val_list]))
        # Remove the final comma and add in a closing bracket
        annotation_string = annotation_string[0: len(annotation_string) - 1] + "]"
        if len(annotation_string) > 2:
            if ":" in self.label:
                self.label = self.label.split(":")[0] + annotation_string + self.label.split(":")[1]
            else:
                self.label = self.label + annotation_string

""" ----------------------------------------------------------------------------------------
    Methods for generating a single tree by clustering, here UPGMA Zvelebil and Baum p. 278
    ----------------------------------------------------------------------------------------"""


def runUPGMA(aln, measure, absoluteDistances=False):
    """ Generate an ultra-metric, bifurcating, rooted tree from an alignment based on pairwise distances.
        Use specified distance metric (see sequence.calcDistances).
        If absoluteDistances is True, the tree will be assigned the total distance from provided species.
        Otherwise, the relative addition at each path will be assigned."""
    D = {}
    N = {}  # The number of sequences in each node
    M = aln.calcDistances(measure)  # determine all pairwise distances
    nodes = [PhyloNode(label=seq.name) for seq in aln.seqs]  # construct all leaf nodes
    """ For each node-pair, assign the distance between them. """
    for i in range(len(nodes)):
        nodes[i].sequence = aln.seqs[i]
        nodes[i].dist = 0.0
        N[nodes[i]] = 1  # each cluster contains a single sequence
        for j in range(0, i):
            D[frozenset([nodes[i], nodes[j]])] = M[i, j]
    """ Now: treat each node as a cluster,
        until there is only one cluster left,
        find the *closest* pair of clusters, and
        merge that pair into a new cluster (to replace the two that merged).
        In each case, the new cluster is represented by the (phylo)node that is formed. """
    while len(N) > 1:  # N will contain all "live" clusters, to be reduced to a signle below
        closest_pair = (None, None)  # The two nodes that are closest to one another according to supplied metric
        closest_dist = None  # The distance between them
        for pair in D:  # check all pairs which should be merged
            dist = D[pair]
            if closest_dist == None or dist < closest_dist:
                closest_dist = dist
                closest_pair = list(pair)
        # So we know the closest, now we need to merge...
        x = closest_pair[0]  # See Zvelebil and Baum p. 278 for notation
        y = closest_pair[1]
        z = PhyloNode()  # create a new node for the cluster z
        z.dist = D.pop(frozenset([x, y])) / 2.0  # assign the absolute distance, travelled so far, note: this will change to relative distance later
        Nx = N.pop(x)  # find number of sequences in x, remove the cluster from list N
        Ny = N.pop(y)  # find number of sequences in y, remove the cluster from list N
        dz = {}  # new distances to cluster z
        for w in N:  # for each node w ...
            # we will merge x and y into a new cluster z, so need to consider w (which is not x or y)
            dxw = D.pop(frozenset([x, w]))  # retrieve and remove distance from D: x to w
            dyw = D.pop(frozenset([y, w]))  # retrieve and remove distance from D: y to w
            dz[w] = (Nx * dxw + Ny * dyw) / (Nx + Ny)  # distance: z to w
        N[z] = Nx + Ny  # total number of sequences in new cluster, insert new cluster in list N
        for w in dz:  # we have to run through the nodes again, now not including the removed x and y
            D[frozenset([z, w])] = dz[w]  # for each "other" cluster, update distance per EQ8.16 (Z&B p. 278)
        x.parent = z
        y.parent = z
        z.children = [x, y]
        nodes.append(z)
    if not absoluteDistances:
        x._propagateDistance(z.dist)  # convert absolute distances to relative by recursing down left path
        y._propagateDistance(z.dist)  # convert absolute distances to relative by recursing down right path
        z.dist = 0.0  # root z is at distance 0 from merged x and y
    return PhyloTree(z)  # make it to tree, return

""" ----------------------------------------------------------------------------------------
    Methods for processing files of trees on the Newick format
    ----------------------------------------------------------------------------------------"""

def _findComma(string, level=0):
    """ Find all commas at specified level of embedding """
    mylevel = 0
    commas = []
    for i in range(len(string)):
        if string[i] == '(':
            mylevel += 1
        elif string[i] == ')':
            mylevel -= 1
        elif string[i] == ',' and mylevel == level:
            commas.append(i)
    return commas

def parseNewickNode(string):
    """ Utility function that recursively parses embedded string using Newick format. """
    first = string.find('(')
    last = string[::-1].find(')')  # look from the back
    if first == -1 and last == -1:  # we are at leaf
        y = string.split(':')
        node = PhyloNode(label=y[0])
        if len(y) >= 2:
            node.dist = float(y[1])
        return node
    elif first >= 0 and last >= 0:
        # remove parentheses
        last = len(string) - last - 1  # correct index to refer from start instead of end of string
        embed = string[first + 1:last]
        tail = string[last + 1:]
        # find where corresp comma is
        commas = _findComma(embed)
        if len(commas) < 1:
            raise RuntimeError('Invalid format: invalid placement of "," in sub-string "' + embed + '"')
        prev_comma = 0
        child_tokens = []
        for comma in commas:
            child_tokens.append(embed[prev_comma:comma].strip())
            prev_comma = comma + 1
        child_tokens.append(embed[prev_comma:].strip())
        y = tail.split(':')
        node = PhyloNode(label=y[0])  # node is an instance of the PhyloNode() class
        if len(y) >= 2:
            node.dist = float(y[1])
        node.children = []
        for tok in child_tokens:
            child = parseNewickNode(tok)
            child.parent = node
            node.children.append(child)
        return node
    else:
        raise RuntimeError('Invalid format: unbalanced parentheses in sub-string "' + string + '"')


def parseNewick(string):
    """ Main method for parsing a Newick string into a (phylogenetic) tree.
        Handles labels (on both leaves and internal nodes), and includes distances (if provided).
        Returns an instance of a PhyloTree. """
    if string.find(';') != -1:
        string = string[:string.find(';')]
    return PhyloTree(parseNewickNode(string))


def readNewick(filename):
    """ Read file on Newick format.
        Returns an instance of a PhyloTree."""
    f = open(filename)
    string = ''.join(f)
    return parseNewick(string)


def writeNewickFile(filename, my_tree):
    with open(filename, 'w') as fh:
        print(my_tree, end="", file=fh)


def read_nexus(filename):
    """
    Read a file in Nexus format
    :param filename:
    :return:
    """
    f = open(filename)
    return parse_nexus(f)

def parse_nexus(string):
    string = string.read()
    lines = string.split("\n")
    annotation_dict = defaultdict(dict)
    for num, line in enumerate(lines):
        # print (line)
        if line.strip().startswith("dimensions ntax="):
                taxon_number = line.strip().split("dimensions ntax=")[1].split(";")[0]
        if line.strip().startswith("taxlabels"):
            taxon_num = num + 1
            while not lines[taxon_num].strip().startswith(";"):
                taxon_name = lines[taxon_num].split("[")[0].strip()
                for annot_line in lines[taxon_num].split("[&")[1].split(","):
                    #TODO: Make these regex calls
                    # print ("Annotation Key is ", annot_line.split("=")[0])
                    annot_key = annot_line.split("=")[0]
                    # print (annot_line.split("=")[1])
                    if '"' in annot_line:
                        annot_val = annot_line.split("=")[1].split("\"")[1]
                    else:
                        annot_val = annot_line.split("=")[1].split("]")[0]

                    annotation_dict[taxon_name][annot_key.strip()] = annot_val
                taxon_num +=1
        if line.strip().startswith("begin trees"):
            tree_num = num + 1
            tree = (lines[tree_num].split("[&R]")[1])
    phylo_tree = parseNewick(tree)
    nexus_annotations = annotations.NexusAnnotations(tree=phylo_tree)
    nexus_annotations.add_annotations(annotation_dict)
    phylo_tree.putAnnotations(nexus_annotations)
    ## Extract all of the annotations from the tree and add them to the NexusAnnotations object
    print ("Number of taxons is %s " % (taxon_number))
    return phylo_tree

""" ----------------------------------------------------------------------------------------
    Method for generating a PhyloTree with unique tip names
    ----------------------------------------------------------------------------------------"""

def get_unique_tree(tree):
    unique_tree = tree
    unique_labels = {}
    for node in unique_tree.getNodes():
        if node.isLeaf() and node.label in unique_labels:
            unique_labels[node.label] = unique_labels[node.label] + 1
            node.label = node.label + str(unique_labels[node.label])
        elif node.isLeaf():
            unique_labels[node.label] = 1
    return unique_tree

def unpack_list(list):
    return (" ".join(["%s"] * len(list)) + "!") % (x for x in list)

if __name__ == '__main__':
    tree = readNewick('/Users/mikael/simhome/ASR/edge1.nwk')
    print(tree)