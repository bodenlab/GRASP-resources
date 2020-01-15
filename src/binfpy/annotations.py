from collections import defaultdict
from phylo import *
import phylo
import matplotlib
import random

matplotlib_colours={'aliceblue':'#F0F8FF','aqua':'#00FFFF','aquamarine':'#7FFFD4','azure':'#F0FFFF','beige':'#F5F5DC','bisque':'#FFE4C4','blanchedalmond':'#FFEBCD','blue':'#0000FF','blueviolet':'#8A2BE2','brown':'#A52A2A','burlywood':'#DEB887','cadetblue':'#5F9EA0','chartreuse':'#7FFF00','chocolate':'#D2691E','coral':'#FF7F50','cornflowerblue':'#6495ED','crimson':'#DC143C','cyan':'#00FFFF','darkblue':'#00008B','darkcyan':'#008B8B','darkgoldenrod':'#B8860B','darkgreen':'#006400','darkkhaki':'#BDB76B','darkmagenta':'#8B008B','darkolivegreen':'#556B2F','darkorange':'#FF8C00','darkorchid':'#9932CC','darkred':'#8B0000','darksalmon':'#E9967A','darkseagreen':'#8FBC8F','darkslateblue':'#483D8B','darkturquoise':'#00CED1','darkviolet':'#9400D3','deeppink':'#FF1493','deepskyblue':'#00BFFF','dodgerblue':'#1E90FF','firebrick':'#B22222','floralwhite':'#FFFAF0','forestgreen':'#228B22','fuchsia':'#FF00FF','gainsboro':'#DCDCDC','ghostwhite':'#F8F8FF','gold':'#FFD700','goldenrod':'#DAA520','green':'#008000','greenyellow':'#ADFF2F','honeydew':'#F0FFF0','hotpink':'#FF69B4','indianred':'#CD5C5C','indigo':'#4B0082','khaki':'#F0E68C','lavender':'#E6E6FA','lavenderblush':'#FFF0F5','lawngreen':'#7CFC00','lemonchiffon':'#FFFACD','lightblue':'#ADD8E6','lightcoral':'#F08080','lightcyan':'#E0FFFF','lightgoldenrodyellow':'#FAFAD2','lightgreen':'#90EE90','lightgray':'#D3D3D3','lightpink':'#FFB6C1','lightsalmon':'#FFA07A','lightseagreen':'#20B2AA','lightskyblue':'#87CEFA','lightslategray':'#778899','lightsteelblue':'#B0C4DE','lightyellow':'#FFFFE0','lime':'#00FF00','limegreen':'#32CD32','magenta':'#FF00FF','maroon':'#800000','mediumaquamarine':'#66CDAA','mediumblue':'#0000CD','mediumorchid':'#BA55D3','mediumpurple':'#9370DB','mediumseagreen':'#3CB371','mediumslateblue':'#7B68EE','mediumspringgreen':'#00FA9A','mediumturquoise':'#48D1CC','mediumvioletred':'#C71585','midnightblue':'#191970','mintcream':'#F5FFFA','mistyrose':'#FFE4E1','moccasin':'#FFE4B5','navajowhite':'#FFDEAD','navy':'#000080','oldlace':'#FDF5E6','olive':'#808000','olivedrab':'#6B8E23','orange':'#FFA500','orangered':'#FF4500','orchid':'#DA70D6','palegoldenrod':'#EEE8AA','palegreen':'#98FB98','paleturquoise':'#AFEEEE','palevioletred':'#DB7093','papayawhip':'#FFEFD5','peachpuff':'#FFDAB9','peru':'#CD853F','pink':'#FFC0CB','plum':'#DDA0DD','powderblue':'#B0E0E6','purple':'#800080','red':'#FF0000','rosybrown':'#BC8F8F','royalblue':'#4169E1','saddlebrown':'#8B4513','salmon':'#FA8072','sandybrown':'#FAA460','seagreen':'#2E8B57','seashell':'#FFF5EE','sienna':'#A0522D','skyblue':'#87CEEB','slateblue':'#6A5ACD','springgreen':'#00FF7F','steelblue':'#4682B4','tan':'#D2B48C','teal':'#008080','thistle':'#D8BFD8','tomato':'#FF6347','turquoise':'#40E0D0','violet':'#EE82EE','wheat':'#F5DEB3','yellow':'#FFFF00','yellowgreen':'#9ACD32'}

twenty_colours ={'dodgerblue':'#1E90FF', 'orangered':'#FF4500', 'greenyellow':'#ADFF2F', 'orchid':'#DA70D6'}

symbols = ["*", "^", "!", "#", "~", "+", ":", "<", ">", "@", "%", "=", "-", "a", "b", "c", "d", "e" "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh", "ii", "jj", "kk", "ll", "mm", "nn", "oo", "pp", "qq", "rr", "ss", "tt", "uu", "vv", "ww", "xx", "yy", "zz"]

class NexusAnnotations():
    """
    Defines a set of annotations for a phylogenetic tree, to be written to NEXUS format
    """
    node_annotations = defaultdict(list)
    leaf_annotations = defaultdict(list)
    used_colours = set()
    annotated_nodes = set()

    def __init__(self, colour_dict=matplotlib_colours, symbol_list = symbols, tree=None, node_annotations=defaultdict(list), leaf_annotations=defaultdict(list), annotation_symbols={}):
        self.tree = tree
        self.node_annotations = node_annotations
        self.leaf_annotations = leaf_annotations
        self.annotation_symbols = annotation_symbols
        self.colour_dict = colour_dict
        self.symbol_list = symbol_list
        self.used_colours = set()
        self.annotated_nodes = set()

        # if type(self.tree) != PhyloTree:
        #     raise RuntimeError("NexusAnnotations need the tree to be a PhyloTree object")

        if tree:
            self.add_original_annotations()

    def add_original_annotations(self):
        """
        Add an entry for the original labels of the tree, so we can map back if we need to
        """
        for node in self.tree.getNodes():
            if node.isLeaf():
                self.leaf_annotations[node] = {"Original" : node.label}

    def add_annotations(self, annotation_dict):
        for node in self.tree.getNodes():
            if node.label in annotation_dict:
                for key, val in annotation_dict[node.label].items():
                    self.leaf_annotations[node][key] = val
        # print (self.leaf_annotations)


    def add_annotation(self, node, key="", annotation=[], annotate_descendants=False):

        self.node_annotations[node] = {key:annotation}
        self.annotated_nodes.add(node)

        for annot in annotation:
            if annot not in self.annotation_symbols:
                symbol = self.generate_annotation_symbol(annot)
                self.add_annotation_symbol(annot, symbol)
            else:
                symbol = self.generate_annotation_symbol(annot)
                self.add_annotation_symbol(annot, symbol)


        # Add in a symbol to represent this annotation if one doesn't already exist



        if annotate_descendants:
            for descendant in node.getDescendants():
                self.add_annotation(descendant, key, annotation)


    def add_colour_annotation(self, node, colour=None, random_colour=False, colour_descendants=True):
        """
        Add a single colour annotation to a set of nodes
        :param labels: The nodes to annotate
        :param colour: The colour to annotate with
        """

        if not colour:
            colour = self.generate_colour()


        if colour_descendants:
            for descendant in node.getDescendants(transitive=True):
                self.node_annotations[descendant] = {"!color": colour}
                self.annotated_nodes.add(descendant)
                self.used_colours.add(colour)


    def add_colour_annotations(self, nodes, colour_list=None, colour_descendants=False):
        """
        Add multiple colour annotations to a list of nodes
        :param nodes: A list of lists of nodes to annotate
        :param colour_list: The colour to annotate each list of nodes with
        """

        if not colour_list:
            colour_list = self.generate_colour_list(len(nodes))

        for index, node_set in enumerate(nodes):
            set_colour = self.colour_dict[colour_list[index]]


            print(set_colour)
            print(index, node_set)
            for node in node_set:
                self.add_colour_annotation(node, set_colour, colour_descendants=colour_descendants)

    def generate_annotation_symbol(self, annotation):

        i = 0
        while i < len(self.symbol_list):
            symbol = random.choice(self.symbol_list)
            if symbol not in self.annotation_symbols.values():
                return symbol
            else:
                return symbol
        i+=1

    def add_annotation_symbol(self, symbol, annotation):

        self.annotation_symbols[symbol] = annotation



    def generate_colour(self):
        """
        Generate a colour that hasn't been used yet in this set of Nexus Annotations
        :return: A unique colour
        """

        # i = 0
        # while i < len(self.colour_dict.values()):
        #     colour = random.choice(list(self.colour_dict.values()))
        #     if colour not in self.used_colours:
        #         return colour
        #     else:
        #         return colour
        # i+=1
        colour = random.choice(list(self.colour_dict.values()))
        return colour
    def generate_colour_list(self, num):
        return num

# tree = phylo.readNewick("/Users/gabefoley/Dropbox/PhD/Projects/Phylo Island/Species trees/species_tree_115.nwk")

