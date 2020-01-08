import unittest
from tree_object import TreeObject
from tree_controller import TreeController

class TestStringMethods(unittest.TestCase):

    def test_same_tree(self):
        t1 = TreeObject('/Users/ariane/Documents/grasp_comparison/gabe/GRASP_recons/CYP2/20191216_CYP2U_Paper_Version_reconstructed-tree_GRASP.nwk') #'/Users/ariane/Documents/grasp_comparison/gabe/input/KARI/bkari1_rr99_sub0r.nwk') #("/Users/ariane/PycharmProjects/SimilarNodeFinder/data/0_10_dhad_28102018.nwk")
        t2 = TreeObject('/Users/ariane/Documents/grasp_comparison/gabe/GRASP_recons/CYP2/20191218_CYP2U-CYP2R-CYP2D_Paper_Version_reconstructed-tree_GRASP.nwk') #'/Users/ariane/Documents/grasp_comparison/gabe/input/KARI/bkari1_rr99_sub2r.nwk') #("/Users/ariane/PycharmProjects/SimilarNodeFinder/data/0_10_dhad_28102018.nwk")

        tc = TreeController()
        node_macthing = tc.get_similar_nodes_between_trees(t1, t2, True)

        for label_original, label_match in node_macthing.items():
            if label_original == 'N0':
                print(label_original, label_match.get_label())

    def test_not_same_tree(self):
        t1 = TreeObject("/Users/ariane/PycharmProjects/SimilarNodeFinder/data/0_10_dhad_28102018.nwk") #CYP2U1_Solo_Met_Rerooted.nwk")#0_10_dhad_28102018.nwk")
        t2 = TreeObject("/Users/ariane/PycharmProjects/SimilarNodeFinder/data/20_40_dhad_28102018.nwk") #CYP2U1_CYP2R1_Realigned_Met_Rerooted_Rotated.nwk")#

        tc = TreeController()
        node_macthing = tc.get_similar_nodes_between_trees(t1, t2, True)

        for label, node in node_macthing.items():
            print(label, node.get_label())

        self.assertEqual(node_macthing["N0"].get_label(), "N0")
        self.assertEqual(node_macthing["N8_0.974"].get_label(), "N11_0.969")
        self.assertEqual(node_macthing["N5_0.985"].get_label(), "N9_0.991")
        self.assertEqual(node_macthing["N3_0.994"].get_label(), "N6_0.033")

    def test_n0_ancestor(self):
        t1 = TreeObject("/Users/ariane/PycharmProjects/SimilarNodeFinder/data/2U_2R_2D_Ancestors.aln")
        t2 = TreeObject("/Users/ariane/PycharmProjects/SimilarNodeFinder/data/20_40_dhad_28102018.nwk")

        tc = TreeController()
        node_macthing = tc.get_similar_nodes_between_trees(t1, t2, True)

    def test_grasp(self):

        processed_dir = '/Users/ariane/Documents/grasp_paper/julian/Julian_GRASP_dump/Processed_results/'
        cyp2_comp = processed_dir + 'CYP2/CYP_tool_group_and_method_distance_comparison.txt'
        output_dir = '/Users/ariane/Documents/grasp_paper/julian/output/'
        # Here we want to do the comparison with a couple of different vis methods
        # Read in the file and convert it into a matrix
        order = []  # We want to have the order kept
        mean_dict = dict()
        stdev_dict = dict()
        order2 = []

        first = True
        with open(cyp2_comp, 'r+') as filein:
            for line in filein:
                if first:
                    first = False
                else:
                    line = line.split('\t')
                    label1 = line[0] + '-' + line[1] + '-' + line[2] + '-' + line[3]
                    label2 = line[4] + '-' + line[5] + '-' + line[6] + '-' + line[7]
                    order.append(label1)
                    order2.append(label2)
                    if not mean_dict.get(label1):
                        mean_dict[label1] = dict()
                        stdev_dict[label1] = dict()
                    mean_dict[label1][label2] = line[8]
                    stdev_dict[label1][label2] = line[9]

        with open(output_dir + 'cpy2_mean.csv', 'w+') as fileout:
            line = "test"
            for label in order2:
                line += "," + label
            fileout.write(line[:-1] + "\n")
            for label1 in order:
                line = label1
                for label2 in order2:
                    try:
                        line += "," + mean_dict[label1][label2]
                    except:
                        line += ", "
                fileout.write(line + "\n")


if __name__ == '__main__':
    unittest.main()