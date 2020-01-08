from tree_node_object import TreeNodeObject

class TreeController:

    def __init__(self):
        self.breakout = False
        self.best_node = None


    def _get_intersection(self, tree_known_ancs, tree_unknown_ancs):
        known_extent_dict = tree_known_ancs.get_extent_dict()
        unknown_extent_dict = tree_unknown_ancs.get_extent_dict()

        # Get the interection of the smaller
        search_list = known_extent_dict.keys()
        search_dict = unknown_extent_dict
        if len(known_extent_dict) > len(unknown_extent_dict):
            search_list = unknown_extent_dict.keys()
            search_dict = known_extent_dict

        intersection = dict()
        for label in search_list:
            if search_dict.get(label):
                intersection[label] = search_dict.get(label)

        return intersection

    """
    Similar to the get_similar_nodes_efficient implementation in GRASP.

    It finds a single matching node for a particular label
    """

    def _get_similar_node_for_node(self, tree_known_ancs, tree_unknown_ancs, node_label, same_tree):
        intersection = []  # The intersection between the known and unknown trees, array of tree node objects

        # Set the IDs in the intersection to be based on the placement in array
        id = 0
        for tree_node_object in intersection:
            tree_node_object.setId(id)
            id += 1

        tree_known_ancs.get_root().build_intersection_label_mapping(intersection)
        tree_unknown_ancs.get_root().build_intersection_label_mapping(intersection)

        # Just look for a single child.


    """
    Similar to the get_similar_nodes_efficient implementation in GRASP.
    
    It matches the nodes between two trees, the tree_known_ancs is used as the reference.
    
    TreeUnknownAncs is used as the one which gets mapped to the tree_known.
    """
    def get_similar_nodes_between_trees(self, tree_known_ancs, tree_unknown_ancs, same_tree):
        intersection = self._get_intersection(tree_known_ancs, tree_unknown_ancs) # The intersection between the known and unknown trees, array of tree node objects

        # Set the IDs in the intersection to be based on the placement in array
        id = 0
        for label, tree_node_object in intersection.items():
            tree_node_object.set_id(id)
            id += 1

        tree_known_ancs.get_root().build_intersection_label_mapping(intersection)
        tree_unknown_ancs.get_root().build_intersection_label_mapping(intersection)

        nodes = tree_known_ancs.get_ancestor_dict()
        tree_unknown_ancs_root = tree_unknown_ancs.get_root()

        similar_node_dict = dict()

        for label, ancestor_node in nodes.items():
            ancestor_intersection = ancestor_node.get_intersect_ids()

            if len(ancestor_intersection) < 1:
                similar_node_dict[ancestor_node.get_label()] = None

            else:
                self._score_nodes(ancestor_intersection, tree_unknown_ancs_root, ancestor_node)
                similar_node_dict[ancestor_node.get_label()] = self.best_node

            # reset for next node
            self.best_node = None

            tree_unknown_ancs.clear_scores()
            tree_known_ancs.clear_scores()

        # We want to go through and check if any of the nodes had a score of 0. This would mean 1 of 2 things:
        # 1. equal number of correctly and incorrectly placed nodes
        # 2. there were no extents in both trees under the node of interest
        # For situation 2 we can place it base on the correctly placed nodes

        return similar_node_dict


    """
    Scores nodes (similar to efficient version in GRASP)
    Determines the score based on how many extents were included under a particular node
    
    extent_id_map = dict of ID's that map to nodes
    node = a tree node object that we're currently inspecting
    
    The extend_id_map includes the extents under a specific node
    """
    def _score_nodes(self, extent_id_list, node, ancs_node_compared_with):
        if self.breakout:
            return 0.0

        value = 1.0
        score = 0.0
        if node.is_extent():
            # If this is not in both trees then consider it neutral and don't include in the score
            if not node.is_in_intersection():
                return 0.0
            # If it's in the list we want, then we return a negative (i.e. good low is best) score
            if node.get_label() in extent_id_list:
                return -value
            # If it's in the intersection and in the
            if not node.get_label() in extent_id_list:
                return value

        for child_tree_node in node.get_children():
            score += self._score_nodes(extent_id_list, child_tree_node, ancs_node_compared_with)
            if self.breakout:
                return 0.0

        node.add_to_score(score)

        if not self.best_node:
            self.best_node = node
        elif node.get_score() < self.best_node.get_score():
            if ancs_node_compared_with.label == 'N0':
                i = 0

            self.best_node = node
        elif node.get_score() == self.best_node.get_score():
            if ancs_node_compared_with.label == 'N0':
                i = 0
            # Here we use the number of other extents included to be the tiebreaker
            # we choose the one which has the fewest extra extents (i.e. should give the most similar reconstruction
            if node.get_other_extent_count() < self.best_node.get_other_extent_count():
                self.best_node = node

            if node.get_other_extent_count() == self.best_node.get_other_extent_count():
                this_dist_diff = abs(node.get_dist_to_root() - ancs_node_compared_with.get_dist_to_root())
                best_node_dist_diff = abs(self.best_node.get_dist_to_root() - ancs_node_compared_with.get_dist_to_root())
                #print(this_dist_diff, best_node_dist_diff)
                if this_dist_diff < best_node_dist_diff:
                    self.best_node = node

        return node.get_score()

