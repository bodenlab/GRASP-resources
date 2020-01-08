

class TreeNodeObject:

    """
    label = str
    parent = TreeNodeObject
    dist = double
    id = int
    extent = boolean
    """
    def __init__(self, label, parent, dist, id, extent):

        self.parent = parent
        self.id = id
        self.extent = extent
        self.leaves_count = 0
        self._format_label(label, extent)
        if not dist:
            dist = 0.0
        self.dist = dist

        # init the arrays needed
        self.children = []
        self.leaves = dict() # dict of leaf node objects indexed on label
        self.intersect_id_list = []
        self.original_label = label # the label param will get formatted
        self.id = id
        self.score = 0.0
        self.in_intersection = False
        self.other_extent_count = 0
        self.dist_to_root = None

    def set_dist_to_root(self):
        if self.dist_to_root:
            return self.dist_to_root
        if self.parent:
            if not self.parent.get_dist_to_root():
                self.parent.set_dist_to_root()
            self.dist_to_root = self.dist + self.parent.get_dist_to_root()
            return self.dist_to_root
        # Parent or all (i.e. root)
        self.dist_to_root = self.dist
        return self.dist

    def get_dist_to_root(self):
        return self.dist_to_root

    def set_id(self, id):
        self.id = id

    def get_label(self):
        return self.label

    def get_id(self):
        return self.id

    def add_child(self, child):
        self.children.append(child)

    def set_leaves_count(self):
        if self.label == 'N0' or self.label == 'N1':
            i = 0
        if self.leaves_count:
            return self.leaves_count

        if self.extent:
            self.leaves_count = 1
            return self.leaves_count

        for child in self.children:
            self.leaves_count += child.set_leaves_count()
        if self.label == 'N0' or self.label == 'N1':
            i = 0
        return self.leaves_count

    def get_leaves_count(self):
        return self.leaves_count

    def _format_label(self, raw_label, extent):
        is_number = False
        try:
            float(raw_label)
            is_number = True
        except:
            is_number = False

        if not is_number and len(raw_label.split('\\|')) > 1:
            splitted = raw_label.split('\\|')
            if len(splitted) == 2:
                self.label = splitted[1]
            else:
                self.label = splitted[0]
        elif not extent and len(raw_label.split('_')) > 1:
            self.label = raw_label.split('_')[0]
        elif extent or raw_label == 'N0' or not is_number:
            self.label = raw_label
        else:
            print('--------- WARNING NO LABEL! ---------------')

    """
     * Make a mapping of the intersection and the labels contained in this node.
     *
     * Here we want to use the placement as an ID. This will be consistent acrocss all mappings.
     *
     * Note: we do this with the extentIntersection containing the nodes from the UNKNOWN tree,
     * this allows us to map the ID's back easily.
     *
     * @param extentIntersection
     """
    def build_intersection_label_mapping(self, extent_intersection_dict):
        if self.extent:
            # Check if it is in the intersection
            if extent_intersection_dict.get(self.label):
                self.in_intersection = True
                self.intersect_id_list.append(self.label)
                return True
            # Add this extent to the extents that weren't in the intersection count
            self.other_extent_count += 1
            return False
        for tree_node_object in self.children:
            tree_node_object.build_intersection_label_mapping(extent_intersection_dict)
            self.intersect_id_list += tree_node_object.get_intersect_ids()

        # Potentially keep track of the count
        if self.label == 'N0' or self.label == 'N1':
            i = 0
        self.other_extent_count = self.leaves_count - len(self.intersect_id_list)

        return False

    def reset_score(self):
        self.score = 0
        self.other_extent_count = 0
        for child in self.children:
            child.reset_score()

    def is_extent(self):
        return self.extent

    def is_in_intersection(self):
        return self.in_intersection

    def get_label(self):
        return self.label

    def add_to_score(self, score):
        self.score += score

    def get_other_extent_count(self):
        self.other_extent_count = self.leaves_count - len(self.intersect_id_list)
        return self.other_extent_count

    def get_intersect_ids(self):
        return self.intersect_id_list

    def get_children(self):
        return self.children

    def get_score(self):
        return self.score