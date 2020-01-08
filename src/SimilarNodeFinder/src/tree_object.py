from tree_node_object import TreeNodeObject


class TreeObject:

    def __init__(self, file_path):
        self.root = None

        # Dictionary with node ids as the key
        self.node_dict = dict()

        # Dictionaries with the node labels as the key
        self.ancestor_node_dict = dict()
        self.extent_node_dict = dict()
        self.file_path = file_path
        self._load()

    def get_extent_dict(self):
        return self.extent_node_dict

    def get_ancestor_dict(self):
        return self.ancestor_node_dict

    def clear_scores(self):
        self.root.reset_score()

    def _load(self):
        self._load_newick(self.file_path)
        # set the distance to the root for all the nodes
        for node in self.node_dict:
            self.node_dict[node].set_dist_to_root()
        for node in self.ancestor_node_dict:
            self.ancestor_node_dict[node].set_dist_to_root()
            self.ancestor_node_dict[node].set_leaves_count()
        for node in self.extent_node_dict:
            self.extent_node_dict[node].set_dist_to_root()
            self.extent_node_dict[node].set_leaves_count()

    def get_root(self):
        if not self.root:
            self._load_newick(self.file_path)
        return self.root

    def get_extent_list(self):
        return self.extent_node_dict.keys()

    def get_ancestor_list(self):
        return self.ancestor_node_dict.keys()

    def _load_newick(self, file_path):
        newick_str = ""
        with open(file_path, 'r+') as newick_file:
            for line in newick_file:
                newick_str += line

        self.root = self._parse_newick(newick_str, self.root, [], 0, 0)

    def _parse_newick(self, newick_str, parent, node_ids, count, next_id):
        node = None
        newick_str = newick_str.replace('\t', '')
        start_idx = self._get_idx(newick_str, '(')
        end_idx = self._get_last_idx(newick_str, ')')
        if start_idx == -1 and end_idx == -1:
            node = self._parse_leaf_newick(newick_str, parent, next_id)
            next_id += 1
        elif start_idx >= 0 and end_idx >= 0:
            embed = newick_str[start_idx + 1: end_idx]
            tail = newick_str[end_idx + 1:]
            node = self._parse_internal_newick(embed, tail, parent, node_ids, count, next_id)
            next_id += 1

        if not self.node_dict.get(node.get_id()):
            self.node_dict[node.get_id()] = node

        return node

    def _get_dist(self, sub_str, idx_start):
        dist = sub_str[idx_start + 1:].replace(";", '')
        try:
            return float(dist)
        except:
            print("There was a non-digit number that was attempted to be parsed: ", sub_str)
            return 0.0

    def _parse_leaf_newick(self, newick_str, parent, next_id):

        split_idx = self._get_idx(newick_str, ':')

        if split_idx == -1:
            node = TreeNodeObject(newick_str, parent, None, next_id, True)
            next_id += 1

        else:
            label = newick_str[: split_idx].strip()
            if label[:2] == 'sp' or label[:2] == 'tr':
                label = label.split('|')[1]
            else:
                label = label.split('|')[0]
            dist = self._get_dist(newick_str, split_idx)
            node = TreeNodeObject(label, parent, dist, next_id, True)
            next_id += 1
            if not self.root:
                self.root = node

        self.extent_node_dict[node.get_label()] = node
        self.node_dict[node.get_id()] = node

        return node

    def _parse_internal_newick(self, embed, tail, parent, node_ids, count, next_id):
        split_idx = self._get_idx(tail, ':')
        if split_idx == -1:
            if len(tail) > 1:
                label = tail[split_idx + 1:].replace(';', '')
                node = TreeNodeObject(label, parent, None, next_id, False)
                next_id += 1
            else:
                node = TreeNodeObject("N" + str(count), parent, None, next_id, False)
                next_id += 1
        else: # There is a distance
            if len(tail) > 1:
                label = tail[: split_idx]
            else:
                print("------------------ ERROR NO INTERNAL NODE LABELS!!! ----------------------------")
                return None
            dist = self._get_dist(tail, split_idx)
            node = TreeNodeObject(label, parent, dist, count, False)
            if not self.root:
                self.root = node

        node_ids.append(count)
        comma = self._get_comma(embed)

        while comma != -1:

            to_process = embed[: comma]

            while count in node_ids:
                count += 1

            node.add_child(self._parse_newick(to_process, node, node_ids, count, next_id))
            next_id += 1

            if comma + 1 > len(embed):
                break

            embed = embed[comma + 1:]
            comma = self._get_comma(embed)

        self.ancestor_node_dict[node.get_label()] = node
        return node

    def _get_comma(self, newick_str):
        if len(newick_str) == 0:
            return -1
        level = 0
        for i in range(0, len(newick_str)):
            if newick_str[i] == '(':
                level += 1
            elif newick_str[i] == ')':
                level -= 1
            elif newick_str[i] == ',' and level == 0:
                return i

        return len(newick_str)

    def _get_idx(self, string, value):
        return string.find(value)

    def _get_last_idx(self, string, value):
        return string.rfind(value)