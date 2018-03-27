#! /usr/bin/env python
from copy import deepcopy
import dendropy
import os
import sys


def get_tree(tree_file_or_str_or_obj, name_space=None, schema="newick", edge_label=False, **kwargs):
    if not name_space:
        name_space = dendropy.TaxonNamespace()
    if type(tree_file_or_str_or_obj) == dendropy.Tree:
        tree_a = deepcopy(tree_file_or_str_or_obj)
        tree_a.taxon_namespace = name_space
    elif type(tree_file_or_str_or_obj) == str:
        if os.path.isfile(tree_file_or_str_or_obj):
            tree_a = dendropy.Tree.get(path=tree_file_or_str_or_obj, schema=schema, taxon_namespace=name_space,
                                       is_assign_internal_labels_to_edges=edge_label, **kwargs)
        else:
            try:
                tree_a = dendropy.Tree.get_from_string(tree_file_or_str_or_obj, schema=schema,
                                                       taxon_namespace=name_space,
                                                       is_assign_internal_labels_to_edges=edge_label, **kwargs)
            except ValueError:
                raise ValueError("Error: " + str(tree_file_or_str_or_obj) +
                                 " is neither a file or a " + schema + "string\n")
    else:
        sys.stdout.write("Error: parsing tree " + str(tree_file_or_str_or_obj))
        sys.exit()
    return tree_a


def node_bipartition_filter_fun(this_node):
    return lambda n: n.bipartition == this_node.bipartition


def get_node_bipartition(node, all_leaves):
    these_leaves = node.leaf_nodes()
    these_leaves = set([l.taxon.label for l in these_leaves]) & all_leaves
    other_leaves = all_leaves - these_leaves
    these_leaves = tuple(sorted(list(these_leaves)))
    other_leaves = tuple(sorted(list(other_leaves)))
    return these_leaves, other_leaves


def get_customized_node_bipartitions(tree, limiting_set=set(), rooted=False):
    if limiting_set:
        all_leaves = set([l.taxon.label for l in tree.leaf_nodes()]) & limiting_set
    else:
        all_leaves = set([l.taxon.label for l in tree.leaf_nodes()])
    if rooted:
        all_leaves.add("__root__")
    bipartitions = {}
    for node in tree.postorder_node_iter():
        these_leaves, other_leaves = get_node_bipartition(node, all_leaves)
        if these_leaves:
            if these_leaves in bipartitions:
                bipartitions[these_leaves].append(node)
                bipartitions[other_leaves].append(node)
            else:
                bipartitions[these_leaves] = [node]
                bipartitions[other_leaves] = [node]
            bipartitions[node] = sorted([these_leaves, other_leaves])
        else:
            pass
    return bipartitions


def get_corresponding(tree_customized_bipartitions, from_bipartition):
    if from_bipartition in tree_customized_bipartitions:
        return tree_customized_bipartitions[from_bipartition]
    else:
        return []


def label_rates(rooted_phylogram, chronogram):
    phylogram_taxa_set = set([l.taxon.label for l in rooted_phylogram.leaf_nodes()])
    chronogram_taxa_set = set([l.taxon.label for l in chronogram.leaf_nodes()])
    if phylogram_taxa_set == chronogram_taxa_set:
        phylo_bipartitions = get_customized_node_bipartitions(rooted_phylogram, rooted=True)
        for node in chronogram.postorder_node_iter():
            if node.parent_node:
                this_bipartition = get_node_bipartition(node, chronogram_taxa_set)[0]
                corresponding_node = get_corresponding(phylo_bipartitions, this_bipartition)
                if corresponding_node:
                    this_rate = float(corresponding_node[0].edge_length)/float(node.edge_length)
                    node.edge.annotations.add_new("rate", str(this_rate))
                else:
                    raise ValueError("Partition not found! Make sure phylogram and chronogram have the same topology!")
    else:
        raise ValueError("Unequal taxa set!")


def tree_from_labels(from_labeled_tree, to_tree, label_names=("from_tre", "to_tre"),
                     extra_joint_label=True, rooted=(False, False), schema=("newick", "newick"), strict=False):
    name_space = dendropy.TaxonNamespace()
    tree_a = get_tree(from_labeled_tree, name_space, schema[0])
    tree_b = get_tree(to_tree, name_space, schema[1])
    tree_a.is_rooted, tree_b.is_rooted = rooted
    tree_a_taxa_set = set([l.taxon.label for l in tree_a.leaf_nodes()])
    tree_b_taxa_set = set([l.taxon.label for l in tree_b.leaf_nodes()])
    if tree_a_taxa_set == tree_b_taxa_set:
        a_custom_bipartitions = get_customized_node_bipartitions(tree_a, rooted=rooted[0])
        if rooted[1]:
            tree_b_taxa_set.add("__root__")
        for node in tree_b.postorder_node_iter():
            node_label = node.label if (node.label and label_names[0]) else ""
            node.label = None
            # find the corresponding node
            this_b_bipartition = get_node_bipartition(node, tree_b_taxa_set)[0]
            corresponding_node = get_corresponding(a_custom_bipartitions, this_b_bipartition)
            corresponding_node = corresponding_node[0] if corresponding_node else None
            # copy label into annotations
            corresponding_label = corresponding_node.label if corresponding_node and label_names[1] else ""
            label_values = [corresponding_label, node_label]
            for i in range(2):
                if label_values[i]:
                    node.annotations.add_new(label_names[i], label_values[i])
            if extra_joint_label and (corresponding_label or node_label):
                node.annotations.add_new("|".join(label_names), "|".join(label_values))
            # copy annotations
            node.copy_annotations_from(corresponding_node)
            node.edge.copy_annotations_from(corresponding_node.edge)
    else:
        sys.stdout.write("Warning: unbalanced taxa, be aware of using strict="+str(strict)+" mode.")
        taxa_sets = tree_a_taxa_set | tree_b_taxa_set if strict else tree_a_taxa_set & tree_b_taxa_set
        a_custom_bipartitions = get_customized_node_bipartitions(tree_a, taxa_sets, rooted[0])
        for node in tree_b.postorder_node_iter():
            node_label = node.label if node.label and label_names[0] else ""
            node.label = None
            this_b_bipartition = get_node_bipartition(node, taxa_sets)[0]
            corresponding_nodes = get_corresponding(a_custom_bipartitions, this_b_bipartition)
            label_values = [c_node.label for c_node in corresponding_nodes] \
                if corresponding_nodes and label_names[1] else [""]
            label_values.append(node_label)
            len_label_values = len(label_values)
            new_label_names = [label_names[i == len_label_values - 1] for i in range(len(label_values))]
            for i in range(len_label_values):
                if label_values[i]:
                    node.annotations.add_new(new_label_names[i], label_values[i])
            if extra_joint_label and "".join(label_values):
                node.annotations.add_new("|".join(new_label_names), "|".join(label_values))
            for c_node in corresponding_nodes:
                node.copy_annotations_from(c_node)
                node.edge.copy_annotations_from(c_node.edge)
    return tree_b


def get_node_label_mapping(key_label_tree, value_label_tree, rooted=(False, False), schema=("newick", "newick")):
    name_space = dendropy.TaxonNamespace()
    tree_a = get_tree(key_label_tree, name_space, schema[0])
    tree_b = get_tree(value_label_tree, name_space, schema[1])
    tree_a.is_rooted, tree_b.is_rooted = rooted
    tree_a.encode_bipartitions(suppress_storage=True)
    tree_b.encode_bipartitions(suppress_storage=True)
    a_leaves = set([l.taxon.label for l in tree_a.leaf_nodes()])
    b_leaves = set([l.taxon.label for l in tree_a.leaf_nodes()])
    if set(a_leaves) == set(b_leaves):
        mapping = {}
        for node_b in tree_b.postorder_internal_node_iter():
            if node_b.label:
                corresponding_node_a = tree_a.find_node(filter_fn=node_bipartition_filter_fun(node_b))
                if corresponding_node_a:
                    mapping[corresponding_node_a.label] = node_b.label
    else:
        sys.stdout.write("Error: unequal taxa labels.")
        sys.exit()
    return mapping