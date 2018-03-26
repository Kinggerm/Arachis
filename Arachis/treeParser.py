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
