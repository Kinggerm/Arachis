#! /usr/bin/env python
from optparse import OptionParser
from arachis.treeClass import *


def get_options():
    parser = OptionParser(usage="map_tree_labels -f from_tree -t to_tree -o out_put.nex")
    parser.add_option("-f", dest="from_tree",
                      help="get \"node labels\" from this tree.")
    parser.add_option("-t", dest="to_tree",
                      help="get \"node labels\" to this tree.")
    parser.add_option("-l", dest="labels", default="tree1,tree2",
                      help="input labels split by comma. Default: tree1,tree2")
    parser.add_option("-o", dest="out_put",
                      help="output nex format file")
    parser.add_option("--ff", dest="from_format", default="newick",
                      help="format of tree with labels. Default: newick.")
    parser.add_option("--tf", dest="to_format", default="newick",
                      help="format of tree as backbone. Default: newick.")
    parser.add_option("--fr", dest="from_rooted", default=False, action="store_true",
                      help="whether the tree with labels is rooted. Default: False.")
    parser.add_option("--tr", dest="to_rooted", default=False, action="store_true",
                      help="whether the tree as backbone is rooted. Default: False.")
    options, argv = parser.parse_args()
    if not (options.from_tree and options.to_tree and options.out_put):
        sys.stdout.write("Insufficient arguments:\n")
        parser.print_help()
        sys.exit()
    return options, argv


def test():
    options, argv = get_options()
    test_b = tree_from_labels(options.from_tree, options.to_tree, options.labels.split(","),
                              schema=(options.from_format, options.to_format),
                              rooted=(options.from_rooted, options.to_rooted))
    if options.from_rooted and options.to_rooted:
        test_b.ladderize(ascending=True)
    test_b.write(path=options.out_put, schema="nexus")


if __name__ == '__main__':
    test()
