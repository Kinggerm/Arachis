#! /usr/bin/env python
__version__ = "0.2"

# Features
# 1. involve gaps in genomes
# 2. allow gene name to be strings (numbers+alphabet+underline+short-line)
# 3. more relaxed way in parsing tree
# 4. no "$" in the end of a line stands for a circular chromosome, but only one chromosome is allowed (in which case,
#    meaning no "$" found in sequence, multi-line context would be regarded as one single chromosome) due to the
#    limitation of applying to tsp. Applying circular genomes together with linear genomes would be difficult for tsp solver.
# 5. multiprocessing

# Future
# involve branch length
# problem with the gap design in the binary matrix: by default each character would have equal prob (? or the average prob estimated from the data) for two states, which is not true for adjacency.
# why not apply the MULTIGAMMA to the adjacencies rather than apply a content coding together with the adjacencies?
# involve credibility interval?

import os
import sys
import time

major_py_version, minor_py_version = sys.version_info[:2]
if major_py_version == 3 and minor_py_version >= 5:
    python_version = "3.5+"
    from subprocess import getstatusoutput
elif major_py_version == 2 and minor_py_version >= 7:
    python_version = "2.7+"
    from commands import getstatusoutput
else:
    sys.stdout.write("Python version have to be 3.5+ or 2.7+\n")
    sys.exit(0)
dead_code = {"2.7+": (32512, "commands"), "3.5+": (127, "subprocess")}[python_version][0]

this_script_file = str(os.path.realpath(__file__))
path_of_this_script = os.path.split(os.path.realpath(__file__))[0]
# sys.path.append(os.path.join(path_of_this_script, ".."))


import subprocess
from multiprocessing import Pool, Lock
from arachis.genomeClass import *
from arachis.treeClass import get_tree, get_node_label_mapping
from math import log
from optparse import OptionParser


# test whether an external binary is executable
def executable(test_this):
    return True if os.access(test_this, os.X_OK) or getstatusoutput(test_this)[0] != dead_code else False


def get_options(print_title=""):
    usage = "python "+this_script_file
    parser = OptionParser(usage=usage, version=__version__, description=print_title)
    parser.add_option("-d", dest="data_file",
                      help="grimm format data file.")
    parser.add_option("-t", dest="tree_file",
                      help="newick format tree file.")
    parser.add_option("-o", dest="out_folder",
                      help="out put directory.")
    parser.add_option("--seed", dest="random_seed", default="12345",
                      help="random seed for raxml and concorde")
    parser.add_option("-p", dest="processing_num", default=4, type=int,
                      help="Number of processing for main program. [default: %default]")
    parser.add_option("--raxml", dest="which_raxml", default="raxmlHPC",
                      help="path to raxml. Default: \"raxmlHPC\" (in the path)")
    parser.add_option("--concorde", dest="which_concorde", default="concorde",
                      help="path to concorde. Default: \"concorde\" (in the path)")
    parser.add_option("--verbose", dest="verbose", default=False, action="store_true",
                      help="verbose mode.")
    parser.add_option("--keep", dest="keep_temp", default=False, action="store_true",
                      help="keep temp files.")
    parser.add_option("--continue", dest="resume", default=False, action="store_true",
                      help="continue mode.")
    options, argv = parser.parse_args()
    if not (options.data_file and options.tree_file and options.out_folder):
        parser.print_help()
        sys.stdout.write("Insufficient arguments!\n")
        sys.exit()
    elif os.path.exists(options.out_folder) and not options.resume:
        sys.stdout.write("Error: "+options.out_folder+" exists!\n")
        sys.exit()
    else:
        for input_file in (options.data_file, options.tree_file):
            if not os.path.isfile(input_file):
                sys.stdout.write("Error: " + input_file + " not found!\n")
                sys.exit()
        if not options.resume or not os.path.exists(options.out_folder):
            os.mkdir(options.out_folder)
        dependency = True
        for test_dependency in [options.which_raxml, options.which_concorde]:
            if not executable(test_dependency):
                dependency = False
                sys.stdout.write("Error: "+test_dependency+" not found!\n")
        if not dependency:
            sys.exit()
        return options, argv


def pool_multiprocessing(target, iter_args, constant_args, num_process):
    # parse args
    if iter_args:
        if type(iter_args) in {list, tuple}:
            if type(iter_args[0]) not in {list, tuple}:
                iter_args = [[each_arg] for each_arg in iter_args]
        else:
            sys.stderr.write("iter_args must be list/tuple!\n")
    else:
        return
    if constant_args:
        if type(constant_args) not in {list, tuple}:
            constant_args = [constant_args]
    else:
        constant_args = []
    # uninterrupted
    pool = Pool(processes=num_process)
    for this_arg in iter_args:
        pool.apply_async(target, tuple(list(this_arg) + list(constant_args)))
    pool.close()
    try:
        # apply_async().get()
        pool.join()
    except KeyboardInterrupt:
        pool.terminate()
        raise KeyboardInterrupt


def run_it_in_shell(command, verbose=False):
    if verbose:
        sys.stdout.write(command+"\n")
    this_pip = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    output, err = this_pip.communicate()
    if "error" in str(output) or "Error" in str(output):
        raise Exception("Error in running "+command.split()[0]+"\n"+output.decode("utf-8"))


def reroot_trees(tree, out_dir, resume, verbose):
    in_tree = get_tree(tree, rooting='force-rooted')
    anc_labels = []
    anc_trees = []
    # run_command = which_reroot+' '+tree+" > "+os.path.join(out_dir, "RerootedTrees")
    # run_it_in_shell(run_command, verbose)
    # with open(os.path.join(out_dir, "RerootedTrees")) as trees_handler:
    #     line = trees_handler.readline()
    #     while line:
    #         label = line.strip()
    #         tree = trees_handler.readline()
    #         anc_labels.append(label)
    #         anc_trees.append(os.path.join(out_dir, "rerooted."+label))
    #         open(anc_trees[-1], "w").write(tree)
    #         line = trees_handler.readline()
    for node in in_tree.postorder_internal_node_iter():
        if node.parent_node:
            anc_labels.append(node.label)
            anc_trees.append(os.path.join(out_dir, "rerooted."+node.label))
            if not (resume and
                    (os.path.exists(anc_trees[-1])
                     or os.path.exists(os.path.join(out_dir, "OutputGeneOrder_" + node.label)))):
                in_tree.reroot_at_edge(node.edge)
                in_tree.write(path=anc_trees[-1], schema="newick")
    return anc_labels, anc_trees


def run_anc_rec_with_raxml(tree, encoding, out_name, out_dir, model, which_raxml, extra_command="",
                           seed="", verbose=False):
    if seed:
        seed = " -p "+seed
    else:
        seed = ""
    run_command = which_raxml + " -f A -c 6 -m " + model + seed + " -t " + tree + " -s " + encoding + \
                  " -n " + out_name + " -w " + out_dir + extra_command  # -K GTR(default)
    run_it_in_shell(run_command, verbose)


def run_concorde(tsp_file, solution_file, directory, which_concorde, seed, verbose=False):
    original_wd = os.getcwd()
    temp_dir = solution_file+".temp"
    os.mkdir(temp_dir)
    os.chdir(temp_dir)
    if seed:
        seed = " -s " + seed
    else:
        seed = ""
    run_command = which_concorde + seed + " -o " + solution_file + " " + tsp_file
    run_it_in_shell(run_command, verbose)
    os.chdir(original_wd)
    os.system("rm -r "+temp_dir)


def read_raxml_probabilities(prob_file, tags):
    target_dict = {tag: [] for tag in tags}
    with open(prob_file) as prob_handler:
        line = prob_handler.readline()
        while line:
            genome_label = line.strip()
            line = prob_handler.readline()
            while line.strip():
                if genome_label in target_dict:
                    target_dict[genome_label].append([float(prob) for prob in line.strip().split()])
                line = prob_handler.readline()
            line = prob_handler.readline()
    return target_dict


def write_tsp_graph(tsp_dict, out_file):
    out_handler = open(out_file, "w")
    out_handler.write("NAME: "+os.path.basename(out_file)+"\n" +
                      "TYPE: TSP\n" +
                      "DIMENSION: "+str(len(tsp_dict))+"\n" +
                      "EDGE_WEIGHT_TYPE: EXPLICIT\n" +
                      "EDGE_WEIGHT_FORMAT: UPPER_ROW\n" +
                      "EDGE_WEIGHT_SECTION\n")
    max_index = max(tsp_dict)
    for node_1 in range(1, max_index + 1):
        for node_2 in range(node_1 + 1, max_index + 1):
            if node_2 in tsp_dict[node_1]:
                out_handler.write(str(tsp_dict[node_1][node_2]) + " ")
            else:
                out_handler.write("10000 ")
        out_handler.write("\n")
    out_handler.write("EOF\n")
    out_handler.close()


def trans_uni_id_to_gene(tsp_graph, unique_id_to_block, middle_to_block_id, this_anc):
    # print(unique_id_to_block)
    this_chromosome = []
    for i in range(0, len(tsp_graph), 3):
        # ?? next if $Index == $  # LTSPSolutionCycle;
        left = tsp_graph[i]
        middle = tsp_graph[i+1]
        right = tsp_graph[i+2]

        def get_info(input_index):
            if input_index % 2 == 0:
                this_info = {"name": unique_id_to_block[int(input_index / 2)], "head": False}
            else:
                this_info = {"name": unique_id_to_block[int((input_index + 1) / 2)], "head": True}
            return this_info

        def report_it(input_node):
            try:
                return get_info(input_node)
            except KeyError:
                return {"name": unique_id_to_block[middle_to_block_id[input_node]], "middle": True}

        try:
            left_info = get_info(left)
            right_info = get_info(right)
            middle_block_id = middle_to_block_id[middle]
        except KeyError:
            middle_info = get_info(middle)
            sys.stdout.write("Error: Dump Wrong Assembly on " + this_anc + "\n")
            sys.stdout.write("Triple : [{} {} {}]\n".format(left, middle, right))
            sys.stdout.write("Triple : [{} {} {}]\n".format(report_it(left), report_it(middle), report_it(right)))
            sys.exit()
        else:
            middle_name = unique_id_to_block[middle_block_id]
            if left_info["name"] == middle_name == right_info["name"] and left_info["head"] != right_info["head"]:
                this_chromosome.append(("-" if left_info["head"] is False else "")+middle_name)
            else:
                sys.stdout.write("Error: Dump Wrong Assembly on " + this_anc + "\n")
                sys.stdout.write("Triple : [{} {} {}]\n".format(left, middle, right))
                sys.stdout.write("Triple : [{} {} {}]\n".format(left_info, middle_name, right_info))
                sys.exit()
    return this_chromosome


def update_graph(node_1, node_2, value, this_graph):
    for key_1, key_2 in [(node_1, node_2), (node_2, node_1)]:
        if key_1 in this_graph:
            this_graph[key_1][key_2] = value
        else:
            this_graph[key_1] = {key_2: value}


def read_tsp_path(solution_file):
    this_path = []
    with open(solution_file) as file_handler:
        file_handler.readline()
        for line in file_handler:
            for node in line.rstrip().split():
                this_path.append(int(node) + 1)
    return this_path


def add_ancestral_node(tree_object):
    count = 0
    for node in tree_object.postorder_internal_node_iter():
        if node.parent_node:
            count += 1
            node.label = "A" + str(count)


def preparing_tree(input_tree, out_dir, resume):
    # remove branch lengths
    this_tre = get_tree(input_tree, suppress_edge_lengths=True, preserve_underscores=True)
    for node in this_tre.preorder_node_iter():
        if len(node.child_nodes()) != 2:
            raise TypeError("Please make sure your input tree is rooted and bifurcated!")
        break
    this_tre.ladderize(ascending=True)
    add_ancestral_node(this_tre)
    new_tre = os.path.join(out_dir, "gene_order.tre")
    if not (resume and os.path.exists(new_tre)):
        this_tre.write(path=new_tre, schema="newick")
    return new_tre


def reconstruct_ancestral_genomes(data_file, input_tree, out_dir, which_raxml, which_concorde,
                                  random_seed, resume, verbose, keep_temp, num_process):
    """"""
    """ parse grimm format genomes from file """
    genomes = GenomeList(data_file)
    """ preparing tree """
    tree_file = preparing_tree(input_tree, out_dir, resume)
    """ check labels """
    tree_labels = {tree_taxon.label for tree_taxon in get_tree(tree_file, preserve_underscores=True).taxon_namespace}
    grim_labels = genomes.label_set()
    if tree_labels != grim_labels:
        for lb in tree_labels - grim_labels:
            sys.stdout.write("grimm not in tree: " + lb + "\n")
        for lb in grim_labels - tree_labels:
            sys.stdout.write("tree not in grimm: " + lb + "\n")
        os.system("rm -r " + out_dir)
        raise ValueError("The species names in the input tree and alignment file may not match, please check!")
    """ reroot trees """
    ancestor_labels, rerooted_trees = reroot_trees(tree_file, out_dir, resume, verbose)
    """"""
    content_chars = [block.name for block in genomes.content_counter.character_list()]
    sys.stdout.write("Block names: ")
    for con in content_chars:
        sys.stdout.write(str(con) + " ")
    sys.stdout.write("\n")
    content_state = genomes.content_counter.state_set()
    sys.stdout.write("Block state: ")
    for sta in content_state:
        sys.stdout.write(str(sta) + " ")
    sys.stdout.write("\n")

    adjacency_chars = [((adj.left.name, adj.left.direction), (adj.right.name, adj.right.direction))
                       for adj in genomes.adjacency_counter.character_list()]
    # for adj in adjacency_chars:
    #     print(block_to_str(adj[0]), block_to_str(adj[1]))
    with_telomere = genomes.contains_telomere()
    # alphabet = genomes[0].get_phylip_alphabet()

    content_encoding = os.path.join(out_dir, "Content.encoding")
    if not (resume and os.path.exists(content_encoding)):
        open(content_encoding, "w").write(genomes.get_content_phylip())
    order_encoding = os.path.join(out_dir, "Order.encoding")
    if not (resume and os.path.exists(order_encoding)):
        open(order_encoding, "w").write(genomes.get_adjacency_phylip(reverse_code=True))

    sys.stdout.write("Ancestors: "+" ".join(ancestor_labels)+"\n")
    iter_arguments = list(range(len(ancestor_labels)))
    if resume:
        go_to = len(ancestor_labels) - 1
        while go_to >= 0:
            if os.path.exists(os.path.join(out_dir, "OutputGeneOrder_" + ancestor_labels[go_to])):
                del iter_arguments[go_to]
            go_to -= 1
        os.system("rm " + os.path.join(out_dir, "RAxML_") + "*")
    # subprocess will not return error under this code
    if num_process > 1:
        pool_multiprocessing(target=sub_reconstruct_anc_genomes, iter_args=iter_arguments,
                             constant_args=(
                                 ancestor_labels, rerooted_trees, content_encoding, content_chars, content_state,
                                 order_encoding, adjacency_chars, with_telomere, out_dir, which_concorde, which_raxml,
                                 random_seed, verbose, keep_temp),
                             num_process=num_process)
    else:
        for go_to in iter_arguments:
            sub_reconstruct_anc_genomes(go_to, ancestor_labels, rerooted_trees, content_encoding, content_chars,
                                        content_state, order_encoding, adjacency_chars, with_telomere, out_dir,
                                        which_concorde, which_raxml, random_seed, verbose, keep_temp)
    gene_order_out = open(os.path.join(out_dir, "OutputGeneOrder"), "w")
    for go_to in range(len(ancestor_labels)):
        temp_f = open(os.path.join(out_dir, "OutputGeneOrder_" + ancestor_labels[go_to]))
        gene_order_out.write(temp_f.read())
        temp_f.close()
        # os.remove(os.path.join(out_dir, "OutputGeneOrder_" + str(go_to)))
    gene_order_out.close()


def sub_reconstruct_anc_genomes(go_to, ancestor_labels, rerooted_trees, content_encoding, content_chars, content_state,
                                order_encoding, adjacency_chars, with_telomere, out_dir,
                                which_concorde, which_raxml, random_seed, verbose, keep_temp):
    time0 = time.time()
    this_anc = ancestor_labels[go_to]
    this_tree = rerooted_trees[go_to]
    # sys.stdout.write("\nAncestor: "+this_anc+"\n\n")
    # stagger the peak
    # time.sleep({118: 20, 119: 40, 120: 60}.get(go_to, 0))

    """ reconstruct gene content """
    # raxml
    out_content = "Content."+this_anc
    run_anc_rec_with_raxml(this_tree, content_encoding, out_content, os.path.abspath(out_dir),
                           "MULTIGAMMA", which_raxml, "", random_seed, verbose)
    # mapping nodes
    raxml_node_file = os.path.join(out_dir, "RAxML_nodeLabelledRootedTree."+out_content)
    mapping = get_node_label_mapping(this_tree, raxml_node_file, (True, True))
    # print(mapping)
    # read probabilities
    raxml_prob_file = os.path.join(out_dir, "RAxML_marginalAncestralProbabilities."+out_content)
    target = read_raxml_probabilities(raxml_prob_file, tags=("ROOT", mapping[this_anc]))
    # determine states
    content_anc_results = {}
    for go_to_char in range(len(content_chars)):
        probs = [sum([target[tag][go_to_char][j] for tag in target])/2. for j in range(max(content_state) + 1)]
        best_prob = max(probs)
        best_stat = probs.index(best_prob)
        content_anc_results[content_chars[go_to_char]] = {"state": best_stat, "prob": best_prob}

    """ reconstruct gene order """
    # why not MULTIGAMMA for adjacency?
    # raxml
    out_order = "Order."+this_anc
    run_anc_rec_with_raxml(this_tree, order_encoding, out_order, os.path.abspath(out_dir),
                           "BINGAMMA", which_raxml, "", random_seed, verbose)
    # mapping nodes
    raxml_node_file = os.path.join(out_dir, "RAxML_nodeLabelledRootedTree." + out_order)
    mapping = get_node_label_mapping(this_tree, raxml_node_file, (True, True))
    # print(mapping)
    # read probabilities
    raxml_prob_file = os.path.join(out_dir, "RAxML_marginalAncestralProbabilities." + out_order)
    target = read_raxml_probabilities(raxml_prob_file, tags=("ROOT", mapping[this_anc]))
    # determine states
    # PMAG use 0 as existence in raxml, follow which (here reverse_code=True) the result_index should be 0, else 1.
    result_index = 0
    adj_anc_results = {}
    for go_to_char in range(len(adjacency_chars)):
        best_prob = max([target[tag][go_to_char][result_index] for tag in target])
        adj_anc_results[adjacency_chars[go_to_char]] = {"prob": best_prob}

    """ assembly """
    # 1. relabel blocks to continuous numbers: rbcL*2, atpB*1, rps12*2 -> rbcL:1, rbcL:2, atpB:3, rps12:4, rps12:5
    uni_id_to_block = {}
    block_to_id_list = {}
    if with_telomere:
        uni_id_to_block[0] = "$"
        block_to_id_list["$"] = [0]
    count_label = 1
    for go_to_char in range(len(content_chars)):
        this_char = content_chars[go_to_char]
        block_to_id_list[this_char] = []
        for k in range(content_anc_results[this_char]["state"]):
            uni_id_to_block[count_label] = this_char
            # uni_id_to_block[-count_label] = this_char
            block_to_id_list[this_char].append(count_label)
            count_label += 1
    new_adjacency = []
    for adj in adjacency_chars:
        (block_l, dire_l), (block_r, dire_r) = adj
        prob = adj_anc_results[adj]["prob"]
        # print(block_to_str((block_l, dire_l)), block_to_str((block_r, dire_r)), prob)
        try:
            for left_id in block_to_id_list[block_l]:
                for right_id in block_to_id_list[block_r]:
                    new_adjacency.append(((left_id, dire_l), (right_id, dire_r), prob))
                    # print(">", ((left_id, dire_l), (right_id, dire_r), prob))
        except KeyError:
            raise KeyError(this_anc+"Bad adjacency " + str(adj))
    # 2. build graph
    dimension = 2*(len(uni_id_to_block) - int(with_telomere))
    sys.stdout.write(this_anc+ " block end nodes: " + str(int(not with_telomere)) + " - " + str(dimension) + "\n")
    # telomere would have extra label: cap labels
    cap_goto = dimension   # Cap Label starts from dimension
    log_base = 2
    caps = set()
    tsp_graph = {}
    for (id_left, dire_left), (id_right, dire_right), this_prob in new_adjacency:
        # left
        if id_left == 0:
            if this_prob > 0.5:
                cap_goto += 1
                id_left = cap_goto
                caps.add(cap_goto)
            else:
                continue
        elif dire_left:
            id_left = id_left*2
        else:
            id_left = id_left*2 - 1
        # right
        if id_right == 0:
            if this_prob > 0.5:
                cap_goto += 1
                id_right = cap_goto
                caps.add(cap_goto)
            else:
                continue
        elif dire_right:
            id_right = id_right * 2 - 1
        else:
            id_right = id_right * 2
        # print(id_left, id_right, this_prob, int(log((1000000 * (1 - this_prob) + 1), log_base) + 0.5))
        update_graph(id_left, id_right, int(log((1000000 * (1 - this_prob) + 1), log_base) + 0.5), tsp_graph)
    # weight 0 between all pairs of caps
    for cap_a in caps:
        for cap_b in caps:
            if cap_a != cap_b:
                update_graph(cap_a, cap_b, 0, tsp_graph)
    if cap_goto >= dimension + 1:
        sys.stdout.write(this_anc + " terminal nodes: " + str(dimension+1) + " - " + str(cap_goto) + "\n")
    # add middle to ensure connect between head and tail of the same gene
    middle_goto = cap_goto
    middle_to_block = {}
    for unique_id in range(1, int(dimension/2)+1):
        middle_goto += 1
        this_head = unique_id*2 - 1
        this_tail = unique_id*2
        middle_to_block[middle_goto] = unique_id
        update_graph(middle_goto, this_head, 0, tsp_graph)
        update_graph(middle_goto, this_tail, 0, tsp_graph)
    sys.stdout.write(this_anc + " in-block nodes: " + str(cap_goto+1) + " - " + str(middle_goto) + "\n")
    # middle is extremely distant to any other nodes
    for middle_id in range(cap_goto+1, middle_goto+1):
        id_of_same_block = middle_to_block[middle_id]
        head_of_same_block = 2*id_of_same_block - 1
        tail_of_same_block = 2*id_of_same_block
        same_block = {middle_id, head_of_same_block, tail_of_same_block}
        for unique_id_candidate in range(1, middle_goto+1):
            if unique_id_candidate not in same_block:
                update_graph(middle_id, unique_id_candidate, 10000000, tsp_graph)
    # 3. write tsp graph file
    tsp_file = os.path.abspath(os.path.join(out_dir, "pymag.TSPLIB."+this_anc))
    solution_file = os.path.abspath(os.path.join(out_dir, "pymag.TSPSolution." + this_anc))
    write_tsp_graph(tsp_graph, tsp_file)
    # sys.exit()
    run_concorde(tsp_file, solution_file, out_dir, which_concorde, random_seed)
    # parse_tsp_solution
    tsp_path = read_tsp_path(solution_file)
    len_path = len(tsp_path)
    if len_path:
        cap_ids = []
        if with_telomere and caps:
            for i in range(len_path):
                if tsp_path[i] in caps:
                    cap_ids.append(i)
        chromosomes = []
        # convert cycle back to gene order
        if not cap_ids:
            # if circular genome, find the correct start node
            go_to_middle = 0
            for go_to_middle in range(1, len_path):
                if tsp_path[go_to_middle] in middle_to_block:
                    break
            if go_to_middle:
                tsp_path = tsp_path[go_to_middle - 1:] + tsp_path[:go_to_middle - 1]
            # translate
            chromosomes.append(trans_uni_id_to_gene(tsp_path, uni_id_to_block, middle_to_block, this_anc))
        else:
            # sys.stdout.write("result contains " + str(len(cap_ids)) + " telomere.\n")
            for go_to_cap in range(len(cap_ids)-1):
                this_path = tsp_path[cap_ids[go_to_cap]+1: cap_ids[go_to_cap+1]]
                if this_path:
                    # sys.stdout.write("this path: " + " ".join([str(node) for node in this_path])+"\n")
                    chromosomes.append(trans_uni_id_to_gene(this_path, uni_id_to_block, middle_to_block, this_anc)
                                       + ["$"])
            this_path = tsp_path[cap_ids[-1]+1:] + tsp_path[:cap_ids[0]]
            if this_path:
                # sys.stdout.write("this path: " + " ".join([str(node) for node in this_path]) + "\n")
                chromosomes.append(trans_uni_id_to_gene(this_path, uni_id_to_block, middle_to_block, this_anc)
                                   + ["$"])
        """ write to file """
        this_gene_order_out = open(os.path.join(out_dir, "OutputGeneOrder_" + ancestor_labels[go_to]), "w")
        this_gene_order_out.write(">"+this_anc+"\n")
        for chromosome in chromosomes:
            this_gene_order_out.write(" ".join(chromosome)+"\n")
        this_gene_order_out.close()
    else:
        sys.stdout.write(this_anc+" Error: no available TSP solution!\n")
    if not keep_temp:
        os.remove(tsp_file)
        os.remove(solution_file)
        os.remove(this_tree)
        os.system("rm " + os.path.join(out_dir, "RAxML_") + "*." + ancestor_labels[go_to])
    sys.stdout.write(this_anc + " cost: " + str(round(time.time() - time0, 2)) + "s\n")
        # RAxML_marginalAncestralProbabilities.Content.A??
        # RAxML_marginalAncestralStates.Content.A??
        # RAxML_nodeLabelledRootedTree.Content.A??
        # RAxML_info.Content.A??
        # RAxML_marginalAncestralProbabilities.Order.A??
        # RAxML_marginalAncestralStates.Order.A??
        # RAxML_nodeLabelledRootedTree.Order.A??
        # RAxML_info.Order.A??


def main():
    time0 = time.time()
    options, argv = get_options()
    reconstruct_ancestral_genomes(options.data_file, options.tree_file, options.out_folder,
                                  options.which_raxml, options.which_concorde, options.random_seed,
                                  options.resume, options.verbose, options.keep_temp, options.processing_num)
    sys.stdout.write("Cost: "+str(round(time.time()-time0, 2))+"s\n")


if __name__ == '__main__':
    main()
