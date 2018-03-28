#! /usr/bin/env python
from optparse import OptionParser
from arachis.genomeClass import *


def get_options():
    usage = "python this_script.py -1 grimm1 -2 grimm2"
    parser = OptionParser(usage=usage)
    parser.add_option('-1', dest='grimm1',
                      help='')
    parser.add_option('-2', dest='grimm2',
                      help='')
    options, argv = parser.parse_args()
    if not (options.grimm1 and options.grimm2):
        parser.print_help()
        sys.stdout.write('\n######################################\nERROR: Invalid arguments!\n\n')
        exit()
    return options, argv


def compare_genomes(genomes_1, genomes_2):
    labels_1 = genomes_1.label_set()
    labels_2 = genomes_2.label_set()
    statistics = {}
    for genome in labels_1 | labels_2:
        exist_1 = genome in genomes_1
        exist_2 = genome in genomes_2
        if exist_1 and exist_2:
            matched = genomes_1[genome] == genomes_2[genome]
        else:
            matched = False
        statistics[genome] = {1: exist_1, 2: exist_2, "equal": matched}
    return statistics


def color_it(print_it, determine):
    if determine:
        return '\033[93m' + print_it + '\033[0m'
    else:
        return '\033[91m' + print_it + '\033[0m'


def main():
    options, argv = get_options()
    grimm1, grimm2 = options.grimm1, options.grimm2
    list1 = GenomeList(grimm1)
    list2 = GenomeList(grimm2)

    stat = compare_genomes(list1, list2)

    # preparing result table
    all_genomes = sorted(list(stat.keys()))
    total = len(all_genomes)
    result = [["genome", "in file1", "in file2", "match"]]
    for genome in all_genomes:
        result.append([genome, str(stat[genome][1]), str(stat[genome][2]), str(stat[genome]["equal"])])
    sys.stdout.write("\n")

    # print nice
    print_width = [max([len(row[column]) for row in result]) + 1 for column in range(len(result[0]))]
    for col in range(len(result[0])):
        sys.stdout.write(result[0][col].rjust(print_width[col]))
    sys.stdout.write("\n")
    difference = 0
    for row in result[1:]:
        match = row[3] == "True"
        difference += int(not match)
        for col in range(len(result[0])):
            sys.stdout.write(color_it(row[col].rjust(print_width[col]), match))
        sys.stdout.write("\n")
    sys.stdout.write("Matched/Unmatched: "+str(total - difference)+"/"+str(difference)+"\n\n")


if __name__ == '__main__':
    main()
