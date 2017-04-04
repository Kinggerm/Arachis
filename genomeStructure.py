#! /usr/bin/env python

try:
    from io import StringIO
except:
    import StringIO
import os
import sys
legal_block_name_char = set("0123456789:;<=>ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz|~[]")
legal_phylip_alphabet = "0123456789ABCDEFGHIJKLMNOPQRSTUV"


def check_validity(check_string):
    for char in list(check_string):
        if char not in legal_block_name_char:
            raise Exception("Unrecognised notation for \""+char+"\" in \""+check_string+"\"")


def block_from_str(block_str):
    if block_str.startswith("-"):
        check_validity(block_str[1:])
        return block_str[1:], False
    elif block_str == "$":
        return "$", None
    elif block_str == "*":
        return "*", None
    else:
        check_validity(block_str)
        return block_str, True


def block_to_str(block_tuple):
    name, direction = block_tuple
    return ("-" if direction is False else "")+name


def reverse_block(block):
    block_name, block_direction = block
    return block_name, None if block_direction is None else not block_direction


def reverse_adjacency(adjacency):
    return reverse_block(adjacency[1]), reverse_block(adjacency[0])


class Genome:
    def __init__(self, label="", seq="", verbose=False):
        self.label = label
        """ converted as grimm format defined in our analysis """
        # # if no "$" in the end of a line stands for circular chromosome, then use following code instead
        # self.sequence = "\n".join([block.strip()
        #                            for block in seq.replace("$", " $\n").split("\n")
        #                            if block.strip()]).strip() if seq else seq
        if "$" in seq:
            """ converted as standardized grimm format described in: http://grimm.ucsd.edu/GRIMM/grimm_instr.html """
            self.sequence = "".join([block+" $\n"
                                     for block in seq.replace("\n", "").split("$")
                                     if block.strip()]).strip() if seq else seq
        else:
            """ if no "$" in the end of a line stands for circular chromosome, but only one chromosome is allowed 
            (multi-line context would be regarded as one single chromosome) due to the limitation of applying to tsp """
            self.sequence = seq.replace("\n", "")
        self.__verbose = verbose
        self.__phylip_alphabet = legal_phylip_alphabet
        self.__contains_telomere = False
        self.__contains_gap = False
        # when the data contains gaps, this variable stores those blocks next to the gap
        self.__gap_adjacency = set()
        self.content = {}
        self.__content_list = []
        self.adjacency = {}
        self.__adjacency_list = []
        self.update()

    def get_grimm(self):
        return ">"+self.label+"\n"+self.sequence

    def check_content(self):
        if max(self.content.values()) > len(self.__phylip_alphabet):
            raise Exception("The total number of states is over " + str(len(self.__phylip_alphabet)))

    def update(self):

        def __add_content(ctt):
            if ctt == "$":
                self.__contains_telomere = True
            elif ctt == "*":
                self.__contains_gap = True
            elif ctt in self.content:
                self.content[ctt] += 1
            else:
                self.content[ctt] = 1
                self.__content_list.append(ctt)

        def __add_adjacency(adj):
            rev_adj = reverse_adjacency(adj)
            combined = sorted([adj, rev_adj])[0]
            if combined[0][0] == "*":
                self.__gap_adjacency.add(adj)
                self.__gap_adjacency.add(rev_adj)
            elif combined in self.adjacency:
                self.adjacency[combined] += 1
            else:
                self.adjacency[combined] = 1
                self.__adjacency_list.append(combined)

        self.content = {}
        self.__content_list = []
        self.adjacency = {}
        self.__adjacency_list = []
        for chromosome in self.sequence.split("\n"):
            if chromosome:
                chromosome = [block_from_str(block) for block in chromosome.split()]
                len_chrome = len(chromosome)
                for i in range(len_chrome):
                    __add_content(chromosome[i][0])  # could be (i - 1) % len_chrome
                    __add_adjacency((chromosome[(i - 1) % len_chrome], chromosome[i]))
        self.check_content()
        if self.__verbose:
            sys.stdout.write("Maximum number of states in "+self.label+" : "+str(max(self.content.values())))

    def __str__(self):
        return self.get_grimm()

    def __repr__(self):
        return "<Genome Object>\n>"+self.get_grimm()

    def __eq__(self, other):
        if type(other) == type(self):
            if self.adjacency == other.adjacency \
                    and self.content == other.content \
                    and self.contains_telomere() == other.contains_telomere() \
                    and self.contains_gap() == other.contains_gap():
                return True
            else:
                return False
        else:
            return False

    def content_list(self):
        return self.__content_list

    def adjacency_list(self):
        return self.__adjacency_list

    def phylip_alphabet(self):
        return self.__phylip_alphabet

    def contains_telomere(self):
        return self.__contains_telomere

    def contains_gap(self):
        return self.__contains_gap

    def unknown_adjacency(self):
        return self.__gap_adjacency


class GenomeList(list):
    def __init__(self, grimm_file_or_string=""):
        list.__init__([])
        self.__pointer = {}
        self.__contains_telomere = False
        self.content_counter = self.__ElementCounter()
        self.adjacency_counter = self.__ElementCounter()
        self.parse_grimm(grimm_file_or_string)

    def parse_grimm(self, grimm_file_or_string, add_mode=False):
        if not add_mode:
            self.clear()
            self.__pointer = {}
            self.content_counter = self.__ElementCounter()
            self.adjacency_counter = self.__ElementCounter()
        if grimm_file_or_string.strip():
            if os.path.isfile(grimm_file_or_string):
                grimm_handler = open(grimm_file_or_string)
            else:
                grimm_handler = StringIO(grimm_file_or_string)
            line = grimm_handler.readline()
            while line:
                line = line.split("#")[0].strip()
                if line.startswith(">"):
                    this_name = line[1:]
                    this_seq = ""
                    line = grimm_handler.readline()
                    while line and not line.startswith(">"):
                        this_seq += line.split("#")[0].strip() + "\n"
                        line = grimm_handler.readline()
                    self.append(Genome(this_name, this_seq))
                    self.content_counter.add_items(self[-1].content_list(), self[-1].content)
                    self.adjacency_counter.add_items(self[-1].adjacency_list(), self[-1].adjacency)
                else:
                    line = grimm_handler.readline()

    def get_content_phylip(self, gap_as_ambiguous=True, binary=False):
        state_set = self.content_counter.state_set()
        block_list = self.content_counter.character_list()
        # fill the jumpy states e.g. 1 2 4, then we fill 3 to make it continuous
        fake_sites = "" if binary else "".join([str(st) for st in range(max(state_set) + 1) if st not in state_set])
        head = " " + str(len(self)) + " " + str(len(block_list) + len(fake_sites)) + "\n"
        seqs = []
        for genome in self:
            # if there's a gap in the genome, add "-"(unknown) rather than "0"(doesn't exist) to missing block
            if gap_as_ambiguous and genome.contains_gap():
                alphabet = "-"+genome.phylip_alphabet()[1:]
            else:
                alphabet = genome.phylip_alphabet()
            max_index = 1 if binary else len(alphabet) - 1
            seq = "".join([alphabet[min(max_index, genome.content.get(b, 0))] for b in block_list]) + fake_sites
            seqs.append(genome.label + " " * 10 + seq + "\n")
        return head + "".join(seqs)

    def get_adjacency_phylip(self, gap_as_ambiguous=True, binary=True, reverse_code=False):
        state_set = self.adjacency_counter.state_set()
        adj_list = self.adjacency_counter.character_list()
        fake_sites = "" if binary else "".join([str(st) for st in range(max(state_set) + 1) if st not in state_set])
        head = " " + str(len(self)) + " " + str(len(adj_list) + len(fake_sites)) + "\n"
        seqs = []
        alphabet = self[0].phylip_alphabet()
        for genome in self:
            # if there's a gap in the genome, add "-" to candidate adjacency that involves the (block next to the gap)
            # or (block that are missing due to the gap)
            treat_gap = gap_as_ambiguous and genome.contains_gap()
            if treat_gap:
                block_set = self.content_counter.character_set()
                unknown_blocks = {b for b in block_set if b not in genome.content}
            else:
                unknown_blocks = {}
            seq = []
            max_index = 1 if binary else len(alphabet)-1
            for adj in adj_list:
                if adj in genome.adjacency:
                    seq.append(alphabet[min(max_index, genome.adjacency[adj])])
                elif treat_gap:
                    start, end = adj
                    missing_start = start[0] in unknown_blocks
                    missing_end = end[0] in unknown_blocks
                    start_t_break = (start, ("*", None)) in genome.unknown_adjacency()
                    break_t_end = (("*", None), end) in genome.unknown_adjacency()
                    if (missing_end or break_t_end) and (missing_start or start_t_break):
                        seq.append("-")
                    else:
                        seq.append(alphabet[0])
                else:
                    seq.append(alphabet[0])
            if binary and reverse_code:
                trans = {"1": "0", "0": "1"}
                seq = [trans[base] for base in seq]
            seqs.append(genome.label + " " * 10 + "".join(seq) + "\n")
        return head + "".join(seqs)

    class __ElementCounter:
        def __init__(self):
            self.__character_list = []
            self.__character_set = set()
            self.__state_set = set()

        def add_items(self, element_list, element_dict):
            for key in element_list:
                value = element_dict[key]
                if key not in self.__character_set:
                    self.__character_list.append(key)
                    self.__character_set.add(key)
                if value not in self.__state_set:
                    self.__state_set.add(value)

        def character_list(self):
            return list(self.__character_list)

        def state_set(self):
            return set(self.__state_set)

        def character_set(self):
            return set(self.__character_set)

    def __getitem__(self, item):
        if type(item) == int:
            return super(GenomeList, self).__getitem__(item)
        elif type(item) == str:
            return super(GenomeList, self).__getitem__(self.__pointer[item])

    def __contains__(self, item):
        if type(item) == str:
            return self.__pointer.__contains__(item)
        else:
            return super(GenomeList, self).__contains__(item)

    def append(self, item):
        super(GenomeList, self).append(item)
        self.__pointer[item.label] = len(self) - 1
        self.__contains_telomere = self.__contains_telomere or item.contains_telomere()

    def update(self):
        self.__contains_telomere = False
        self.__pointer = {}
        for i in range(len(self)):
            self.__pointer[self[i].label] = i
            self.__contains_telomere = self.__contains_telomere or self[i].contains_telomere()

    def contains_telomere(self):
        return self.__contains_telomere

    def label_set(self):
        return set(self.__pointer.keys())

    def label_list(self):
        return [genome.label for genome in self]

    def __str__(self):
        return "".join([str(genome)+"\n" for genome in self])

    def __repr__(self):
        return "<GenomeList Object>\n"+"".join([str(genome)+"\n" for genome in self])


def test():
    test_genome = GenomeList(">species1\na b -c d e$\n>species2\na -b c d e$\n>species3\na -b c * d e$\n")
    # test_genome = Genome("test", "-a $ b\n\nc $\nd $\ne * -f$")
    # test_genome = Genome("test", "e * -f$")
    print(test_genome)
    # print(test_genome[0].adjacency)
    # print(repr(test_genome["species1"]))
    print(test_genome.content_counter.character_list())
    print(test_genome.get_content_phylip())
    print([">".join([block_to_str(b) for b in site]) for site in test_genome.adjacency_counter.character_list()])
    print(test_genome.get_adjacency_phylip())

if __name__ == '__main__':
    test()
