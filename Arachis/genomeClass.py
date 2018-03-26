#! /usr/bin/env python
import sys
import os
import time
from copy import deepcopy
from itertools import combinations
from io import StringIO
try:
    import psutil
except ImportError:
    m_process = None
else:
    m_process = psutil.Process(os.getpid())

LEGAL_BLOCK_NAME_CHAR = set(".0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz|~[]")
LEGAL_PHYLIP_ALPHABET = "0123456789ABCDEFGHIJKLMNOPQRSTUV"
NO_DIRECTION_NAME = {"*", "$"}
SB_TABLE = {}


def check_char_validity(check_string):
    for c_id in range(len(check_string)):
        if check_string[c_id] not in LEGAL_BLOCK_NAME_CHAR:
            if check_string[c_id] in {"*", "$"}:
                raise ValueError(check_string[c_id] + " does not accept \"-\" operation!")
            elif check_string[c_id] == "-" and c_id != 0:
                pass
            else:
                raise ValueError("Invalid notation for \"" + check_string[c_id] + "\" in \"" + check_string + "\"\n"
                                 "Avail characters for 1st letter: " + "".join(list(LEGAL_BLOCK_NAME_CHAR)) + "\n"
                                 "Avail characters for others:" + "".join(list(LEGAL_BLOCK_NAME_CHAR)) + "-")


class SignedBlock(object):
    __slots__ = ["name", "direction"]

    def __init__(self, name, check_char=True):
        if type(name) == str:
            if name[0] == "-":
                if check_char:
                    check_char_validity(name[1:])
                here_name, direction = name[1:], False
            elif name in NO_DIRECTION_NAME:
                here_name, direction = name, None
            else:
                if check_char:
                    check_char_validity(name)
                here_name, direction = name, True
            object.__setattr__(self, "name", here_name)
            object.__setattr__(self, "direction", direction)
            if name not in SB_TABLE:
                SB_TABLE[name] = self
        else:
            raise ValueError("name must be a str type!")

    def __neg__(self):
        this_name = ("-" if self.direction else "") + self.name
        if this_name in SB_TABLE:
            return SB_TABLE[this_name]
        else:
            return SignedBlock(this_name, check_char=False)

    def __abs__(self):
        if self.name in SB_TABLE:
            return SB_TABLE[self.name]
        else:
            return SignedBlock(self.name, check_char=False)

    def __str__(self):
        return ("-" if self.direction is False else "") + self.name

    def __repr__(self):
        return str(type(self)) + " " + str(self)

    def __eq__(self, other):
        return True if type(
            other) == SignedBlock and self.name == other.name and self.direction == other.direction else False

    def __ne__(self, other):
        return True if self.name != other.name or self.direction != other.direction else False

    def __lt__(self, other):
        return (self.name, self.direction).__lt__((other.name, other.direction))

    def __le__(self, other):
        return (self.name, self.direction).__le__((other.name, other.direction))

    def __gt__(self, other):
        return (self.name, self.direction).__gt__((other.name, other.direction))

    def __ge__(self, other):
        return (self.name, self.direction).__ge__((other.name, other.direction))

    def __hash__(self):
        return hash((self.name, self.direction))

    def __deepcopy__(self, memodict={}):
        return self

    def __mul__(self, other):
        if other == 1:
            return self
        elif other == -1:
            return -self
        else:
            raise TypeError("unsupported operand for *: " + str(other) + " 'SignedBlock'")

    def __rmul__(self, other):
        return self.__mul__(other)


GapBlock = SignedBlock("*")
TerminalBlock = SignedBlock("$")


class Adjacency:
    def __init__(self, left, right):
        assert type(left) == SignedBlock
        assert type(right) == SignedBlock
        candidate = sorted([(left, right), (-right, -left)])
        self.left, self.right = candidate[0]
        self.__block_set = {left, right, -left, -right}

    def __eq__(self, other):
        return True if self.left == other.left and self.right == other.right else False

    def __lt__(self, other):
        return (self.left, self.right).__lt__((other.left, other.right))

    def __le__(self, other):
        return (self.left, self.right).__le__((other.left, other.right))

    def __gt__(self, other):
        return (self.left, self.right).__gt__((other.left, other.right))

    def __ge__(self, other):
        return (self.left, self.right).__ge__((other.left, other.right))

    def __hash__(self):
        return hash((self.left, self.right))

    def __str__(self):
        return str(self.left) + ">>" + str(self.right)

    def __repr__(self):
        return str(type(self)) + " " + str(self)

    def __contains__(self, item):
        return item in self.__block_set

    def __getitem__(self, item):
        if type(item) == SignedBlock:
            if item == self.left:
                return 1
            elif item == -self.left:
                return -1
            elif item == self.right:
                return 2
            elif item == -self.right:
                return -2
            else:
                raise KeyError(str(item) + " not Found!")
        elif item == 1:
            return self.left
        elif item == 2:
            return self.right
        elif item == -1:
            return -self.left
        elif item == -2:
            return -self.right
        else:
            raise KeyError("The key must be SignedBlock type or integer: 1, -1, 2, -2!")

    def __sub__(self, other):
        if type(other) == SignedBlock:
            return self[[None, 2, 1, -1, -2][self[other]]]
        else:
            raise TypeError("Only \"Adjacency - SignedBlock\" was allowed!")

    def is_gap_adj(self):
        return True if self.left == GapBlock or self.right == GapBlock else False

    def is_terminal_adj(self):
        return True if self.left == TerminalBlock or self.right == TerminalBlock else False


class Chromosome:
    def __init__(self, chromosome_str_or_list, build_bl_hash=True, build_adj_hash=True, build_adjacent_blocks=True,
                 check_redundant_gaps=True, check_content=True, mirror=False):
        self.__list = []
        self.circular = True
        self.__build_bl_hash = build_bl_hash
        self.__build_adj_hash = build_adj_hash
        self.__build_adjacent_blocks = build_adjacent_blocks
        self.block_pointer = {}
        self.adjacency_pointer = {}
        self.__uniq_adjacency_list = []
        self.__gap_adjacency = {}
        self.block_downstream = {}
        self.__hash = None
        self.__branching_points = None
        # optional
        self.__gaps_filled = None
        self.__block_set = None
        self.__block_list = None
        further_update = True
        if chromosome_str_or_list:
            if type(chromosome_str_or_list) == list:
                if check_content:
                    for block in chromosome_str_or_list:
                        type_b = type(block)
                        if type_b == str:
                            if block in SB_TABLE:
                                self.__list.append(SB_TABLE[block])
                            else:
                                self.__list.append(SignedBlock(block))
                        elif type_b == SignedBlock:
                            self.__list.append(block)
                        else:
                            raise ValueError(
                                "Invalid input chromosome: list of " + str(block) + "(Type " + str(type_b) + ")!")
                else:
                    self.__list = list(chromosome_str_or_list)
            elif type(chromosome_str_or_list) == str:
                for block in chromosome_str_or_list.split():
                    if block in SB_TABLE:
                        self.__list.append(SB_TABLE[block])
                    else:
                        self.__list.append(SignedBlock(block))
            elif type(chromosome_str_or_list) == Chromosome \
                    or issubclass(type(chromosome_str_or_list), Chromosome):
                self.__list = list(chromosome_str_or_list)
                self.__len = len(chromosome_str_or_list)
                self.__contains_gap = chromosome_str_or_list.contains_gap()
                self.circular = chromosome_str_or_list.circular
                further_update = False
                if chromosome_str_or_list.is_hashed():
                    self.__hash = hash(chromosome_str_or_list)
                if chromosome_str_or_list.is_branching_points_detected():
                    self.__branching_points = tuple(chromosome_str_or_list.get_branching_points())
                if chromosome_str_or_list.is_block_list_detected():
                    self.__block_list = chromosome_str_or_list.block_list()
                    self.__block_set = chromosome_str_or_list.block_set()
                if chromosome_str_or_list.is_gap_filled():
                    self.__gaps_filled = chromosome_str_or_list.gap_filled()
                if mirror and chromosome_str_or_list.block_pointer:
                    self.block_pointer = chromosome_str_or_list.block_pointer
                else:
                    self.build_bl_hash()
                if mirror and chromosome_str_or_list.adjacency_pointer:
                    self.adjacency_pointer = chromosome_str_or_list.adjacency_pointer
                    self.__uniq_adjacency_list = chromosome_str_or_list.adjacency_unique_list()
                    self.__gap_adjacency = chromosome_str_or_list.unknown_adjacency()
                else:
                    self.build_adj_hash()
                if mirror and chromosome_str_or_list.block_downstream:
                    self.block_downstream = chromosome_str_or_list.block_downstream
                else:
                    self.build_adjacent_blocks()
            else:
                raise ValueError("Invalid input chromosome type:" + str(type(chromosome_str_or_list)) + "!")
            if further_update:
                if self.__list[-1].name == "$":
                    self.circular = False
                else:
                    self.circular = True

                # remove redundant gaps
                if check_redundant_gaps:
                    len_p = len(self.__list)
                    i = 0
                    while i < len_p:
                        if self.__list[i] == self.__list[(i + 1) % len_p] == GapBlock:
                            del self.__list[i]
                            len_p = len(self.__list)
                        else:
                            i += 1
                self.__len = len(self.__list)
                # make comparison faster when it is circular, make find_all faster
                if self.__build_bl_hash:
                    self.build_bl_hash()
                if self.__build_adj_hash:
                    self.build_adj_hash()
                if self.__build_adjacent_blocks:
                    self.build_adjacent_blocks()
                if GapBlock in self:
                    self.__contains_gap = True
                else:
                    self.__contains_gap = False

    def __str__(self):
        return " ".join([str(signed_b) for signed_b in self])

    def __repr__(self):
        return str(type(self)) + str(self)

    def __eq__(self, other):
        if type(other) == type(self) and self.circular == other.circular:
            if self.circular:
                for match_start_id in self.find_all_block(other[0]):
                    order1 = self.__list[match_start_id:] + self.__list[:match_start_id]
                    order2 = [-bl for bl in (self.__list[match_start_id + 1:] + self.__list[:match_start_id + 1])[::-1]]
                    if list(other) == order1:
                        return True
                    elif list(other) == order2:
                        return True
                else:
                    return False
            else:
                if self.__list == list(other) or self.reverse_list() == list(other):
                    return True
                else:
                    return False
        return False

    def __len__(self):
        return self.__len

    def __iter__(self):
        for block in self.__list:
            yield block

    def __getitem__(self, item):
        if type(item) == int:
            return self.__list[item]
        elif type(item) == slice:
            here_step = 1 if item.step is None else item.step
            if item.start is not None and item.stop is not None:
                if item.start >= item.stop and here_step > 0:
                    if self.circular:
                        return self.__list[item.start::here_step] + self.__list[:item.stop:here_step]
                    else:
                        raise IndexError("Invalid index form:" + str(item) + " for a non-circular chromosome!")
                elif item.stop >= item.start and here_step < 0:
                    if self.circular:
                        return [-block
                                for block in self.__list[:item.stop:here_step] + self.__list[item.start::here_step]]
                    else:
                        raise IndexError("Invalid index form:" + str(item) + " for a non-circular chromosome!")
            if here_step < 0:
                return [-block for block in self.__list[item]]
            else:
                return self.__list[item]
        else:
            raise IndexError("Invalid index form: " + str(type(item)) + "!")

    def __contains__(self, item):
        if self.__build_bl_hash:
            return item in self.block_pointer
        else:
            return item in self.__list or -item in self.__list

    def __bool__(self):
        return bool(self.__list)

    def __hash__(self, update=False):
        if self.__hash and not update:
            return self.__hash
        else:
            comp_rev = self.reverse_list()
            if self.circular:
                start_block = sorted(self, key=lambda x: x.name, reverse=True)[0]
                candidates = []
                start_block_ids_1 = self.find_all_block(start_block)
                for start in start_block_ids_1:
                    candidates.append(tuple(self[start:] + self[:start]))
                # reverse
                for b_id, candidate_b in enumerate(comp_rev):
                    if candidate_b.name == start_block.name:
                        candidates.append(tuple(comp_rev[b_id:] + comp_rev[:b_id]))
                candidates.sort()
                self.__hash = hash(candidates[0])
            else:
                self.__hash = hash(sorted([tuple(list(self)), tuple(list(comp_rev))])[0])
            return self.__hash

    def __list__(self):
        return list(self.__list)

    def list(self):
        return list(self.__list)

    def len(self):
        return self.__len

    def is_hashed(self):
        return False if self.__hash is None else True

    def block_set(self, signed=True, update=False):
        if update:
            self.__block_set = None
            self.__block_list = None
        if self.__block_set is None:
            self.__block_set = set()
            self.__block_list = []
            for block in self:
                if block not in self.__block_set and -block not in self.__block_set:
                    self.__block_set.add(block)
                    self.__block_list.append(block)
        if signed:
            return set(self.__block_set)
        else:
            return {abs(b) for b in self.__block_set}

    def block_list(self, signed=True, update=False):
        if update:
            self.__block_set = None
            self.__block_list = None
        if self.__block_list is None:
            self.__block_set = set()
            self.__block_list = []
            for block in self:
                if block not in self.__block_set and -block not in self.__block_set:
                    self.__block_set.add(block)
                    self.__block_list.append(block)
        if signed:
            return list(self.__block_list)
        else:
            return [abs(b) for b in self.__block_list]

    def is_block_list_detected(self):
        return False if self.__block_set is None else True

    def contains_gap(self):
        return self.__contains_gap

    def find_all_block(self, item):
        if type(item) == str:
            item = SB_TABLE.get(item, SignedBlock(item))
        if self.__build_bl_hash:
            if item in self.block_pointer:
                return set(self.block_pointer[item])
            else:
                return set()
        else:
            return set([go_to for go_to, candidate in enumerate(self) if candidate == item or candidate == -item])

    def reverse_list(self, bl_container=None):
        if not bl_container:
            bl_container = self
        if type(bl_container) == list:
            return [-bl for bl in bl_container[::-1]]
        else:
            return bl_container[::-1]

    def complementary_reverse(self):
        if self.circular:
            self.__list = self.reverse_list()
            if self.__build_bl_hash:
                self.build_bl_hash()
        else:
            self.__list[:-1] = self.reverse_list()[1:]
            if self.__build_bl_hash:
                self.build_bl_hash()

    def complementary_reversed(self):
        if self.circular:
            return Chromosome(self.reverse_list(),
                              build_bl_hash=self.__build_bl_hash, build_adj_hash=self.__build_adj_hash,
                              build_adjacent_blocks=self.__build_adjacent_blocks, check_redundant_gaps=False,
                              check_content=False)
        else:
            return Chromosome(self.reverse_list()[1:] + [TerminalBlock],
                              build_bl_hash=self.__build_bl_hash, build_adj_hash=self.__build_adj_hash,
                              build_adjacent_blocks=self.__build_adjacent_blocks, check_redundant_gaps=False,
                              check_content=False)

    def inverse(self, start_b, end_b):
        start_b, end_b = min(start_b, end_b) % len(self), max(start_b, end_b) % len(self)
        self.__list[start_b:end_b + 1] = self.reverse_list(self.__list[start_b:end_b + 1])
        if self.__build_bl_hash:
            self.build_bl_hash()
        if self.__build_adj_hash:
            self.build_adj_hash()
        if self.__build_adjacent_blocks:
            self.build_adjacent_blocks()

    def inverted(self, start_b, end_b, **build):
        start_b, end_b = min(start_b, end_b) % len(self), max(start_b, end_b) % len(self)
        build_bl_hash = self.__build_bl_hash if "build_bl_hash" not in build else build["build_bl_hash"]
        build_adj_hash = self.__build_adj_hash if "build_adj_hash" not in build else build["build_adj_hash"]
        build_adj_b = self.__build_adjacent_blocks if "build_adjacent_blocks" not in build \
            else build["build_adjacent_blocks"]
        b_list = self.__list[:start_b] + self.reverse_list(self.__list[start_b:end_b + 1]) + self.__list[end_b + 1:]
        return Chromosome(b_list, build_bl_hash=build_bl_hash, build_adj_hash=build_adj_hash,
                          build_adjacent_blocks=build_adj_b, check_redundant_gaps=False, check_content=False)

    def insert(self, at, item):
        if type(item) == Chromosome:
            self.__list = self.__list[:at] + list(item) + self.__list[at:]
            if self.__block_set:
                for b in item:
                    if b not in self.__block_set:
                        self.__block_set.add(b)
                        self.__block_list.append(b)
        elif type(item) == SignedBlock:
            self.__list.insert(at, item)
            if self.__block_set:
                if item not in self.__block_set:
                    self.__block_set.add(item)
                    self.__block_list.append(item)
        else:
            raise ValueError("Invalid item type to insert")
        self.__len = len(self.__list)
        if self.__build_bl_hash:
            self.build_bl_hash()
        if self.__build_adj_hash:
            self.build_adj_hash()
        if self.__build_adjacent_blocks:
            self.build_adjacent_blocks()

    def inserted(self, at, item, **build):
        if type(item) == Chromosome or type(item) == list:
            new_list = self.__list[:at] + list(item) + self.__list[at:]
        elif type(item) == SignedBlock:
            new_list = self.__list[:at] + [item] + self.__list[at:]
        else:
            raise ValueError("Invalid item type to insert")
        build_bl_hash = self.__build_bl_hash if "build_bl_hash" not in build else build["build_bl_hash"]
        build_adj_hash = self.__build_adj_hash if "build_adj_hash" not in build else build["build_adj_hash"]
        build_adj_b = self.__build_adjacent_blocks if "build_adjacent_blocks" not in build \
            else build["build_adjacent_blocks"]
        return Chromosome(new_list, build_bl_hash=build_bl_hash, build_adj_hash=build_adj_hash,
                          build_adjacent_blocks=build_adj_b, check_redundant_gaps=False, check_content=False)

    def delete(self, start, stop=None):
        stop = start + 1 if not stop else stop
        if self.__block_set:
            to_del = list(self.__list[start: stop])
        del self.__list[start: stop]
        self.__len = len(self.__list)
        if self.__build_bl_hash:
            self.build_bl_hash()
        if self.__build_adj_hash:
            self.build_adj_hash()
        if self.__build_adjacent_blocks:
            self.build_adjacent_blocks()
        if self.__block_set:
            for b in to_del:
                if b not in self.block_pointer:
                    self.__block_set.remove(b)
                    self.__block_list.remove(b)

    def deleted(self, start, stop=None, **build):
        stop = start + 1 if not stop else stop
        build_bl_hash = self.__build_bl_hash if "build_bl_hash" not in build else build["build_bl_hash"]
        build_adj_hash = self.__build_adj_hash if "build_adj_hash" not in build else build["build_adj_hash"]
        build_adj_b = self.__build_adjacent_blocks if "build_adjacent_blocks" not in build \
            else build["build_adjacent_blocks"]
        return Chromosome(self.__list[:start] + self.__list[stop:],
                          build_bl_hash=build_bl_hash, build_adj_hash=build_adj_hash,
                          build_adjacent_blocks=build_adj_b, check_redundant_gaps=False, check_content=False)

    def scj_broken(self, fission_before_bl):
        if self.circular:
            return [Chromosome(self[fission_before_bl:fission_before_bl] + [TerminalBlock],
                               build_bl_hash=self.__build_bl_hash, build_adj_hash=self.__build_adj_hash,
                               build_adjacent_blocks=self.__build_adjacent_blocks, check_redundant_gaps=False,
                               check_content=False)]
        else:
            return [Chromosome(self.__list[:fission_before_bl] + [TerminalBlock],
                               build_bl_hash=self.__build_bl_hash, build_adj_hash=self.__build_adj_hash,
                               build_adjacent_blocks=self.__build_adjacent_blocks, check_redundant_gaps=False,
                               check_content=False),
                    Chromosome(self.__list[fission_before_bl:-1] + [TerminalBlock],
                               build_bl_hash=self.__build_bl_hash, build_adj_hash=self.__build_adj_hash,
                               build_adjacent_blocks=self.__build_adjacent_blocks, check_redundant_gaps=False,
                               check_content=False)]

    def scj_join(self, other_chromosome=None):
        new_list = self[:-1] if other_chromosome is None else self[:-1] + list(other_chromosome)
        this_circular = self.circular if other_chromosome is None else self.circular or other_chromosome.circular
        if this_circular:
            raise TypeError("scj fussion operation does not apply to circular chromosomes!")
        else:
            return Chromosome(new_list, build_bl_hash=self.__build_bl_hash, build_adj_hash=self.__build_adj_hash,
                              build_adjacent_blocks=self.__build_adjacent_blocks, check_redundant_gaps=False,
                              check_content=False)

    def build_bl_hash(self):
        self.block_pointer = {}
        for block_id in range(len(self)):
            this_block = self[block_id]
            if this_block not in self.block_pointer:
                self.block_pointer[this_block] = {block_id}
            else:
                self.block_pointer[this_block].add(block_id)
            neg_block = -this_block
            if neg_block not in self.block_pointer:
                self.block_pointer[neg_block] = {block_id}
            else:
                self.block_pointer[neg_block].add(block_id)

    def build_adj_hash(self):
        self.adjacency_pointer = {}
        for go_to in range(len(self)):
            go_to_minus_1 = (go_to - 1) % self.__len
            this_adj = Adjacency(self[go_to_minus_1], self[go_to])
            if this_adj.is_gap_adj():
                if this_adj in self.__gap_adjacency:
                    self.__gap_adjacency[this_adj].append(go_to)
                else:
                    self.__gap_adjacency[this_adj] = [go_to]
            elif this_adj in self.adjacency_pointer:
                self.adjacency_pointer[this_adj].append(go_to)
            else:
                self.adjacency_pointer[this_adj] = [go_to]
                self.__uniq_adjacency_list.append(this_adj)

    def build_adjacent_blocks(self):
        self.block_downstream = {}
        for go_to in range(len(self)):
            go_to_minus_1 = (go_to - 1) % self.__len
            this_block = self[go_to]
            this_name = this_block.name
            this_direction = this_block.direction
            if this_block != GapBlock:
                if this_name not in self.block_downstream:
                    self.block_downstream[this_name] = {"head": {}, "tail": {}}
                go_to_add_1 = (go_to + 1) % self.__len
                if this_direction:
                    if -self[go_to_minus_1] not in self.block_downstream[this_name]["head"]:
                        self.block_downstream[this_name]["head"][-self[go_to_minus_1]] = {(go_to, go_to_minus_1)}
                    else:
                        self.block_downstream[this_name]["head"][-self[go_to_minus_1]].add((go_to, go_to_minus_1))
                    if self[go_to_add_1] not in self.block_downstream[this_name]["tail"]:
                        self.block_downstream[this_name]["tail"][self[go_to_add_1]] = {(go_to, go_to_add_1)}
                    else:
                        self.block_downstream[this_name]["tail"][self[go_to_add_1]].add((go_to, go_to_add_1))
                else:
                    if self[go_to_add_1] not in self.block_downstream[this_name]["head"]:
                        self.block_downstream[this_name]["head"][self[go_to_add_1]] = {(go_to, go_to_add_1)}
                    else:
                        self.block_downstream[this_name]["head"][self[go_to_add_1]].add((go_to, go_to_add_1))
                    if -self[go_to_minus_1] not in self.block_downstream[this_name]["tail"]:
                        self.block_downstream[this_name]["tail"][-self[go_to_minus_1]] = {(go_to, go_to_minus_1)}
                    else:
                        self.block_downstream[this_name]["tail"][-self[go_to_minus_1]].add((go_to, go_to_minus_1))

    def adjacency_unique_list(self):
        if not self.__uniq_adjacency_list:
            self.build_adj_hash()
        return list(self.__uniq_adjacency_list)

    def unknown_adjacency(self):
        if not self.__uniq_adjacency_list:
            self.build_adj_hash()
        return deepcopy(self.__gap_adjacency)

    def adjacency_differs_from(self, other, gap_contains_gene=True):
        differences = {"gap": True, "record": {}, "count": 0}
        # self_adj = deepcopy(self.adjacency_pointer)
        # other_adj = deepcopy(other.adjacency_pointer)
        if gap_contains_gene and self.contains_gap() and other.contains_gap():
            pass
        elif gap_contains_gene and self.contains_gap():
            for adj in self.__uniq_adjacency_list:
                this_count = len(self.adjacency_pointer[adj])
                that_count = len(other.adjacency_pointer.get(adj, []))
                if this_count > that_count:
                    differences["record"][adj] = [this_count, that_count]
                    differences["count"] += this_count - that_count
        elif gap_contains_gene and other.contains_gap():
            for adj in other.adjacency_unique_list():
                this_count = len(self.adjacency_pointer.get(adj, []))
                that_count = len(other.adjacency_pointer[adj])
                if this_count < that_count:
                    differences["record"][adj] = [this_count, that_count]
                    differences["count"] += that_count - this_count
        else:
            if not self.contains_gap() and not other.contains_gap():
                differences["gap"] = False
            total_adj = self.adjacency_unique_list()
            for adj in other.adjacency_unique_list():
                if adj not in self.adjacency_pointer:
                    total_adj.append(adj)
            for adj in total_adj:
                this_count = len(self.adjacency_pointer.get(adj, []))
                that_count = len(other.adjacency_pointer.get(adj, []))
                if this_count != that_count:
                    differences["record"][adj] = [this_count, that_count]
                    differences["count"] += abs(this_count - that_count)
        return differences

    def content_differs_from(self, other, gap_contains_gene=True):
        differences = {"gap": True, "record": {}, "count": 0}
        if gap_contains_gene and self.contains_gap() and other.contains_gap():
            pass
        elif gap_contains_gene and self.contains_gap():
            for block in self.block_list():
                if block != GapBlock:
                    this_count = len(self.find_all_block(block))
                    that_count = len(other.find_all_block(block))
                    if this_count > that_count:
                        differences["record"][block] = [this_count, that_count]
                        differences["count"] += this_count - that_count
        elif gap_contains_gene and other.contains_gap():
            for block in other.block_list():
                if block != GapBlock:
                    this_count = len(self.find_all_block(block))
                    that_count = len(other.find_all_block(block))
                    if this_count < that_count:
                        differences["record"][block] = [this_count, that_count]
                        differences["count"] += that_count - this_count
        else:
            if not (self.contains_gap() or other.contains_gap()):
                differences["gap"] = False
            # should not use total_content = set(self) | set(other), because different direction means different block.
            total_content = self.block_list()
            for block in other.block_list():
                if block not in self.__block_set:
                    total_content.append(block)
            for block in total_content:
                if block != GapBlock:
                    this_count = len(self.find_all_block(block))
                    that_count = len(other.find_all_block(block))
                    if this_count != that_count:
                        differences["record"][block] = [this_count, that_count]
                        differences["count"] += abs(this_count - that_count)
        return differences

    def get_branching_points(self):
        if not self.__branching_points:
            if not self.block_downstream:
                self.build_adj_hash()
            # using adjacency shifts, more efficient than detecting repeats
            branching_points = []
            for check_block_name in self.block_downstream:
                # trifurcate
                for this_end in ["head", "tail"]:
                    if len(self.block_downstream[check_block_name][this_end]) > 1:
                        all_adj = sorted(set([pair
                                              for pair_set in self.block_downstream[check_block_name][this_end].values()
                                              for pair in pair_set]))
                        forward = []
                        reverse = []
                        for pair in all_adj:
                            if self[pair[0]].direction:
                                forward.append(pair)
                            else:
                                reverse.append(pair)
                        # combinations
                        for pair_i in forward:
                            for pair_j in reverse:
                                branching_points.append((pair_i, pair_j))
            self.__branching_points = tuple(branching_points)
        return list(self.__branching_points)

    def is_branching_points_detected(self):
        return False if self.__branching_points is None else True

    def get_isomers(self, **kwargs):
        go_to = -1
        new_isomers = []
        new_changes = []
        frozen_points = []
        next_isomer = self
        branching_points = self.get_branching_points()
        inherited_change = []
        new_iso_set = kwargs.get("recorded_isomers", set())
        new_iso_set.add(hash(self))
        while True:
            for br_id in range(len(branching_points)):
                (start, start_next), (end, end_next) = branching_points[br_id]
                candidate_isomer = next_isomer.inverted(start, end, build_bl_hash=False,
                                                        build_adj_hash=False, build_adjacent_blocks=False)
                iso_hash = hash(candidate_isomer)
                if iso_hash not in new_iso_set:
                    candidate_isomer.build_bl_hash()
                    candidate_isomer.build_adj_hash()
                    candidate_isomer.build_adjacent_blocks()
                    new_iso_set.add(iso_hash)
                    new_isomers.append(candidate_isomer)
                    new_changes.append(inherited_change + [(start, end)])
                    frozen_points.append(branching_points[br_id])
            go_to += 1
            if go_to >= len(new_isomers):
                break
            next_isomer = new_isomers[go_to]
            branching_points = next_isomer.get_branching_points()
            branching_points.remove(frozen_points[go_to])
            inherited_change = deepcopy(new_changes[go_to])
        return new_isomers, new_changes

    def event_from(self, other, fill_gap=True):
        changes = []
        next_other, block_changed = self.content_event_from(other, changes)
        if fill_gap and block_changed is None:
            gap_filled = self.gap_filling(next_other)
            while gap_filled:
                next_other, block_changed = self.content_event_from(next_other, changes)
                if block_changed:
                    gap_filled = self.gap_filling(next_other)
                else:
                    break
        changes += self.inversion_event_from(next_other, check_content=(block_changed is None))["changes"]
        return changes

    def is_gap_filled(self):
        return False if self.__gaps_filled is None else False

    def gap_filled(self):
        return self.__gaps_filled

    def gap_filling(self, reference):
        if self.__gaps_filled is None:
            self.__gaps_filled = []
        if self.contains_gap():
            changed = False
            update = True
            while update:
                update = False
                adj_dif_record = self.adjacency_differs_from(reference, gap_contains_gene=False)["record"]
                all_gap_loc = sorted(list(self.find_all_block(GapBlock)))
                len_gap_loc = len(all_gap_loc)
                con_dif = self.content_differs_from(reference, False)
                lack_block = False
                if con_dif["count"]:
                    for b, (this_count, that_count) in con_dif["record"].items():
                        if this_count == 0:  # this_count < that_count:
                            lack_block = True
                            break
                # 1. join contigs
                if len_gap_loc == 1:
                    left = (all_gap_loc[0] - 1) % self.__len
                    right = (all_gap_loc[0] + 1) % self.__len
                    new_adj = Adjacency(self[left], self[right])
                    if new_adj in adj_dif_record:
                        this_c, ref_c = adj_dif_record[new_adj]
                        if this_c < ref_c:
                            self.__gaps_filled.append({all_gap_loc[0] - 1: -1, all_gap_loc[0] + 1: 1})
                            self.delete(all_gap_loc[0])
                            self.__contains_gap = False
                            return True
                        # minimize change
                    if not lack_block:
                        self.__gaps_filled.append({all_gap_loc[0] - 1: -1, all_gap_loc[0] + 1: 1})
                        self.delete(all_gap_loc[0])
                        self.__contains_gap = False
                        return True
                else:
                    for go_gap in range(len_gap_loc):
                        for left_is_head in (-1, 1):
                            for go_next in range(go_gap, len_gap_loc):
                                for right_is_head in (-1, 1):
                                    if (go_next == go_gap and left_is_head == right_is_head) or (go_next == (
                                            go_gap + left_is_head) % len_gap_loc and left_is_head != right_is_head):
                                        pass
                                    else:
                                        left = (all_gap_loc[go_gap] + left_is_head) % self.__len
                                        right = (all_gap_loc[go_next] + right_is_head) % self.__len
                                        new_adj = Adjacency((-left_is_head) * self[left], right_is_head * self[right])
                                        if new_adj in adj_dif_record:
                                            this_c, ref_c = adj_dif_record[new_adj]
                                            if this_c < ref_c:
                                                self.__gaps_filled.append({left: left_is_head, right: right_is_head})
                                                if go_gap == go_next:
                                                    self.delete((left + 1) % self.__len)
                                                    if GapBlock in self:
                                                        self.__contains_gap = True
                                                    else:
                                                        self.__contains_gap = False
                                                elif left_is_head == -1 and right_is_head == -1:
                                                    self.inverse(left + 2, right)
                                                    self.delete((left + 1) % self.__len)
                                                    if GapBlock in self:
                                                        self.__contains_gap = True
                                                    else:
                                                        self.__contains_gap = False
                                                elif left_is_head == 1 and right_is_head == 1:
                                                    self.inverse(left, right - 2)
                                                    self.delete((right - 1) % self.__len)
                                                    if GapBlock in self:
                                                        self.__contains_gap = True
                                                    else:
                                                        self.__contains_gap = False
                                                else:
                                                    if left_is_head == 1 and right_is_head == -1:
                                                        intermediate = all_gap_loc[go_gap + 1]
                                                        self.inverse(left, intermediate - 1)
                                                        self.inverse(right, intermediate + 1)
                                                        self.delete(intermediate % self.__len)
                                                        if GapBlock in self:
                                                            self.__contains_gap = True
                                                        else:
                                                            self.__contains_gap = False
                                                    else:
                                                        intermediate = all_gap_loc[go_next + 1]
                                                        self.inverse(left + 2, intermediate - 1)
                                                        self.inverse(left + 2, left + 1 + intermediate - right)
                                                        self.delete((left + 1) % self.__len)
                                                        if GapBlock in self:
                                                            self.__contains_gap = True
                                                        else:
                                                            self.__contains_gap = False
                                                update = True
                                                changed = True
                                                break
                                if update:
                                    break
                            if update:
                                break
                        if update:
                            break
                if update:
                    continue
                # 2. fill blocks into gaps
                if con_dif["count"]:
                    for block, (this_count, that_count) in con_dif["record"].items():
                        if this_count < that_count:
                            for go_gap in range(len_gap_loc):
                                gap_loc = all_gap_loc[go_gap]
                                # insertion before the gap
                                gap_loc_minus_1 = (gap_loc - 1) % self.__len
                                left_b = self[gap_loc_minus_1]
                                for be_forward_block in (-1, 1):
                                    right_b = be_forward_block * block
                                    new_adj = Adjacency(left_b, right_b)
                                    if new_adj in adj_dif_record:
                                        this_c, ref_c = adj_dif_record[new_adj]
                                        if this_c < ref_c:
                                            self.__gaps_filled.append((gap_loc, right_b))
                                            self.insert(gap_loc, right_b)
                                            update = True
                                            changed = True
                                            break
                                # insertion after the gap
                                gap_loc_add_1 = (gap_loc + 1) % self.__len
                                right_b = self[gap_loc_add_1]
                                for be_forward_block in (-1, 1):
                                    left_b = be_forward_block * block
                                    new_adj = Adjacency(left_b, right_b)
                                    if new_adj in adj_dif_record:
                                        this_c, ref_c = adj_dif_record[new_adj]
                                        if this_c < ref_c:
                                            self.__gaps_filled.append((gap_loc_add_1, left_b))
                                            self.insert(gap_loc_add_1, left_b)
                                            update = True
                                            changed = True
                                            break
                                if update:
                                    break
                            if update:
                                break
            return changed
        else:
            return False

    def content_event_from(self, other, changes):
        con_events = self.content_differs_from(other)
        adj_events = self.adjacency_differs_from(other)
        if con_events["gap"]:
            return other, None
        else:
            changed = False
            while con_events["count"]:
                changed = False
                for block, (this_count, that_count) in con_events["record"].items():
                    if this_count > that_count:
                        associated_adj = []
                        for adj, (here_count, there_count) in adj_events["record"].items():
                            if block in adj and here_count > there_count:
                                associated_adj.append(adj)
                        cp_candidates = []
                        for adj in associated_adj:
                            candidate_loc = []
                            for loc in other.find_all_block(adj - block):
                                candidate_loc.append(loc)
                            for loc in candidate_loc:
                                directed_b = adj - other[loc]
                                temp = other.inserted(loc, directed_b)
                                cp_candidates.append([temp,
                                                      self.adjacency_differs_from(temp)["count"],
                                                      ("gain", ("insert", loc, directed_b))])
                                temp = other.inserted(loc + 1, directed_b)
                                cp_candidates.append([temp,
                                                      self.adjacency_differs_from(temp)["count"],
                                                      ("gain", ("insert", loc + 1, directed_b))])
                        if cp_candidates:
                            cp_candidates.sort(key=lambda x: x[1])
                            other, con_dif_count, single_event = cp_candidates[0]
                            changes.append(single_event)
                            changed = True
                    else:
                        for loc in other.find_all_block(block):
                            new_adj = Adjacency(other[(loc - 1) % len(other)], other[(loc + 1) % len(other)])
                            if new_adj in adj_events["record"]:
                                self_count, other_count = adj_events["record"][new_adj]
                                if self_count > other_count:
                                    changed = True
                                    changes.append(("loss", ("delete", loc, other[loc])))
                                    other = other.deleted(loc)
                                    break
                        if not changed:
                            candidates = [[other.deleted(test_loc), test_loc]
                                          for test_loc in other.find_all_block(block)]
                            candidates.sort(key=lambda x: self.adjacency_differs_from(x[0])["count"])
                            if candidates:
                                temp, loc = candidates[0]
                                changes.append(("loss", ("delete", loc, other[loc])))
                                other = temp
                                changed = True
                    if changed:
                        break
                if not changed:
                    sys.stdout.write("Content dif (" + str(con_events["count"]) + "): " +
                                     str(con_events["record"].items()) + "\n")
                    sys.stdout.write("from: " + str(other) + "\nto:  " + str(self) + "\n")
                    raise InterruptedError("Cannot find insertion position!")
                con_events = self.content_differs_from(other)
                adj_events = self.adjacency_differs_from(other)
            return other, changed

    # untested
    def breakpoints_from(self, other_chromosome, check_content=True):
        if check_content and self.content_differs_from(other_chromosome)["count"] != 0:
            raise NotImplementedError("Calculating breakpoints between unequal-content ones is not allowed!")
        else:
            # if equals
            if hash(self) == hash(other_chromosome):
                return 0
            else:
                adjacency_dif = self.adjacency_differs_from(other_chromosome)
                if adjacency_dif["gap"]:
                    sys.stdout.write("\tWarning: Calculating breakpoints between Gap-containing may be inaccurate!\n")
                return int(adjacency_dif["count"] / 2)

    def adjacency_event_from(self, other_chromosome, changes, check_content=True):
        if check_content and self.content_differs_from(other_chromosome)["count"] != 0:
            raise NotImplementedError("Calculating breakpoints between unequal-content ones is not allowed!")
        else:
            # if equals
            if hash(self) == hash(other_chromosome):
                return other_chromosome, False
            adjacency_dif = self.adjacency_differs_from(other_chromosome)
            for adj in adjacency_dif["record"].keys():
                for count_adj in range(abs(adjacency_dif["record"][adj][0] - adjacency_dif["record"][adj][1])):
                    changes.append(("adjacency_shifts", adj))
            if adjacency_dif["gap"]:
                sys.stdout.write("\tWarning: Enumerating adjacency event between Gap-containing may be inaccurate!\n")
            for adj in adjacency_dif["record"].keys():
                self_count, other_count = adjacency_dif["record"][adj]
                if self_count > other_count:
                    for k in range(self_count - other_count):
                        changes.append(["emerging_adjacency", [adj, None]])
            return other_chromosome, None

    # exhausted search for one shortest path
    def inversion_event_from(self, other_chromosome, check_content=True):
        time0 = time.time()
        changes = []
        if check_content and self.content_differs_from(other_chromosome)["count"] != 0:
            raise NotImplementedError("Inferring inversion events between unequal-content ones is not allowed!")
        else:
            # if equals
            if hash(self) == hash(other_chromosome):
                return other_chromosome, False
            adjacency_dif = self.adjacency_differs_from(other_chromosome)
            for adj in adjacency_dif["record"].keys():
                for count_adj in range(abs(adjacency_dif["record"][adj][0] - adjacency_dif["record"][adj][1])):
                    changes.append(("adjacency_shifts", adj))
            sys.stdout.write("\tBreakpoints: " + str(int(adjacency_dif["count"] / 2)) + "\n")
            if adjacency_dif["gap"]:
                sys.stdout.write("\tWarning: "
                                 "Inferring inversion events between Gap-containing ones is not allowed!\n"
                                 "\tReturning adjacency shifts instead (may also be inaccurate)!\n")
                # if gap, add adjacency to emerging adjacency
                for adj in adjacency_dif["record"].keys():
                    self_count, other_count = adjacency_dif["record"][adj]
                    if self_count > other_count:
                        for k in range(self_count - other_count):
                            changes.append(["emerging_adjacency", [adj, None]])
                sys.stdout.write("\tTotal inversion time: " + str(round(time.time() - time0, 4)) + "s\n")
                return {"changes": changes, "changed": None}

            other_isomers, iso_inv_changes = other_chromosome.get_isomers()
            for iso_id in range(len(other_isomers)):
                isomer = other_isomers[iso_id]
                if hash(self) == hash(isomer):
                    for start, end in iso_inv_changes[iso_id]:
                        changes.append(["inversion_with_iso",
                                        [Adjacency(isomer[start - 1], isomer[start]),
                                         Adjacency(isomer[end], isomer[end + 1]),
                                         start, end]])
                    sys.stdout.write("\n\tTotal inversion time: " + str(round(time.time() - time0, 4)) + "s\n")
                    return {"changes": changes, "changed": True}
            else:
                recorded_chromosome = set()
                recorded_chromosome.add(hash(other_chromosome))
                inv_changes = []
                next_step = [[(other_chromosome,
                               inv_changes,
                               (None, None),
                               self.adjacency_differs_from(other_chromosome))]]

                def adding_isomers(in_isomers, in_iso_inv_changes, in_inv_changes, in_adj_dif, in_pre_step):
                    for i_id in range(len(in_isomers)):
                        in_isomer = in_isomers[i_id]
                        if hash(in_isomer) not in recorded_chromosome:
                            h_inv_changes = deepcopy(in_inv_changes)
                            len_in = len(in_isomer)
                            for in_head, in_tail in in_iso_inv_changes[i_id]:
                                in_head_m1 = (in_head - 1) % len_in
                                in_tail_a1 = (in_tail + 1) % len_in
                                h_inv_changes.append(
                                    ["inversion_with_iso",
                                     [Adjacency(in_isomer[in_head_m1], in_isomer[in_head]),
                                      Adjacency(in_isomer[in_tail], in_isomer[in_tail_a1]),
                                      in_head_m1, in_tail_a1]])
                            in_pre_step.append([(in_isomer, h_inv_changes, (None, None), in_adj_dif)])

                adding_isomers(other_isomers, iso_inv_changes, inv_changes, adjacency_dif, next_step)
                for isomer in other_isomers:
                    recorded_chromosome.add(hash(isomer))
                while True:
                    time_r = time.time()
                    for step_id in range(len(next_step)):
                        next_step[step_id].sort(key=lambda x: (len(x[1]), x[3]["count"]))
                        for sub_step_id in range(len(next_step[step_id])):
                            recorded_chromosome.add(hash(next_step[step_id][sub_step_id][0]))
                    next_step.sort(key=lambda x: (len(x[0][1]), x[0][3]["count"]))
                    previous_step = next_step
                    next_step = []
                    sys.stdout.write("\tinherited combinations: " + str(sum([len(sub_s) for sub_s in previous_step]))
                                     + "; inversion sites: ")
                    sys.stdout.flush()
                    for sub_step in previous_step:
                        for prev_other_p, prev_inv_changes, prev_sites, prev_adj_change in sub_step:
                            # create a new subset
                            next_step.append([])
                            new_step_id = len(next_step) - 1
                            # find inversion sites
                            inversion_sites = set()
                            for adj in prev_adj_change["record"]:
                                self_count, other_count = prev_adj_change["record"][adj]
                                if other_count > self_count:
                                    for adj_loc in prev_other_p.adjacency_pointer[adj]:
                                        inversion_sites.add(adj_loc)
                            sys.stdout.write(" " + str(len(inversion_sites)))
                            sys.stdout.flush()

                            inv_site_list = sorted(list(inversion_sites))
                            for go_1, go_2 in combinations(list(range(len(inversion_sites))), 2):
                                site_1, site_2 = inv_site_list[go_1], inv_site_list[go_2]
                                if (site_1, site_2) != prev_sites:
                                    new_p = prev_other_p.inverted(site_1, site_2 - 1)
                                    hash_new = hash(new_p)
                                    if hash_new in recorded_chromosome:
                                        continue
                                    else:
                                        recorded_chromosome.add(hash_new)
                                        new_inv_changes = deepcopy(prev_inv_changes)
                                        len_pre = len(prev_other_p)
                                        site_1_m1 = (site_1 - 1) % len_pre
                                        site_2_m1 = (site_2 - 1) % len_pre
                                        adj_1 = Adjacency(prev_other_p[site_1_m1], prev_other_p[site_1])
                                        adj_2 = Adjacency(prev_other_p[site_2_m1], prev_other_p[site_2])
                                        inv_record = [adj_1, adj_2, site_1, site_2_m1]
                                        new_inv_changes.append(["inversion", inv_record])
                                        new_inv_changes.append(["inversion_with_iso", inv_record])
                                        new_inv_changes.append(
                                            ["emerging_adjacency",
                                             [Adjacency(prev_other_p[site_1_m1], -prev_other_p[site_2_m1]), site_1_m1]])
                                        new_inv_changes.append(
                                            ["emerging_adjacency",
                                             [Adjacency(-prev_other_p[site_1], prev_other_p[site_2]), site_2_m1]])
                                        if hash_new == hash(self):
                                            sys.stdout.write("; time: " + str(round(time.time() - time_r, 4)) + "s")
                                            if m_process:
                                                this_m = \
                                                    "; memory: " + \
                                                    str(round(m_process.memory_info().rss / 1024.0 / 1024 / 1024, 2)) \
                                                    + "G\n"
                                            else:
                                                this_m = "\n"
                                            sys.stdout.write(this_m)
                                            changes += new_inv_changes
                                            count_inv = len([x for x in new_inv_changes if x[0] == "inversion"])
                                            count_iso = len(
                                                [x for x in new_inv_changes if x[0] == "inversion_with_iso"])
                                            sys.stdout.write(
                                                "\tInversions: " + str(count_inv) + " + " +
                                                str(count_iso - count_inv) + "(iso)" + "\n"
                                                                                       "\tTotal inversion time: " + str(
                                                    round(time.time() - time0, 4)) + "s\n")
                                            sys.stdout.flush()
                                            return {"changes": changes, "changed": True}
                                        new_isomers, iso_inv_changes = new_p.get_isomers(recorded_isomers
                                                                                         =recorded_chromosome)
                                        for iso_id in range(len(new_isomers)):
                                            new_isomer = new_isomers[iso_id]
                                            len_io = len(new_isomer)
                                            if hash(new_isomer) == hash(self):
                                                for start, end in iso_inv_changes[iso_id]:
                                                    start_m1 = (start - 1) % len_io
                                                    end_a1 = (end + 1) % len_io
                                                    new_inv_changes.append(
                                                        ["inversion_with_iso",
                                                         [Adjacency(new_isomer[start_m1], new_isomer[start]),
                                                          Adjacency(new_isomer[end], new_isomer[end_a1]),
                                                          start, end]])
                                                changes += new_inv_changes
                                                sys.stdout.write(
                                                    "; time: " + str(round(time.time() - time_r, 4)) + "s")
                                                if m_process:
                                                    this_m = "; memory: " + str(round(
                                                        m_process.memory_info().rss / 1024.0 / 1024 / 1024, 2)) + "G\n"
                                                else:
                                                    this_m = "\n"
                                                sys.stdout.write(this_m)
                                                count_inv = len([x for x in new_inv_changes if x[0] == "inversion"])
                                                count_iso = len(
                                                    [x for x in new_inv_changes if x[0] == "inversion_with_iso"])
                                                sys.stdout.write(
                                                    "\tInversions: " + str(count_inv) + " + " +
                                                    str(count_iso - count_inv) + "(iso)" + "\n\tTotal inversion time: "
                                                    + str(round(time.time() - time0, 4)) + "s\n")
                                                sys.stdout.flush()
                                                return {"changes": changes, "changed": True}

                                        new_adj_dif = self.adjacency_differs_from(new_p)
                                        next_step[new_step_id].append(
                                            (new_p, new_inv_changes, (site_1, site_2), new_adj_dif))
                                        adding_isomers(new_isomers, iso_inv_changes, new_inv_changes, new_adj_dif,
                                                       next_step)
                            if not next_step[-1]:
                                del next_step[-1]
                    if m_process:
                        this_m = "; memory: " + \
                                 str(round(m_process.memory_info().rss / 1024.0 / 1024 / 1024, 2)) \
                                 + "G"
                    else:
                        this_m = ""
                    sys.stdout.write(this_m)
                    sys.stdout.flush()
                    sys.stdout.write("; time: " + str(round(time.time() - time_r, 4)) + "s\n")


class Genome:
    def __init__(self, label="", seq="", verbose=False):
        self.label = label
        self.__chromosomes = []
        """ converted as grimm format defined in our analysis """
        if "$" in seq:
            """ converted as standardized grimm format described in: http://grimm.ucsd.edu/GRIMM/grimm_instr.html """
            self.grimm_seq = ""
            for block in seq.replace("\n", "").split("$"):
                if block.strip():
                    self.grimm_seq += block + " $\n"
                    self.__chromosomes.append(Chromosome(block + " $"))
            self.grimm_seq = self.grimm_seq.strip()
        else:
            """ that no "$" in the end of a line stands for a single circular chromosome
            (multi-line context would be regarded as one single chromosome due to the limitation of applying to tsp) """
            self.grimm_seq = seq.replace("\n", " ")
            self.__chromosomes.append(Chromosome(self.grimm_seq))
        self.__chromosomes.sort()
        self.__verbose = verbose
        self.__phylip_alphabet = LEGAL_PHYLIP_ALPHABET
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
        return ">" + self.label + "\n" + self.grimm_seq

    def check_num_state(self):
        if max(list(self.content.values()) + [0]) > len(self.__phylip_alphabet):
            raise Exception("The total number of states is over " + str(len(self.__phylip_alphabet)))

    def update(self):
        self.content = {}
        self.__content_list = []
        self.adjacency = {}
        self.__adjacency_list = []
        self.__contains_gap = False
        self.__contains_telomere = False
        for chromosome in self.__chromosomes:
            for block in chromosome.block_list(signed=False):
                if block == GapBlock:
                    self.__contains_gap = True
                elif block == TerminalBlock:
                    self.__contains_telomere = True
                else:
                    self.content[block] = len(chromosome.find_all_block(block))
                    self.__content_list.append(block)
            for adj in chromosome.adjacency_unique_list():
                self.adjacency[adj] = len(chromosome.adjacency_pointer[adj])
                self.__adjacency_list.append(adj)
            for adj in sorted(chromosome.unknown_adjacency()):
                self.__gap_adjacency.add(adj)
        self.check_num_state()
        if self.__verbose:
            sys.stdout.write("Maximum number of states in " + self.label + " : " + str(max(self.content.values())))

    def __str__(self):
        return self.get_grimm()

    def __repr__(self):
        return str(type(self)) + "\n>" + self.get_grimm()

    def __eq__(self, other):
        if type(other) == type(self) and self.__chromosomes == other.chromosomes():
            return True
        else:
            return False

    def chromosomes(self):
        return list(self.__chromosomes)

    def content_list(self):
        return list(self.__content_list)

    def adjacency_list(self):
        return list(self.__adjacency_list)

    def phylip_alphabet(self):
        return self.__phylip_alphabet

    def contains_telomere(self):
        return bool(self.__contains_telomere)

    def contains_gap(self):
        return bool(self.__contains_gap)

    def unknown_adjacency(self):
        return set(self.__gap_adjacency)


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
            del self[:]
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
        fake_sites = "" if binary else "".join(
            [LEGAL_PHYLIP_ALPHABET[st] for st in range(max(state_set) + 1) if st not in state_set])
        head = " " + str(len(self)) + " " + str(len(block_list) + len(fake_sites)) + "\n"
        seqs = []
        for genome in self:
            # if there's a gap in the genome, add "-"(unknown) rather than "0"(doesn't exist) to missing block
            if gap_as_ambiguous and genome.contains_gap():
                alphabet = "-" + genome.phylip_alphabet()[1:]
            else:
                alphabet = genome.phylip_alphabet()
            max_index = 1 if binary else len(alphabet) - 1
            seq = "".join([alphabet[min(max_index, genome.content.get(b, 0))] for b in block_list]) + fake_sites
            seqs.append(genome.label + " " * 10 + seq + "\n")
        return head + "".join(seqs)

    def get_adjacency_phylip(self, gap_as_ambiguous=True, binary=True, reverse_code=False):
        state_set = self.adjacency_counter.state_set()
        adj_list = self.adjacency_counter.character_list()
        fake_sites = "" if binary else "".join(
            [LEGAL_PHYLIP_ALPHABET[st] for st in range(max(state_set) + 1) if st not in state_set])
        head = " " + str(len(self)) + " " + str(len(adj_list) + len(fake_sites)) + "\n"
        seqs = []
        alphabet = LEGAL_PHYLIP_ALPHABET
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
            max_index = 1 if binary else len(alphabet) - 1
            for adj in adj_list:
                if adj in genome.adjacency:
                    seq.append(alphabet[min(max_index, genome.adjacency[adj])])
                elif treat_gap:
                    start, end = adj.left, adj.right
                    missing_start = abs(start) in unknown_blocks
                    missing_end = abs(end) in unknown_blocks
                    start_t_break = Adjacency(start, GapBlock) in genome.unknown_adjacency()
                    break_t_end = Adjacency(GapBlock, end) in genome.unknown_adjacency()
                    if (missing_end or break_t_end) and (missing_start or start_t_break):
                        seq.append("-")
                    else:
                        seq.append(alphabet[0])
                else:
                    seq.append(alphabet[0])
            # pmag use the reverse code
            if binary and reverse_code:
                trans = {"1": "0", "0": "1", "-": "-"}
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
        return "".join([str(genome) + "\n" for genome in self])

    def __repr__(self):
        return str(type(self))+"\n" + "".join([str(genome) + "\n" for genome in self])
