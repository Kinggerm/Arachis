#! /usr/bin/env python
from arachis.genomeClass import *


class SubChromosome(Chromosome):
    def __init__(self, chromosome_list, build_bl_hash=True, build_adj_hash=True, build_adjacent_blocks=True, 
                 check_redundant_gaps=False, check_content=False, mirror=False):
        Chromosome.__init__(self, chromosome_str_or_list=chromosome_list, build_bl_hash=build_bl_hash,
                            build_adj_hash=build_adj_hash, build_adjacent_blocks=build_adjacent_blocks,
                            check_redundant_gaps=check_redundant_gaps, check_content=check_content,
                            mirror=mirror)
        self.circular = False


class Plastome(Chromosome):
    def __init__(self, chromosome_str_or_list, build_bl_hash=True, build_adj_hash=True, build_adjacent_blocks=True,
                 check_redundant_gaps=True, check_content=True, detect_ir=True, mirror=False):
        self.__build_bl_hash = build_bl_hash
        self.__build_adj_hash = build_adj_hash
        self.__build_adjacent_blocks = build_adjacent_blocks
        self.__get_ir = detect_ir
        self.__hash = None
        Chromosome.__init__(self, chromosome_str_or_list, build_bl_hash=build_bl_hash,
                            build_adj_hash=build_adj_hash, build_adjacent_blocks=build_adjacent_blocks,
                            check_redundant_gaps=check_redundant_gaps, check_content=check_content,
                            mirror=mirror)
        if check_content and not self.circular:
            raise TypeError("\n" + str(type(self)) + " is designed as a circular signed permutation!\n"
                            "Please use: arachis.genomeClass.Chromosome(\"" + str(chromosome_str_or_list) + "\")\n"
                            "Forbidden : arachis.plastomeClass.Plastome(\"" + str(chromosome_str_or_list) + "\")\n")
        self.ir = None
        self.lsc = None
        self.ssc = None
        self.ir_locations = None
        self.lsc_location = None
        self.ssc_location = None
        if detect_ir:
            self.update_ir()

    def __eq__(self, other):
        if (type(other) == type(self) or type(other) == Chromosome) and self.circular == other.circular:
            if self.circular:
                for match_start_id in self.find_all_block(other[0]):
                    order1 = self[match_start_id:] + self[:match_start_id]
                    new_match_start = (match_start_id + 1) % len(self)
                    order2 = [-bl for bl in (self[new_match_start:] + self[:new_match_start])[::-1]]
                    if list(other) == order1:
                        return True
                    elif list(other) == order2:
                        return True
            else:
                other_list = list(other)
                if list(self) == other_list or self.reverse_list() == other_list:
                    return True
        return False

    def __hash__(self, update=False):
        if self.__hash and not update:
            return self.__hash
        else:
            # direction does not matter here
            start_block = sorted(self, key=lambda x: x.name, reverse=True)[0]
            candidates = []
            # forward
            start_block_ids_1 = self.find_all_block(start_block)
            for start in start_block_ids_1:
                candidates.append(tuple(self[start:] + self[:start]))
            # reverse
            comp_rev = self.reverse_list()
            for b_id, candidate_b in enumerate(comp_rev):
                if candidate_b.name == start_block.name:
                    candidates.append(tuple(comp_rev[b_id:] + comp_rev[:b_id]))
            candidates.sort()
            self.__hash = hash(candidates[0])
        return self.__hash

    def super_class(self):
        return super(Plastome, self)

    def get_isomers(self, **kwargs):
        isomers, changes = super(Plastome, self).get_isomers(**kwargs)
        detect_ir = False if "detect_ir" not in kwargs else kwargs["detect_ir"]
        return [Plastome(isomer, build_bl_hash=self.__build_bl_hash,
                            build_adj_hash=self.__build_adj_hash,
                            build_adjacent_blocks=self.__build_adjacent_blocks, check_redundant_gaps=False,
                            check_content=False, detect_ir=detect_ir, mirror=True)
                for isomer in isomers], changes

    def gap_filling(self, reference):
        return super(Plastome, self).gap_filling(reference)

    def event_from(self, other, fill_gap=True):
        changes = []
        next_other, ir_changed, is_interrupted = self.naive_ir_event_from(other, changes)
        next_other, block_changed = self.content_event_from(next_other, changes)
        if fill_gap:
            gap_filled = self.gap_filling(next_other)
            while gap_filled:
                next_other, ir_changed, is_interrupted = self.naive_ir_event_from(next_other, changes)
                next_other, block_changed = self.content_event_from(next_other, changes)
                if ir_changed or block_changed:
                    gap_filled = self.gap_filling(next_other)
                else:
                    break
        next_other = Chromosome(next_other, build_bl_hash=True, build_adj_hash=True, build_adjacent_blocks=True,
                                check_redundant_gaps=False, check_content=False, mirror=True)
        changes += self.inversion_event_from(next_other, check_content=False)["changes"]
        return changes

    def inversion_event_from(self, other_cp_chromosome, check_content=True, detect_ir=False):
        return super(Plastome, self).inversion_event_from(other_cp_chromosome, check_content)

    def content_event_from(self, other, changes, detect_ir=False):
        new_cp, block_changed = super(Plastome, self).content_event_from(other, changes)
        new_cp = Plastome(new_cp, build_bl_hash=self.__build_bl_hash,
                             build_adj_hash=self.__build_adj_hash,
                             build_adjacent_blocks=self.__build_adjacent_blocks,
                             check_redundant_gaps=False, detect_ir=detect_ir, mirror=True)
        return new_cp, block_changed

    def naive_ir_event_from(self, other_cp, changes, verbose=False):
        len_self = len(self)
        interrupted = False
        initial_events = len(changes)
        if self.ir and not other_cp.ir:
            insert_adjacency = set()
            for loc in self.ir_locations:
                start, end, direction = loc["start"], loc["end"], loc["direction"]
                insert_adjacency.add(Adjacency(self[(start - direction)%len_self], self[(end + direction)%len_self]))
            for block_id in range(len(other_cp)):
                block_id_minus_1 = (block_id-1) % len(other_cp)
                if Adjacency(other_cp[block_id_minus_1], other_cp[block_id]) in insert_adjacency:
                    if Adjacency(other_cp[block_id_minus_1], self.ir[0]) in self.adjacency_pointer:
                        insert_part = list(self.ir)
                        other_cp = other_cp.inserted(block_id, insert_part, detect_ir=True)
                        changes.append(
                            ["lsc_contraction", ["insert", block_id, insert_part]])
                    elif Adjacency(other_cp[block_id_minus_1], self.ir[-2]) in self.adjacency_pointer:
                        insert_part = list(self.ir.reverse_list())
                        other_cp = other_cp.inserted(block_id, insert_part, detect_ir=True)
                        changes.append(
                            ["lsc_contraction", ["insert", block_id, insert_part]])
                    else:
                        interrupted = True
                        if verbose:
                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                            sys.stdout.write("\n1: Rearrangement around IR ... cannot simply interpret results!\n")
                    break
            else:
                interrupted = True
                if verbose:
                    sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                    sys.stdout.write("\n2: Rearrangement around IR ... cannot simply interpret results!\n")
        elif not self.ir and other_cp.ir:
            ir_location_1, ir_location_2 = other_cp.ir_locations
            ir_l1_s = ir_location_1["start"]
            ir_l2_s = ir_location_2["start"]
            ir_l1_e_a1 = (ir_location_1["end"] + 1) % len(other_cp)
            ir_l2_e_a1 = (ir_location_2["end"] + 1) % len(other_cp)
            next_other_1 = other_cp.deleted(ir_l1_s, ir_l1_e_a1, detect_ir=True)
            next_other_2 = other_cp.deleted(ir_l2_s, ir_l2_e_a1, detect_ir=True)
            adj_dif_1 = self.adjacency_differs_from(next_other_1)
            adj_dif_2 = self.adjacency_differs_from(next_other_2)
            if adj_dif_1["count"] <= adj_dif_2["count"]:
                changes.append(
                    ["lsc_expansion", ["delete", (ir_l1_s, ir_l1_e_a1), other_cp[ir_l1_s:ir_l1_e_a1]]])
                other_cp = next_other_1
            else:

                changes.append(
                    ["lsc_expansion", ["delete", (ir_l2_s, ir_l2_e_a1), other_cp[ir_l2_s:ir_l2_e_a1]]])
                other_cp = next_other_2
        elif self.ir and other_cp.ir:
            # detect IR movement according to boundary shifts
            # case 1/4: lsc/ssc expansion
            # if start side of other.ir should be moved into lsc/ssc
            for block_id in range(len(other_cp.ir)):
                if other_cp.ir[block_id] in self.ir:
                    # block_id > 0: some blocks at the start of other.ir was not in self.ir
                    if block_id > 0:
                        count_neither = 0
                        not_assigned = True
                        to_lsc = None
                        for test_id in range(block_id):
                            if other_cp.ir[test_id] in self.lsc:
                                if not_assigned:
                                    to_lsc = True
                                    not_assigned = False
                                else:
                                    if not to_lsc:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n3: Rearrangement between SSC and LSC ... cannot "
                                                             "simply interpret results!\n")
                            elif other_cp.ir[test_id] in self.ssc:
                                if not_assigned:
                                    to_lsc = False
                                    not_assigned = False
                                else:
                                    if to_lsc:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n4: Rearrangement between SSC and LSC ... cannot "
                                                             "simply interpret results!\n")
                            else:
                                if not_assigned:
                                    count_neither += 1
                        block_id -= count_neither
                        if block_id > 0:
                            if to_lsc:
                                lsc_start, lsc_end = other_cp.lsc_location["start"], other_cp.lsc_location["end"]
                                lsc_s_m1_mb = (lsc_start - 1 - block_id) % len(other_cp)
                                lsc_s_mb = (lsc_start - block_id) % len(other_cp)
                                lsc_e_a1_ab = (lsc_end + 1 + block_id) % len(other_cp)
                                lsc_e_a1 = (lsc_end + 1) % len(other_cp)
                                if Adjacency(other_cp[lsc_s_m1_mb], other_cp[lsc_start]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[lsc_s_m1_mb], other_cp[(lsc_start + 1) % len(other_cp)]) \
                                        in self.adjacency_pointer:
                                    changes.append(
                                        ["lsc_expansion",
                                         ["delete", (lsc_s_mb, lsc_start), other_cp[lsc_s_mb:lsc_start]]])
                                    other_cp = other_cp.deleted(lsc_s_mb, lsc_start, detect_ir=True)
                                elif Adjacency(other_cp[lsc_end], other_cp[lsc_e_a1_ab]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[(lsc_end - 1) % len(other_cp)], other_cp[lsc_e_a1_ab]) \
                                        in self.adjacency_pointer:
                                    changes.append(
                                        ["lsc_expansion",
                                         ["delete", (lsc_e_a1, lsc_e_a1_ab), other_cp[lsc_e_a1:lsc_e_a1_ab]]])
                                    other_cp = other_cp.deleted(lsc_e_a1, lsc_e_a1_ab, detect_ir=True)
                                else:
                                    other_cp_1 = other_cp.deleted(lsc_s_mb, lsc_start, detect_ir=True)
                                    other_cp_2 = other_cp.deleted(lsc_e_a1, lsc_e_a1_ab, detect_ir=True)
                                    old_dist = other_cp.adjacency_differs_from(self, False)["count"]
                                    event_dist_1 = other_cp_1.adjacency_differs_from(self, False)["count"]
                                    event_dist_2 = other_cp_2.adjacency_differs_from(self, False)["count"]
                                    if min(event_dist_1, event_dist_2) < old_dist:
                                        if event_dist_1 <= event_dist_2:
                                            changes.append(
                                                ["lsc_expansion",
                                                 ["delete", (lsc_s_mb, lsc_start), other_cp[lsc_s_mb:lsc_start]]])
                                            other_cp = other_cp_1
                                        else:
                                            changes.append(
                                                ["lsc_expansion",
                                                 ["delete", (lsc_e_a1, lsc_e_a1_ab), other_cp[lsc_e_a1:lsc_e_a1_ab]]])
                                            other_cp = other_cp_2
                                    else:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n5: Rearrangement around IR ... cannot "
                                                             "simply interpret results!\n")
                            else:
                                ssc_start, ssc_end = other_cp.ssc_location["start"], other_cp.ssc_location["end"]
                                ssc_s_m1_mb = (ssc_start - 1 - block_id) % len(other_cp)
                                ssc_s_mb = (ssc_start - block_id) % len(other_cp)
                                ssc_e_a1_ab = (ssc_end + 1 + block_id) % len(other_cp)
                                ssc_e_a1 = (ssc_end + 1) % len(other_cp)
                                if Adjacency(other_cp[ssc_s_m1_mb], other_cp[ssc_start]) in self.adjacency_pointer \
                                        or Adjacency(other_cp[ssc_s_m1_mb], other_cp[(ssc_start + 1) % len(other_cp)])\
                                        in self.adjacency_pointer:
                                    changes.append(["ssc_expansion",
                                                    ["delete", (ssc_s_mb, ssc_start), other_cp[ssc_s_mb:ssc_start]]])
                                    other_cp = other_cp.deleted(ssc_s_mb, ssc_start, detect_ir=True)
                                elif Adjacency(other_cp[ssc_end], other_cp[ssc_e_a1_ab]) in self.adjacency_pointer \
                                        or Adjacency(other_cp[(ssc_end - 1) % len(other_cp)], other_cp[ssc_e_a1_ab])\
                                        in self.adjacency_pointer:
                                    changes.append(["ssc_expansion",
                                                    ["delete", (ssc_e_a1, ssc_e_a1_ab), other_cp[ssc_e_a1:ssc_e_a1_ab]]])
                                    other_cp = other_cp.deleted(ssc_e_a1, ssc_e_a1_ab, detect_ir=True)
                                else:
                                    other_cp_1 = other_cp.deleted(ssc_s_mb, ssc_start, detect_ir=True)
                                    other_cp_2 = other_cp.deleted(ssc_e_a1, ssc_e_a1_ab, detect_ir=True)
                                    old_dist = other_cp.adjacency_differs_from(self, False)["count"]
                                    event_dist_1 = other_cp_1.adjacency_differs_from(self, False)["count"]
                                    event_dist_2 = other_cp_2.adjacency_differs_from(self, False)["count"]
                                    if min(event_dist_1, event_dist_2) < old_dist:
                                        if event_dist_1 <= event_dist_2:
                                            changes.append(
                                                ["ssc_expansion",
                                                 ["delete", (ssc_s_mb, ssc_start), other_cp[ssc_s_mb:ssc_start]]])
                                            other_cp = other_cp_1
                                        else:
                                            changes.append(
                                                ["ssc_expansion",
                                                 ["delete", (ssc_e_a1, ssc_e_a1_ab), other_cp[ssc_e_a1:ssc_e_a1_ab]]])
                                            other_cp = other_cp_2
                                    else:
                                        if verbose:
                                            interrupted = True
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n6: Rearrangement around IR ... cannot "
                                                             "simply interpret results!\n")
                    break
            # case 2/4: lsc/ssc expansion
            # if end side of other.ir should be moved into lsc/ssc
            if not other_cp.ir:
                if verbose:
                    sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                    sys.stdout.write("\n19: Intermediate structure have no/too-short ir!\n")
                return other_cp, bool(len(changes) - initial_events), True
            for block_id in range(1, len(other_cp.ir) + 1):
                if other_cp.ir[-block_id] in self.ir:
                    # block_id > 1: some blocks at the end of other.ir was not in self.ir
                    if block_id > 1:
                        count_neither = 0
                        not_assigned = True
                        to_lsc = None
                        for test_id in range(1, block_id):
                            if other_cp.ir[-test_id] in self.lsc:
                                if not_assigned:
                                    to_lsc = True
                                    not_assigned = False
                                else:
                                    if not to_lsc:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n7: Rearrangement between SSC and LSC ... cannot "
                                                             "simply interpret results!\n")
                            elif other_cp.ir[-test_id] in self.ssc:
                                if not_assigned:
                                    to_lsc = False
                                    not_assigned = False
                                else:
                                    if to_lsc:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n8: Rearrangement between SSC and LSC ... cannot "
                                                             "simply interpret results!\n")
                            else:
                                if not_assigned:
                                    count_neither += 1
                        block_id -= count_neither
                        if block_id > 1:
                            if to_lsc:
                                lsc_start, lsc_end = other_cp.lsc_location["start"], other_cp.lsc_location["end"]
                                lsc_s_mb_a1 = (lsc_start - block_id + 1) % len(other_cp)
                                lsc_s_mb = (lsc_start - block_id) % len(other_cp)
                                lsc_e_a1 = (lsc_end + 1) % len(other_cp)
                                lsc_e_ab = (lsc_end + block_id) % len(other_cp)
                                if Adjacency(other_cp[lsc_s_mb], other_cp[lsc_start]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[lsc_s_mb], other_cp[(lsc_start + 1) % len(other_cp)]) \
                                        in self.adjacency_pointer:
                                    changes.append(
                                        ["lsc_expansion",
                                         ["delete", (lsc_s_mb_a1, lsc_start), other_cp[lsc_s_mb_a1:lsc_start]]])
                                    other_cp = other_cp.deleted(lsc_s_mb_a1, lsc_start, detect_ir=True)
                                elif Adjacency(other_cp[lsc_end], other_cp[lsc_e_ab]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[(lsc_end - 1) % len(other_cp)], other_cp[lsc_e_ab])\
                                        in self.adjacency_pointer:
                                    changes.append(
                                        ["lsc_expansion",
                                         ["delete", (lsc_e_a1, lsc_e_ab), other_cp[lsc_e_a1:lsc_e_ab]]])
                                    other_cp = other_cp.deleted(lsc_e_a1, lsc_e_ab, detect_ir=True)
                                else:
                                    other_cp_1 = other_cp.deleted(lsc_s_mb_a1, lsc_start, detect_ir=True)
                                    other_cp_2 = other_cp.deleted(lsc_e_a1, lsc_e_ab, detect_ir=True)
                                    old_dist = other_cp.adjacency_differs_from(self, False)["count"]
                                    event_dist_1 = other_cp_1.adjacency_differs_from(self, False)["count"]
                                    event_dist_2 = other_cp_2.adjacency_differs_from(self, False)["count"]
                                    if min(event_dist_1, event_dist_2) < old_dist:
                                        if event_dist_1 <= event_dist_2:
                                            changes.append(
                                                ["lsc_expansion",
                                                 ["delete", (lsc_s_mb_a1, lsc_start), other_cp[lsc_s_mb_a1:lsc_start]]])
                                            other_cp = other_cp_1
                                        else:
                                            changes.append(
                                                ["lsc_expansion",
                                                 ["delete", (lsc_e_a1, lsc_e_ab), other_cp[lsc_e_a1:lsc_e_ab]]])
                                            other_cp = other_cp_2
                                    else:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n9: Rearrangement around IR ... cannot "
                                                             "simply interpret results!\n")
                            else:
                                ssc_start, ssc_end = other_cp.ssc_location["start"], other_cp.ssc_location["end"]
                                # if delete the extra ir, the new adjacency would match self.adjacency
                                ssc_s_mb = (ssc_start - block_id) % len(other_cp)
                                ssc_s_mb_a1 = (ssc_start - block_id + 1) % len(other_cp)
                                ssc_e_ab = (ssc_end + block_id) % len(other_cp)
                                ssc_e_a1 = (ssc_end + 1) % len(other_cp)
                                if Adjacency(other_cp[ssc_s_mb], other_cp[ssc_start]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[ssc_s_mb], other_cp[(ssc_start + 1) % len(other_cp)])\
                                        in self.adjacency_pointer:
                                    changes.append(
                                        ["ssc_expansion",
                                         ["delete", (ssc_s_mb_a1, ssc_start), other_cp[ssc_s_mb_a1:ssc_start]]])
                                    other_cp = other_cp.deleted(ssc_s_mb_a1, ssc_start, detect_ir=True)
                                elif Adjacency(other_cp[ssc_end], other_cp[ssc_e_ab]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[(ssc_end - 1) % len(other_cp)], other_cp[ssc_e_ab])\
                                        in self.adjacency_pointer:
                                    changes.append(
                                        ["ssc_expansion",
                                        ["delete", (ssc_e_a1, ssc_e_ab), other_cp[ssc_e_a1:ssc_e_ab]]])
                                    other_cp = other_cp.deleted(ssc_e_a1, ssc_e_ab, detect_ir=True)
                                else:
                                    other_cp_1 = other_cp.deleted(ssc_s_mb_a1, ssc_start, detect_ir=True)
                                    other_cp_2 = other_cp.deleted(ssc_e_a1, ssc_e_ab, detect_ir=True)
                                    old_dist = other_cp.adjacency_differs_from(self, False)["count"]
                                    event_dist_1 = other_cp_1.adjacency_differs_from(self, False)["count"]
                                    event_dist_2 = other_cp_2.adjacency_differs_from(self, False)["count"]
                                    if min(event_dist_1, event_dist_2) < old_dist:
                                        if event_dist_1 <= event_dist_2:
                                            changes.append(
                                                ["ssc_expansion",
                                                 ["delete", (ssc_s_mb_a1, ssc_start), other_cp[ssc_s_mb_a1:ssc_start]]])
                                            other_cp = other_cp_1
                                        else:
                                            changes.append(
                                                ["ssc_expansion",
                                                 ["delete", (ssc_e_a1, ssc_e_ab), other_cp[ssc_e_a1:ssc_e_ab]]])
                                            other_cp = other_cp_2
                                    else:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n10: Rearrangement around IR ... cannot "
                                                             "simply interpret results!\n")
                    break
            # case 3/4: lsc/ssc contraction
            if not other_cp.ir:
                if verbose:
                    sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                    sys.stdout.write("\n20: Intermediate structure have no/too-short ir!\n")
                return other_cp, bool(len(changes) - initial_events), True
            for block_id in range(len(self.ir)):
                if self.ir[block_id] in other_cp.ir:
                    if block_id > 0:
                        count_neither = 0
                        not_assigned = True
                        to_lsc = None
                        for test_id in range(block_id):
                            if self.ir[test_id] in other_cp.lsc:
                                if not_assigned:
                                    to_lsc = True
                                    not_assigned = False
                                else:
                                    if not to_lsc:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n11: Rearrangement between SSC and LSC ... cannot "
                                                             "simply interpret results!\n")
                            elif self.ir[test_id] in other_cp.ssc:
                                if not_assigned:
                                    to_lsc = False
                                    not_assigned = False
                                else:
                                    if to_lsc:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n12: Rearrangement between SSC and LSC ... cannot "
                                                             "simply interpret results!\n")
                            else:
                                if not_assigned:
                                    count_neither += 1
                        block_id -= count_neither
                        if block_id > 0:
                            if to_lsc:
                                lsc_start, lsc_end = other_cp.lsc_location["start"], other_cp.lsc_location["end"]
                                lsc_s_m1 = (lsc_start - 1) % len(other_cp)
                                lsc_s_ab = (lsc_start + block_id) % len(other_cp)
                                lsc_e_a1 = (lsc_end + 1) % len(other_cp)
                                lsc_e_mb_a1 = (lsc_end - block_id + 1) % len(other_cp)
                                if Adjacency(other_cp[lsc_s_m1], other_cp[lsc_start]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[lsc_s_m1], other_cp[(lsc_start + 1) % len(other_cp)])\
                                        in self.adjacency_pointer:
                                    insert_part = self.reverse_list(other_cp[lsc_start: lsc_s_ab])
                                    changes.append(["lsc_contraction", ["insert", lsc_e_a1, insert_part]])
                                    other_cp = other_cp.inserted(lsc_e_a1, insert_part, detect_ir=True)
                                elif Adjacency(other_cp[lsc_end], other_cp[lsc_e_a1]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[(lsc_end - 1) % len(other_cp)], other_cp[lsc_e_a1])\
                                        in self.adjacency_pointer:
                                    insert_part = self.reverse_list(other_cp[lsc_e_mb_a1: lsc_e_a1])
                                    changes.append(["lsc_contraction", ["insert", lsc_start, insert_part]])
                                    other_cp = other_cp.inserted(lsc_start, insert_part, detect_ir=True)
                                else:
                                    insert_part1 = self.reverse_list(other_cp[lsc_start: lsc_s_ab])
                                    insert_part2 = self.reverse_list(other_cp[lsc_e_mb_a1: lsc_e_a1])
                                    other_cp_1 = other_cp.inserted(lsc_e_a1, insert_part1, detect_ir=True)
                                    other_cp_2 = other_cp.inserted(lsc_start, insert_part2, detect_ir=True)
                                    old_dist = other_cp.adjacency_differs_from(self, False)["count"]
                                    event_dist_1 = other_cp_1.adjacency_differs_from(self, False)["count"]
                                    event_dist_2 = other_cp_2.adjacency_differs_from(self, False)["count"]
                                    if min(event_dist_1, event_dist_2) < old_dist:
                                        if event_dist_1 <= event_dist_2:
                                            changes.append(["lsc_contraction", ["insert", lsc_e_a1, insert_part1]])
                                            other_cp = other_cp_1
                                        else:
                                            changes.append(["lsc_contraction", ["insert", lsc_start, insert_part2]])
                                            other_cp = other_cp_2
                                    else:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n13: Rearrangement around IR ... cannot "
                                                             "simply interpret results!\n")
                            else:
                                ssc_start, ssc_end = other_cp.ssc_location["start"], other_cp.ssc_location["end"]
                                ssc_s_m1 = (ssc_start - 1) % len(other_cp)
                                ssc_s_ab = (ssc_start + block_id) % len(other_cp)
                                ssc_e_mb_a1 = (ssc_end - block_id + 1) % len(other_cp)
                                ssc_e_a1 = (ssc_end + 1) % len(other_cp)
                                if Adjacency(other_cp[ssc_s_m1], other_cp[ssc_start]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[ssc_s_m1], other_cp[(ssc_start + 1) % len(self)])\
                                        in self.adjacency_pointer:
                                    insert_part = self.reverse_list(other_cp[ssc_start: ssc_s_ab])
                                    changes.append(["ssc_contraction", ["insert", ssc_e_a1, insert_part]])
                                    other_cp = other_cp.inserted(ssc_e_a1, insert_part, detect_ir=True)
                                elif Adjacency(other_cp[ssc_end], other_cp[ssc_e_a1]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[(ssc_end - 1) % len(other_cp)], other_cp[ssc_e_a1])\
                                        in self.adjacency_pointer:
                                    insert_part = self.reverse_list(other_cp[ssc_e_mb_a1: ssc_e_a1])
                                    changes.append(["ssc_contraction", ["insert", ssc_start, insert_part]])
                                    other_cp = other_cp.inserted(ssc_start, insert_part, detect_ir=True)
                                else:
                                    insert_part1 = self.reverse_list(other_cp[ssc_start: ssc_s_ab])
                                    insert_part2 = self.reverse_list(other_cp[ssc_e_mb_a1: ssc_e_a1])
                                    other_cp_1 = other_cp.inserted(ssc_e_a1, insert_part1, detect_ir=True)
                                    other_cp_2 = other_cp.inserted(ssc_start, insert_part2, detect_ir=True)
                                    old_dist = other_cp.adjacency_differs_from(self, False)["count"]
                                    event_dist_1 = other_cp_1.adjacency_differs_from(self, False)["count"]
                                    event_dist_2 = other_cp_2.adjacency_differs_from(self, False)["count"]
                                    if min(event_dist_1, event_dist_2) < old_dist:
                                        if event_dist_1 <= event_dist_2:
                                            changes.append(["ssc_contraction", ["insert", ssc_e_a1, insert_part1]])
                                            other_cp = other_cp_1
                                        else:
                                            changes.append(["ssc_contraction", ["insert", ssc_start, insert_part2]])
                                            other_cp = other_cp_2
                                    else:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n14: Rearrangement around IR ... cannot "
                                                             "simply interpret results!\n")
                    break
            # case 4/4: lsc/ssc contraction
            if not other_cp.ir:
                if verbose:
                    sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                    sys.stdout.write("\n21: Intermediate structure have no/too-short ir!\n")
                return other_cp, bool(len(changes) - initial_events), True
            for block_id in range(1, len(self.ir) + 1):
                if self.ir[-block_id] in other_cp.ir:
                    if block_id > 1:
                        count_neither = 0
                        not_assigned = True
                        to_lsc = None
                        for test_id in range(1, block_id):
                            if self.ir[-test_id] in other_cp.lsc:
                                if not_assigned:
                                    to_lsc = True
                                    not_assigned = False
                                else:
                                    if not to_lsc:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n15: Rearrangement between SSC and LSC ... cannot "
                                                             "simply interpret results!\n")
                            elif self.ir[-test_id] in other_cp.ssc:
                                if not_assigned:
                                    to_lsc = False
                                    not_assigned = False
                                else:
                                    if to_lsc:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n16: Rearrangement between SSC and LSC ... cannot "
                                                             "simply interpret results!\n")
                            else:
                                if not_assigned:
                                    count_neither += 1
                        block_id -= count_neither
                        if block_id > 1:
                            if to_lsc:
                                lsc_start, lsc_end = other_cp.lsc_location["start"], other_cp.lsc_location["end"]
                                lsc_s_m1 = (lsc_start - 1) % len(other_cp)
                                lsc_s_ab_m1 = (lsc_start + block_id - 1) % len(other_cp)
                                lsc_e_a1 = (lsc_end + 1) % len(other_cp)
                                lsc_e_mb_a2 = (lsc_end - block_id + 2) % len(other_cp)
                                if Adjacency(other_cp[lsc_s_m1], other_cp[lsc_start]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[lsc_s_m1], other_cp[(lsc_start + 1) % len(other_cp)])\
                                        in self.adjacency_pointer:
                                    insert_part = self.reverse_list(other_cp[lsc_start: lsc_s_ab_m1])
                                    other_cp = other_cp.inserted(lsc_e_a1, insert_part, detect_ir=True)
                                    changes.append(["lsc_contraction", ["insert", lsc_e_a1, insert_part]])
                                elif Adjacency(other_cp[lsc_end], other_cp[lsc_e_a1]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[(lsc_end - 1) % len(other_cp)], other_cp[lsc_e_a1])\
                                        in self.adjacency_pointer:
                                    insert_part = self.reverse_list(other_cp[lsc_e_mb_a2: lsc_e_a1])
                                    changes.append(["lsc_contraction", ["insert", lsc_start, insert_part]])
                                    other_cp = other_cp.inserted(lsc_start, insert_part, detect_ir=True)
                                else:
                                    insert_part1 = self.reverse_list(other_cp[lsc_start: lsc_s_ab_m1])
                                    insert_part2 = self.reverse_list(other_cp[lsc_e_mb_a2: lsc_e_a1])
                                    other_cp_1 = other_cp.inserted(lsc_e_a1, insert_part1, detect_ir=True)
                                    other_cp_2 = other_cp.inserted(lsc_start, insert_part2, detect_ir=True)
                                    old_dist = other_cp.adjacency_differs_from(self, False)["count"]
                                    event_dist_1 = other_cp_1.adjacency_differs_from(self, False)["count"]
                                    event_dist_2 = other_cp_2.adjacency_differs_from(self, False)["count"]
                                    if min(event_dist_1, event_dist_2) < old_dist:
                                        if event_dist_1 <= event_dist_2:
                                            changes.append(["lsc_contraction", ["insert", lsc_e_a1, insert_part1]])
                                            other_cp = other_cp_1
                                        else:
                                            changes.append(["lsc_contraction", ["insert", lsc_start, insert_part2]])
                                            other_cp = other_cp_2
                                    else:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n17: Rearrangement around IR ... cannot "
                                                             "simply interpret results!\n")
                            else:
                                ssc_start, ssc_end = other_cp.ssc_location["start"], other_cp.ssc_location["end"]
                                ssc_s_m1 = (ssc_start - 1) % len(other_cp)
                                ssc_s_ab_m1 = (ssc_start + block_id - 1) % len(other_cp)
                                ssc_e_a1 = (ssc_end + 1) % len(other_cp)
                                ssc_e_mb_a2 = (ssc_end - block_id + 2) % len(other_cp)
                                if Adjacency(other_cp[ssc_s_m1], other_cp[ssc_start]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[ssc_s_m1], other_cp[(ssc_start + 1) % len(other_cp)])\
                                        in self.adjacency_pointer:
                                    insert_part = self.reverse_list(other_cp[ssc_start: ssc_s_ab_m1])
                                    changes.append(["ssc_contraction", ["insert", ssc_e_a1, insert_part]])
                                    other_cp = other_cp.inserted(ssc_e_a1, insert_part, detect_ir=True)
                                elif Adjacency(other_cp[ssc_end], other_cp[ssc_e_a1]) in self.adjacency_pointer\
                                        or Adjacency(other_cp[(ssc_end - 1) % len(other_cp)], other_cp[ssc_e_a1])\
                                        in self.adjacency_pointer:
                                    insert_part = self.reverse_list(other_cp[ssc_e_mb_a2: ssc_e_a1])
                                    changes.append(["ssc_contraction", ["insert", ssc_start, insert_part]])
                                    other_cp = other_cp.inserted(ssc_start, insert_part, detect_ir=True)
                                else:
                                    insert_part1 = self.reverse_list(other_cp[ssc_start: ssc_s_ab_m1])
                                    insert_part2 = self.reverse_list(other_cp[ssc_e_mb_a2: ssc_e_a1])
                                    other_cp_1 = other_cp.inserted(ssc_e_a1, insert_part1, detect_ir=True)
                                    other_cp_2 = other_cp.inserted(ssc_start, insert_part2, detect_ir=True)
                                    old_dist = other_cp.adjacency_differs_from(self, False)["count"]
                                    event_dist_1 = other_cp_1.adjacency_differs_from(self, False)["count"]
                                    event_dist_2 = other_cp_2.adjacency_differs_from(self, False)["count"]
                                    if min(event_dist_1, event_dist_2) < old_dist:
                                        if event_dist_1 <= event_dist_2:
                                            changes.append(["ssc_contraction", ["insert", ssc_e_a1, insert_part1]])
                                            other_cp = other_cp_1
                                        else:
                                            changes.append(["ssc_contraction", ["insert", ssc_start, insert_part2]])
                                            other_cp = other_cp_2
                                    else:
                                        interrupted = True
                                        if verbose:
                                            sys.stdout.write("\nfrom: " + str(other_cp) + "\nto:  " + str(self))
                                            sys.stdout.write("\n18: Rearrangement around IR ... cannot "
                                                             "simply interpret results!\n")
                    break
        return other_cp, bool(len(changes) - initial_events), interrupted

    def inverse(self, start_b, end_b, update_ir=True):
        super(Plastome, self).inverse(start_b, end_b)
        if update_ir:
            self.update_ir()

    def inverted(self, start_b, end_b, **build):
        build_bl_hash = self.__build_bl_hash if "build_bl_hash" not in build else build["build_bl_hash"]
        build_adj_hash = self.__build_adj_hash if "build_adj_hash" not in build else build["build_adj_hash"]
        build_adj_b = self.__build_adjacent_blocks if "build_adjacent_blocks" not in build \
            else build["build_adjacent_blocks"]
        detect_ir = False if "detect_ir" not in build else build["detect_ir"]
        return Plastome(super(Plastome, self).inverted(start_b, end_b, **{"build_bl_hash": False,
                                                                          "build_adj_hash": False,
                                                                          "build_adjacent_blocks": False}),
                        build_bl_hash=build_bl_hash, build_adj_hash=build_adj_hash, build_adjacent_blocks=build_adj_b,
                        check_redundant_gaps=False, detect_ir=detect_ir, mirror=True)

    def insert(self, at, item, update_ir=True):
        super(Plastome, self).insert(at, item)
        if update_ir:
            self.update_ir()

    def inserted(self, at, item, **build):
        build_bl_hash = self.__build_bl_hash if "build_bl_hash" not in build else build["build_bl_hash"]
        build_adj_hash = self.__build_adj_hash if "build_adj_hash" not in build else build["build_adj_hash"]
        build_adj_b = self.__build_adjacent_blocks if "build_adjacent_blocks" not in build \
            else build["build_adjacent_blocks"]
        detect_ir = False if "detect_ir" not in build else build["detect_ir"]
        return Plastome(super(Plastome, self).inserted(at, item, **{"build_bl_hash": False, "build_adj_hash": False,
                                                                    "build_adjacent_blocks": False}),
                        build_bl_hash=build_bl_hash, build_adj_hash=build_adj_hash, build_adjacent_blocks=build_adj_b,
                        check_redundant_gaps=False, detect_ir=detect_ir)

    def delete(self, start, stop=None, update_ir=True):
        super(Plastome, self).delete(start, stop)
        if update_ir:
            self.update_ir()

    def deleted(self, start, stop=None, **build):
        stop = start + 1 if not stop else stop
        build_bl_hash = self.__build_bl_hash if "build_bl_hash" not in build else build["build_bl_hash"]
        build_adj_hash = self.__build_adj_hash if "build_adj_hash" not in build else build["build_adj_hash"]
        build_adj_b = self.__build_adjacent_blocks if "build_adjacent_blocks" not in build \
            else build["build_adjacent_blocks"]
        detect_ir = False if "detect_ir" not in build else build["detect_ir"]
        return Plastome(list(self)[:start] + list(self)[stop:], build_bl_hash=build_bl_hash,
                        build_adj_hash=build_adj_hash, build_adjacent_blocks=build_adj_b,
                        check_redundant_gaps=False, detect_ir=detect_ir)

    def scj_broken(self, fission_before_bl):
        sys.stdout.write("Warning: Plastome does not support fission operation! "
                         "Recreate result as Chromosome!\n")
        return super(Plastome, self).scj_broken(fission_before_bl)

    def scj_join(self, other_chromosome=None):
        sys.stdout.write("Warning: Plastome does not support fussion operation! "
                         "Recreate result as Chromosome!\n")
        return super(Plastome, self).scj_join(other_chromosome)

    def update_ir(self, min_repeat_length=3):
        # assume the longest is ir
        all_repeats = self.detect_list_repeats(min_repeat_length=min_repeat_length)
        len_self = len(self)
        if all_repeats:
            ir_locations = sorted(all_repeats[0], key=lambda x: -x["direction"])
            if len(ir_locations) != 2 or ir_locations[0]["direction"] == ir_locations[1]["direction"]:
                pass
            else:
                self.ir_locations = ir_locations
                ir_start, ir_end, ir_j = ir_locations[0]["start"], ir_locations[0]["end"], ir_locations[0]["direction"]
                self.ir = SubChromosome(self[ir_start: (ir_end + ir_j) % len_self: ir_j],
                                        build_bl_hash=self.__build_bl_hash,
                                        build_adj_hash=False, build_adjacent_blocks=False, check_redundant_gaps=False)
                # as sorted, direct1==1 and direct2==-1
                start1, end1, direct1 = ir_locations[0]["start"], ir_locations[0]["end"], ir_locations[0]["direction"]
                start2, end2, direct2 = ir_locations[1]["start"], ir_locations[1]["end"], ir_locations[1]["direction"]
                sc = [[SubChromosome(self[(end1 + 1) % len_self: end2], build_bl_hash=self.__build_bl_hash,
                                     build_adj_hash=False, build_adjacent_blocks=False, check_redundant_gaps=False),
                       ((end1 + 1) % len(self), (end2 - 1) % len(self))],
                      [SubChromosome(self[(start2 + 1) % len_self: start1], build_bl_hash=self.__build_bl_hash,
                                     build_adj_hash=False, build_adjacent_blocks=False, check_redundant_gaps=False),
                       ((start2 + 1) % len(self), (start1 - 1) % len(self))]]
                sc.sort(key=lambda x: len(x[0]))
                [self.ssc, ssc_location], [self.lsc, lsc_location] = sc
                self.ssc_location = {"start": ssc_location[0], "end": ssc_location[1]}
                self.lsc_location = {"start": lsc_location[0], "end": lsc_location[1]}
    # more robust than "branching point algorithm" when gaps involved
    def detect_list_repeats(self, outer_list=None, min_repeat_length=3, circular=True, word_size=3):
        word_size = min(word_size, min_repeat_length)
        check_len = len(self) if outer_list is None else len(outer_list)
        if check_len < min_repeat_length:
            return []
        if circular:
            here_seq = list(self) + self[:word_size-1] if outer_list is None \
                else outer_list + outer_list[:word_size-1]
        else:
            here_seq = list(self) if outer_list is None else outer_list
        here_seq_r = self.reverse_list(here_seq)
        here_seq_length = len(here_seq)
        """initialization"""
        words_to_index = {}
        index_to_words = {}

        def add_to_words(add_index, this_forward, this_reverse):
            if this_forward in words_to_index:
                words_to_index[this_forward].add((add_index, 1))
                words_to_index[this_reverse].add((add_index, -1))
            else:
                words_to_index[this_forward] = {(add_index, 1)}
                if this_reverse == this_forward:
                    words_to_index[this_reverse].add((add_index, -1))
                else:
                    words_to_index[this_reverse] = {(add_index, -1)}

        for i in range(0, here_seq_length):
            forward_s = tuple(here_seq[i:i + word_size])
            reverse_s = tuple(here_seq_r[here_seq_length - i - word_size: here_seq_length - i])
            add_to_words(i, forward_s, reverse_s)
            index_to_words[i] = forward_s

        """find repeats"""
        repeat_indices = set()
        for repeat_tuples in words_to_index.values():
            if len(repeat_tuples) >= 2:
                for repeat_index in repeat_tuples:
                    repeat_indices.add(repeat_index[0])
        repeat_indices = sorted(list(repeat_indices))
        repeats = []
        points_to_repeats = {}
        last_connection = set()
        len_indices = len(repeat_indices)
        for i in range(len_indices):
            this_index = repeat_indices[i]
            this_word = index_to_words[this_index]
            this_connection = words_to_index[this_word]

            """part 1: dealing with old connection"""
            # Loop 1: find repeats_to_stop
            # Loop 2: delete the pointers pointing to the stopped repeats
            # Loop 3: update the pointers
            # Loop 4: add new pointers if shorter repeats should be continued with less alias
            # Loop 5: update the repeats according to pointers
            repeats_to_stop = {}
            kinds_del_from_points = set()
            for one_connection in last_connection:
                here_id, direction_trans = one_connection
                candidate_new_connect = (here_id + direction_trans, direction_trans)
                if candidate_new_connect not in this_connection:
                    # if one_connection in points_to_repeats:
                    for group, repeat_num in points_to_repeats[one_connection]:
                        if group in repeats_to_stop:
                            repeats_to_stop[group][repeat_num] = one_connection
                        else:
                            repeats_to_stop[group] = {repeat_num: one_connection}
                        kinds_del_from_points.add(group)

            for group in kinds_del_from_points:
                for now_start, now_go_to, now_direction in repeats[group]:
                    connection_del_from_points = (now_go_to - (word_size - 1) * (now_direction == 1),
                                                  now_direction)
                    if connection_del_from_points in points_to_repeats:
                        count_this_group = 0
                        while count_this_group < len(points_to_repeats[connection_del_from_points]):
                            if points_to_repeats[connection_del_from_points][count_this_group][0] == group:
                                del points_to_repeats[connection_del_from_points][count_this_group]
                            else:
                                count_this_group += 1
                        if not len(points_to_repeats[connection_del_from_points]):
                            del points_to_repeats[connection_del_from_points]

            for one_connection in last_connection:
                here_id, direction_trans = one_connection
                candidate_new_connect = (here_id + direction_trans, direction_trans)
                if candidate_new_connect in this_connection:
                    if one_connection in points_to_repeats:
                        points_to_repeats[candidate_new_connect] = []
                        for one_repeat_id in points_to_repeats[one_connection]:
                            points_to_repeats[candidate_new_connect].append(one_repeat_id)
                        del points_to_repeats[one_connection]

            for group in repeats_to_stop:
                if len(repeats[group]) - len(repeats_to_stop[group]) >= 2:
                    repeat_to_be_continued = False
                    for repeat_num in range(len(repeats[group])):
                        if repeat_num not in repeats_to_stop[group]:
                            start_id, go_to_id, go_to_direction = repeats[group][repeat_num]
                            new_connect = (go_to_id - (word_size - 1) * (go_to_direction == 1) + go_to_direction,
                                           go_to_direction)
                            if new_connect in this_connection and new_connect not in points_to_repeats:
                                repeat_to_be_continued = True
                                break
                    if repeat_to_be_continued:
                        repeats.append([])
                        for inside_repeat_num in range(len(repeats[group])):
                            if inside_repeat_num not in repeats_to_stop[group]:
                                start_id, go_to_id, go_to_direction = repeats[group][inside_repeat_num]
                                repeats[-1].append([start_id, go_to_id, go_to_direction])
                                new_connect = (
                                    go_to_id - (word_size - 1) * (go_to_direction == 1) + go_to_direction,
                                    go_to_direction)
                                if new_connect in points_to_repeats:
                                    points_to_repeats[new_connect].append((len(repeats) - 1, len(repeats[-1]) - 1))
                                else:
                                    points_to_repeats[new_connect] = [(len(repeats) - 1, len(repeats[-1]) - 1)]
            for one_connection in points_to_repeats:
                for group, repeat_num in points_to_repeats[one_connection]:
                    start_id, previous_id, this_direction = repeats[group][repeat_num]
                    repeats[group][repeat_num][1] += this_direction

            """part 2: dealing with new connection"""
            for one_connection in this_connection:
                here_id, direction_trans = one_connection
                candidate_last_connect = (here_id - direction_trans, direction_trans)
                if candidate_last_connect not in last_connection:
                    repeats.append([])
                    for inside_connection in this_connection:
                        inside_id, inside_direction = inside_connection
                        repeats[-1].append([inside_id + (word_size - 1) * (inside_direction == -1),
                                            inside_id + (word_size - 1) * (inside_direction == 1),
                                            inside_direction])
                        if (inside_id, inside_direction) in points_to_repeats:
                            points_to_repeats[inside_connection].append((len(repeats) - 1, len(repeats[-1]) - 1))
                        else:
                            points_to_repeats[inside_connection] = [(len(repeats) - 1, len(repeats[-1]) - 1)]
                    break

            if i + 1 < len_indices:
                next_index = repeat_indices[i + 1]
            else:
                next_index = None
            if next_index != this_index + 1:
                points_to_repeats = {}
                last_connection = set()
            else:
                last_connection = this_connection

        """aftertreatment"""
        # delete repeated repeats
        final_repeat = []
        repeat_dicts = set()
        for repeat_group in repeats:
            if tuple(repeat_group[0]) not in repeat_dicts:
                for single_repeat in repeat_group:
                    here_start, here_end, here_direction = single_repeat
                    repeat_dicts.add((here_start, here_end, here_direction))
                    repeat_dicts.add((here_end, here_start, -here_direction))
                final_repeat.append(repeat_group)
            else:
                continue

        if circular:
            # merge repeats of two ends
            def from_start_or_not(repeat_group_inside):
                for c_from, c_to, c_direction in repeat_group_inside:
                    if c_from == 0:
                        return True
                return False

            def to_end_or_not(repeat_group_inside):
                for c_from, c_to, c_direction in repeat_group_inside:
                    if (c_from, c_to)[c_direction == 1] == here_seq_length - 1:
                        return True
                return False

            def be_merging_or_not(candidate_repeat_1, candidate_repeat_2):
                inside_1_from, inside_1_to, inside_1_direction = candidate_repeat_1
                inside_2_from, inside_2_to, inside_2_direction = candidate_repeat_2
                if inside_1_direction != inside_2_direction:
                    return False
                else:
                    if inside_1_to == here_seq_length - 1 and (inside_2_from, inside_2_to)[
                                inside_2_direction == -1] == 0:
                        return True
                    else:
                        if (inside_1_to - inside_2_from + inside_1_direction) * inside_1_direction == word_size - 1:
                            return True
                        else:
                            return False

            count_from_start = 0
            while count_from_start < len(final_repeat) and from_start_or_not(final_repeat[count_from_start]):
                count_from_start += 1
            possible_ends = set()
            for candidate_group_id in range(len(final_repeat)):
                if to_end_or_not(final_repeat[candidate_group_id]):
                    possible_ends.add(candidate_group_id)
            # transform the direction of end to positive
            for i in possible_ends:
                e_direction = None
                for e_from, e_to, e_direction in final_repeat[i]:
                    if (e_from, e_to)[e_direction == 1] == here_seq_length - 1:
                        break
                if final_repeat[i] and e_direction == -1:
                    for j in range(len(final_repeat[i])):
                        e_from, e_to, e_direction = final_repeat[i][j]
                        final_repeat[i][j] = [e_to, e_from, -e_direction]

            repeat_group_to_del = set()
            for i in range(count_from_start):
                for j in possible_ends:
                    matches = {}
                    for candidate_id_1 in range(len(final_repeat[j])):
                        for candidate_id_2 in range(len(final_repeat[i])):
                            if be_merging_or_not(final_repeat[j][candidate_id_1], final_repeat[i][candidate_id_2]):
                                matches[candidate_id_1] = candidate_id_2
                                break
                    if len(matches) < 2:
                        # no effective matches
                        continue
                    elif len(matches) == len(final_repeat[j]) == len(final_repeat[i]):
                        # all matched, merge
                        repeat_group_to_del.add(i)
                        for merge_1_id, merge_2_id in matches.items():
                            final_repeat[j][merge_1_id][1] = final_repeat[i][merge_2_id][1]
                    elif len(matches) == len(final_repeat[j]):
                        # all of one side matches, extend this side
                        for merge_1_id, merge_2_id in matches.items():
                            final_repeat[j][merge_1_id][1] = final_repeat[i][merge_2_id][1]
                    elif len(matches) == len(final_repeat[i]):
                        # all of next side matches, extend next side
                        for merge_1_id, merge_2_id in matches.items():
                            final_repeat[i][merge_2_id][0] = final_repeat[j][merge_1_id][0]
                    else:
                        # some of the repeat matches, create a new repeat group
                        final_repeat.append([])
                        for merge_1_id, merge_2_id in matches.items():
                            final_repeat[-1].append([final_repeat[j][merge_1_id][0],
                                                     final_repeat[i][merge_2_id][1],
                                                     final_repeat[i][merge_2_id][2]])
            repeat_group_to_del = sorted(list(repeat_group_to_del), reverse=True)
            for id_to_del in repeat_group_to_del:
                del final_repeat[id_to_del]

            # circularize number
            for count_group_ in range(len(final_repeat)):
                for j in range(len(final_repeat[count_group_])):
                    here_start, here_end, here_direction = final_repeat[count_group_][j]
                    final_repeat[count_group_][j] = [here_start % check_len,
                                                     here_end % check_len,
                                                     here_direction]

        # delete small repeats
        count_group__ = 0
        while count_group__ < len(final_repeat):
            here_start, here_end, here_direction = final_repeat[count_group__][0]
            if ((here_end - here_start + here_direction) * here_direction) % check_len < min_repeat_length:
                del final_repeat[count_group__]
            else:
                count_group__ += 1
        # sort according to occurrence
        for group_to_sort in range(len(final_repeat)):
            # sort by direction
            # calculate repeat length
            start, end, direction = final_repeat[group_to_sort][0]
            if start == end:
                this_len = 1
            else:
                if (end - start) * direction > 0:
                    this_len = (end - start) * direction + 1
                else:
                    this_len = (start, end)[direction == 1] + check_len - (start, end)[direction != 1] + 1
            # transform into dict
            final_repeat[group_to_sort] = [{"start": start, "end": end, "direction": direction, "length": this_len}
                                           for start, end, direction in final_repeat[group_to_sort]]
        # sort according to length: from longest to shortest
        final_repeat.sort(key=lambda x: -x[0]["length"])
        return final_repeat
