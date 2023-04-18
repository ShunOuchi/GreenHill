#!/usr/bin/env python

import sys
import re
from math import ceil, floor


def main(file_name):
    contig_re = re.compile(r"[nN]+")
    record_lists = []
    with open(file_name) as fin:
        record_i = 1
        record_list = []
        for (name, seq) in read_fasta(fin):
            frag_i = 1
            pre_gap_start = 0
            pre_gap_end = 0
            record_list = []
            for m in contig_re.finditer(seq):
                ctg_start = pre_gap_end - ((pre_gap_end - pre_gap_start) // 2 + (pre_gap_end - pre_gap_start) % 2)
                ctg_end = m.start() + ((m.end() - m.start()) // 2)
                print(f">{name}:::fragment_{frag_i} {record_i} {ctg_end - ctg_start}")
                record_list.append(record_i)
                record_i += 1
                frag_i += 1
                pre_gap_start = m.start()
                pre_gap_end = m.end()
            ctg_start = pre_gap_end - ((pre_gap_end - pre_gap_start) // 2 + (pre_gap_end - pre_gap_start) % 2)
            ctg_end = len(seq)
            print(f">{name}:::fragment_{frag_i} {record_i} {ctg_end - ctg_start}")
            record_list.append(record_i)
            record_i += 1
            frag_i += 1
            record_lists.append(record_list)

        for record_list in record_lists:
            print(" ".join([str(_) for _ in record_list]))


def read_fasta(file):
    name_re = re.compile(r"^>(\S+)")
    name = ""
    seq = ""
    for line in file:
        match = name_re.match(line)
        if match:
            name = match.group(1)
            break

    for line in file:
        match = name_re.match(line)
        if match:
            yield (name, seq)
            name = match.group(1)
            seq = ""
        else:
            seq += line.rstrip("\n")

    yield (name, seq)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage:", sys.argv[0], "in.fasta >for_juicebox.assembly", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1])
