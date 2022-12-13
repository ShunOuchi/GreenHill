#!/usr/bin/env python

import sys
import re


def main(file_name):
    contig_re = re.compile(r"[acgtACGT]+")
    record_lists = []
    with open(file_name) as fin:
        record_i = 1
        record_list = []
        for (name, seq) in read_fasta(fin):
            frag_i = 1
            pre_end = 0
            record_list = []
            for m in contig_re.finditer(seq):
                if pre_end != 0:
                    print(f">{name}:::gap_{frag_i} {record_i} {m.start() - pre_end}")
                    record_list.append(record_i)
                    record_i += 1
                    frag_i += 1
                print(f">{name}:::ctg_{frag_i} {record_i} {m.end() - m.start()}")
                record_list.append(record_i)
                record_i += 1
                frag_i += 1
                pre_end = m.end()
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
