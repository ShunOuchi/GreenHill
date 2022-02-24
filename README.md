# Platanus-allee README.md

## Description
Platanus-allee is a de novo haplotype assembler (phasing tool), which assembles each haplotype
sequence in a diploid genome. Compared to the read mapping-based haplotype phasing tools,
Platanus-allee is especially useful to analyze highly divergent (heterozygous) regions in
which haplotypes extremely differ. This tool requires at least one Illumina library and can
accepts Illumina mate-pairs, PacBio/Oxford-Nanopore long reads and 10X linked-reads. In addition
to the haplotype phasing function, Platanus-allee can construct consensus sequences (pseudo
haploid genome), which have mosaic structures of haplotypes (i.e., maternal and paternal haplotypes
are mixed) and can be used as the conventional "draft genome". Note that Platanus-allee was
previously called "Platanus2" and renamed according to emphasize the difference from Platanus,
which is the tool to assemble consensus sequences of haplotypes.

## Version
v1.0.0

## Web site
<http://platanus.bio.titech.ac.jp/>

## Author
Shun Ouchi and Rei Kajitani at Tokyo Institute of Technology wrote key source codes.
Address for this tool: <platanus@bio.titech.ac.jp>


## Requirements
* GCC
    - <https://gcc.gnu.org/>
    - version >= 4.4, with OpenMP
    - To compile the source code.

* Minimap2
    - <https://github.com/lh3/minimap2>
    - Only required to use PacBio/Oxford-Nanopore long reads.

## Installation
```sh
cd src
make
cp platanus_3D <installation_path>
```


## Synopsis
### Inputs
* Input assembly (required)
    * Haplotype-aware style input: primaryBubble.fa secondaryBubble.fa nonBubble.fa
    * Pseudo-haplotype style input: contigs.fa
* Input reads
    * Hi-C reads: HIC_1.fq HIC_2.fq (required)
    * Illumina paired-end: PE_1.fq PE_2.fq (optional)
    * Illumina mate-pair : MP_1.fq MP_2.fq (optional)
    * PacBio/Oxford-Nanopore long reads  : longread.fq (optional)

### Commands
```
#for Haplotype-aware style input
platanus_3D solveDBG -3D \
-c nonBubble.fa
-b primaryBubble.fa secondaryBubble.fa\
-IP1 PE_1.fq PE_2.fq \
-OP2 MP_1.fq MP_2.fq \
-p longread.fq \
-HIC HIC_1.fq HIC_2.fq \
2>3D.log

#for Pseudo-haplotype style input
platanus_3D solveDBG -3D \
-cph contigs.fa
-IP1 PE_1.fq PE_2.fq \
-OP2 MP_1.fq MP_2.fq \
-p longread.fq \
-HIC HIC_1.fq HIC_2.fq \
2>3D.log
```

### Final output
    afterPhase.fa (phased diploid scaffolds)

---
## Usage
### Command
```sh
platanus_3D solveDBG -3D [OPTIONS] 2>log
```
### Options
    -o STR                             : prefix of output file and directory (do not use "/", default out, length <= 200)
    -c FILE1 [FILE2 ...]               : contig (or scaffold) file (fasta format; for Haplotype-aware style input)
    -b FILE1 [FILE2 ...]               : bubble seq file (fasta format; for Haplotype-aware style input)
    -cph FILE1 [FILE2 ...]             : contig (or scaffold) file (fasta format; for Pseudo-haplotype style input; only effective without -c, -b option)
    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (interleaved file, fasta or fastq)
    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (separate forward and reverse files, fasta or fastq)
    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (interleaved, fasta or fastq)
    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (separate forward and reverse files, fasta or fastq)
    -p PAIR1 [PAIR2 ...]               : long-read file (PacBio, Nanopore) (fasta or fastq)
    -hic PAIR1 [PAIR2 ...]             : HiC_pair_files (reads in 1 file, fasta or fastq)
    -HIC FWD1 REV1 [FWD2 REV2 ...]     : HiC_pair_files (reads in 2 files, fasta or fastq)
    -t INT                             : number of threads (default 1)
    -tmp DIR                           : directory for temporary files (default .)
    -l INT                             : minimum number of links to scaffold (default 3)
    -k INT                             : minimum number of links to phase variants (default 1)
    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default 32 64 96)
    -mapper FILE                       : path of mapper executable file (default minimap2, only effective with -p option)
    -minimap2_sensitive                : sensitive mode for minimap2 (default, off; only effective with -p option)

### Input format:
   Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -ip, -IP, -op, -OP, -p, -hic and -HIC option.

### Outputs:
   afterPhase.fa

---
## Notes
* Compressed input files
Both uncompressed and compressed (gzip or bzip2) FASTA/FASTQ files are accepted.
Formats are auto-detected. Internally, "file -bL", "gzip -cd" and "bzip2 -cd" commands, which can be
used in most of the UNIX OSs, are utilized.

* Minimap2
This tool is used to align PacBio/Oxford-Nanopore long reads and included in Platanus-allee package
as the directory of "minimap2-2.0-r191". When long reads are input through the -p option of "phase"
and "consensus" commands, please check Minimap2 is installed as "minimap2" command or specify the
path of Minimap2 using the -mapper option.

* Paired-end (mate-pair) input
Paired libraries are classified into "inward-pair" and "outward-pair" according to the sequence direction.
For file formats, separate and interleaved files can be input through -IP (-OP) and -ip (-op)
options, respectively.

Inward-pair (usually called "paired-end", accepted in options "-IP" or "-ip"):

    FWD --->
        5' -------------------- 3'
        3' -------------------- 5'
                        <--- REV

Outward-pair (usually called "mate-pair", accepted in options "-OP" or "-op"):

                        ---> REV
        5' -------------------- 3'
        3' -------------------- 5'
    FWD <---

Example inputs:

    Inward-pair (separate, insert=300)   : PE300_1.fq PE300_2.fq
    Inward-pair (interleaved, insert=500): PE500_pair.fq
    Outward-pair (separate, insert=2k)   : MP2k_1.fa MP2k_2.fq

Corresponding options:

    -IP1 PE300_1_pair.fq PE300_2.fq \
    -ip2 PE500_pair.fq \
    -OP3 MP2k_1.fq MP2k_2.fq

* Compressed input files
Compressed files can be input using the process substitution function ("<(command)") of bash/zsh
without temporary uncompressed files.
e.g.
```sh
# Files compressed by gzip
platanus_allee assemble -f <(zcat PE_1.fq.gz) <(zcat PE_2.fq.gz)
# Files compressed by bzip2
platanus_allee assemble -f <(bzcat PE_1.fq.bz2) <(bzcat PE_2.fq.bz2)
```
