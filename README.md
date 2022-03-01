# Platanus-3D README.md

## Description
Platanus-3D is a de novo chromosome-level scaffolding and phasing tool using Hi-C.
Platanus-3D generates chromosome-level haplotypes by scaffolding and phasing
the input contigs using a combination of information from Hi-C and other reads (PE, MP, LongRead).

## Version
v1.0.0

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
git clone https://github.com/ShunOuchi/Platanus-3D.git
cd src
make
cp platanus_3D <installation_path>
```


## Synopsis
### Inputs
* Input assembly (required)
    * Haplotype-aware style input (such as [Platanus-allee](http://platanus.bio.titech.ac.jp/platanus2)):<br>
      primaryBubble.fa secondaryBubble.fa nonBubble.fa
    * Pseudo-haplotype style input (such as [FALCON-Unzip](https://github.com/PacificBiosciences/FALCON_unzip)) or Mixed-haplotype style input (such as [Canu](https://github.com/marbl/canu)):<br>
      contigs.fa
* Input reads
    * Hi-C reads: HIC_1.fq HIC_2.fq (required)
    * Illumina paired-end: PE_1.fq PE_2.fq (optional)
    * Illumina mate-pair : MP_1.fq MP_2.fq (optional)
    * PacBio/Oxford-Nanopore long reads  : longread.fq (optional)

### Commands

* for Haplotype-aware style input
```
platanus_3D \
-c nonBubble.fa \
-b primaryBubble.fa secondaryBubble.fa \
-IP1 PE_1.fq PE_2.fq \
-OP2 MP_1.fq MP_2.fq \
-p longread.fq \
-HIC HIC_1.fq HIC_2.fq \
2>3D.log
```

* for Pseudo-haplotype or Mixed-haplotype style input
```
platanus_3D \
-cph contigs.fa \
-IP1 PE_1.fq PE_2.fq \
-OP2 MP_1.fq MP_2.fq \
-p longread.fq \
-HIC HIC_1.fq HIC_2.fq \
2>3D.log
```

### Final output
    out_afterPhase.fa (phased diploid scaffolds)


---
## Example
Below is showing examples how to run Platanus-3D using test dataset.
The test dataset is the simulated diploid dataset of _Caenorhabditis elegans_ chr1.
### Example 1. I have Platanus-allee assembly (Haplotype-aware style input)
```
platanus_3D \
-c Platanus-allee_result/out_nonBubbleOther.fa \
-b Platanus-allee_result/out_primaryBubble.fa Platanus-allee_result/out_secondaryBubble.fa \
-IP1 reads/PE_*.fq.gz \
-OP2 reads/MP5k_*.fq.gz \
-OP3 reads/MP9k_*.fq.gz \
-p reads/longread.fq.gz \
-HIC reads/HIC_1.fq.gz reads/HIC_2.fq.gz
```
### Example 2. I have FALCON-Unzip assembly (Psuedo-haplotype style input)
```
platanus_3D \
-cph FALCON-Unzip_result/cns_p_ctg.fa FALCON-Unzip_result/cns_h_ctg.fa \
-p reads/longread.fq,gz \
-HIC reads/HIC_1.fq.gz reads/HIC_2.fq.gz
```
### Example 3. I have Canu assembly (Mixed-haplotype style input)
```
platanus_3D \
-cph Canu_result/asm.contigs.fa \
-p reads/longread.fq.gz \
-HIC reads/HIC_1.fq.gz reads/HIC_2.fq.gz
```

---
## Usage
### Command
```sh
platanus_3D [OPTIONS] 2>log
```
### Options
    -o STR                             : prefix of output file and directory (do not use "/", default out, length <= 200)
    -c FILE1 [FILE2 ...]               : contig (or scaffold) file (fasta format; for Haplotype-aware style input)
    -b FILE1 [FILE2 ...]               : bubble seq file (fasta format; for Haplotype-aware style input)
    -cph FILE1 [FILE2 ...]             : contig (or scaffold) file (fasta format; for Pseudo-haplotype or Mixed-haplotype style input; only effective without -c, -b option)
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
   PREFIX_ConsensusOutput.fa

   PREFIX_afterPhase.fa

PREFIX is specified by -o


---
## Notes
* Compressed input files

Both uncompressed and compressed (gzip or bzip2) FASTA/FASTQ files are accepted.
Formats are auto-detected. Internally, "file -bL", "gzip -cd" and "bzip2 -cd" commands, which can be
used in most of the UNIX OSs, are utilized.

* Minimap2

This tool is used to align PacBio/Oxford-Nanopore long reads.
When long reads are input through the -p option, please check Minimap2 is installed as "minimap2" command
or specify the path of Minimap2 using the -mapper option.

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

