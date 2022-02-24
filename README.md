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
v2.2.2

## Web site
<http://platanus.bio.titech.ac.jp/>

## Author
Rei Kajitani at Tokyo Institute of Technology wrote key source codes.  
Address for this tool: <platanus@bio.titech.ac.jp>


## Requirements
* GCC 
    - <https://gcc.gnu.org/>
    - version >= 4.4, with OpenMP
    - To compile the source code.

* Minimap2
    - <https://github.com/lh3/minimap2>
    - Only required to use PacBio/Oxford-Nanopore long reads.

* Long Ranger
    - <https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger>
    - Only required to prepare 10X linked-reads input (barcoded.fastq).
   
   
## Installation
```sh
make
cp platanus_allee <installation_path>
```


## Synopsis
### Inputs
* Illumina paired-end: PE_1.fq PE_2.fq (mandatory)
* Illumina mate-pair : MP_1.fq MP_2.fq (optional)
* PacBio long reads  : PacBio.fq (optional)
* 10X linked-reads   : 10X_barcoded.fq (optional)

### Commands (normal mode)
```
platanus_allee assemble -f PE_1.fq PE_2.fq 2>assemble.log

platanus_allee phase \
-c out_contig.fa \
-IP1 PE_1.fq PE_2.fq \
-OP2 MP_1.fq MP_2.fq \
-p PacBio.fq \
-x 10X_barcoded.fq \
2>phase.log

platanus_allee consensus \
-c out_consensusInput.fa \
-IP1 PE_1.fq PE_2.fq \
-OP2 MP_1.fq MP_2.fq \
-p PacBio.fq -x 10X_barcoded.fq \
-x 10X_barcoded.fq \
2>consensus.log
```

### Final output
    out_allPhasedScaffold.fa (phased diploid scaffolds)
    out_consensusScaffolds.fa (consensus haploid scaffolds)


---
## Contig assembly usage
### Command
```sh
platanus_allee assemble [OPTIONS] 2>log
```
### Options
    -o STR               : prefix of output files (default out, length <= 200)
    -f FILE1 [FILE2 ...] : reads file (fasta or fastq, number <= 100)
    -t INT               : number of threads (<= 100, default 1)
    -m INT               : memory limit for making kmer distribution (GB, >=1, default 16)
    -tmp DIR             : directory for temporary files (default .)
    -k INT               : initial k-mer size (default 32)
    -K FLOAT             : maximum-k-mer factor (maximum-k = FLOAT*read-length, default  0.5)
    -s INT               : step size of k-mer extension (>= 1, default 20)
    -n INT               : initial k-mer coverage cutoff (default 0, 0 means auto)
    -c INT               : minimum k-mer coverage (default 2)
    -a FLOAT             : k-mer extension safety level (default 10.0)
    -u FLOAT             : maximum difference for bubble crush (identity, default 0)
    -d FLOAT             : maximum difference for branch cutting (coverage ratio, default 0.5)
    -e FLOAT             : k-mer coverage depth (k = initial k-mer size specified by -k) of homozygous region (default auto)

### Input format:
   Uncompressed and compressed (gzip or bzip2) files are accepted for -f option.

### Outputs:
    PREFIX_contig.fa
    PREFIX_kmerFrq.tsv

PREFIX is specified by -o
  
  
## Haplotype phasing usage
### Command
```sh
platanus_allee phase [OPTIONS] 2>log
```
### Options
    -o STR                             : prefix of output file and directory (do not use "/", default out, length <= 200)
    -c FILE1 [FILE2 ...]               : contig (or scaffold) file (fasta format)
    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (interleaved file, fasta or fastq)
    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (separate forward and reverse files, fasta or fastq)
    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (interleaved, fasta or fastq)
    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (separate forward and reverse files, fasta or fastq)
    -p PAIR1 [PAIR2 ...]               : long-read file (PacBio, Nanopore) (fasta or fastq)
    -x PAIR1 [PAIR2 ...]               : linked-reads files (paired-ends, 10x Genomics) (interleaved, fasta or fastq)
    -X FWD1 REV1 [FWD2 REV2 ...]       : linked-reads files (paired-ends, 10x Genomics) (separate forward and reverse files, fasta or fastq)
    -t INT                             : number of threads (default 1)
    -tmp DIR                           : directory for temporary files (default .)
    -i INT                             : number of iterations (default 2)
    -l INT                             : minimum number of links to scaffold (default 3)
    -k INT                             : minimum number of links to phase variants (default 1)
    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default 32 64 96)
    -mapper FILE                       : path of mapper executable file (default minimap2, only effective with -p option)
    -minimap2_sensitive                : sensitive mode for minimap2 (default, off; only effective with -p option)

### Input format:
   Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -ip, -IP, -op, -OP, -p, -x and -X option.

### Outputs:
    PREFIX_allPhasedScaffold.fa (including sequences below)
    PREFIX_primaryBubble.fa
    PREFIX_secondaryBubble.fa
    PREFIX_nonBubbleHetero.fa
    PREFIX_nonBubbleOther.fa
    PREFIX_consensusInput.fa (for "consensus" command (-c))


PREFIX is specified by -o
  
  
## Consensus scaffold construction usage
### Command
```sh
platanus_allee consensus [OPTIONS] 2>log
```

### Options
    -o STR                             : prefix of output file (default out, length <= 200)
    -c FILE1 [FILE2 ...]               : input_scaffolds (fasta format)
    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (interleaved file, fasta or fastq)
    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (separate forward and reverse files, fasta or fastq)
    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (interleaved, fasta or fastq)
    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (separate forward and reverse files, fasta or fastq)
    -p PAIR1 [PAIR2 ...]               : long-read file (PacBio, Nanopore) (fasta or fastq)
    -x PAIR1 [PAIR2 ...]               : linked-reads files (paired-ends, 10x Genomics) (interleaved, fasta or fastq)
    -X FWD1 REV1 [FWD2 REV2 ...]       : linked-reads files (paired-ends, 10x Genomics) (separate forward and reverse files, fasta or fastq)
    -t INT                             : number of threads (default 1)
    -tmp DIR                           : directory for temporary files (default .)
    -e FLOAT                           : coverage depth of homozygous region (default auto)
    -L INT                             : maximum fragment length of tag (10x Genomics) (default 200000)
    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default 32 64 96)
    -l INT                             : minimum number of links to scaffold (default 3)
    -mapper FILE                       : path of mapper executable file (default, minimap, only effective with -p option)
    -minimap2_sensitive                : sensitive mode for minimap2 (default, off; only effective with -p option)
    -no_partial                        : not close gaps partially, i.e. only close ones completely (default, off)

### Input format:
   Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -ip, -IP, -op, -OP, -p, -x and -X option.

### Outputs
    PREFIX_consensusScaffold.fa (final)


PREFIX is specified by -o
  
  
## Dividing erroneous sequences usage
### Command
```sh
platanus_allee divide [OPTIONS] 2>log
```

### Options
    -o STR                             : prefix of output file (default out, length <= 200)
    -c FILE1 [FILE2 ...]               : contig_file (fasta format)
    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (interleaved file, fasta or fastq)
    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (separate forward and reverse files, fasta or fastq)
    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (interleaved, fasta or fastq)
    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (separate forward and reverse files, fasta or fastq)
    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default 32 64 96)
    -t INT                             : number of threads (<= 1, default 1)
    -tmp DIR                           : directory for temporary files (default .)

### Outputs
    PREFIX_divided.fa
    PREFIX_dividedComponent.bed

PREFIX is specified by -o


---
## Notes
* Options related to run time
Although -t (number of threads) of all commands and -m (memory amount) of the "assemble" command are 
not mandatory to run, it is recommended to set the values adjusting your machine-environment.
These options may severely effect the run time.  
e.g.,  
Available number of threads and memory amount are 12 and 256GB, respectively.  
->  -t 12 -m 256

* Compressed input files
Both uncompressed and compressed (gzip or bzip2) FASTA/FASTQ files are accepted.
Formats are auto-detected. Internally, "file -bL", "gzip -cd" and "bzip2 -cd" commands, which can be
used in most of the UNIX OSs, are utilized.

* Minimap2
This tool is used to align PacBio/Oxford-Nanopore long reads and included in Platanus-allee package 
as the directory of "minimap2-2.0-r191". When long reads are input through the -p option of "phase" 
and "consensus" commands, please check Minimap2 is installed as "minimap2" command or specify the 
path of Minimap2 using the -mapper option.

* 10X linked-reads input
Platanus-allee accepts barcoded fastq files, which have the barcode information in name lines 
(@...) as the "BX:Z:" tags. These files, "barcoded.fastq", can be generated using the 
"longranger basic" command of Long Ranger (https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger).  
e.g.
```sh
longranger basic --fastqs=/raw/fastq/directory --id=your_id --sample=your_sample_name
# Result: your_id/outs/barcoded.fastq.gz
```

* Paired-end (mate-pair) input  
The "phase" and "consensus" accept paired-end and/or mate-pair libraries. Paired libraries are 
classified into "inward-pair" and "outward-pair" according to the sequence direction. 
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
