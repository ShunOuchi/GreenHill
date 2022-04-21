# Platanus-allee ChangeLog.md

## Version 2.2.2 (17 May 2019)
* To reduce false-detection of heterozygous regions, the fork-anchor function is deprecated.
* To reduce mis-assemblies, the partial-gap-close function in the phasing step is deprecated.
* Add -no_partial option of the consensus command to turn off the partial-gap-close function.

## Version 2.2.1 (8 May 2019)
* Improve performance to utilize long-reads.
* Add aggressive mode for phase and consensus commands (-aggressive). It is experimental and not recommended. 
* Simplify output-files of consensus command.

## Version 2.2.0 (25 Apr. 2019)
* Rename and simplify output file-names.
* Add gap-closing function in the "consensus" command.
* Improve ability to detect homologous haplotype-pairs.
* Reduce mis-assemblies for many cases.
* Reduce gap-rates using partial-gap-close function.
* Improve the speed of read-mapping.

## Version 2.1.0 (20 Nov. 2018)
* Add "divide" command and functions to divide erroneous sequences in a phasing step.

## Version 2.0.2 (16 Sep. 2018)
* Fixed a bug in the untangling function for cross-structures in graphs (phase, solve_DBG commands).
* Compressed (gzip or bzip2) FASTA/FASTQ files can be accepted.

## Version 2.0.1 (17 Aug. 2018)
* Fixed a bug causing the segmentation fault in gap_close.
* The result is the same as that of Version 2.0.0 if it finishes a process.

## Version 2.0.0 (9 Jul. 2018)
* The version used to write the manuscript (under review).
* Initial release as an exutable file for Linux x86_64.
