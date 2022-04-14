# <img src="assets/nucleotide_essentials_logo.svg" width="30%" align="right" /> Nucleotide_Essentials.jl

[![CI](https://github.com/phorve/Nucleotide_Essentials.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/phorve/Nucleotide_Essentials.jl/actions/workflows/CI.yml)
[![Runtests](https://github.com/phorve/Nucleotide_Essentials.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/phorve/Nucleotide_Essentials.jl/actions/workflows/Runtests.yml)
[![codecov](https://codecov.io/gh/phorve/Nucleotide_Essentials.jl/branch/main/graph/badge.svg?token=ET2C1TKZGP)](https://codecov.io/gh/phorve/Nucleotide_Essentials.jl)

**Nucleotide_Essentials.jl** is a collection of tools for working with Illumina next-generation sequencing reads currently under development and testing.

## Installation:

Active development is still underway but the most current version can be downloaded from the REPL using: 

The package is registered in the General registry and so can be installed at the REPL with: 
```julia
] add Nucleotide_Essentials
```

## Current Functions: 
* `readFastq` - Import a .fastq file into Julia
* `readFasta` - Import a .fasta file into Julia
* `writeFasta` - Write a .fasta or .fasta.gz to file
* `FilterQuality_se` - Filter single-ended Illumina reads
* `FilterQuality_pe` - Filter paired-end Illumina reads
* `FastqtoFasta` - Convert a .fastq file or FastqRecord to a FastaRecord
* `potential_mismatches` - Create a list of potential mismatch barcodes (mutations + deletions)
* `demultiplex_se` - Demultiplex single-end Illumina reads
* `demultiplex_pe` - Demultiplex paired-end Illumina reads
* `PlotQuality` - Visualize quality of .fastq file

For bug fixes, new features, and future suggestions, open an issue! 

Full documentation can be found at [here](https://www.patrickfhorve.com/Nucleotide_Essentials.jl/v0.2/)