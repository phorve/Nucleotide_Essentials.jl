# Nucleotide_Essentials.jl

[![Build status](https://github.com/phorve/Nucleotide_Essentials.jl/workflows/CI/badge.svg)](https://github.com/phorve/Nucleotide_Essentials.jl/actions)

[![codecov.io](http://codecov.io/github/phorve/Nucleotide_Essentials.jl/coverage.svg?branch=main)](http://codecov.io/github/phorve/Nucleotide_Essentials.jl?branch=main)

**Nucleotide_Essentials.jl** is a collection of tools for working with next-generation sequencing reads currently under development and testing.

## Installation:

Active development is still underway but the most current version can be downloaded from the REPL using: 
```julia 
using Pkg; Pkg.add("https://github.com/phorve/Nucleotide_Essentials.jl")
```
## Current Functions: 
* `readFastq` - Import a .fastq file into Julia
* `potential_mismatches` - Create a list of potential mismatch barcodes (mutations + deletions)
* `demultiplex_se` - Demultiplex single-end Illumina reads
* `demultiplex_pe` - Demultiplex paired-end Illumina reads
* `PlotQuality` - Visualize quality of .fastq file

Full documentation can be found at [here](https://www.patrickfhorve.com/Nucleotide_Essentials.jl/dev/)