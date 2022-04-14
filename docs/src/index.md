# Nucleotide_Essentials.jl

```@contents
```

## Data Types
```@docs
FastqRecord
FastaRecord
```

## Functions

```@docs
readFastq
readFasta
writeFasta
FastqtoFasta
FilterQuality_se
FilterQuality_pe
PlotQuality
potential_mismatches
reverse_complement
demultiplex_se
demultiplex_pe
```

## Index

```@index
```
## Change Log 

#### Nucleotide_Essentials v0.2.0

* Added support for quality filtering of .fastq reads 
* Added support for Gzip compressed files 
* Performance improvements in ```PlotQuality()``` and added support for exporting quality plots
* Added support for automatic quality profile encoding detection (Phred+64 and Phred+33 encoding)
* Minor documentation updates 
