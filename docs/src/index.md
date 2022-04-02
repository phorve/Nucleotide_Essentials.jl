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

* Add support for quality filtering .fastq reads 
* Add support for automatic quality profile encoding detection 
* Minor documentation updates 
