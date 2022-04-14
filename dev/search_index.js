var documenterSearchIndex = {"docs":
[{"location":"#Nucleotide_Essentials.jl","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.jl","text":"","category":"section"},{"location":"","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.jl","text":"","category":"page"},{"location":"#Data-Types","page":"Nucleotide_Essentials.jl","title":"Data Types","text":"","category":"section"},{"location":"","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.jl","text":"FastqRecord\nFastaRecord","category":"page"},{"location":"#Nucleotide_Essentials.FastqRecord","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.FastqRecord","text":"Nucleotide_Essentials.FastqRecord\n\nComponents\n\nID: The unique sequence identifier associated with that entry\nsequence: The nucleotide sequence of that entry\nquality: The quality scores of that entry \nfilename: The original file name \n\n\n\n\n\n","category":"type"},{"location":"#Nucleotide_Essentials.FastaRecord","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.FastaRecord","text":"Nucleotide_Essentials.FastaRecord\n\nComponents\n\nID: The unique sequence identifier associated with that entry\nsequence: The nucleotide sequence of that entry\nfilename: The original file name \n\n\n\n\n\n","category":"type"},{"location":"#Functions","page":"Nucleotide_Essentials.jl","title":"Functions","text":"","category":"section"},{"location":"","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.jl","text":"readFastq\nreadFasta\nwriteFasta\nFastqtoFasta\nFilterQuality_se\nFilterQuality_pe\nPlotQuality\npotential_mismatches\nreverse_complement\ndemultiplex_se\ndemultiplex_pe","category":"page"},{"location":"#Nucleotide_Essentials.readFastq","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.readFastq","text":"Nucleotide_Essentials.readFastq\nreadFastq(Path::String)\n\n.fastq file => readFastq(Path) => FastqRecord(ID, sequence, quality, filename)\n\nsupported keyword arguments include: \n\nPath::String: The full or relative path to a .fastq file\n\nExample:\n\n# Supply the path to a .Fastq file that you would like to import\nmyfastq = readFastq(\"myfastq.fastq\")\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.readFasta","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.readFasta","text":"Nucleotide_Essentials.readFasta\n\nImports a .fasta file into julia\n\nreadFasta(Path::String)\n\n.fasta file => readFasta(Path) => FastaRecord(ID, sequence, filename)\n\nSupported keyword arguments include: \n\nPath::String: The full or relative path to a .fasta file\n\nExample:\n\n# Supply the path to a .fasta file that you would like to import - it is recommended to include `;` in your command to prevent printing potentially large .fasta files in the REPL\nmyfasta = readFasta(\"myfasta.fasta\");\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.writeFasta","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.writeFasta","text":"Nucleotide_Essentials.writeFasta\n\nreadFasta(Path::String) FastaRecord => writefasta(inputfasta, out, compressed) => .fasta file/.fasta.gz file \n\nCreates a single or multiple entry FastaRecord and outputs either a .fasta or compressed .fasta.gz file to the desired directory\n\nSupported keyword arguments include: \n\ninput_fasta::FastaRecord: A FastaRecord with either a single entry or multiple entries \nout::String: The full or relative path to the directory where files should be written to\ncompressed::Bool: Whether or not to write the .fasta files as compressed files or not. \nIf true, files will written as .fasta.gz files\nIf false, files will written as .fasta files\n\nExample:\n\n# .fasta files can be written as from an already imported FastaRecord in Julia \nmyfasta = readFasta(\"myfasta.fasta\");\nwriteFasta(input_fasta, \"example/output/directory, false)\n\n# .fasta files can be written as a .fasta.gz from an already imported FastaRecord in Julia \nmyfasta = readFasta(\"myfasta.fasta\");\nwriteFasta(input_fasta, \"example/output/directory\", true)\n\n# .fasta files with multiple sequences can be read and written as individual .fasta or .fasta.gz in the same step\nmyfasta = readFasta(\"myfasta.fasta\");\nwriteFasta(readFasta(\"/myfasta.fasta\"), \"example/output/directory\", true);\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.FastqtoFasta","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.FastqtoFasta","text":"Nucleotide_Essentials.FastqtoFasta\n\nConverts a FastqRecord to a FastaRecord. Can also input and convert a .fastq file to a FastaRecord in the same function. \n\nFastqtoFasta(Fastq::Union{String, FastqRecord}) .fastq file => FastqtoFasta(Fastq) => FastaRecord(ID, sequence, filename) FastqRecord(ID, sequence, quality, filename) => FastqtoFasta(Fastq) => FastaRecord(ID, sequence, filename)\n\nSupported keyword arguments include: \n\nFastq::Union{String, FastqRecord}: \nThe full or relative path to a .fastq file\nA FastqRecord\n\nExample:\n\n# Supply the path to a .fastq file that you would like to convert to a FastRecord\nmyfasta = FastqtoFasta(\"myfastq.fastq\");\n\n# Alternatively, a FastqRecord can be used as the input \nmyfasta = FastqtoFasta(myFastqRecord);\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.FilterQuality_se","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.FilterQuality_se","text":"Nucleotide_Essentials.FilterQuality_se\n\nFilters an input .fastq file based upon the encoded Phred+33 or Phred+64 quality scores. The encoding of the reads is automatically deteremined by looking for unique encoding in Phred+33 and Phred+64. Phred+64 encoding is identified by searching for ^, a, ], and f.\n\nReads are filtered based upon the number of expected errors (mathrmE) based on the error rate based on quality score and the sum of error probabilities, following the equation: \n\nmathrmE = sum_ip_i = sum_i10^frac-Q_i10\n\nStringent filtering (maxEE = 1) is used by default but can be adjusted by the user. \n\nReads that pass the filtering parameters are output to a file ending in _FilteredReads.fastq in the user-determined directory, as indicated by out. \n\nSupported keyword arguments include: \n\nread1::String: Path to the reads to undergo quality filtering \nout::String: Path to the directory where reads that pass the quality filtering should be written \nmaxEE::Int64 (optional): The max number of expected errors a read can include as the filtering parameter (default: maxEE = 1)\nverbose::Bool (optional): Whether or not to show some intermediary feedback on the progress of the function (default = false)\n\nExample:\n\nFilterQuality_se(\"forward_R1.fasta\", \"/outdirectory\")\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.FilterQuality_pe","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.FilterQuality_pe","text":"Nucleotide_Essentials.FilterQuality_pe\n\nFilters an input .fastq file based upon the encoded Phred+33 or Phred+64 quality scores. The encoding of the reads is automatically deteremined by looking for unique encoding in Phred+33 and Phred+64. Phred+64 encoding is identified by searching for ^, a, ], and f.\n\nReads are filtered based upon the number of expected errors (mathrmE) based on the error rate based on quality score and the sum of error probabilities, following the equation: \n\nmathrmE = sum_ip_i = sum_i10^frac-Q_i10\n\nStringent filtering (maxEE = 1) is used by default but can be adjusted by the user. \n\nOutput Files: \n\nIf both the forward and reverse reads pass the filtering parameters: \nForward reads are output to a file ending in R1_Paired_filtered.fastq in the user-determined directory, as indicated by out\nReverse reads are output to a file ending in R2_Paired_filtered.fastq in the user-determined directory, as indicated by out\nIf only the forward read passes the filtering parameters: \nForward reads are output to a file ending in R1_Unpaired_filtered.fastq in the user-determined directory, as indicated by out\nReverse reads are not written to a file \nIf only the reverse read passes the filtering parameters: \nReverse reads are output to a file ending in R2_Unpaired_filtered.fastq in the user-determined directory, as indicated by out\nForward reads are not written to a file \n\nSupported keyword arguments include: \n\nread1::String: Path to the forward reads to undergo quality filtering \nread2::String: Path to the reverse reads to undergo quality filtering \nout::String: Path to the directory where reads that pass the quality filtering should be written \nmaxEE::Int64 (optional): The max number of expected errors a read can include as the filtering parameter (default: maxEE = 1)\nverbose::Bool (optional): Whether or not to show some intermediary feedback on the progress of the function (default = false)\n\nExample:\n\nFilterQuality_pe(\"forward_R1.fasta\", \"reverse_R2.fasta\", \"/outdirectory\")\n\n# changing the filtering parameters \nFilterQuality_pe(\"forward_R1.fasta\", \"reverse_R2.fasta\", \"/outdirectory\", 2, true)\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.PlotQuality","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.PlotQuality","text":"Nucleotide_Essentials.PlotQuality\n\nReturns a plot of the quality profile of a .fastq or .fastq.gz file \n\nThis function plots a visual summary of the distribution of quality scores (automatically detects Phred+33 or Phred+64 encoding) as a function of sequence position for the input fastq file(s).\n\nThe plotted lines show summary statistics at each sequence position:\n\ngreen is the mean \ndashed red lines are the 25th and 75th quantiles\n\nSupported keyword arguments include:\n\nInput::FastqRecord: The name of a FastqRecord for plotting\nverbose::Bool (optional): Whether or not to show some intermediary feedback on the progress of the function (default = false)\noutputfigure::Bool (optional): Whether or not to output a .png file with the created QualityPlot (default = false)\nfigurepath::String (optional): If outputting a .png figure to file, specify the path to a directory where the file should be written to (default = pwd())\n\nExample:\n\n# A quality profile can be created by supply the path to a .fastq or .fastq.gz file \nPlotQuality(\"path/to/my/file.fastq\")\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.potential_mismatches","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.potential_mismatches","text":"Nucleotide_Essentials.potential_mismatches\n\nReturns an Vector{Any} of potential barcodes with a single nucleotide change, including both deletions and substitutions\n\nSupported keyword arguments include:\n\nPath::String: The full or relative path to a .fastq file\nmismatch::Int64: The number of altered nucleotides to include (1 is only supported at this time)\n\nExample:\n\npotential_mismatches(\"GCGT\", 1)\n17-element Vector{Any}:\n\"GCGT\"\n\"CCGT\" \n\"ACGT\" \n\"TCGT\" \n\"GGGT\" \n\"GAGT\" \n\"GTGT\" \n\"GCCT\" \n\"GCAT\" \n\"GCTT\" \n\"GCGG\" \n\"GCGC\" \n\"GCGA\" \n\"CGT\"\n\"GGT\"\n\"GCT\" \n\"GCG\"\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.reverse_complement","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.reverse_complement","text":"Nucleotide_Essentials.reverse_complement\n\nTakes a string of nucleotide bases and returns the reverse complement of that string. Accepts inputs of String and SubString{String} (input from a FastqRecord)\n\nSupported keyword arguments include:\n\nsequence::Union{String, SubString{String}}: A string sequence of nucleotide bases or sequence entry from a FastqRecord\n\nExample:\n\nreverse_complement(\"ATCGT\")\n\"ACGAT\"\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.demultiplex_se","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.demultiplex_se","text":"Nucleotide_Essentials.demultiplex_se\n\nCompares a list of provided barcodes with the provided multiplexed reads and separates the reads into individual .fastq files. If a barcode is found within the read, the barcode is removed from the sequence. The quality data of the reads is preserved and written to the outputted .fastq file. If a barcode is not found, the sequnce and quality is written to the unassigned .fastq file unchanged. \n\nThe mapping file must be either a .csv or .txt file with two columns. The first column heading must be SampleID and the second column heading must be BarcodeSequence. \n\nEXAMPLE MAPPING FILE: \n\nSampleID BarcodeSequence\nSample1 Barcode1\nSample2 Barcode2\nSample3 Barcode3\nSample4 Barcode4\nSample5 Barcode5\nSample6 Barcode6\nSample7 Barcode7\nSample8 Barcode8\n\nSupported keyword arguments include:\n\nR1::String: Path to multiplexed reads   \nMap::String: Path to the mapping file \nmismatch::Int64=0 (optional): Number of allowed mismatches in barcode. Potential options include 0 or 1. If 1 mismatch, computation time will significantly increase. Default is to allow for 0 mismatches (exact matches only). \ndebug::Bool=false (optional): If true, a log file will be created and debugging data will be printed while the function is running (default is false).\n\nExample:\n\ndemultiplex_se(\"multiplexreads.fastq\", \"mapping_file.fastq\")\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.demultiplex_pe","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.demultiplex_pe","text":"Nucleotide_Essentials.demultiplex_pe\n\nCompares a list of provided barcodes with the provided paired-end multiplexed reads and separates the reads into individual .fastq files. If a barcode is found within R1 reads, the barcode is removed from the sequence. The quality data of the reads is preserved and written to the outputted .fastq file. If a barcode is not found, the sequnce and quality is written to the R1 unassigned .fastq file unchanged. If a barcode is found within R2 reads, the barcode is removed from the sequence. The quality data of the reads is preserved and written to the outputted .fastq file. If a barcode is not found, the sequnce and quality is written to the R2 unassigned .fastq file unchanged. \n\nDual-indexed reads are not yet supported\n\nThe mapping file must be either a .csv or .txt file with two columns. The first column heading must be SampleID and the second column heading must be BarcodeSequence. \n\nEXAMPLE MAPPING FILE: \n\nSampleID BarcodeSequence\nSample1 Barcode1\nSample2 Barcode2\nSample3 Barcode3\nSample4 Barcode4\nSample5 Barcode5\nSample6 Barcode6\nSample7 Barcode7\nSample8 Barcode8\n\nSupported keyword arguments include:\n\nR1::String: Path to forward multiplexed reads   \nR2::String: Path to reverse multiplexed reads   \nMap::String: Path to the mapping file \nmismatch::Int64=0 (optional): Number of allowed mismatches in barcode. Potential options include 0 or 1. If 1 mismatch, computation time will significantly increase. Default is to allow for 0 mismatches (exact matches only). \ndebug::Bool=false (optional): If true, a log file will be created and debugging data will be printed while the function is running (default is false).\n\nExample:\n\ndemultiplex_pe(\"forward_multiplexreads.fastq\", \"reverse_multiplexreads.fastq\", \"mapping_file.fastq\")\n\n\n\n\n\n","category":"function"},{"location":"#Index","page":"Nucleotide_Essentials.jl","title":"Index","text":"","category":"section"},{"location":"","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.jl","text":"","category":"page"},{"location":"#Change-Log","page":"Nucleotide_Essentials.jl","title":"Change Log","text":"","category":"section"},{"location":"#Nucleotide_Essentials-v0.2.0","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials v0.2.0","text":"","category":"section"},{"location":"","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.jl","text":"Added support for quality filtering of .fastq reads \nAdded support for Gzip compressed files \nPerformance improvements in PlotQuality() and added support for exporting quality plots\nAdded support for automatic quality profile encoding detection (Phred+64 and Phred+33 encoding)\nMinor documentation updates ","category":"page"}]
}
