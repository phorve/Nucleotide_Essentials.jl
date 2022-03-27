var documenterSearchIndex = {"docs":
[{"location":"#Nucleotide_Essentials.jl","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.jl","text":"","category":"section"},{"location":"","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.jl","text":"","category":"page"},{"location":"#Data-Types","page":"Nucleotide_Essentials.jl","title":"Data Types","text":"","category":"section"},{"location":"","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.jl","text":"FastqRecord\nFastaRecord","category":"page"},{"location":"#Nucleotide_Essentials.FastqRecord","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.FastqRecord","text":"Nucleotide_Essentials.FastqRecord\n\nComponents\n\nID: The unique sequence identifier associated with that entry\nsequence: The nucleotide sequence of that entry\nquality: The quality scores of that entry \nfilename: The original file name \n\n\n\n\n\n","category":"type"},{"location":"#Nucleotide_Essentials.FastaRecord","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.FastaRecord","text":"Nucleotide_Essentials.FastaRecord\n\nComponents\n\nID: The unique sequence identifier associated with that entry\nsequence: The nucleotide sequence of that entry\nfilename: The original file name \n\n\n\n\n\n","category":"type"},{"location":"#Functions","page":"Nucleotide_Essentials.jl","title":"Functions","text":"","category":"section"},{"location":"","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.jl","text":"readFastq\nreadFasta\nPlotQuality\npotential_mismatches\nreverse_complement\ndemultiplex_se\ndemultiplex_pe","category":"page"},{"location":"#Nucleotide_Essentials.readFastq","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.readFastq","text":"Nucleotide_Essentials.readFastq\n\nreadFastq(Path::String)\n\n.fastq file => readFastq(Path) => FastqRecord(ID, sequence, quality, filename)\n\nsupported keyword arguments include: \n\n'Path::String': The full or relative path to a .fastq file\n\nExample:\n\n# Supply the path to a .Fastq file that you would like to import\nmyfastq = readFastq(\"myfastq.fastq\")\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.readFasta","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.readFasta","text":"Nucleotide_Essentials.readFasta\n\nreadFasta(Path::String)\n\n.fasta file => readFasta(Path) => FastaRecord(ID, sequence, filename)\n\nsupported keyword arguments include: \n\n'Path::String': The full or relative path to a .fastq file\n\nExample:\n\n# Supply the path to a .Fastq file that you would like to import - it is recommended to include `;` in your command to prevent printing potentially large .fasta files in the REPL\nmyfasta = readFasta(\"myfasta.fasta\");\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.PlotQuality","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.PlotQuality","text":"Nucleotide_Essentials.PlotQuality\n\nReturns a plot of the quality profile of a FastqRecord\n\nThis function plots a visual summary of the distribution of quality scores as a function of sequence position for the input fastq file(s).\n\nThe plotted lines show summary statistics at each sequence position:\n\ngreen is the mean \ndashed red lines are the 25th and 75th quantiles\n\nsupported keyword arguments include:\n\n'Input::FastqRecord': The name of a FastqRecord for plotting\n'QualityType::String' (optional): The system to use for parsing the fastq quality information. Potential quality encodings include: \nPhred64 (defualt)\nAscii (not yet supported)\nPhred33 (not yet supported)    \n'verbose::Bool' (optional): Whether or not to show some intermediary feedback on the progress of the function (default = false)\n\nExample:\n\n# A quality profile can be created by directly calling an already read FastqRecord\nPlotQuality(myfastq)\n\n# Alternatively, a quality profile can be created by nesting readFastq() inside of PlotQuality() \nPlotQuality(readFastq(\"myfastq.fastq\"))\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.potential_mismatches","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.potential_mismatches","text":"Nucleotide_Essentials.potential_mismatches\n\nReturns an Vector{Any} of potential barcodes with a single nucleotide change, including both deletions and substitutions\n\nsupported keyword arguments include:\n\n'Path::String': The full or relative path to a .fastq file\n'mismatch::Int64': The number of altered nucleotides to include (1 is only supported at this time)\n\nExample:\n\npotential_mismatches(\"GCGT\", 1)\n17-element Vector{Any}:\n\"GCGT\"\n\"CCGT\" \n\"ACGT\" \n\"TCGT\" \n\"GGGT\" \n\"GAGT\" \n\"GTGT\" \n\"GCCT\" \n\"GCAT\" \n\"GCTT\" \n\"GCGG\" \n\"GCGC\" \n\"GCGA\" \n\"CGT\"\n\"GGT\"\n\"GCT\" \n\"GCG\"\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.reverse_complement","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.reverse_complement","text":"Nucleotide_Essentials.reverse_complement\n\nTakes a string of nucleotide bases and returns the reverse complement of that string. Accepts inputs of String and SubString{String} (input from a FastqRecord)\n\nsupported keyword arguments include:\n\n'sequence::Union{String, SubString{String}}': A string sequence of nucleotide bases or sequence entry from a FastqRecord\n\nExample:\n\nreverse_complement(\"ATCGT\")\n\"ACGAT\"\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.demultiplex_se","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.demultiplex_se","text":"Nucleotide_Essentials.demultiplex_se\n\nCompares a list of provided barcodes with the provided multiplexed reads and separates the reads into individual .fastq files. If a barcode is found within the read, the barcode is removed from the sequence. The quality data of the reads is preserved and written to the outputted .fastq file. If a barcode is not found, the sequnce and quality is written to the unassigned .fastq file unchanged. \n\nThe mapping file must be either a .csv or .txt file with two columns. The first column heading must be SampleID and the second column heading must be BarcodeSequence. \n\nEXAMPLE MAPPING FILE: \n\nSampleID BarcodeSequence\nSample1 Barcode1\nSample2 Barcode2\nSample3 Barcode3\nSample4 Barcode4\nSample5 Barcode5\nSample6 Barcode6\nSample7 Barcode7\nSample8 Barcode8\n\nsupported keyword arguments include:\n\n'R1::String': Path to multiplexed reads   \n'Map::String': Path to the mapping file \n\n'mismatch::Int64=0' (optional): Number of allowed mismatches in barcode. Potential options include 0 or 1. If 1 mismatch, computation time will significantly increase. Default is to allow for 0 mismatches (exact matches only).\n\n'debug::Bool=false' (optional): If true, a log file will be created and debugging data will be printed while the function is running (default is false).\n\nExample:\n\ndemultiplex_se(\"multiplexreads.fastq\", \"mapping_file.fastq\")\n\n\n\n\n\n","category":"function"},{"location":"#Nucleotide_Essentials.demultiplex_pe","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.demultiplex_pe","text":"Nucleotide_Essentials.demultiplex_pe\n\nCompares a list of provided barcodes with the provided paired-end multiplexed reads and separates the reads into individual .fastq files. If a barcode is found within R1 reads, the barcode is removed from the sequence. The quality data of the reads is preserved and written to the outputted .fastq file. If a barcode is not found, the sequnce and quality is written to the R1 unassigned .fastq file unchanged. If a barcode is found within R2 reads, the barcode is removed from the sequence. The quality data of the reads is preserved and written to the outputted .fastq file. If a barcode is not found, the sequnce and quality is written to the R2 unassigned .fastq file unchanged. \n\nDual-indexed reads are not yet supported\n\nThe mapping file must be either a .csv or .txt file with two columns. The first column heading must be SampleID and the second column heading must be BarcodeSequence. \n\nEXAMPLE MAPPING FILE: \n\nSampleID BarcodeSequence\nSample1 Barcode1\nSample2 Barcode2\nSample3 Barcode3\nSample4 Barcode4\nSample5 Barcode5\nSample6 Barcode6\nSample7 Barcode7\nSample8 Barcode8\n\nsupported keyword arguments include:\n\n'R1::String': Path to forward multiplexed reads   \n'R2::String': Path to reverse multiplexed reads   \n'Map::String': Path to the mapping file \n'mismatch::Int64=0' (optional): Number of allowed mismatches in barcode. Potential options include 0 or 1. If 1 mismatch, computation time will significantly increase. Default is to allow for 0 mismatches (exact matches only). \n'debug::Bool=false' (optional): If true, a log file will be created and debugging data will be printed while the function is running (default is false).\n\nExample:\n\ndemultiplex_pe(\"forward_multiplexreads.fastq\", \"reverse_multiplexreads.fastq\", \"mapping_file.fastq\")\n\n\n\n\n\n","category":"function"},{"location":"#Index","page":"Nucleotide_Essentials.jl","title":"Index","text":"","category":"section"},{"location":"","page":"Nucleotide_Essentials.jl","title":"Nucleotide_Essentials.jl","text":"","category":"page"}]
}
