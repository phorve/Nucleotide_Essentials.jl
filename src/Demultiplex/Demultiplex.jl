using DataFrames, CSV, Glob, Logging, StatsBase

"""
    Nucleotide_Essentials.potential_mismatches

Returns an ```Vector{Any}``` of potential barcodes with a single nucleotide change, including both deletions and substitutions

supported keyword arguments include:

* 'Path::String': The full or relative path to a .fastq file
* 'mismatch::Int64': The number of altered nucleotides to include (1 is only supported at this time)

# Example:
```julia
potential_mismatches("GCGT", 1)
17-element Vector{Any}:
"GCGT"
"CCGT" 
"ACGT" 
"TCGT" 
"GGGT" 
"GAGT" 
"GTGT" 
"GCCT" 
"GCAT" 
"GCTT" 
"GCGG" 
"GCGC" 
"GCGA" 
"CGT"
"GGT"
"GCT" 
"GCG"
```
"""
function potential_mismatches(og_barcode::Union{String15, String}, mismatch::Int64)
    potentials = []
    barcode = collect(og_barcode)
    if mismatch == 1
        for spot in 1:length(barcode)
            barcode = collect(og_barcode)
            barcode[spot] = 'G'
            push!(potentials, String(barcode))
            barcode = collect(og_barcode)
            barcode[spot] = 'C'
            push!(potentials, String(barcode))
            barcode = collect(og_barcode)
            barcode[spot] = 'A'
            push!(potentials, String(barcode))
            barcode = collect(og_barcode)
            barcode[spot] = 'T'
            push!(potentials, String(barcode))
        end 
        for spot in 1:length(barcode)    
            # Now handle potential deletions in the barcode 
            barcode = collect(og_barcode)
            deleteat!(barcode, spot)
            del_barcode = join(barcode)
            push!(potentials, String(barcode))
        end 
        potentials = unique(potentials)
    end
    return potentials 
end 

"""
    Nucleotide_Essentials.reverse_complement

Takes a string of nucleotide bases and returns the reverse complement of that string. Accepts inputs of ```String``` and ```SubString{String}``` (input from a FastqRecord)

supported keyword arguments include:

* 'sequence::Union{String, SubString{String}}': A string sequence of nucleotide bases or sequence entry from a FastqRecord

# Example: 
```julia
reverse_complement("ATCGT")
"ACGAT"
```
"""
function reverse_complement(sequence::Union{String, SubString{String}})
    rc = reverse(replace(sequence, 
                "C" => "G",
                "G" => "C",
                "A" => "T",
                "T" => "A"))
end 

"""
    Nucleotide_Essentials.demultiplex_se

Compares a list of provided barcodes with the provided multiplexed reads and separates the reads into individual .fastq files. If a barcode is found within the read, the barcode is removed from the sequence. The quality data of the reads is preserved and written to the outputted .fastq file. If a barcode is not found, the sequnce and quality is written to the unassigned .fastq file unchanged. 

The mapping file must be either a .csv or .txt file with two columns. The first column heading must be `SampleID` and the second column heading must be `BarcodeSequence`. 

EXAMPLE MAPPING FILE: 

| SampleID | BarcodeSequence |
|----------|-----------------|
| Sample1  | Barcode1        |
| Sample2  | Barcode2        |
| Sample3  | Barcode3        |
| Sample4  | Barcode4        |
| Sample5  | Barcode5        |
| Sample6  | Barcode6        |
| Sample7  | Barcode7        |
| Sample8  | Barcode8        |

supported keyword arguments include:

* 'R1::String': Path to multiplexed reads   
* 'Map::String': Path to the mapping file 
* 'mismatch::Int64=0' (optional): Number of allowed mismatches in barcode. Potential options include 0 or 1. If 1 mismatch, computation time will significantly increase. Default is to allow for 0 mismatches (exact matches only). 
* 'debug::Bool=false' (optional): If true, a log file will be created and debugging data will be printed while the function is running (default is false).

# Example: 
```julia
demultiplex_se("multiplexreads.fastq", "mapping_file.fastq")
```
"""
function demultiplex_se(
    R1::String,
    Map::String,
    mismatch::Int64=0,
    debug::Bool=false
    )

    @debug("Making output directory")
    Dir = string(rsplit(abspath(R1), "/", limit = 2)[1])
    Out = "/demultiplex_output"
    mkdir(string(Dir, Out))

    if debug === true
        # Set up location for the log 
        ENV["JULIA_DEBUG"] = Main
        if isdir(string(Dir, "/logs")) === false
            mkdir(string(Dir, "/logs"))
        end
        io = open(string(Dir, "/logs/log_demultiplex.txt"), "w+")
        logger = ConsoleLogger(io)
        global_logger(logger)
    else
        Logging.disable_logging(Logging.Info)
    end

    pre_dmx = readFastq(R1) # these are our forward multiplexed reads
    seqs = pre_dmx.sequence
    qual = pre_dmx.quality
    og_name = pre_dmx.ID

    # set up our mapping file - this is our mapping file in txt format 
    if endswith(Map, ".csv") == true
        @debug("Reading in our mapping file - it is a .csv file")
        mapping = CSV.read(Map, DataFrame)
        barcodes = Array(mapping[!, :BarcodeSequence])
        barcode_lengths = []
        for b in barcodes
            hold_length = length(b)
            push!(barcode_lengths, hold_length)
        end
        barcode_length = median(barcode_lengths)
        IDs = Array(mapping[!, :SampleID])

        # Make sure the number of barcodes matches the number of sample IDs
        if length(unique(barcodes)) != length(unique(IDs))
            error("The number of barcodes and unique samples should match, check mapping file for errors")
        end
    end

    # set up our mapping file - this is our mapping file in txt format 
    if endswith(Map, ".txt") == true
        @debug("Reading in our mapping file - it is a .txt file")
        mapping = CSV.read(Map, DataFrame)
        barcodes = Array(mapping[!, :BarcodeSequence])
        barcode_lengths = []
        for b in barcodes
            hold_length = length(b)
            push!(barcode_lengths, hold_length)
        end
        barcode_length = median(barcode_lengths)
        IDs = Array(mapping[!, :SampleID])

        # Make sure the number of barcodes matches the number of sample IDs
        if length(unique(barcodes)) != length(unique(IDs))
            error("The number of barcodes and unique samples should match, check mapping file for errors")
        end
    end

    # dmuxing
    @debug("Beginning demultiplexing algorithm")
    for s in 1:length(seqs)
        # demultiplexing not allowing for any mismatches 
        if mismatch == 0
            @debug(string("Beginning demultiplexing for sequence #", s, " of ", length(seqs)))
            lookat = Integer(barcode_length + 5)
            this_seq = first(seqs[s], lookat)
            for b in 1:length(barcodes)
                @debug(string("Querying against barcode #", b))
                if occursin(barcodes[b], this_seq)
                    sample = string(IDs[b])
                    file_end = string("_R1.fastq")
                    if isfile(string(Dir, "/", Out, "/", sample, file_end)) == false
                        @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                        out_string = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                        write(string(Dir, "/", Out, "/", sample, file_end), out_string)
                        break
                    else
                        @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                        out_string2 = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                        iostream = open(string(Dir, "/", Out, "/", sample, file_end), "a")
                        write(iostream, out_string2)
                        close(iostream)
                        break
                    end
                elseif b == length(barcodes)
                    @debug("Beginning last barcode")
                    if occursin(barcodes[b], this_seq)
                        sample = string(IDs[b])
                        file_end = string("_R1.fastq")
                        if isfile(string(Dir, "/", Out, "/", sample, file_end)) == false
                            @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                            out_string = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                            write(string(Dir, "/", Out, "/", sample, file_end), out_string)
                            break
                        else
                            @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                            out_string2 = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                            iostream = open(string(Dir, "/", Out, "/", sample, file_end), "a")
                            write(iostream, out_string2)
                            close(iostream)
                            break
                        end
                    else
                        if isfile(string(Dir, "/", Out, "/unassigned_R1.fastq")) == false
                            @debug("No unassigned Fastq file detected - Creating file")
                            @debug("Reached last barcode - writing to unassigned Fastq file.")
                            out_string3 = string("@", og_name[s], "\n", seqs[s], "\n+\n", qual[s])
                            write(string(Dir, "/", Out, "/unassigned_R1.fastq"), out_string3)
                        else
                            @debug("Reached last barcode - writing to unassigned Fastq file.")
                            out_string4 = string("\n@", og_name[s], "\n", seqs[s], "\n+\n", qual[s])
                            iostream = open(string(Dir, "/", Out, "/unassigned_R1.fastq"), "a")
                            write(iostream, out_string4)
                            close(iostream)
                        end
                    end
                else
                    continue
                end
            end

            # demultiplexing allowing for mismatches 
        else
            @debug(string("Beginning demultiplexing for sequence #", s, " of ", length(seqs)))
            lookat = Integer(barcode_length + 5)
            this_seq = first(seqs[s], lookat)
            for b in 1:length(barcodes)
                @debug(string("Building barcode pool for barcode #", b, ": ", barcodes[b]))
                temp_barcode_pool = potential_mismatches(barcodes[b], mismatch)
                @debug("Finished building potential barcode pool")
                for drop in 1:length(temp_barcode_pool)
                    @debug(string("Querying for variant #", drop, " from barcode #", b))
                    if occursin(temp_barcode_pool[drop], this_seq)
                        sample = string(IDs[b])
                        file_end = string("_R1.fastq")
                        if isfile(string(Dir, "/", Out, "/", sample, file_end)) == false
                            @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                            out_string = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                            write(string(Dir, "/", Out, "/", sample, file_end), out_string)
                            break
                        else
                            @debug("Else after if 1.1")
                            @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                            out_string2 = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                            iostream = open(string(Dir, "/", Out, "/", sample, file_end), "a")
                            write(iostream, out_string2)
                            close(iostream)
                            break
                        end
                    elseif drop == length(temp_barcode_pool)
                        @debug("Beginning last variant")
                        if occursin(temp_barcode_pool[drop], this_seq)
                            sample = string(IDs[b])
                            file_end = string("_R1.fastq")
                            if isfile(string(Dir, "/", Out, "/", sample, file_end)) == false
                                @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                                out_string = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                                write(string(Dir, "/", Out, "/", sample, file_end), out_string)
                                break
                            else
                                @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                                out_string2 = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                                iostream = open(string(Dir, "/", Out, "/", sample, file_end), "a")
                                write(iostream, out_string2)
                                close(iostream)
                                break
                            end
                        else
                            if b == length(barcodes)
                                @debug("Reached last barcode - writing to unassigned Fastq file.")
                                if isfile(string(Dir, "/", Out, "/unassigned_R1.fastq")) == false
                                    @debug("No unassigned Fastq file detected - Creating file")
                                    out_string3 = string("@", og_name[s], "\n", seqs[s], "\n+\n", qual[s])
                                    write(string(Dir, "/", Out, "/unassigned_R1.fastq"), out_string3)
                                else
                                    @debug("Reached last barcode - writing to unassigned Fastq file.")
                                    out_string4 = string("\n@", og_name[s], "\n", seqs[s], "\n+\n", qual[s])
                                    iostream = open(string(Dir, "/", Out, "/unassigned_R1.fastq"), "a")
                                    write(iostream, out_string4)
                                    close(iostream)
                                end
                            else
                                @debug("Reached last variant in barcode pool - moving on to the next pool")
                                continue
                            end
                        end
                    else
                        @debug("I'm not sure?")
                        continue
                    end
                end
            end
        end
    end

    # output statistics on demultiplexed samples 
    sampleholds = []
    reads = []
    files = glob("*.fastq", string(Dir, "/", Out))
    for file in files
        push!(sampleholds, rsplit(file, "/", limit=2)[2])
        push!(reads, (countlines(file) / 4))
    end
    stats = DataFrame(Filename=sampleholds, AssignedReads=reads)
    CSV.write(string(Dir, "/", Out, "/demultiplex_stats.csv"), stats)

    if debug === true
        # stop logger
        close(io)
    end
end # end function demultiplex_se

"""
    Nucleotide_Essentials.demultiplex_pe

Compares a list of provided barcodes with the provided paired-end multiplexed reads and separates the reads into individual .fastq files. If a barcode is found within R1 reads, the barcode is removed from the sequence. The quality data of the reads is preserved and written to the outputted .fastq file. If a barcode is not found, the sequnce and quality is written to the R1 unassigned .fastq file unchanged. If a barcode is found within R2 reads, the barcode is removed from the sequence. The quality data of the reads is preserved and written to the outputted .fastq file. If a barcode is not found, the sequnce and quality is written to the R2 unassigned .fastq file unchanged. 

*Dual-indexed reads are not _yet_ supported*

The mapping file must be either a .csv or .txt file with two columns. The first column heading must be `SampleID` and the second column heading must be `BarcodeSequence`. 

EXAMPLE MAPPING FILE: 

| SampleID | BarcodeSequence |
|----------|-----------------|
| Sample1  | Barcode1        |
| Sample2  | Barcode2        |
| Sample3  | Barcode3        |
| Sample4  | Barcode4        |
| Sample5  | Barcode5        |
| Sample6  | Barcode6        |
| Sample7  | Barcode7        |
| Sample8  | Barcode8        |

supported keyword arguments include:

* 'R1::String': Path to forward multiplexed reads   
* 'R2::String': Path to reverse multiplexed reads   
* 'Map::String': Path to the mapping file 
* 'mismatch::Int64=0' (optional): Number of allowed mismatches in barcode. Potential options include 0 or 1. If 1 mismatch, computation time will significantly increase. Default is to allow for 0 mismatches (exact matches only). 
* 'debug::Bool=false' (optional): If true, a log file will be created and debugging data will be printed while the function is running (default is false).

# Example: 
```julia
demultiplex_pe("forward_multiplexreads.fastq", "reverse_multiplexreads.fastq", "mapping_file.fastq")
```
"""
function demultiplex_pe(
    R1::String,
    R2::String,
    Map::String,
    mismatch::Int64=0,
    debug::Bool=false
    )

    @debug("Making output directory")
    Dir = string(rsplit(abspath(R1), "/", limit = 2)[1])
    Out = "/demultiplex_output"
    mkdir(string(Dir, Out))

    if debug === true
        # Set up location for the log 
        ENV["JULIA_DEBUG"] = Main
        if isdir(string(Dir, "/logs")) === false
            mkdir(string(Dir, "/logs"))
        end
        io = open(string(Dir, "/logs/log_demultiplex.txt"), "w+")
        logger = ConsoleLogger(io)
        global_logger(logger)
    else
        Logging.disable_logging(Logging.Info)
    end

    pre_dmx = readFastq(R1) # these are our forward multiplexed reads
    pre_dmx2 = readFastq(R2) # these are our forward multiplexed reads
    seqs = pre_dmx.sequence
    seqs2 = pre_dmx2.sequence
    qual = pre_dmx.quality
    qual2 = pre_dmx2.quality
    og_name = pre_dmx.ID
    og_name2 = pre_dmx2.ID

    # set up our mapping file - this is our mapping file in txt format 
    if endswith(Map, ".csv") == true
        @debug("Reading in our mapping file - it is a .csv file")
        mapping = CSV.read(Map, DataFrame)
        barcodes = Array(mapping[!, :BarcodeSequence])
        barcode_lengths = []
        for b in barcodes
            hold_length = length(b)
            push!(barcode_lengths, hold_length)
        end
        barcode_length = median(barcode_lengths)
        IDs = Array(mapping[!, :SampleID])

        # Make sure the number of barcodes matches the number of sample IDs
        if length(unique(barcodes)) != length(unique(IDs))
            error("The number of barcodes and unique samples should match, check mapping file for errors")
        end
    end

    # set up our mapping file - this is our mapping file in txt format 
    if endswith(Map, ".txt") == true
        @debug("Reading in our mapping file - it is a .txt file")
        mapping = CSV.read(Map, DataFrame)
        barcodes = Array(mapping[!, :BarcodeSequence])
        barcode_lengths = []
        for b in barcodes
            hold_length = length(b)
            push!(barcode_lengths, hold_length)
        end
        barcode_length = median(barcode_lengths)
        IDs = Array(mapping[!, :SampleID])

        # Make sure the number of barcodes matches the number of sample IDs
        if length(unique(barcodes)) != length(unique(IDs))
            error("The number of barcodes and unique samples should match, check mapping file for errors")
        end
    end

    # dmuxing forward reads
    @debug("Beginning demultiplexing algorithm for forward reads")
    for s in 1:length(seqs)
        # demultiplexing not allowing for any mismatches 
        if mismatch == 0
            @debug(string("Beginning demultiplexing for sequence #", s, " of ", length(seqs)))
            lookat = Integer(barcode_length + 5)
            this_seq = first(seqs[s], lookat)
            for b in 1:length(barcodes)
                @debug(string("Querying against barcode #", b))
                if occursin(barcodes[b], this_seq)
                    sample = string(IDs[b])
                    file_end = string("_R1.fastq")
                    if isfile(string(Dir, "/", Out, "/", sample, file_end)) == false
                        @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                        out_string = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                        write(string(Dir, "/", Out, "/", sample, file_end), out_string)
                        break
                    else
                        @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                        out_string2 = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                        iostream = open(string(Dir, "/", Out, "/", sample, file_end), "a")
                        write(iostream, out_string2)
                        close(iostream)
                        break
                    end
                elseif b == length(barcodes)
                    @debug("Beginning last barcode")
                    if occursin(barcodes[b], this_seq)
                        sample = string(IDs[b])
                        file_end = string("_R1.fastq")
                        if isfile(string(Dir, "/", Out, "/", sample, file_end)) == false
                            @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                            out_string = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                            write(string(Dir, "/", Out, "/", sample, file_end), out_string)
                            break
                        else
                            @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                            out_string2 = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                            iostream = open(string(Dir, "/", Out, "/", sample, file_end), "a")
                            write(iostream, out_string2)
                            close(iostream)
                            break
                        end
                    else
                        if isfile(string(Dir, "/", Out, "/unassigned_R1.fastq")) == false
                            @debug("No unassigned Fastq file detected - Creating file")
                            @debug("Reached last barcode - writing to unassigned Fastq file.")
                            out_string3 = string("@", og_name[s], "\n", seqs[s], "\n+\n", qual[s])
                            write(string(Dir, "/", Out, "/unassigned_R1.fastq"), out_string3)
                        else
                            @debug("Reached last barcode - writing to unassigned Fastq file.")
                            out_string4 = string("\n@", og_name[s], "\n", seqs[s], "\n+\n", qual[s])
                            iostream = open(string(Dir, "/", Out, "/unassigned_R1.fastq"), "a")
                            write(iostream, out_string4)
                            close(iostream)
                        end
                    end
                else
                    continue
                end
            end

            # demultiplexing allowing for mismatches 
        else
            @debug(string("Beginning demultiplexing for sequence #", s, " of ", length(seqs)))
            lookat = Integer(barcode_length + 5)
            this_seq = first(seqs[s], lookat)
            for b in 1:length(barcodes)
                @debug(string("Building barcode pool for barcode #", b, ": ", barcodes[b]))
                temp_barcode_pool = potential_mismatches(barcodes[b], mismatch)
                @debug("Finished building potential barcode pool")
                for drop in 1:length(temp_barcode_pool)
                    @debug(string("Querying for variant #", drop, " from barcode #", b))
                    if occursin(temp_barcode_pool[drop], this_seq)
                        sample = string(IDs[b])
                        file_end = string("_R1.fastq")
                        if isfile(string(Dir, "/", Out, "/", sample, file_end)) == false
                            @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                            out_string = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                            write(string(Dir, "/", Out, "/", sample, file_end), out_string)
                            break
                        else
                            @debug("Else after if 1.1")
                            @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                            out_string2 = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                            iostream = open(string(Dir, "/", Out, "/", sample, file_end), "a")
                            write(iostream, out_string2)
                            close(iostream)
                            break
                        end
                    elseif drop == length(temp_barcode_pool)
                        @debug("Beginning last variant")
                        if occursin(temp_barcode_pool[drop], this_seq)
                            sample = string(IDs[b])
                            file_end = string("_R1.fastq")
                            if isfile(string(Dir, "/", Out, "/", sample, file_end)) == false
                                @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                                out_string = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                                write(string(Dir, "/", Out, "/", sample, file_end), out_string)
                                break
                            else
                                @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                                out_string2 = string("@", IDs[b], "\n", chop(seqs[s], head = Integer(barcode_length)), "\n+\n", chop(qual[s], head = Integer(barcode_length)))
                                iostream = open(string(Dir, "/", Out, "/", sample, file_end), "a")
                                write(iostream, out_string2)
                                close(iostream)
                                break
                            end
                        else
                            if b == length(barcodes)
                                @debug("Reached last barcode - writing to unassigned Fastq file.")
                                if isfile(string(Dir, "/", Out, "/unassigned_R1.fastq")) == false
                                    @debug("No unassigned Fastq file detected - Creating file")
                                    out_string3 = string("@", og_name[s], "\n", seqs[s], "\n+\n", qual[s])
                                    write(string(Dir, "/", Out, "/unassigned_R1.fastq"), out_string3)
                                else
                                    @debug("Reached last barcode - writing to unassigned Fastq file.")
                                    out_string4 = string("\n@", og_name[s], "\n", seqs[s], "\n+\n", qual[s])
                                    iostream = open(string(Dir, "/", Out, "/unassigned_R1.fastq"), "a")
                                    write(iostream, out_string4)
                                    close(iostream)
                                end
                            else
                                @debug("Reached last variant in barcode pool - moving on to the next pool")
                                continue
                            end
                        end
                    else
                        @debug("I'm not sure?")
                        continue
                    end
                end
            end
        end
    end

    # dmuxing reverse reads
    @debug("Beginning demultiplexing algorithm for reverse reads")
    for s in 1:length(seqs2)
        # demultiplexing not allowing for any mismatches 
        if mismatch == 0
            @debug(string("Beginning demultiplexing for sequence #", s, " of ", length(seqs2)))
            lookat = Integer(barcode_length + 5)
            this_line = reverse_complement(seqs2[s])
            this_seq = (first(this_line, lookat))
            for b in 1:length(barcodes)
                @debug(string("Querying against barcode #", b))
                if occursin(barcodes[b], this_seq)
                    sample = string(IDs[b])
                    file_end = string("_R2.fastq")
                    if isfile(string(Dir, "/", Out, "/", sample, file_end)) == false
                        @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                        out_string = string("@", IDs[b], "\n", (chop(seqs2[s], tail = Integer(barcode_length))), "\n+\n", (chop(qual2[s], tail = Integer(barcode_length))))
                        write(string(Dir, "/", Out, "/", sample, file_end), out_string)
                        break
                    else
                        @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                        out_string2 = string("@", IDs[b], "\n", (chop(seqs2[s], tail = Integer(barcode_length))), "\n+\n", (chop(qual2[s], tail = Integer(barcode_length))))
                        iostream = open(string(Dir, "/", Out, "/", sample, file_end), "a")
                        write(iostream, out_string2)
                        close(iostream)
                        break
                    end
                elseif b == length(barcodes)
                    @debug("Beginning last barcode")
                    if occursin(barcodes[b], this_seq)
                        sample = string(IDs[b])
                        file_end = string("_R2.fastq")
                        if isfile(string(Dir, "/", Out, "/", sample, file_end)) == false
                            @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                            out_string = string("@", IDs[b], "\n", (chop(seqs2[s], tail = Integer(barcode_length))), "\n+\n", (chop(qual2[s], tail = Integer(barcode_length))))
                            write(string(Dir, "/", Out, "/", sample, file_end), out_string)
                            break
                        else
                            @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                            out_string2 = string("@", IDs[b], "\n", (chop(seqs2[s], tail = Integer(barcode_length))), "\n+\n", (chop(qual2[s], tail = Integer(barcode_length))))
                            iostream = open(string(Dir, "/", Out, "/", sample, file_end), "a")
                            write(iostream, out_string2)
                            close(iostream)
                            break
                        end
                    else
                        if isfile(string(Dir, "/", Out, "/unassigned_R2.fastq")) == false
                            @debug("No unassigned Fastq file detected - Creating file")
                            @debug("Reached last barcode - writing to unassigned Fastq file.")
                            out_string3 = string("@", og_name2[s], "\n", seqs2[s], "\n+\n", qual2[s])
                            write(string(Dir, "/", Out, "/unassigned_R2.fastq"), out_string3)
                        else
                            @debug("Reached last barcode - writing to unassigned Fastq file.")
                            out_string4 = string("\n@", og_name2[s], "\n", seqs2[s], "\n+\n", qual2[s])
                            iostream = open(string(Dir, "/", Out, "/unassigned_R2.fastq"), "a")
                            write(iostream, out_string4)
                            close(iostream)
                        end
                    end
                else
                    continue
                end
            end

            # demultiplexing allowing for mismatches 
        else
            @debug(string("Beginning demultiplexing for sequence #", s, " of ", length(seqs2)))
            lookat = Integer(barcode_length + 5)
            this_line = reverse_complement(seqs2[s])
            this_seq = (first(this_line, lookat))
            for b in 1:length(barcodes)
                @debug(string("Building barcode pool for barcode #", b, ": ", barcodes[b]))
                temp_barcode_pool = potential_mismatches(barcodes[b], mismatch)
                @debug("Finished building potential barcode pool")
                for drop in 1:length(temp_barcode_pool)
                    @debug(string("Querying for variant #", drop, " from barcode #", b))
                    if occursin(temp_barcode_pool[drop], this_seq)
                        sample = string(IDs[b])
                        file_end = string("_R2.fastq")
                        if isfile(string(Dir, "/", Out, "/", sample, file_end)) == false
                            @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                            out_string = string("@", IDs[b], "\n", (chop(seqs2[s], tail = Integer(barcode_length))), "\n+\n", (chop(qual2[s], tail = Integer(barcode_length))))
                            write(string(Dir, "/", Out, "/", sample, file_end), out_string)
                            break
                        else
                            @debug("Else after if 1.1")
                            @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                            out_string2 = string("@", IDs[b], "\n", (chop(seqs2[s], tail = Integer(barcode_length))), "\n+\n", (chop(qual2[s], tail = Integer(barcode_length))))
                            iostream = open(string(Dir, "/", Out, "/", sample, file_end), "a")
                            write(iostream, out_string2)
                            close(iostream)
                            break
                        end
                    elseif drop == length(temp_barcode_pool)
                        @debug("Beginning last variant")
                        if occursin(temp_barcode_pool[drop], this_seq)
                            sample = string(IDs[b])
                            file_end = string("_R2.fastq")
                            if isfile(string(Dir, "/", Out, "/", sample, file_end)) == false
                                @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                                out_string = string("@", IDs[b], "\n", (chop(seqs2[s], tail = Integer(barcode_length))), "\n+\n", (chop(qual2[s], tail = Integer(barcode_length))))
                                write(string(Dir, "/", Out, "/", sample, file_end), out_string)
                                break
                            else
                                @debug(string("Barcode found - writing to fastq file for ", IDs[b]))
                                out_string2 = string("@", IDs[b], "\n", (chop(seqs2[s], tail = Integer(barcode_length))), "\n+\n", (chop(qual2[s], tail = Integer(barcode_length))))
                                iostream = open(string(Dir, "/", Out, "/", sample, file_end), "a")
                                write(iostream, out_string2)
                                close(iostream)
                                break
                            end
                        else
                            if b == length(barcodes)
                                @debug("Reached last barcode - writing to unassigned Fastq file.")
                                if isfile(string(Dir, "/", Out, "/unassigned_R2.fastq")) == false
                                    @debug("No unassigned Fastq file detected - Creating file")
                                    out_string3 = string("@", og_name2[s], "\n", seqs2[s], "\n+\n", qual2[s])
                                    write(string(Dir, "/", Out, "/unassigned_R2.fastq"), out_string3)
                                else
                                    @debug("Reached last barcode - writing to unassigned Fastq file.")
                                    out_string4 = string("\n@", og_name2[s], "\n", seqs2[s], "\n+\n", qual2[s])
                                    iostream = open(string(Dir, "/", Out, "/unassigned_R2.fastq"), "a")
                                    write(iostream, out_string4)
                                    close(iostream)
                                end
                            else
                                @debug("Reached last variant in barcode pool - moving on to the next pool")
                                continue
                            end
                        end
                    else
                        @debug("I'm not sure?")
                        continue
                    end
                end
            end
        end
    end

    # output statistics on demultiplexed samples 
    sampleholds = []
    reads = []
    files = glob("*.fastq", string(Dir, "/", Out))
    for file in files
        push!(sampleholds, rsplit(file, "/", limit=2)[2])
        push!(reads, (countlines(file) / 4))
    end
    stats = DataFrame(Filename=sampleholds, AssignedReads=reads)
    CSV.write(string(Dir, "/", Out, "/demultiplex_stats.csv"), stats)

    if debug === true
        # stop logger
        close(io)
    end
end # end function demultiplex_pe