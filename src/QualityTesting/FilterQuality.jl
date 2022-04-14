using DataFrames, DataFramesMeta
import Base.Threads.@threads
"""
    Nucleotide_Essentials.FilterQuality_se

Filters an input .fastq file based upon the encoded Phred+33 or Phred+64 quality scores. The encoding of the reads is automatically deteremined by looking for unique encoding in Phred+33 and Phred+64. Phred+64 encoding is identified by searching for `^`, `a`, `]`, and `f`.

Reads are filtered based upon the number of [expected errors](https://academic.oup.com/bioinformatics/article/31/21/3476/194979?login=true) (``\\mathrm{E}``) based on the error rate based on quality score and the sum of error probabilities, following the equation: 

``\\mathrm{E} = \\sum{_ip_i} = \\sum{_i}10^{\\frac{-Q_i}{10}}``

Stringent filtering (```maxEE = 1```) is used by default but can be adjusted by the user. 

Reads that pass the filtering parameters are output to a file ending in `_FilteredReads.fastq` in the user-determined directory, as indicated by ```out```. 

Supported keyword arguments include: 
    
* ```read1::String```: Path to the reads to undergo quality filtering 
* ```out::String```: Path to the directory where reads that pass the quality filtering should be written 
* ```maxEE::Int64``` (optional): The max number of expected errors a read can include as the filtering parameter (default: ```maxEE = 1```)
* ```verbose::Bool``` (optional): Whether or not to show some intermediary feedback on the progress of the function (default = false)

# Example: 
```julia
FilterQuality_se("forward_R1.fasta", "/outdirectory")
```
"""
function FilterQuality_se(read1::String, out::String, maxEE::Int64 = 1, verbose::Bool = false)

    if endswith(read1, ".gz") === true
        if verbose == true
            println(string("Unzipping ", read1))
        end 
        run(`gunzip $read1`); #this works
        unzipped = true
    else
        unzipped = false
    end 

    R1 = readFastq(read1)
    quals = R1.quality
    seqs = R1.sequence
    IDs = R1.ID
    Passing = 0
    Failing = 0

    if occursin("^", join(quals)) == true && occursin("a", join(quals)) == true && occursin("]", join(quals)) == true && occursin("f", join(quals)) == true
        encoding = 64;
        println("Detected Phred+64 quality encoding")
    end 
    if occursin("^", join(quals)) == false && occursin("a", join(quals)) == false && occursin("]", join(quals)) == false && occursin("f", join(quals)) == false
        println("Detected Phred+33 quality encoding")
        encoding = 33;
    end 

    if encoding == 33
        scores_33 = Dict(
            ('!') => "1.0000", 
            ('\"') => "0.79433",
            ('#') => "0.63096",
            ('\$') => "0.50119",
            ('%') => "0.39811",
            ('&') => "0.31623",
            ('\'') => "0.25119",
            ('(') => "0.19953",
            (')') => "0.15849",
            ('*') => "0.12589",
            ('+') => "0.10000",
            (',') => "0.07943",
            ('-') => "0.06310",
            ('.') => "0.05012",
            ('/') => "0.03981",
            ('0') => "0.03162",
            ('1') => "0.02512",
            ('2') => "0.01995",
            ('3') => "0.01585",
            ('4') => "0.01259",
            ('5') => "0.01000",
            ('6') => "0.00794",
            ('7') => "0.00631",
            ('8') => "0.00501",
            ('9') => "0.00398",
            (':') => "0.00316",
            (';') => "0.00251",
            ('<') => "0.00200",
            ('=') => "0.00158",
            ('>') => "0.00126",
            ('?') => "0.00100",
            ('@') => "0.00079",
            ('A') => "0.00063",
            ('B') => "0.00050",
            ('C') => "0.00040",
            ('D') => "0.00032",
            ('E') => "0.00025",
            ('F') => "0.00020",
            ('G') => "0.00016",
            ('H') => "0.00013",
            ('I') => "0.00010", 
            ('J') => "0.00008",
            ('K') => "0.00006"
        )
    end 
    if encoding == 64
        scores_64 = Dict(
            ('@') => "1.0000", 
            ('A') => "0.79433",
            ('B') => "0.63096",
            ('C') => "0.50119",
            ('D') => "0.39811",
            ('E') => "0.31623",
            ('F') => "0.25119",
            ('G') => "0.19953",
            ('H') => "0.15849",
            ('I') => "0.12589",
            ('J') => "0.10000",
            ('K') => "0.07943",
            ('L') => "0.06310",
            ('M') => "0.05012",
            ('N') => "0.03981",
            ('O') => "0.03162",
            ('P') => "0.02512",
            ('Q') => "0.01995",
            ('R') => "0.01585",
            ('S') => "0.01259",
            ('T') => "0.01000",
            ('U') => "0.00794",
            ('V') => "0.00631",
            ('W') => "0.00501",
            ('X') => "0.00398",
            ('Y') => "0.00316",
            ('Z') => "0.00251",
            ('[') => "0.00200",
            ('\\') => "0.00158",
            (']') => "0.00126",
            ('^') => "0.00100",
            ('_') => "0.00079",
            ('`') => "0.00063",
            ('a') => "0.00050",
            ('b') => "0.00040",
            ('c') => "0.00032",
            ('d') => "0.00025",
            ('e') => "0.00020",
            ('f') => "0.00016",
            ('g') => "0.00013",
            ('h') => "0.00010", 
            ('i') => "0.00008",
            ('j') => "0.00006"
        )
    end 

    alreadyprinted = false
    alreadyprinted1 = false
    alreadyprinted2 = false
    alreadyprinted3 = false
    alreadyprinted4 = false
    alreadyprinted5 = false
    alreadyprinted6 = false
    alreadyprinted7 = false
    alreadyprinted8 = false
    alreadyprinted9 = false
    alreadyprinted10 = false
    alreadyprinted11= false
    alreadyprinted12 = false
    alreadyprinted13= false
    alreadyprinted14 = false
    alreadyprinted15 = false
    alreadyprinted16 = false
    alreadyprinted17 = false
    alreadyprinted18 = false

    if verbose == true
        println("Prepping the reads for quality filtering")
    end
    
    function convert_phred(c)
        if encoding == 33
            scores_33[c]
        else 
            scores_64[c]
        end 
    end

    for i in 1:length(quals)
        ee1 = sum(parse.(Float64, (map(convert_phred, collect(quals[i])))))
        if ee1 >= maxEE
            Failing = Failing + 1
        else 
            Passing = Passing + 1
            if isfile(string(out, "/",split(R1.filename, ".")[1], "_filtered.fastq")) == false
                out_string = string("@", IDs[i], "\n", seqs[i], "\n+\n", quals[i])
                write(string(out, "/",split(R1.filename, ".")[1], "_filtered.fastq"), out_string)
            else 
                out_string2 = string("\n@", IDs[i], "\n", seqs[i], "\n+\n", quals[i])
                iostream = open(string(out, "/",split(R1.filename, ".")[1], "_filtered.fastq"), "a")
                write(iostream, out_string2)
                close(iostream)
            end 
        end
        if verbose == true
            if (round((i/length(seqs))*100)) > 5 && (round((i/length(seqs))*100)) < 10 && alreadyprinted == false
                alreadyprinted = true
                println("Filtering 5% Complete")
            end 
            if (round((i/length(seqs))*100)) > 10 && (round((i/length(seqs))*100)) < 15 && alreadyprinted1 == false
                alreadyprinted1 = true
                println("Filtering 10% Complete")
            end 
            if (round((i/length(seqs))*100)) > 15 && (round((i/length(seqs))*100)) < 20 && alreadyprinted2 == false
                alreadyprinted2 = true
                println("Filtering 15% Complete")
            end 
            if (round((i/length(seqs))*100)) > 20 && (round((i/length(seqs))*100)) < 25 && alreadyprinted3 == false
                alreadyprinted3 = true
                println("Filtering 20% Complete")
            end 
            if (round((i/length(seqs))*100)) > 25 && (round((i/length(seqs))*100)) < 30 && alreadyprinted4 == false
                alreadyprinted4 = true
                println("Filtering 25% Complete")
            end 
            if (round((i/length(seqs))*100)) > 30 && (round((i/length(seqs))*100)) < 35 && alreadyprinted5 == false
                alreadyprinted5 = true
                println("Filtering 30% Complete")
            end 
            if (round((i/length(seqs))*100)) > 35 && (round((i/length(seqs))*100)) < 40 && alreadyprinted6 == false
                alreadyprinted6 = true
                println("Filtering 35% Complete")
            end 
            if (round((i/length(seqs))*100)) > 40 && (round((i/length(seqs))*100)) < 45 && alreadyprinted7 == false
                alreadyprinted7 = true
                println("Filtering 40% Complete")
            end 
            if (round((i/length(seqs))*100)) > 45 && (round((i/length(seqs))*100)) < 50 && alreadyprinted8 == false
                alreadyprinted8 = true
                println("Filtering 45% Complete")
            end 
            if (round((i/length(seqs))*100)) > 50 && (round((i/length(seqs))*100)) < 55 && alreadyprinted9 == false
                alreadyprinted9 = true
                println("Filtering 50% Complete")
            end 
            if (round((i/length(seqs))*100)) > 55 && (round((i/length(seqs))*100)) < 60 && alreadyprinted10 == false
                alreadyprinted10 = true
                println("Filtering 55% Complete")
            end 
            if (round((i/length(seqs))*100)) > 60 && (round((i/length(seqs))*100)) < 65 && alreadyprinted11 == false
                alreadyprinted11 = true
                println("Filtering 60% Complete")
            end 
            if (round((i/length(seqs))*100)) > 65 && (round((i/length(seqs))*100)) < 70 && alreadyprinted12 == false
                alreadyprinted12 = true
                println("Filtering 65% Complete")
            end 
            if (round((i/length(seqs))*100)) > 70 && (round((i/length(seqs))*100)) < 75 && alreadyprinted13 == false
                alreadyprinted13 = true
                println("Filtering 70% Complete")
            end 
            if (round((i/length(seqs))*100)) > 75 && (round((i/length(seqs))*100)) < 80 && alreadyprinted14 == false
                alreadyprinted14 = true
                println("Filtering 75% Complete")
            end 
            if (round((i/length(seqs))*100)) > 80 && (round((i/length(seqs))*100)) < 85 && alreadyprinted15 == false
                alreadyprinted15 = true
                println("Filtering 80% Complete")
            end 
            if (round((i/length(seqs))*100)) > 85 && (round((i/length(seqs))*100)) < 90 && alreadyprinted16 == false
                alreadyprinted16 = true
                println("Filtering 85% Complete")
            end 
            if (round((i/length(seqs))*100)) > 90 && (round((i/length(seqs))*100)) < 95 && alreadyprinted17 == false
                alreadyprinted17 = true
                println("Filtering 90% Complete")
            end 
            if (round((i/length(seqs))*100)) > 95 && (round((i/length(seqs))*100)) < 100 && alreadyprinted18 == false
                alreadyprinted18 = true
                println("Filtering 95% Complete")
            end 
        end 
    end 

    if verbose == true
        println("Cleaning things up ")
    end 
    
    # Zip the outputted files 
    outfile1 = string(out, "/",split(R1.filename, ".")[1], "_filtered.fastq")
    run(`gzip $outfile1`); #this works

    if unzipped == true
        if verbose == true
            println(string("Zipping ", read1))
        end 
        run(`gzip $read1`); #this works
    end 

    println("")
    println("")
    printstyled(string(Passing, " (", (round(Passing/length(seqs)*100, digits = 2)), "%) reads passed quality filtering\n"), color = :green)
    println("")
    println("")
    printstyled(string(Failing, " (", (round(Failing/length(seqs)*100, digits = 2)), "%) reads failed quality filtering\n"), color = :red)
end

"""
    Nucleotide_Essentials.FilterQuality_pe

Filters an input .fastq file based upon the encoded Phred+33 or Phred+64 quality scores. The encoding of the reads is automatically deteremined by looking for unique encoding in Phred+33 and Phred+64. Phred+64 encoding is identified by searching for `^`, `a`, `]`, and `f`.

Reads are filtered based upon the number of [expected errors](https://academic.oup.com/bioinformatics/article/31/21/3476/194979?login=true) (``\\mathrm{E}``) based on the error rate based on quality score and the sum of error probabilities, following the equation: 

``\\mathrm{E} = \\sum{_ip_i} = \\sum{_i}10^{\\frac{-Q_i}{10}}``

Stringent filtering (```maxEE = 1```) is used by default but can be adjusted by the user. 

Output Files: 

* If both the forward and reverse reads pass the filtering parameters: 
    * Forward reads are output to a file ending in `R1_Paired_filtered.fastq` in the user-determined directory, as indicated by ```out```
    * Reverse reads are output to a file ending in `R2_Paired_filtered.fastq` in the user-determined directory, as indicated by ```out```
* If only the forward read passes the filtering parameters: 
    * Forward reads are output to a file ending in `R1_Unpaired_filtered.fastq` in the user-determined directory, as indicated by ```out```
    * Reverse reads are not written to a file 
* If only the reverse read passes the filtering parameters: 
    * Reverse reads are output to a file ending in `R2_Unpaired_filtered.fastq` in the user-determined directory, as indicated by ```out```
    * Forward reads are not written to a file 
    

Supported keyword arguments include: 
    
* ```read1::String```: Path to the forward reads to undergo quality filtering 
* ```read2::String```: Path to the reverse reads to undergo quality filtering 
* ```out::String```: Path to the directory where reads that pass the quality filtering should be written 
* ```maxEE::Int64``` (optional): The max number of expected errors a read can include as the filtering parameter (default: ```maxEE = 1```)
* ```verbose::Bool``` (optional): Whether or not to show some intermediary feedback on the progress of the function (default = false)

# Example: 
```julia
FilterQuality_pe("forward_R1.fasta", "reverse_R2.fasta", "/outdirectory")

# changing the filtering parameters 
FilterQuality_pe("forward_R1.fasta", "reverse_R2.fasta", "/outdirectory", 2, true)
```
"""
function FilterQuality_pe(read1::String, read2::String, out::String, maxEE::Int64 = 1, verbose::Bool = false)
    
    if endswith(read1, ".gz") === true
        if verbose == true
            println(string("Unzipping ", read1))
        end 
        run(`gunzip $read1`); #this works
        unzipped1 = true
    else
        unzipped1 = false
    end
    if endswith(read2, ".gz") === true
        if verbose == true
            println(string("Unzipping ", read2))
        end 
        run(`gunzip $read2`); #this works
        unzipped2 = true
    else
        unzipped2 = false
    end


    # Forward reads 
    if endswith(read1, ".gz") === true
        R1_fixed = replace(read1, ".fastq.gz"=>".fastq")
        R1 = readFastq(R1_fixed)
    else 
        R1 = readFastq(read1)
    end 
    quals = R1.quality
    seqs = R1.sequence
    IDs = R1.ID

    # Reverse reads 
    if endswith(read2, ".gz") === true
        R2_fixed = replace(read2, ".fastq.gz"=>".fastq")
        R2 = readFastq(R2_fixed)
    else 
        R2 = readFastq(read2)
    end 
    quals2 = R2.quality
    seqs2 = R2.sequence
    IDs2 = R2.ID

    # Empty objects for tracking passing/failing reads
    Passing = 0
    UnpairedR1Succes = 0
    UnpairedR2Succes = 0
    bothfailed = 0

    alreadyprinted = false
    alreadyprinted1 = false
    alreadyprinted2 = false
    alreadyprinted3 = false
    alreadyprinted4 = false
    alreadyprinted5 = false
    alreadyprinted6 = false
    alreadyprinted7 = false
    alreadyprinted8 = false
    alreadyprinted9 = false
    alreadyprinted10 = false
    alreadyprinted11= false
    alreadyprinted12 = false
    alreadyprinted13= false
    alreadyprinted14 = false
    alreadyprinted15 = false
    alreadyprinted16 = false
    alreadyprinted17 = false
    alreadyprinted18 = false
 
    if occursin("^", join(quals)) == true && occursin("a", join(quals)) == true && occursin("]", join(quals)) == true && occursin("f", join(quals)) == true
    encoding = 64;
    println("Detected Phred+64 quality encoding")
    end 
    if occursin("^", join(quals)) == false && occursin("a", join(quals)) == false && occursin("]", join(quals)) == false && occursin("f", join(quals)) == false
        println("Detected Phred+33 quality encoding")
        encoding = 33;
    end 

    if encoding == 33
        scores_33 = Dict(
            ('!') => "1.0000", 
            ('\"') => "0.79433",
            ('#') => "0.63096",
            ('\$') => "0.50119",
            ('%') => "0.39811",
            ('&') => "0.31623",
            ('\'') => "0.25119",
            ('(') => "0.19953",
            (')') => "0.15849",
            ('*') => "0.12589",
            ('+') => "0.10000",
            (',') => "0.07943",
            ('-') => "0.06310",
            ('.') => "0.05012",
            ('/') => "0.03981",
            ('0') => "0.03162",
            ('1') => "0.02512",
            ('2') => "0.01995",
            ('3') => "0.01585",
            ('4') => "0.01259",
            ('5') => "0.01000",
            ('6') => "0.00794",
            ('7') => "0.00631",
            ('8') => "0.00501",
            ('9') => "0.00398",
            (':') => "0.00316",
            (';') => "0.00251",
            ('<') => "0.00200",
            ('=') => "0.00158",
            ('>') => "0.00126",
            ('?') => "0.00100",
            ('@') => "0.00079",
            ('A') => "0.00063",
            ('B') => "0.00050",
            ('C') => "0.00040",
            ('D') => "0.00032",
            ('E') => "0.00025",
            ('F') => "0.00020",
            ('G') => "0.00016",
            ('H') => "0.00013",
            ('I') => "0.00010", 
            ('J') => "0.00008",
            ('K') => "0.00006"
        )
    end 
    if encoding == 64
        scores_64 = Dict(
            ('@') => "1.0000", 
            ('A') => "0.79433",
            ('B') => "0.63096",
            ('C') => "0.50119",
            ('D') => "0.39811",
            ('E') => "0.31623",
            ('F') => "0.25119",
            ('G') => "0.19953",
            ('H') => "0.15849",
            ('I') => "0.12589",
            ('J') => "0.10000",
            ('K') => "0.07943",
            ('L') => "0.06310",
            ('M') => "0.05012",
            ('N') => "0.03981",
            ('O') => "0.03162",
            ('P') => "0.02512",
            ('Q') => "0.01995",
            ('R') => "0.01585",
            ('S') => "0.01259",
            ('T') => "0.01000",
            ('U') => "0.00794",
            ('V') => "0.00631",
            ('W') => "0.00501",
            ('X') => "0.00398",
            ('Y') => "0.00316",
            ('Z') => "0.00251",
            ('[') => "0.00200",
            ('\\') => "0.00158",
            (']') => "0.00126",
            ('^') => "0.00100",
            ('_') => "0.00079",
            ('`') => "0.00063",
            ('a') => "0.00050",
            ('b') => "0.00040",
            ('c') => "0.00032",
            ('d') => "0.00025",
            ('e') => "0.00020",
            ('f') => "0.00016",
            ('g') => "0.00013",
            ('h') => "0.00010", 
            ('i') => "0.00008",
            ('j') => "0.00006"
        )
    end 

    function convert_phred(c)
        if encoding == 33
            scores_33[c]
        else 
            scores_64[c]
        end 
    end

    if verbose == true
        println("Prepping reads for quality filtering")
    end 

    for i in 1:length(quals)
        ee1 = sum(parse.(Float64, (map(convert_phred, collect(quals[i])))))
        ee2 = sum(parse.(Float64, (map(convert_phred, collect(quals2[i])))))
        if ee1 < maxEE && ee2 < maxEE
            Passing = Passing + 1
            if isfile(string(out, "/",split(R1.filename, ".")[1], "_Paired_filtered.fastq")) == false
                out_string1 = string("@", IDs[i], "\n", seqs[i], "\n+\n", quals[i])
                write(string(out, "/",split(R1.filename, ".")[1], "_Paired_filtered.fastq"), out_string1)
            else 
                out_string2 = string("\n@", IDs[i], "\n", seqs[i], "\n+\n", quals[i])
                iostream = open(string(out, "/",split(R1.filename, ".")[1], "_Paired_filtered.fastq"), "a")
                write(iostream, out_string2)
                close(iostream)
            end
            if isfile(string(out, "/",split(R2.filename, ".")[1], "_Paired_filtered.fastq")) == false
                out_string3 = string("@", IDs2[i], "\n", seqs2[i], "\n+\n", quals2[i])
                write(string(out, "/",split(R2.filename, ".")[1], "_Paired_filtered.fastq"), out_string3)
            else 
                out_string4 = string("\n@", IDs2[i], "\n", seqs2[i], "\n+\n", quals2[i])
                iostream = open(string(out, "/",split(R2.filename, ".")[1], "_Paired_filtered.fastq"), "a")
                write(iostream, out_string4)
                close(iostream)
            end
        elseif ee1 < maxEE && ee2 > maxEE  
            UnpairedR1Succes = UnpairedR1Succes + 1
            if isfile(string(out, "/",split(R1.filename, ".")[1], "_Unpaired_filtered.fastq")) == false
                out_string5 = string("@", IDs[i], "\n", seqs[i], "\n+\n", quals[i])
                write(string(out, "/",split(R1.filename, ".")[1], "_Unpaired_filtered.fastq"), out_string5)
            else 
                out_string6 = string("\n@", IDs[i], "\n", seqs[i], "\n+\n", quals[i])
                iostream = open(string(out, "/",split(R1.filename, ".")[1], "_Unpaired_filtered.fastq"), "a")
                write(iostream, out_string6)
                close(iostream)
            end
        elseif ee1 >= maxEE && ee2 < maxEE  
            UnpairedR2Succes = UnpairedR2Succes + 1
            if isfile(string(out, "/",split(R2.filename, ".")[1], "_Unpaired_filtered.fastq")) == false
                out_string7 = string("@", IDs2[i], "\n", seqs2[i], "\n+\n", quals2[i])
                write(string(out, "/",split(R2.filename, ".")[1], "_Unpaired_filtered.fastq"), out_string7)
            else 
                out_string8 = string("\n@", IDs2[i], "\n", seqs2[i], "\n+\n", quals2[i])
                iostream = open(string(out, "/",split(R2.filename, ".")[1], "_Unpaired_filtered.fastq"), "a")
                write(iostream, out_string8)
                close(iostream)
            end 
        elseif ee1 >= maxEE && ee2 >= maxEE  
            bothfailed = bothfailed + 1
        end
        if verbose == true
            if (round((i/length(seqs))*100)) > 5 && (round((i/length(seqs))*100)) < 10 && alreadyprinted == false
                alreadyprinted = true
                println("Filtering 5% Complete")
            end 
            if (round((i/length(seqs))*100)) > 10 && (round((i/length(seqs))*100)) < 15 && alreadyprinted1 == false
                alreadyprinted1 = true
                println("Filtering 10% Complete")
            end 
            if (round((i/length(seqs))*100)) > 15 && (round((i/length(seqs))*100)) < 20 && alreadyprinted2 == false
                alreadyprinted2 = true
                println("Filtering 15% Complete")
            end 
            if (round((i/length(seqs))*100)) > 20 && (round((i/length(seqs))*100)) < 25 && alreadyprinted3 == false
                alreadyprinted3 = true
                println("Filtering 20% Complete")
            end 
            if (round((i/length(seqs))*100)) > 25 && (round((i/length(seqs))*100)) < 30 && alreadyprinted4 == false
                alreadyprinted4 = true
                println("Filtering 25% Complete")
            end 
            if (round((i/length(seqs))*100)) > 30 && (round((i/length(seqs))*100)) < 35 && alreadyprinted5 == false
                alreadyprinted5 = true
                println("Filtering 30% Complete")
            end 
            if (round((i/length(seqs))*100)) > 35 && (round((i/length(seqs))*100)) < 40 && alreadyprinted6 == false
                alreadyprinted6 = true
                println("Filtering 35% Complete")
            end 
            if (round((i/length(seqs))*100)) > 40 && (round((i/length(seqs))*100)) < 45 && alreadyprinted7 == false
                alreadyprinted7 = true
                println("Filtering 40% Complete")
            end 
            if (round((i/length(seqs))*100)) > 45 && (round((i/length(seqs))*100)) < 50 && alreadyprinted8 == false
                alreadyprinted8 = true
                println("Filtering 45% Complete")
            end 
            if (round((i/length(seqs))*100)) > 50 && (round((i/length(seqs))*100)) < 55 && alreadyprinted9 == false
                alreadyprinted9 = true
                println("Filtering 50% Complete")
            end 
            if (round((i/length(seqs))*100)) > 55 && (round((i/length(seqs))*100)) < 60 && alreadyprinted10 == false
                alreadyprinted10 = true
                println("Filtering 55% Complete")
            end 
            if (round((i/length(seqs))*100)) > 60 && (round((i/length(seqs))*100)) < 65 && alreadyprinted11 == false
                alreadyprinted11 = true
                println("Filtering 60% Complete")
            end 
            if (round((i/length(seqs))*100)) > 65 && (round((i/length(seqs))*100)) < 70 && alreadyprinted12 == false
                alreadyprinted12 = true
                println("Filtering 65% Complete")
            end 
            if (round((i/length(seqs))*100)) > 70 && (round((i/length(seqs))*100)) < 75 && alreadyprinted13 == false
                alreadyprinted13 = true
                println("Filtering 70% Complete")
            end 
            if (round((i/length(seqs))*100)) > 75 && (round((i/length(seqs))*100)) < 80 && alreadyprinted14 == false
                alreadyprinted14 = true
                println("Filtering 75% Complete")
            end 
            if (round((i/length(seqs))*100)) > 80 && (round((i/length(seqs))*100)) < 85 && alreadyprinted15 == false
                alreadyprinted15 = true
                println("Filtering 80% Complete")
            end 
            if (round((i/length(seqs))*100)) > 85 && (round((i/length(seqs))*100)) < 90 && alreadyprinted16 == false
                alreadyprinted16 = true
                println("Filtering 85% Complete")
            end 
            if (round((i/length(seqs))*100)) > 90 && (round((i/length(seqs))*100)) < 95 && alreadyprinted17 == false
                alreadyprinted17 = true
                println("Filtering 90% Complete")
            end 
            if (round((i/length(seqs))*100)) > 95 && (round((i/length(seqs))*100)) < 100 && alreadyprinted18 == false
                alreadyprinted18 = true
                println("Filtering 95% Complete")
            end 
        end  
    end   

    # Zip the outputted files 
    if verbose == true
        println("Cleaning things up ")
    end 
    outfile1 = string(out, "/",split(R1.filename, ".")[1], "_Paired_filtered.fastq")
    outfile2 = string(out, "/",split(R2.filename, ".")[1], "_Paired_filtered.fastq")
    outfile3 = string(out, "/",split(R1.filename, ".")[1], "_Unpaired_filtered.fastq")
    outfile4 = string(out, "/",split(R2.filename, ".")[1], "_Unpaired_filtered.fastq")
    if isfile(outfile1) == true
        if verbose == true
            println(string("""|
                        |----Zipping """, outfile1))
        end 
        run(`gzip $outfile1`);
    end 
    if isfile(outfile2) == true
        if verbose == true
            println(string("""|
                        |----Zipping """, outfile2))
        end
        run(`gzip $outfile2`);
    end 
    if isfile(outfile3) == true
        if verbose == true
            println(string("""|
                        |----Zipping """, outfile3))
        end
        run(`gzip $outfile3`);
    end 
    if isfile(outfile4) == true
        if verbose == true
            println(string("""|
                        |----Zipping """, outfile4))
        end
        run(`gzip $outfile4`);
    end 

    if unzipped1 == true
        if verbose == true
            println(string("""|
                        |----Zipping """, R1_fixed))
        end
        run(`gzip $R1_fixed`); #this works
    end 
    if unzipped2 == true
        if verbose == true
            println(string("""|
                        |----Zipping """, R2_fixed))
        end
        run(`gzip $R2_fixed`); #this works
    end 

    println("")
    printstyled(string(Passing, " (", (round(Passing/length(seqs)*100, digits = 2)), "%) read pairs passed quality filtering\n"), color = :green)
    println("")
    printstyled(string(bothfailed, " (", (round(bothfailed/length(seqs)*100, digits = 2)), "%) read pairs failed quality filtering\n"), color = :red)
    println("")
    printstyled(string(UnpairedR1Succes, " (", (round(UnpairedR1Succes/length(seqs)*100, digits = 2)), "%) forward reads passed without a matching reverse read\n"), color = :yellow)
    println("")
    printstyled(string(UnpairedR2Succes, " (", (round(UnpairedR2Succes/length(seqs)*100, digits = 2)), "%) reverse reads passed without a matching forward read\n"), color = :yellow)
end