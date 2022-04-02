# filtering based on https://academic.oup.com/bioinformatics/article/31/21/3476/194979?login=true


"""
    Nucleotide_Essentials.FilterQuality_se

Filters an input .fastq file based upon the encoded Phred+33 or Phred+64 quality scores. The encoding of the reads is automatically deteremined by looking for unique encoding in Phred+33 and Phred+64. Phred+64 encoding is identified by searching for `^`, `a`, `]`, and `f`.

Reads are filtered based upon the number of [expected errors](https://academic.oup.com/bioinformatics/article/31/21/3476/194979?login=true) (``\\mathrm{E}``) based on the error rate based on quality score and the sum of error probabilities, following the equation: 

``\\mathrm{E} = \\sum{_ip_i} = \\sum{_i}10^{\\frac{-Q_i}{10}}``

Stringent filtering (```maxEE = 1```) is used by default but can be adjusted by the user. 

Reads that pass the filtering parameters are output to a file named `R1_FilteredReads.fastq` in the user-determined directory, as indicated by ```out```. 

Supported keyword arguments include: 
    
* `R1::String`: Path to the reads to undergo quality filtering 
* `out::String`: Path to the directory where reads that pass the quality filtering should be written 
* 'maxEE::Int64' (optional): The max number of expected errors a read can include as the filtering parameter (default: ```maxEE = 1```)
* 'verbose::Bool' (optional): Whether or not to show some intermediary feedback on the progress of the function (default = false)

# Example: 
```julia
FilterQuality_se("forward_R1.fasta", "/outdirectory")
```
"""
function FilterQuality_se(R1::String, out::String, maxEE::Int64 = 1, verbose::Bool = false)
    
    R1 = readFastq(R1)

    qualities = R1.quality
    seqs = R1.sequence
    IDs = R1.ID
    Passing = 0
    Failing = 0

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

    if occursin("^", join(qualities)) == true && occursin("a", join(qualities)) == true && occursin("]", join(qualities)) == true && occursin("f", join(qualities)) == true
        encoding = 64
        println("Detected Phred64 quality encoding")
    end 

    if occursin("^", join(qualities)) == false && occursin("a", join(qualities)) == false && occursin("]", join(qualities)) == false && occursin("f", join(qualities)) == false
        println("Detected Phred33 quality encoding")
        encoding = 33
    end 

    for i in 1:length(qualities)
        temp_array = (qualities[i])
        if encoding == 64
            temp_array_2 = replace(temp_array,   
                "@" => "1.0000?", 
                "A" => "0.79433?",
                "B" => "0.63096?",
                "C" => "0.50119?",
                "D" => "0.39811?",
                "E" => "0.31623?",
                "F" => "0.25119?",
                "G" => "0.19953?",
                "H" => "0.15849?",
                "I" => "0.12589?",
                "J" => "0.10000?",
                "K" => "0.07943?",
                "L" => "0.06310?",
                "M" => "0.05012?",
                "N" => "0.03981?",
                "O" => "0.03162?",
                "P" => "0.02512?",
                "Q" => "0.01995?",
                "R" => "0.01585?",
                "S" => "0.01259?",
                "T" => "0.01000?",
                "U" => "0.00794?",
                "V" => "0.00631?",
                "W" => "0.00501?",
                "X" => "0.00398?",
                "Y" => "0.00316?",
                "Z" => "0.00251?",
                "[" => "0.00200?",
                "\\" => "0.00158?",
                "]" => "0.00126?",
                "^" => "0.00100?",
                "_" => "0.00079?",
                "`" => "0.00063?",
                "a" => "0.00050?",
                "b" => "0.00040?",
                "c" => "0.00032?",
                "d" => "0.00025?",
                "e" => "0.00020?",
                "f" => "0.00016?",
                "g" => "0.00013?",
                "h" => "0.00010?", 
                "i" => "0.00008?",
                "j" => "0.00006?"
            )
        end
        if encoding == 33 
            temp_array_2 = replace(temp_array,   
                "!" => "1.0000?", 
                "\"" => "0.79433?",
                "#" => "0.63096?",
                "\$" => "0.50119?",
                "%" => "0.39811?",
                "&" => "0.31623?",
                "'" => "0.25119?",
                "(" => "0.19953?",
                ")" => "0.15849?",
                "*" => "0.12589?",
                "+" => "0.10000?",
                "," => "0.07943?",
                "-" => "0.06310?",
                "." => "0.05012?",
                "/" => "0.03981?",
                "0" => "0.03162?",
                "1" => "0.02512?",
                "2" => "0.01995?",
                "3" => "0.01585?",
                "4" => "0.01259?",
                "5" => "0.01000?",
                "6" => "0.00794?",
                "7" => "0.00631?",
                "8" => "0.00501?",
                "9" => "0.00398?",
                ":" => "0.00316?",
                ";" => "0.00251?",
                "<" => "0.00200?",
                "=" => "0.00158?",
                ">" => "0.00126?",
                "?" => "0.00100?",
                "@" => "0.00079?",
                "A" => "0.00063?",
                "B" => "0.00050?",
                "C" => "0.00040?",
                "D" => "0.00032?",
                "E" => "0.00025?",
                "F" => "0.00020?",
                "G" => "0.00016?",
                "H" => "0.00013?",
                "I" => "0.00010?", 
                "J" => "0.00008?",
                "K" => "0.00006?"
            )
        end 
        out_vector = split(temp_array_2, "?", limit = (length(temp_array_2)))
        out_vector = deleteat!(out_vector, (length(out_vector))) 
        ee1 = sum(Array(parse.(Float64, out_vector)))

        if ee1 >= maxEE
            Failing = Failing + 1
        else 
            Passing = Passing + 1
            if isfile(string(out, "/R1_FilteredReads.fastq")) == false
                out_string = string("@", IDs[i], "\n", seqs[i], "\n+\n", qualities[i])
                write(string(out, "/R1_FilteredReads.fastq"), out_string)
            else 
                out_string2 = string("\n@", IDs[i], "\n", seqs[i], "\n+\n", qualities[i])
                iostream = open(string(out, "/R1_FilteredReads.fastq"), "a")
                write(iostream, out_string2)
                close(iostream)
            end 
        end
        if verbose == true
            if (round((i/length(qualities))*100)) > 5 && (round((i/length(qualities))*100)) < 10 && alreadyprinted == false
                alreadyprinted = true
                println("Filtering 5% Complete")
            end 
            if (round((i/length(qualities))*100)) > 10 && (round((i/length(qualities))*100)) < 15 && alreadyprinted1 == false
                alreadyprinted1 = true
                println("Filtering 10% Complete")
            end 
            if (round((i/length(qualities))*100)) > 15 && (round((i/length(qualities))*100)) < 20 && alreadyprinted2 == false
                alreadyprinted2 = true
                println("Filtering 15% Complete")
            end 
            if (round((i/length(qualities))*100)) > 20 && (round((i/length(qualities))*100)) < 25 && alreadyprinted3 == false
                alreadyprinted3 = true
                println("Filtering 20% Complete")
            end 
            if (round((i/length(qualities))*100)) > 25 && (round((i/length(qualities))*100)) < 30 && alreadyprinted4 == false
                alreadyprinted4 = true
                println("Filtering 25% Complete")
            end 
            if (round((i/length(qualities))*100)) > 30 && (round((i/length(qualities))*100)) < 35 && alreadyprinted5 == false
                alreadyprinted5 = true
                println("Filtering 30% Complete")
            end 
            if (round((i/length(qualities))*100)) > 35 && (round((i/length(qualities))*100)) < 40 && alreadyprinted6 == false
                alreadyprinted6 = true
                println("Filtering 35% Complete")
            end 
            if (round((i/length(qualities))*100)) > 40 && (round((i/length(qualities))*100)) < 45 && alreadyprinted7 == false
                alreadyprinted7 = true
                println("Filtering 40% Complete")
            end 
            if (round((i/length(qualities))*100)) > 45 && (round((i/length(qualities))*100)) < 50 && alreadyprinted8 == false
                alreadyprinted8 = true
                println("Filtering 45% Complete")
            end 
            if (round((i/length(qualities))*100)) > 50 && (round((i/length(qualities))*100)) < 55 && alreadyprinted9 == false
                alreadyprinted9 = true
                println("Filtering 50% Complete")
            end 
            if (round((i/length(qualities))*100)) > 55 && (round((i/length(qualities))*100)) < 60 && alreadyprinted10 == false
                alreadyprinted10 = true
                println("Filtering 55% Complete")
            end 
            if (round((i/length(qualities))*100)) > 60 && (round((i/length(qualities))*100)) < 65 && alreadyprinted11 == false
                alreadyprinted11 = true
                println("Filtering 60% Complete")
            end 
            if (round((i/length(qualities))*100)) > 65 && (round((i/length(qualities))*100)) < 70 && alreadyprinted12 == false
                alreadyprinted12 = true
                println("Filtering 65% Complete")
            end 
            if (round((i/length(qualities))*100)) > 70 && (round((i/length(qualities))*100)) < 75 && alreadyprinted13 == false
                alreadyprinted13 = true
                println("Filtering 70% Complete")
            end 
            if (round((i/length(qualities))*100)) > 75 && (round((i/length(qualities))*100)) < 80 && alreadyprinted14 == false
                alreadyprinted14 = true
                println("Filtering 75% Complete")
            end 
            if (round((i/length(qualities))*100)) > 80 && (round((i/length(qualities))*100)) < 85 && alreadyprinted15 == false
                alreadyprinted15 = true
                println("Filtering 80% Complete")
            end 
            if (round((i/length(qualities))*100)) > 85 && (round((i/length(qualities))*100)) < 90 && alreadyprinted16 == false
                alreadyprinted16 = true
                println("Filtering 85% Complete")
            end 
            if (round((i/length(qualities))*100)) > 90 && (round((i/length(qualities))*100)) < 95 && alreadyprinted17 == false
                alreadyprinted17 = true
                println("Filtering 90% Complete")
            end 
            if (round((i/length(qualities))*100)) > 95 && (round((i/length(qualities))*100)) < 100 && alreadyprinted18 == false
                alreadyprinted18 = true
                println("Filtering 95% Complete")
            end 
        end  
    end     
    
    # Zip the outputted files 
    outfile1 = string(out, "/R1_FilteredReads.fastq")
    run(`gzip $outfile1`); #this works

    println("")
    println("")
    printstyled(string(Passing, " (", (round(Passing/length(qualities)*100, digits = 2)), "%) reads passed quality filtering\n"), color = :green)
    println("")
    println("")
    printstyled(string(Failing, " (", (round(Failing/length(qualities)*100, digits = 2)), "%) reads failed quality filtering\n"), color = :red)
end

"""
    Nucleotide_Essentials.FilterQuality_pe

Filters an input .fastq file based upon the encoded Phred+33 or Phred+64 quality scores. The encoding of the reads is automatically deteremined by looking for unique encoding in Phred+33 and Phred+64. Phred+64 encoding is identified by searching for `^`, `a`, `]`, and `f`.

Reads are filtered based upon the number of [expected errors](https://academic.oup.com/bioinformatics/article/31/21/3476/194979?login=true) (``\\mathrm{E}``) based on the error rate based on quality score and the sum of error probabilities, following the equation: 

``\\mathrm{E} = \\sum{_ip_i} = \\sum{_i}10^{\\frac{-Q_i}{10}}``

Stringent filtering (```maxEE = 1```) is used by default but can be adjusted by the user. 

Output Files: 

* If both the forward and reverse reads pass the filtering parameters: 
    * Forward reads are output to a file named `R1_Paired.fastq` in the user-determined directory, as indicated by ```out```
    * Reverse reads are output to a file named `R2_Paired.fastq` in the user-determined directory, as indicated by ```out```
* If only the forward read passes the filtering parameters: 
    * Forward reads are output to a file named `R1_Unpaired.fastq` in the user-determined directory, as indicated by ```out```
    * Reverse reads are not written to a file 
* If only the reverse read passes the filtering parameters: 
    * Reverse reads are output to a file named `R2_Unpaired.fastq` in the user-determined directory, as indicated by ```out```
    * Forward reads are not written to a file 
    

Supported keyword arguments include: 
    
* `R1::String`: Path to the forward reads to undergo quality filtering 
* `R2::String`: Path to the reverse reads to undergo quality filtering 
* `out::String`: Path to the directory where reads that pass the quality filtering should be written 
* 'maxEE::Int64' (optional): The max number of expected errors a read can include as the filtering parameter (default: ```maxEE = 1```)
* 'verbose::Bool' (optional): Whether or not to show some intermediary feedback on the progress of the function (default = false)

# Example: 
```julia
FilterQuality_pe("forward_R1.fasta", "reverse_R2.fasta", "/outdirectory")
```
"""
function FilterQuality_pe(R1::String, R2::String, out::String, maxEE::Int64 = 1, verbose::Bool = false)
    
    R1 = readFastq(R1)
    R2 = readFastq(R2)

    qualities = R1.quality
    seqs = R1.sequence
    IDs = R1.ID
    qualities2 = R2.quality
    seqs2 = R2.sequence
    IDs2 = R2.ID
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

    if occursin("^", join(qualities)) == true && occursin("a", join(qualities)) == true && occursin("]", join(qualities)) == true && occursin("f", join(qualities)) == true
        encoding = 64
        println("Detected Phred64 quality encoding")
    end 

    if occursin("^", join(qualities)) == false && occursin("a", join(qualities)) == false && occursin("]", join(qualities)) == false && occursin("f", join(qualities)) == false
        println("Detected Phred33 quality encoding")
        encoding = 33
    end 
 
    for i in 1:length(qualities)
        temp_array = (qualities[i])
        if encoding == 33
            temp_array_2 = replace(temp_array,   
                "!" => "1.0000?", 
                "\"" => "0.79433?",
                "#" => "0.63096?",
                "\$" => "0.50119?",
                "%" => "0.39811?",
                "&" => "0.31623?",
                "'" => "0.25119?",
                "(" => "0.19953?",
                ")" => "0.15849?",
                "*" => "0.12589?",
                "+" => "0.10000?",
                "," => "0.07943?",
                "-" => "0.06310?",
                "." => "0.05012?",
                "/" => "0.03981?",
                "0" => "0.03162?",
                "1" => "0.02512?",
                "2" => "0.01995?",
                "3" => "0.01585?",
                "4" => "0.01259?",
                "5" => "0.01000?",
                "6" => "0.00794?",
                "7" => "0.00631?",
                "8" => "0.00501?",
                "9" => "0.00398?",
                ":" => "0.00316?",
                ";" => "0.00251?",
                "<" => "0.00200?",
                "=" => "0.00158?",
                ">" => "0.00126?",
                "?" => "0.00100?",
                "@" => "0.00079?",
                "A" => "0.00063?",
                "B" => "0.00050?",
                "C" => "0.00040?",
                "D" => "0.00032?",
                "E" => "0.00025?",
                "F" => "0.00020?",
                "G" => "0.00016?",
                "H" => "0.00013?",
                "I" => "0.00010?", 
                "J" => "0.00008?",
                "K" => "0.00006?"
            )
        end 
        if encoding == 64
            temp_array_2 = replace(temp_array,   
                "@" => "1.0000?", 
                "A" => "0.79433?",
                "B" => "0.63096?",
                "C" => "0.50119?",
                "D" => "0.39811?",
                "E" => "0.31623?",
                "F" => "0.25119?",
                "G" => "0.19953?",
                "H" => "0.15849?",
                "I" => "0.12589?",
                "J" => "0.10000?",
                "K" => "0.07943?",
                "L" => "0.06310?",
                "M" => "0.05012?",
                "N" => "0.03981?",
                "O" => "0.03162?",
                "P" => "0.02512?",
                "Q" => "0.01995?",
                "R" => "0.01585?",
                "S" => "0.01259?",
                "T" => "0.01000?",
                "U" => "0.00794?",
                "V" => "0.00631?",
                "W" => "0.00501?",
                "X" => "0.00398?",
                "Y" => "0.00316?",
                "Z" => "0.00251?",
                "[" => "0.00200?",
                "\\" => "0.00158?",
                "]" => "0.00126?",
                "^" => "0.00100?",
                "_" => "0.00079?",
                "`" => "0.00063?",
                "a" => "0.00050?",
                "b" => "0.00040?",
                "c" => "0.00032?",
                "d" => "0.00025?",
                "e" => "0.00020?",
                "f" => "0.00016?",
                "g" => "0.00013?",
                "h" => "0.00010?", 
                "i" => "0.00008?",
                "j" => "0.00006?"
            )
        end
        out_vector = split(temp_array_2, "?", limit = (length(temp_array_2)))
        out_vector = deleteat!(out_vector, (length(out_vector))) 
        ee1 = sum(Array(parse.(Float64, out_vector)))

        temp_array2 = (qualities2[i])
        if encoding == 33
            temp_array2_2 = replace(temp_array2,   
                "!" => "1.0000?", 
                "\"" => "0.79433?",
                "#" => "0.63096?",
                "\$" => "0.50119?",
                "%" => "0.39811?",
                "&" => "0.31623?",
                "'" => "0.25119?",
                "(" => "0.19953?",
                ")" => "0.15849?",
                "*" => "0.12589?",
                "+" => "0.10000?",
                "," => "0.07943?",
                "-" => "0.06310?",
                "." => "0.05012?",
                "/" => "0.03981?",
                "0" => "0.03162?",
                "1" => "0.02512?",
                "2" => "0.01995?",
                "3" => "0.01585?",
                "4" => "0.01259?",
                "5" => "0.01000?",
                "6" => "0.00794?",
                "7" => "0.00631?",
                "8" => "0.00501?",
                "9" => "0.00398?",
                ":" => "0.00316?",
                ";" => "0.00251?",
                "<" => "0.00200?",
                "=" => "0.00158?",
                ">" => "0.00126?",
                "?" => "0.00100?",
                "@" => "0.00079?",
                "A" => "0.00063?",
                "B" => "0.00050?",
                "C" => "0.00040?",
                "D" => "0.00032?",
                "E" => "0.00025?",
                "F" => "0.00020?",
                "G" => "0.00016?",
                "H" => "0.00013?",
                "I" => "0.00010?", 
                "J" => "0.00008?",
                "K" => "0.00006?"
            )
        end
        if encoding == 64
            temp_array2_2 = replace(temp_array2,   
                "@" => "1.0000?", 
                "A" => "0.79433?",
                "B" => "0.63096?",
                "C" => "0.50119?",
                "D" => "0.39811?",
                "E" => "0.31623?",
                "F" => "0.25119?",
                "G" => "0.19953?",
                "H" => "0.15849?",
                "I" => "0.12589?",
                "J" => "0.10000?",
                "K" => "0.07943?",
                "L" => "0.06310?",
                "M" => "0.05012?",
                "N" => "0.03981?",
                "O" => "0.03162?",
                "P" => "0.02512?",
                "Q" => "0.01995?",
                "R" => "0.01585?",
                "S" => "0.01259?",
                "T" => "0.01000?",
                "U" => "0.00794?",
                "V" => "0.00631?",
                "W" => "0.00501?",
                "X" => "0.00398?",
                "Y" => "0.00316?",
                "Z" => "0.00251?",
                "[" => "0.00200?",
                "\\" => "0.00158?",
                "]" => "0.00126?",
                "^" => "0.00100?",
                "_" => "0.00079?",
                "`" => "0.00063?",
                "a" => "0.00050?",
                "b" => "0.00040?",
                "c" => "0.00032?",
                "d" => "0.00025?",
                "e" => "0.00020?",
                "f" => "0.00016?",
                "g" => "0.00013?",
                "h" => "0.00010?", 
                "i" => "0.00008?",
                "j" => "0.00006?"
            )
        end 
        out_vector2 = split(temp_array2_2, "?", limit = (length(temp_array2_2)))
        out_vector2 = deleteat!(out_vector2, (length(out_vector2))) 
        ee2 = sum(Array(parse.(Float64, out_vector2)))

        if ee1 < maxEE && ee2 < maxEE
            Passing = Passing + 1
            if isfile(string(out, "/R1_Paired.fastq")) == false
                out_string1 = string("@", IDs[i], "\n", seqs[i], "\n+\n", qualities[i])
                write(string(out, "/R1_Paired.fastq"), out_string1)
            else 
                out_string2 = string("\n@", IDs[i], "\n", seqs[i], "\n+\n", qualities[i])
                iostream = open(string(out, "/R1_Paired.fastq"), "a")
                write(iostream, out_string2)
                close(iostream)
            end
            if isfile(string(out, "/R2_Paired.fastq")) == false
                out_string3 = string("@", IDs2[i], "\n", seqs2[i], "\n+\n", qualities2[i])
                write(string(out, "/R2_Paired.fastq"), out_string3)
            else 
                out_string4 = string("\n@", IDs2[i], "\n", seqs2[i], "\n+\n", qualities2[i])
                iostream = open(string(out, "/R2_Paired.fastq"), "a")
                write(iostream, out_string4)
                close(iostream)
            end
        elseif ee1 < maxEE && ee2 > maxEE  
            UnpairedR1Succes = UnpairedR1Succes + 1
            if isfile(string(out, "/R1_Unpaired.fastq")) == false
                out_string5 = string("@", IDs[i], "\n", seqs[i], "\n+\n", qualities[i])
                write(string(out, "/R1_Unpaired.fastq"), out_string5)
            else 
                out_string6 = string("\n@", IDs[i], "\n", seqs[i], "\n+\n", qualities[i])
                iostream = open(string(out, "/R1_Unpaired.fastq"), "a")
                write(iostream, out_string6)
                close(iostream)
            end
        elseif ee1 >= maxEE && ee2 < maxEE  
            UnpairedR2Succes = UnpairedR2Succes + 1
            if isfile(string(out, "/R2_Unpaired.fastq")) == false
                out_string7 = string("@", IDs2[i], "\n", seqs2[i], "\n+\n", qualities2[i])
                write(string(out, "/R2_Unpaired.fastq"), out_string7)
            else 
                out_string8 = string("\n@", IDs2[i], "\n", seqs2[i], "\n+\n", qualities2[i])
                iostream = open(string(out, "/R2_Unpaired.fastq"), "a")
                write(iostream, out_string8)
                close(iostream)
            end 
        elseif ee1 >= maxEE && ee2 >= maxEE  
            bothfailed = bothfailed + 1
        end
        if verbose == true
            if (round((i/length(qualities))*100)) > 5 && (round((i/length(qualities))*100)) < 10 && alreadyprinted == false
                alreadyprinted = true
                println("Filtering 5% Complete")
            end 
            if (round((i/length(qualities))*100)) > 10 && (round((i/length(qualities))*100)) < 15 && alreadyprinted1 == false
                alreadyprinted1 = true
                println("Filtering 10% Complete")
            end 
            if (round((i/length(qualities))*100)) > 15 && (round((i/length(qualities))*100)) < 20 && alreadyprinted2 == false
                alreadyprinted2 = true
                println("Filtering 15% Complete")
            end 
            if (round((i/length(qualities))*100)) > 20 && (round((i/length(qualities))*100)) < 25 && alreadyprinted3 == false
                alreadyprinted3 = true
                println("Filtering 20% Complete")
            end 
            if (round((i/length(qualities))*100)) > 25 && (round((i/length(qualities))*100)) < 30 && alreadyprinted4 == false
                alreadyprinted4 = true
                println("Filtering 25% Complete")
            end 
            if (round((i/length(qualities))*100)) > 30 && (round((i/length(qualities))*100)) < 35 && alreadyprinted5 == false
                alreadyprinted5 = true
                println("Filtering 30% Complete")
            end 
            if (round((i/length(qualities))*100)) > 35 && (round((i/length(qualities))*100)) < 40 && alreadyprinted6 == false
                alreadyprinted6 = true
                println("Filtering 35% Complete")
            end 
            if (round((i/length(qualities))*100)) > 40 && (round((i/length(qualities))*100)) < 45 && alreadyprinted7 == false
                alreadyprinted7 = true
                println("Filtering 40% Complete")
            end 
            if (round((i/length(qualities))*100)) > 45 && (round((i/length(qualities))*100)) < 50 && alreadyprinted8 == false
                alreadyprinted8 = true
                println("Filtering 45% Complete")
            end 
            if (round((i/length(qualities))*100)) > 50 && (round((i/length(qualities))*100)) < 55 && alreadyprinted9 == false
                alreadyprinted9 = true
                println("Filtering 50% Complete")
            end 
            if (round((i/length(qualities))*100)) > 55 && (round((i/length(qualities))*100)) < 60 && alreadyprinted10 == false
                alreadyprinted10 = true
                println("Filtering 55% Complete")
            end 
            if (round((i/length(qualities))*100)) > 60 && (round((i/length(qualities))*100)) < 65 && alreadyprinted11 == false
                alreadyprinted11 = true
                println("Filtering 60% Complete")
            end 
            if (round((i/length(qualities))*100)) > 65 && (round((i/length(qualities))*100)) < 70 && alreadyprinted12 == false
                alreadyprinted12 = true
                println("Filtering 65% Complete")
            end 
            if (round((i/length(qualities))*100)) > 70 && (round((i/length(qualities))*100)) < 75 && alreadyprinted13 == false
                alreadyprinted13 = true
                println("Filtering 70% Complete")
            end 
            if (round((i/length(qualities))*100)) > 75 && (round((i/length(qualities))*100)) < 80 && alreadyprinted14 == false
                alreadyprinted14 = true
                println("Filtering 75% Complete")
            end 
            if (round((i/length(qualities))*100)) > 80 && (round((i/length(qualities))*100)) < 85 && alreadyprinted15 == false
                alreadyprinted15 = true
                println("Filtering 80% Complete")
            end 
            if (round((i/length(qualities))*100)) > 85 && (round((i/length(qualities))*100)) < 90 && alreadyprinted16 == false
                alreadyprinted16 = true
                println("Filtering 85% Complete")
            end 
            if (round((i/length(qualities))*100)) > 90 && (round((i/length(qualities))*100)) < 95 && alreadyprinted17 == false
                alreadyprinted17 = true
                println("Filtering 90% Complete")
            end 
            if (round((i/length(qualities))*100)) > 95 && (round((i/length(qualities))*100)) < 100 && alreadyprinted18 == false
                alreadyprinted18 = true
                println("Filtering 95% Complete")
            end 
        end  
    end     

    # Zip the outputted files 
    outfile1 = string(out, "/R1_Paired.fastq")
    outfile2 = string(out, "/R2_Paired.fastq")
    outfile3 = string(out, "/R1_Unpaired.fastq")
    outfile4 = string(out, "/R2_Unpaired.fastq")
    run(`gzip $outfile1`);
    run(`gzip $outfile2`);
    run(`gzip $outfile3`);
    run(`gzip $outfile4`);

    println("")
    println("")
    printstyled(string(Passing, " (", (round(Passing/length(qualities)*100, digits = 2)), "%) read pairs passed quality filtering\n"), color = :green)
    println("")
    println("")
    printstyled(string(bothfailed, " (", (round(bothfailed/length(qualities)*100, digits = 2)), "%) read pairs failed quality filtering\n"), color = :red)
    println("")
    println("")
    printstyled(string(UnpairedR1Succes, " (", (round(UnpairedR1Succes/length(qualities)*100, digits = 2)), "%) forward reads passed without a matching reverse read\n"), color = :yellow)
    println("")
    println("")
    printstyled(string(UnpairedR2Succes, " (", (round(UnpairedR2Succes/length(qualities)*100, digits = 2)), "%) reverse reads passed without a matching forward read\n"), color = :yellow)
end

FilterQuality_pe("/Users/patrick/Desktop/demultiplex_ex/Sam78-125_S3_L001_R1_001.fastq", "/Users/patrick/Desktop/demultiplex_ex/Sam78-125_S3_L001_R2_001.fastq", "/Users/patrick/Desktop/demultiplex_ex")