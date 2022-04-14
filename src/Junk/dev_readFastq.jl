import Base.Threads.@threads
using DataFrames
"""
    Nucleotide_Essentials.FastqRecord

# Components

* ID: The unique sequence identifier associated with that entry
* sequence: The nucleotide sequence of that entry
* quality: The quality scores of that entry 
* filename: The original file name 
"""
struct FastqRecord
    ID::Vector{SubString{String}}
    sequence::Vector{SubString{String}}
    quality::Vector{SubString{String}}
    dfquality::DataFrame
    #filename::SubString{String}
end

"""
    Nucleotide_Essentials.readFastq

    readFastq(Path::String)
.fastq file => readFastq(Path) => FastqRecord(ID, sequence, quality, filename)

Supported keyword arguments include: 

* 'Path::String': The full or relative path to a .fastq file

# Example: 

```julia
# Supply the path to a .Fastq file that you would like to import
myfastq = readFastq("myfastq.fastq")
```
"""
function readFastq(Path::String)
    size = filesize(Path)
    #file = rsplit(Path, "/", limit = 2)[2];
    input = read(Path, String);
    input2 = replace(input, "\n@" => "Z");
    input3 = replace(input2, "@" => "Z", count = 1);
    input4 = replace(input3, "Z" => "\n");
    split_fastq = split(input4, "\n");
    og_quals = split_fastq[5:4:end];
    string_quals = og_quals
    max_length = maximum(length.(og_quals));
    if occursin("^", join(og_quals)) == true && occursin("a", join(og_quals)) == true && occursin("]", join(og_quals)) == true && occursin("f", join(og_quals)) == true
        encoding = 64;
        println("Detected Phred+64 quality encoding")
    end 
    if occursin("^", join(og_quals)) == false && occursin("a", join(og_quals)) == false && occursin("]", join(og_quals)) == false && occursin("f", join(og_quals)) == false
        println("Detected Phred+33 quality encoding")
        encoding = 33;
    end 

    if size > 50000000
        merged_df = DataFrame()
        println("Large file detected, beginning staged analysis")
        splitamount = Int64(round(length(og_quals)/50))
        og_quals_sub1 = og_quals[1:splitamount]
        og_quals_sub2 = og_quals[(splitamount+1):(2*splitamount)]
        og_quals_sub3 = og_quals[(2*splitamount+1):(3*splitamount)]
        og_quals_sub4 = og_quals[(3*splitamount+1):(4*splitamount)]
        og_quals_sub5 = og_quals[(4*splitamount+1):(5*splitamount)]
        og_quals_sub6 = og_quals[(5*splitamount+1):(6*splitamount)]
        og_quals_sub7 = og_quals[(6*splitamount+1):(7*splitamount)]
        og_quals_sub8 = og_quals[(7*splitamount+1):(8*splitamount)]
        og_quals_sub9 = og_quals[(8*splitamount+1):(9*splitamount)]
        og_quals_sub10 = og_quals[(9*splitamount+1):(10*splitamount)]
        og_quals_sub11 = og_quals[(10*splitamount+1):(11*splitamount)]
        og_quals_sub12 = og_quals[(11*splitamount+1):(12*splitamount)]
        og_quals_sub13 = og_quals[(12*splitamount+1):(13*splitamount)]
        og_quals_sub14 = og_quals[(13*splitamount+1):(14*splitamount)]
        og_quals_sub15 = og_quals[(14*splitamount+1):(15*splitamount)]
        og_quals_sub16 = og_quals[(15*splitamount+1):(16*splitamount)]
        og_quals_sub17 = og_quals[(16*splitamount+1):(17*splitamount)]
        og_quals_sub18 = og_quals[(17*splitamount+1):(18*splitamount)]
        og_quals_sub19 = og_quals[(18*splitamount+1):(19*splitamount)]
        og_quals_sub20 = og_quals[(19*splitamount+1):(20*splitamount)]
        og_quals_sub21 = og_quals[(20*splitamount+1):(21*splitamount)]
        og_quals_sub22 = og_quals[(21*splitamount+1):(22*splitamount)]
        og_quals_sub23 = og_quals[(22*splitamount+1):(23*splitamount)]
        og_quals_sub24 = og_quals[(23*splitamount+1):(24*splitamount)]
        og_quals_sub25 = og_quals[(24*splitamount+1):(25*splitamount)]
        og_quals_sub26 = og_quals[(25*splitamount+1):(26*splitamount)]
        og_quals_sub27 = og_quals[(26*splitamount+1):(27*splitamount)]
        og_quals_sub28 = og_quals[(27*splitamount+1):(28*splitamount)]
        og_quals_sub29 = og_quals[(28*splitamount+1):(29*splitamount)]
        og_quals_sub30 = og_quals[(29*splitamount+1):(30*splitamount)]
        og_quals_sub31 = og_quals[(30*splitamount+1):(31*splitamount)]
        og_quals_sub32 = og_quals[(31*splitamount+1):(32*splitamount)]
        og_quals_sub33 = og_quals[(32*splitamount+1):(33*splitamount)]
        og_quals_sub34 = og_quals[(33*splitamount+1):(34*splitamount)]
        og_quals_sub35 = og_quals[(34*splitamount+1):(35*splitamount)]
        og_quals_sub36 = og_quals[(35*splitamount+1):(36*splitamount)]
        og_quals_sub37 = og_quals[(36*splitamount+1):(37*splitamount)]
        og_quals_sub38 = og_quals[(37*splitamount+1):(38*splitamount)]
        og_quals_sub39 = og_quals[(38*splitamount+1):(39*splitamount)]
        og_quals_sub40 = og_quals[(39*splitamount+1):(40*splitamount)]
        og_quals_sub41 = og_quals[(40*splitamount+1):(41*splitamount)]
        og_quals_sub42 = og_quals[(41*splitamount+1):(42*splitamount)]
        og_quals_sub43 = og_quals[(42*splitamount+1):(43*splitamount)]
        og_quals_sub44 = og_quals[(43*splitamount+1):(44*splitamount)]
        og_quals_sub45 = og_quals[(44*splitamount+1):(45*splitamount)]
        og_quals_sub46 = og_quals[(45*splitamount+1):(46*splitamount)]
        og_quals_sub47 = og_quals[(46*splitamount+1):(47*splitamount)]
        og_quals_sub48 = og_quals[(47*splitamount+1):(48*splitamount)]
        og_quals_sub49 = og_quals[(48*splitamount+1):(49*splitamount)]
        og_quals_sub50 = og_quals[(49*splitamount+1):end]

        subgroups = [og_quals_sub1, 
        og_quals_sub2, 
        og_quals_sub3, 
        og_quals_sub4, 
        og_quals_sub5, 
        og_quals_sub6, 
        og_quals_sub7, 
        og_quals_sub8, 
        og_quals_sub9, 
        og_quals_sub10,
        og_quals_sub11,
        og_quals_sub12,
        og_quals_sub13,
        og_quals_sub14,
        og_quals_sub15,
        og_quals_sub16,
        og_quals_sub17,
        og_quals_sub18,
        og_quals_sub19,
        og_quals_sub20,
        og_quals_sub21,
        og_quals_sub22,
        og_quals_sub23,
        og_quals_sub24,
        og_quals_sub25,
        og_quals_sub26,
        og_quals_sub27,
        og_quals_sub28,
        og_quals_sub29,
        og_quals_sub30,
        og_quals_sub31,
        og_quals_sub32,
        og_quals_sub33,
        og_quals_sub34,
        og_quals_sub35,
        og_quals_sub36,
        og_quals_sub37,
        og_quals_sub38,
        og_quals_sub39,
        og_quals_sub40,
        og_quals_sub41,
        og_quals_sub42,
        og_quals_sub43,
        og_quals_sub44,
        og_quals_sub45,
        og_quals_sub46,
        og_quals_sub47,
        og_quals_sub48,
        og_quals_sub49,
        og_quals_sub50]
        
        for s in 1:length(subgroups)
            println(string("Beginning batched reading stage ", s, " of ", length(subgroups)))
            if encoding == 33
                og_quals = replace.(subgroups[s],   
                            "!" => "00?", 
                            "\"" => "01?",
                            "#" => "02?",
                            "\$" => "03?",
                            "%" => "04?",
                            "&" => "05?",
                            "'" => "06?",
                            "(" => "07?",
                            ")" => "08?",
                            "*" => "09?",
                            "+" => "10?",
                            "," => "11?",
                            "-" => "12?",
                            "." => "13?",
                            "/" => "14?",
                            "0" => "15?",
                            "1" => "16?",
                            "2" => "17?",
                            "3" => "18?",
                            "4" => "19?",
                            "5" => "20?",
                            "6" => "21?",
                            "7" => "22?",
                            "8" => "23?",
                            "9" => "24?",
                            ":" => "25?",
                            ";" => "26?",
                            "<" => "27?",
                            "=" => "28?",
                            ">" => "29?",
                            "?" => "30?",
                            "@" => "31?",
                            "A" => "32?",
                            "B" => "33?",
                            "C" => "34?",
                            "D" => "35?",
                            "E" => "36?",
                            "F" => "37?",
                            "G" => "38?",
                            "H" => "39?",
                            "I" => "40?",
                            "J" => "41?",
                            "K" => "42?")
            end 
            if encoding == 64
                og_quals = replace.(subgroups[s],   
                            "@" => "00?", 
                            "A" => "01?",
                            "B" => "02?",
                            "C" => "03?",
                            "D" => "04?",
                            "E" => "05?",
                            "F" => "06?",
                            "G" => "07?",
                            "H" => "08?",
                            "I" => "09?",
                            "J" => "10?",
                            "K" => "11?",
                            "L" => "12?",
                            "M" => "13?",
                            "N" => "14?",
                            "O" => "15?",
                            "P" => "16?",
                            "Q" => "17?",
                            "R" => "18?",
                            "S" => "19?",
                            "T" => "20?",
                            "U" => "21?",
                            "V" => "22?",
                            "W" => "23?",
                            "X" => "24?",
                            "Y" => "25?",
                            "Z" => "26?",
                            "[" => "27?",
                            "\\" => "28?",
                            "]" => "29?",
                            "^" => "30?",
                            "_" => "31?",
                            "`" => "32?",
                            "a" => "33?",
                            "b" => "34?",
                            "c" => "35?",
                            "d" => "36?",
                            "e" => "37?",
                            "f" => "38?",
                            "g" => "39?",
                            "h" => "40?", 
                            "i" => "41?",
                            "j" => "42?")            
            end 
            og_quals = Array(split.(og_quals, "?"))
            pop!.(og_quals)
            df = DataFrame([String[] for i in 1:max_length], :auto)
            for i in 1:length(og_quals)
                #println(string("starting loop #", i, " of ", length(og_quals)))
                temp = collect(og_quals[i])
                if length(temp) < max_length
                    difference = max_length - length(temp)
                    for r in 1:difference 
                        append!(temp, ["NaN"])
                    end 
                end 
                parse.(Float64, temp)
                push!(df, temp)
            end 
            df = passmissing(parse).(Float64, df);
            append!(merged_df,df)
        end 
        for col in eachcol(merged_df)
            replace!(col,NaN => Float64(100))
        end
        out = FastqRecord(split_fastq[2:4:end], split_fastq[3:4:end], string_quals, merged_df)#, file)
    else 
        if encoding == 33
            og_quals = replace.(og_quals,   
                        "!" => "00?", 
                        "\"" => "01?",
                        "#" => "02?",
                        "\$" => "03?",
                        "%" => "04?",
                        "&" => "05?",
                        "'" => "06?",
                        "(" => "07?",
                        ")" => "08?",
                        "*" => "09?",
                        "+" => "10?",
                        "," => "11?",
                        "-" => "12?",
                        "." => "13?",
                        "/" => "14?",
                        "0" => "15?",
                        "1" => "16?",
                        "2" => "17?",
                        "3" => "18?",
                        "4" => "19?",
                        "5" => "20?",
                        "6" => "21?",
                        "7" => "22?",
                        "8" => "23?",
                        "9" => "24?",
                        ":" => "25?",
                        ";" => "26?",
                        "<" => "27?",
                        "=" => "28?",
                        ">" => "29?",
                        "?" => "30?",
                        "@" => "31?",
                        "A" => "32?",
                        "B" => "33?",
                        "C" => "34?",
                        "D" => "35?",
                        "E" => "36?",
                        "F" => "37?",
                        "G" => "38?",
                        "H" => "39?",
                        "I" => "40?",
                        "J" => "41?",
                        "K" => "42?")
        end 
        if encoding == 64
            og_quals = replace.(og_quals,   
                        "@" => "00?", 
                        "A" => "01?",
                        "B" => "02?",
                        "C" => "03?",
                        "D" => "04?",
                        "E" => "05?",
                        "F" => "06?",
                        "G" => "07?",
                        "H" => "08?",
                        "I" => "09?",
                        "J" => "10?",
                        "K" => "11?",
                        "L" => "12?",
                        "M" => "13?",
                        "N" => "14?",
                        "O" => "15?",
                        "P" => "16?",
                        "Q" => "17?",
                        "R" => "18?",
                        "S" => "19?",
                        "T" => "20?",
                        "U" => "21?",
                        "V" => "22?",
                        "W" => "23?",
                        "X" => "24?",
                        "Y" => "25?",
                        "Z" => "26?",
                        "[" => "27?",
                        "\\" => "28?",
                        "]" => "29?",
                        "^" => "30?",
                        "_" => "31?",
                        "`" => "32?",
                        "a" => "33?",
                        "b" => "34?",
                        "c" => "35?",
                        "d" => "36?",
                        "e" => "37?",
                        "f" => "38?",
                        "g" => "39?",
                        "h" => "40?", 
                        "i" => "41?",
                        "j" => "42?")            
        end 
        og_quals = Array(split.(og_quals, "?"))
        pop!.(og_quals)
        df = DataFrame([String[] for i in 1:max_length], :auto)
        for i in 1:length(og_quals)
            println(string("starting loop #", i, " of ", length(og_quals)))
            temp = collect(og_quals[i])
            if length(temp) < max_length
                difference = max_length - length(temp)
                for r in 1:difference 
                    append!(temp, ["NaN"])
                end 
            end 
            parse.(Float64, temp)
            push!(df, temp)
        end 
        df = passmissing(parse).(Float64, df);
        for col in eachcol(df)
            replace!(col,NaN => Float64(100))
        end
        out = FastqRecord(split_fastq[2:4:end], split_fastq[3:4:end], string_quals, df)#, file)
    end  
    return out
end # end function readFastq

# Some benchmarking timing stats based on .fastq size 
    # 2 MB - 4.049943 seconds (54.18 M allocations: 3.418 GiB, 17.70% gc time, 0.11% compilation time)
    # 34.8 MB - 61.989768 seconds (900.48 M allocations: 56.452 GiB, 10.68% gc time)
    # 33.2 MB - 53.795583 seconds (825.51 M allocations: 51.120 GiB, 8.06% gc time)
    # 1.41 GB - 2646.207447 seconds (32.56 G allocations: 1.949 TiB, 21.35% gc time)