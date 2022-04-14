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
    filename::String
end

"""
    Nucleotide_Essentials.readFastq
    readFastq(Path::String)
.fastq file => readFastq(Path) => FastqRecord(ID, sequence, quality, filename)

supported keyword arguments include: 

* ```Path::String```: The full or relative path to a .fastq file

# Example: 
```julia
# Supply the path to a .Fastq file that you would like to import
myfastq = readFastq("myfastq.fastq")
```
"""
function readFastq(Path::String)
    try
        global file = rsplit(Path, "/", limit = 2)[2]
    catch err
        if isa(err, BoundsError)
            global file = rsplit(Path, "/", limit = 2)[1]
        end 
    end 
    input = read(Path, String)
    input2 = replace(input, "\n@" => "Z")
    input3 = replace(input2, "@" => "Z", count = 1)
    input4 = replace(input3, "Z" => "\n")
    split_fastq = split(input4, "\n")
    out = FastqRecord(split_fastq[2:4:end], split_fastq[3:4:end], split_fastq[5:4:end], file) 
    return out
end # end function readFastq