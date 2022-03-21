"""
readFastQ.jl


"""
struct FastqRecord
    ID::Vector{SubString{String}}
    sequence::Vector{SubString{String}}
    quality::Vector{SubString{String}}
end

function readFastq(Path::String, quality::Number=1)
    input = read(Path, String)
    input2 = replace(input, "\n@" => "Z")
    input3 = replace(input2, "@" => "Z", count = 1)
    input4 = replace(input3, "Z" => "\n")
    split_fastq = split(input4, "\n")
    out = FastqRecord(split_fastq[2:4:end], split_fastq[3:4:end], split_fastq[5:4:end]) 
    return out
end # end function readFastq