"""
    Nucleotide_Essentials.FastaRecord

# Components

* ID: The unique sequence identifier associated with that entry
* sequence: The nucleotide sequence of that entry
* filename: The original file name 
"""
struct FastaRecord
    ID::Vector{SubString{String}}
    sequence::Vector{SubString{String}}
    filename::String
end

"""
    Nucleotide_Essentials.readFasta

Imports a .fasta file into julia

    readFasta(Path::String)
.fasta file => readFasta(Path) => FastaRecord(ID, sequence, filename)

supported keyword arguments include: 

* 'Path::String': The full or relative path to a .fasta file

# Example: 

```julia
# Supply the path to a .fasta file that you would like to import - it is recommended to include `;` in your command to prevent printing potentially large .fasta files in the REPL
myfasta = readFasta("myfasta.fasta");
```
"""
function readFasta(path::String)
    file = rsplit(path, "/", limit = 2)[2]
    input = read(path, String)
    split_fastq = split(input, ">", keepempty = false)
    hold_temp = []
    for s in 1:length(split_fastq)
        temp = split(split_fastq[s], "\n", limit = 2)
        hold_temp = vcat(hold_temp, temp)
    end 
out = FastaRecord(hold_temp[1:2:end], hold_temp[2:2:end], rsplit(path, "/", limit = 2)[2]);
displaysize(30,30)
return(out)
end 