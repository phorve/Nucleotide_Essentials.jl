"""
    Nucleotide_Essentials.FastqtoFasta

Converts a FastqRecord to a FastaRecord. Can also input and convert a .fastq file to a FastaRecord in the same function. 

FastqtoFasta(Fastq::Union{String, FastqRecord})
.fastq file => FastqtoFasta(Fastq) => FastaRecord(ID, sequence, filename)
FastqRecord(ID, sequence, quality, filename) => FastqtoFasta(Fastq) => FastaRecord(ID, sequence, filename)

Supported keyword arguments include: 

* 'Fastq::Union{String, FastqRecord}': 
    * The full or relative path to a .fastq file
    * A FastqRecord

# Example: 

```julia
# Supply the path to a .fastq file that you would like to convert to a FastRecord
myfasta = FastqtoFasta("myfastq.fastq");

# Alternatively, a FastqRecord can be used as the input 
myfasta = FastqtoFasta(myFastqRecord);
```
"""
function FastqtoFasta(Fastq::Union{String, FastqRecord})
    if typeof(Fastq) == String
        input = readFastq(Fastq);
        newfasta = FastaRecord(input.ID, input.sequence, input.filename); 
    else
        newfasta = FastaRecord(Fastq.ID, Fastq.sequence, Fastq.filename);   
    end 
    return(newfasta)
end