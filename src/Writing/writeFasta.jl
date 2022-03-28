using CodecZlib

"""
    Nucleotide_Essentials.writeFasta

    readFasta(Path::String)
    FastaRecord => write_fasta(input_fasta, out, compressed) => .fasta file/.fasta.gz file 

Creates a single or multiple entry FastaRecord and outputs either a .fasta or compressed .fasta.gz file to the desired directory

supported keyword arguments include: 

* 'input_fasta::FastaRecord': A FastaRecord with either a single entry or multiple entries 
* 'out::String': The full or relative path to the directory where files should be written to
* 'compressed::Bool': Whether or not to write the .fasta files as compressed files or not. 
    * If `true`, files will written as .fasta.gz files
    * If `false`, files will written as .fasta files

# Example: 

```julia
# .fasta files can be written as from an already imported FastaRecord in Julia 
myfasta = readFasta("myfasta.fasta");
writeFasta(input_fasta, "example/output/directory, false)

# .fasta files can be written as a .fasta.gz from an already imported FastaRecord in Julia 
myfasta = readFasta("myfasta.fasta");
writeFasta(input_fasta, "example/output/directory", true)

# .fasta files with multiple sequences can be read and written as individual .fasta or .fasta.gz in the same step
myfasta = readFasta("myfasta.fasta");
writeFasta(readFasta("/myfasta.fasta"), "example/output/directory", true);
```
"""
function writeFasta(input_fasta::FastaRecord, out::String, compressed::Bool)
    if compressed == true
        for i in 1:length(input_fasta.ID)
            temp_record = string(">", input_fasta.ID[i], "\n", input_fasta.sequence[i])        
            out_file = string(out, "/", input_fasta.ID[i], ".fasta")
            write(string(out, "/", input_fasta.ID[i], ".fasta"), temp_record)
            run(`gzip $out_file`); #this works
        end 
    else 
        for i in 1:length(input_fasta.ID)
            temp_record = string(">", input_fasta.ID[i], "\n", input_fasta.sequence[i])
            write(string(out, "/", input_fasta.ID[i], ".fasta"), temp_record)
        end 
    end     
end