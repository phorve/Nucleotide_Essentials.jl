module Nucleotide_Essentials

import  CSV, 
        DataFrames, 
        Glob

export  readFastq,
        FastqRecord,
        potential_mismatches,
        demultiplex 

include("Reading/readFastq.jl")
include("Demultiplex/Demultiplex.jl")

end  # end module 