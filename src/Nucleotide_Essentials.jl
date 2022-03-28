module Nucleotide_Essentials

import  CSV, 
        DataFrames, 
        Glob,
        StatsBase,
        DataFramesMeta,
        Plots,
        StatsPlots,
        Logging

export  readFastq,
        readFasta,
        FastqRecord,
        FastaRecord,
        writeFasta,
        potential_mismatches,
        demultiplex_se,
        demultiplex_pe,
        reverse_complement, 
        PlotQuality

include("Reading/readFastq.jl")
include("Reading/readFasta.jl")
include("Writing/writeFasta.jl")
include("Demultiplex/Demultiplex.jl")
include("Visualization/PlotQuality.jl")

end  # end module 