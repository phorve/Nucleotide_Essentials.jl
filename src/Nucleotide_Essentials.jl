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
        FastqRecord,
        potential_mismatches,
        demultiplex_se,
        demultiplex_pe,
        reverse_complement, 
        PlotQuality

include("Reading/readFastq.jl")
include("Demultiplex/Demultiplex.jl")
include("Visualization/PlotQuality.jl")

end  # end module 