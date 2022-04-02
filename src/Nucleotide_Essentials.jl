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
        FastqtoFasta,
        potential_mismatches,
        demultiplex_se,
        demultiplex_pe,
        reverse_complement, 
        PlotQuality,
        FilterQuality_se,
        FilterQuality_pe

include("Reading/readFastq.jl")
include("Reading/readFasta.jl")
include("Writing/writeFasta.jl")
include("Writing/FastqtoFasta.jl")
include("Demultiplex/Demultiplex.jl")
include("QualityTesting/PlotQuality.jl")
include("QualityTesting/FilterQuality.jl")

end  # end module 