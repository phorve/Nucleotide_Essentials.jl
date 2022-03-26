using   DataFrames,
        StatsBase,
        DataFramesMeta,
        Plots,
        StatsPlots

"""
    Nucleotide_Essentials.PlotQuality

Returns a plot of the quality profile of a FastqRecord

This function plots a visual summary of the distribution of quality scores as a function of sequence position for the input fastq file(s).

The plotted lines show summary statistics at each sequence position:

* green is the mean 
* dashed red lines are the 25th and 75th quantiles

supported keyword arguments include:

* 'Input::FastqRecord': The name of a FastqRecord for plotting
* 'QualityType::String' (optional): The system to use for parsing the fastq quality information. Potential quality encodings include: 
    * Phred64 (defualt)
    * Ascii (not yet supported)
    * Phred33 (not yet supported)    
* 'verbose::Bool' (optional): Whether or not to show some intermediary feedback on the progress of the function (default = false)

# Example:
```julia
# A quality profile can be created by directly calling an already read FastqRecord
PlotQuality(myfastq)

# Alternatively, a quality profile can be created by nesting readFastq() inside of PlotQuality() 
PlotQuality(readFastq("myfastq.fastq"))
```
"""
function PlotQuality(Input::FastqRecord, QualityType::String = "phred64", verbose::Bool = false)
    # make some empty arrays for data to be stored 
    quality_data = []
    lengths = []
    
    # information on the number of reads coming in for our plot 
    reads_thisone = length(Input.sequence)
    annotation = "Reads: $reads_thisone"
    
    # read in the sequence quality information 
    qualities = Input.quality
    
    # loop through the quality information and convert to something we can plot 
    for q in 1:length(qualities)
        if verbose == true
            println(round((q/length(qualities))*100; digits = 3), "% complete")
        end 
        q_vector = (String(collect(qualities[q])))
        q_string = replace(q_vector,   
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
            "I" => "40?")
        
        out_vector = Array(split(q_string, "?", limit = (length(q_vector)+1)))
        out_vector = deleteat!(out_vector, (length(q_vector)+1))      
        temp_data = Array(parse.(Float64, out_vector))
        this_length = length(temp_data)
        quality_data = push!(quality_data,temp_data)
        lengths = push!(lengths, this_length)
    end 
    if verbose == true
        println("=====Calculating quality statistics=====")
    end 
    output = resize!.(quality_data, maximum(length, quality_data))
    readlength = Int64(median(lengths))
    output = DataFrame(output, :auto)
    output = first(output, readlength)
    output1 = ifelse.(output .< 0.1, NaN, output)
    # Prep the data for plotting 
    output2 = DataFrame([[names(output1)]; collect.(eachrow(output1))], [:column; Symbol.(axes(output1, 1))])
    output3 = select!(output2, Not(:1))
    output4 = filter(row -> all(x -> !(x isa Number && isnan(x)), row), output3)
    output_long = stack(output4)
    summarized = describe(output4)
    summarized[!, :Cycle] = (1:nrow(summarized))
    summarized.Cycle = 1:nrow(summarized)
    cycles = Vector(output_long.variable)
    cycles = Array(parse.(Float64, cycles))
    output_long.Cycle = cycles
    res = @transform groupby(output_long, :variable) begin
       :quantile_1 = quantile(:value, 0.25)
       :quantile_2 = quantile(:value, 0.5)
       :quantile_3 = quantile(:value, 0.75)
    end
    res.Cycle = cycles

    # make our actual plot 
    if verbose == true
        println("=====Creating Quality Plot=====")
    end 
    gr()
    qualityplot = begin
        @df summarized Plots.plot(
            :Cycle, 
            :mean, 
            legend = false, 
            title = Input.filename,
            linecolor = :green, 
            titlefontsize = 8,
            xlabel = "Cycle", 
            ylabel = "Quality Score")
        @df res Plots.plot!(:Cycle, :quantile_1, linealpha = 0.7, linestyle = :dot, linecolor = :red)
        @df res Plots.plot!(:Cycle, :quantile_3, linealpha = 0.7, linestyle = :dot, linecolor = :red)
    ylims!((0,40))
    annotate!(50, 2, Plots.text(annotation))
    end 

    return qualityplot
end 