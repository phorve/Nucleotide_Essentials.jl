using   DataFrames,
        StatsBase,
        DataFramesMeta,
        Plots,
        StatsPlots

"""
    Nucleotide_Essentials.PlotQuality

Returns a plot of the quality profile of a FastqRecord

This function plots a visual summary of the distribution of quality scores (automatically detects Phred+33 or Phred+64 encoding) as a function of sequence position for the input fastq file(s).

The plotted lines show summary statistics at each sequence position:

* green is the mean 
* dashed red lines are the 25th and 75th quantiles

Supported keyword arguments include:

* 'Input::FastqRecord': The name of a FastqRecord for plotting
* 'verbose::Bool' (optional): Whether or not to show some intermediary feedback on the progress of the function (default = false)

# Example:
```julia
# A quality profile can be created by directly calling an already read FastqRecord
PlotQuality(myfastq)

# Alternatively, a quality profile can be created by nesting readFastq() inside of PlotQuality() 
PlotQuality(readFastq("myfastq.fastq"))
```
"""
function PlotQuality(Input::FastqRecord, verbose::Bool = false)
    # make some empty arrays for data to be stored 
    quality_data = []
    lengths = []
    
    if occursin("^", join(qualities)) == true && occursin("a", join(qualities)) == true && occursin("]", join(qualities)) == true && occursin("f", join(qualities)) == true
    encoding = 64
    println("Detected Phred64 quality encoding")
    end 

    if occursin("^", join(qualities)) == false && occursin("a", join(qualities)) == false && occursin("]", join(qualities)) == false && occursin("f", join(qualities)) == false
        println("Detected Phred33 quality encoding")
        encoding = 33
    end 

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
        if encoding == 33
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
        end 
        if encoding == 64
            q_string = replace(q_vector,   
                "@" => "00?", 
                "A" => "01?",
                "B" => "02?",
                "C" => "03?",
                "D" => "04?",
                "E" => "05?",
                "F" => "06?",
                "G" => "07?",
                "H" => "08?",
                "I" => "09?",
                "J" => "10?",
                "K" => "11?",
                "L" => "12?",
                "M" => "13?",
                "N" => "14?",
                "O" => "15?",
                "P" => "16?",
                "Q" => "17?",
                "R" => "18?",
                "S" => "19?",
                "T" => "20?",
                "U" => "21?",
                "V" => "22?",
                "W" => "23?",
                "X" => "24?",
                "Y" => "25?",
                "Z" => "26?",
                "[" => "27?",
                "\\" => "28?",
                "]" => "29?",
                "^" => "30?",
                "_" => "31?",
                "`" => "32?",
                "a" => "33?",
                "b" => "34?",
                "c" => "35?",
                "d" => "36?",
                "e" => "37?",
                "f" => "38?",
                "g" => "39?",
                "h" => "40?", 
                "i" => "41?",
                "j" => "42?")
        end 
        
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