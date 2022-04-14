using   DataFrames,
        StatsBase,
        DataFramesMeta,
        Plots,
        StatsPlots

"""
    Nucleotide_Essentials.PlotQuality

Returns a plot of the quality profile of a .fastq or .fastq.gz file 

This function plots a visual summary of the distribution of quality scores (automatically detects Phred+33 or Phred+64 encoding) as a function of sequence position for the input fastq file(s).

The plotted lines show summary statistics at each sequence position:

* green is the mean 
* dashed red lines are the 25th and 75th quantiles

Supported keyword arguments include:

* ```Input::FastqRecord```: The name of a FastqRecord for plotting
* ```verbose::Bool``` (optional): Whether or not to show some intermediary feedback on the progress of the function (default = false)
* ```outputfigure::Bool``` (optional): Whether or not to output a .png file with the created QualityPlot (default = false)
* ```figurepath::String``` (optional): If outputting a .png figure to file, specify the path to a directory where the file should be written to (default = ```pwd()```)

# Example:
```julia
# A quality profile can be created by supply the path to a .fastq or .fastq.gz file 
PlotQuality("path/to/my/file.fastq")
```
"""
function PlotQuality(Input::String, verbose::Bool = false, outputfigure::Bool = false, figurepath::String = pwd())
    if verbose == true
        println(string("Reading file: ", Input))
    end 

    if endswith(Input, ".gz") === true
        if verbose == true
            println(string("""|
                        |----Unzipping """, Input))
        end 
        run(`gunzip $Input`); #this works
        unzipped = true
    else
        unzipped = false
    end 

    if endswith(Input, ".gz") === true
        R1_fixed = replace(Input, ".fastq.gz"=>".fastq")
        R1 = readFastq(R1_fixed)
    else 
        R1 = readFastq(Input)
    end 
    quals = R1.quality

    if occursin("^", join(quals)) == true && occursin("a", join(quals)) == true && occursin("]", join(quals)) == true && occursin("f", join(quals)) == true
        encoding = 64;
        if verbose == true
            println("""|
            |----Detected Phred+64 quality encoding""")
        end 
    end 
    if occursin("^", join(quals)) == false && occursin("a", join(quals)) == false && occursin("]", join(quals)) == false && occursin("f", join(quals)) == false
        if verbose == true
            println("""|
            |----Detected Phred+33 quality encoding""")
        end 
        encoding = 33;
    end 

    if encoding == 33
        scores_33 = Dict(
            ('!') => "0", 
            ('\"') => "1",
            ('#') => "2",
            ('\$') => "3",
            ('%') => "4",
            ('&') => "5",
            ('\'') => "6",
            ('(') => "7",
            (')') => "8",
            ('*') => "9",
            ('+') => "10",
            (',') => "11",
            ('-') => "12",
            ('.') => "13",
            ('/') => "14",
            ('0') => "15",
            ('1') => "16",
            ('2') => "17",
            ('3') => "18",
            ('4') => "19",
            ('5') => "20",
            ('6') => "21",
            ('7') => "22",
            ('8') => "23",
            ('9') => "24",
            (':') => "25",
            (';') => "26",
            ('<') => "27",
            ('=') => "28",
            ('>') => "29",
            ('?') => "30",
            ('@') => "31",
            ('A') => "32",
            ('B') => "33",
            ('C') => "34",
            ('D') => "35",
            ('E') => "36",
            ('F') => "37",
            ('G') => "38",
            ('H') => "39",
            ('I') => "40", 
            ('J') => "41",
            ('K') => "42"
        )
    end 
    if encoding == 64
        scores_64 = Dict(
            ('@') => "0", 
            ('A') => "1",
            ('B') => "2",
            ('C') => "3",
            ('D') => "4",
            ('E') => "5",
            ('F') => "6",
            ('G') => "7",
            ('H') => "8",
            ('I') => "9",
            ('J') => "10",
            ('K') => "11",
            ('L') => "12",
            ('M') => "13",
            ('N') => "14",
            ('O') => "15",
            ('P') => "16",
            ('Q') => "17",
            ('R') => "18",
            ('S') => "19",
            ('T') => "20",
            ('U') => "21",
            ('V') => "22",
            ('W') => "23",
            ('X') => "24",
            ('Y') => "25",
            ('Z') => "26",
            ('[') => "27",
            ('\\') => "28",
            (']') => "29",
            ('^') => "30",
            ('_') => "31",
            ('`') => "32",
            ('a') => "33",
            ('b') => "34",
            ('c') => "35",
            ('d') => "36",
            ('e') => "37",
            ('f') => "38",
            ('g') => "39",
            ('h') => "40", 
            ('i') => "41",
            ('j') => "42"
        )
    end 
    
    function convert_phred(c)
        if encoding == 33
            scores_33[c]
        else 
            scores_64[c]
        end 
    end
    
    max_length = maximum(length.(quals))
    converted_quals = DataFrame([Float64[] for i in 1:max_length], :auto)
    for i in 1:length(quals)
        temp = parse.(Float64, (map(convert_phred, collect(quals[i]))))
        if length(temp) < max_length
            difference = max_length - length(temp)
            for r in 1:difference 
                append!(temp, [0.00])
            end 
        end 
        converted_quals = push!(converted_quals, temp)
    end 
    
    summarized = describe(converted_quals)
    summarized[!, :Cycle] = (1:nrow(summarized))
    summarized.Cycle = 1:nrow(summarized)

    output_long = stack(converted_quals)
    output_long.variable .= replace.(output_long.variable, "x" => "")
    cycles = Vector(output_long.variable)
    cycles = Array(parse.(Float64, cycles))
    output_long.Cycle = cycles
    res = @transform groupby(output_long, :variable) begin
       :quantile_1 = quantile(:value, 0.25)
       :quantile_2 = quantile(:value, 0.5)
       :quantile_3 = quantile(:value, 0.75)
    end
    res.Cycle = cycles
    annotation = nrow(converted_quals)
    # make our actual plot 
    if verbose == true
        println("""|
        |----Creating Quality Plot""")
    end 
    gr()
    qualityplot = begin
        @df summarized Plots.plot(
            :Cycle, 
            :mean, 
            legend = false, 
            title = R1.filename,
            linecolor = :green, 
            titlefontsize = 8,
            xlabel = "Cycle", 
            ylabel = "Quality Score")
        @df res Plots.plot!(:Cycle, :quantile_1, linealpha = 0.7, linestyle = :dot, linecolor = :red)
        @df res Plots.plot!(:Cycle, :quantile_3, linealpha = 0.7, linestyle = :dot, linecolor = :red)
    ylims!((0,42))
    #vline!([median(length.(quals))]) # show the median length of the reads 
    annotate!(50, 2, Plots.text(string("Reads: ", annotation), 12))
    end 

    if outputfigure == true
        fig_out = string(figurepath, "/", R1.filename, "_qualityplot.png")
        if verbose == true
            println(string("""|
            |----Saving output file to """, fig_out))
        end 
        savefig(qualityplot, fig_out)
    end 

    if unzipped == true
        if verbose == true
            println(string("""|
            |----Cleaning things up"""))
        end 
        run(`gzip $R1_fixed`);
    end 
    return qualityplot
end 