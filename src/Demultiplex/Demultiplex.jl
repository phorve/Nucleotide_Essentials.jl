"""
Demultiplex.jl


"""

using DataFrames, CSV, Glob

function potential_mismatches(og_barcode::String, mismatch::Int64 = 1)
    potentials = []
    barcode = collect(og_barcode)
    if mismatch == 1
        for spot in 1:length(barcode)
            barcode = collect(og_barcode)
            barcode[spot] = 'G'
            push!(potentials, String(barcode))
            barcode[spot] = 'C'
            push!(potentials, String(barcode))
            barcode[spot] = 'A'
            push!(potentials, String(barcode))
            barcode[spot] = 'T'
            push!(potentials, String(barcode))
        end 
        potentials = unique(potentials)
    end
    if mismatch == 2
        for spot in 1:length(barcode)
            if spot == length(barcode)
                continue
            else
                barcode = collect(og_barcode)
                barcode[spot] = 'G'
                barcode[spot+1] = 'G'
                push!(potentials, String(barcode))
                barcode[spot] = 'G'
                barcode[spot+1] = 'C'
                push!(potentials, String(barcode))
                barcode[spot] = 'G'
                barcode[spot+1] = 'A'
                push!(potentials, String(barcode))
                barcode[spot] = 'G'
                barcode[spot+1] = 'T'
                push!(potentials, String(barcode))

                barcode[spot] = 'C'
                barcode[spot+1] = 'G'
                push!(potentials, String(barcode))
                barcode[spot] = 'C'
                barcode[spot+1] = 'C'
                push!(potentials, String(barcode))
                barcode[spot] = 'C'
                barcode[spot+1] = 'A'
                push!(potentials, String(barcode))
                barcode[spot] = 'C'
                barcode[spot+1] = 'T'
                push!(potentials, String(barcode))

                barcode[spot] = 'A'
                barcode[spot+1] = 'G'
                push!(potentials, String(barcode))
                barcode[spot] = 'A'
                barcode[spot+1] = 'C'
                push!(potentials, String(barcode))
                barcode[spot] = 'A'
                barcode[spot+1] = 'A'
                push!(potentials, String(barcode))
                barcode[spot] = 'A'
                barcode[spot+1] = 'T'
                push!(potentials, String(barcode))

                barcode[spot] = 'T'
                barcode[spot+1] = 'G'
                push!(potentials, String(barcode))
                barcode[spot] = 'T'
                barcode[spot+1] = 'C'
                push!(potentials, String(barcode))
                barcode[spot] = 'T'
                barcode[spot+1] = 'A'
                push!(potentials, String(barcode))
                barcode[spot] = 'T'
                barcode[spot+1] = 'T'
                push!(potentials, String(barcode))
            end 
        end 
        potentials = unique(potentials)
    end 
    return potentials 
end 


function demultiplex(R1::String, Map::String, Dir::String, Reads::String="Single", R2::String="", Out::String="Demultiplexing_Output", verbose::Bool=false)
    mkdir(string(Dir,"/",Out))
    pre_dmx = readFastq(R1) # this is our forward multiplexed reads
    if Reads == "Paired"
        pre_dmx2 = readFastq(R2) # this is our forward multiplexed reads
    end 
    seqs = pre_dmx.sequence
    qual = pre_dmx.quality
    og_name = pre_dmx.ID
    if Reads == "Paired"
        seqs2 = pre_dmx2.sequence
        qual2 = pre_dmx2.quality
        og_name2 = pre_dmx2.ID
    end 

    # set up our mapping file - this is our mapping file in txt format 
    if endswith(Map, ".csv") == true
        mapping = CSV.read(Map, DataFrame) 
        barcodes = Array(mapping[!,:BarcodeSequence])
        IDs = Array(mapping[!,:"#SampleID"])
        
        # Make sure the number of barcodes matches the number of sample IDs
        if length(unique(barcodes)) != length(unique(IDs))
            error("The number of barcodes and unique samples should match, check mapping file for errors")
        end 
    end 

    # set up our mapping file - this is our mapping file in txt format 
    if endswith(Map, ".txt") == true
        mapping = CSV.read(map, DataFrame) 
        barcodes = Array(mapping[!,:BarcodeSequence])
        IDs = Array(mapping[!,:"#SampleID"])

        # Make sure the number of barcodes matches the number of sample IDs
        if length(unique(barcodes)) != length(unique(IDs))
            error("The number of barcodes and unique samples should match, check mapping file for errors")
        end 
    end 

    # dmuxing
    for b in 1:length(barcodes)
        for s in 1:length(seqs)
            if occursin(barcodes[b], seqs[s])
                sample = string(IDs[b])
                file_end = string("_R1.fastq")
                if isfile(string(Dir,"/",Out, "/",sample, file_end)) == false
                    out_string = string("@",IDs[b],"\n",seqs[s],"\n+\n",qual[s])
                    write(string(Dir,"/",Out, "/",sample, file_end), out_string)
                else 
                    out_string2 = string("\n@",IDs[b],"\n",seqs[s],"\n+\n",qual[s])
                    iostream =  open(string(Dir,"/",Out, "/",sample, file_end), "a")
                    write(iostream, out_string2);
                    close(iostream)
                end
            end 
        end 
    end 
    for s in 1:length(seqs)
        count = 0
        for b in barcodes
            if occursin(b, seqs[s]) == false
                count = count+1
            else 
                count = count-1
            end 
        end
        if count == length(barcodes)
            if isfile(string(Dir,"/",Out, "/unassigned_R1.fastq")) == false
                out_string3 = string("@",og_name[s],"\n",seqs[s],"\n+\n",qual[s])
                write(string(Dir,"/",Out, "/unassigned_R1.fastq"), out_string3)
            else 
                out_string4 = string("\n@",og_name[s],"\n",seqs[s],"\n+\n",qual[s])
                iostream =  open(string(Dir,"/",Out, "/unassigned_R1.fastq"), "a")
                write(iostream, out_string4);
                close(iostream)
            end
        end 
    end  
    # We still need to work on this paired portion
    if Reads == "Paired"
        for b in 1:length(barcodes)
            for s in 1:length(seqs2)
                if occursin(barcodes[b], seqs2[s])
                    sample = string(IDs[b])
                    file_end = string("_R2.fastq")
                    if isfile(string(Dir,"/",Out, "/",sample, file_end)) == false
                        out_string5 = string("@",IDs[b],"\n",seqs2[s],"\n+\n",qual2[s])
                        write(string(Dir,"/",Out, "/",sample, file_end), out_string5)
                    else 
                        out_string6 = string("\n@",IDs[b],"\n",seqs2[s],"\n+\n",qual2[s])
                        iostream =  open(string(Dir,"/",Out, "/",sample, file_end), "a")
                        write(iostream, out_string6);
                        close(iostream)
                    end
                end 
            end 
        end 
        for s in 1:length(seqs2)
            count = 0
            for b in barcodes
                if occursin(b, seqs2[s]) == false
                    count = count+1
                else 
                    count = count-1
                end 
            end
            if count == length(barcodes)
                if isfile(string(Dir,"/",Out, "/unassigned_R2.fastq")) == false
                    out_string7 = string("@",og_name2[s],"\n",seqs2[s],"\n+\n",qual2[s])
                    write(string(Dir,"/",Out, "/unassigned_R2.fastq"), out_string7)
                else 
                    out_string8 = string("\n@",og_name2[s],"\n",seqs2[s],"\n+\n",qual2[s])
                    iostream =  open(string(Dir,"/",Out, "/unassigned_R2.fastq"), "a")
                    write(iostream, out_string8);
                    close(iostream)
                end
            end 
        end   
    end 

    # output statistics on demultiplexed samples 
    sampleholds = []
    reads = []
    files = glob("*.fastq", string(Dir,"/",Out))
    for file in files 
        push!(sampleholds, rsplit(file, "/", limit = 2)[2])
        push!(reads, (countlines(file)/4))
    end 
    stats = DataFrame(Filename=sampleholds, AssignedReads=reads)
    CSV.write(string(Dir,"/",Out,"/demultiplex_stats.csv"), stats)
end # end function demultiplex