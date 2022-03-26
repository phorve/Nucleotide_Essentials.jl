# demultiplex_test.jl
@testset "demultiplex test" begin 
    @test reverse_complement("ATCGT") == "ACGAT"
    @test reverse_complement("ATCGTG") == "CACGAT"
    @test reverse_complement("TGCCGTGACGT") == "ACGTCACGGCA"
    @test reverse_complement("ATCG") == "CGAT"
    @test potential_mismatches("GCGT", 1) == Any["GCGT", 
                                                 "CCGT", 
                                                 "ACGT", 
                                                 "TCGT", 
                                                 "GGGT", 
                                                 "GAGT", 
                                                 "GTGT", 
                                                 "GCCT", 
                                                 "GCAT", 
                                                 "GCTT", 
                                                 "GCGG", 
                                                 "GCGC", 
                                                 "GCGA", 
                                                 "CGT", 
                                                 "GGT", 
                                                 "GCT", 
                                                 "GCG"]
end 