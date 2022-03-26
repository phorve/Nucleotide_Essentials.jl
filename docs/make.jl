using Documenter
using Nucleotide_Essentials

makedocs(
    sitename = "Nucleotide_Essentials",
    format = Documenter.HTML(),
    modules = [Nucleotide_Essentials]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(;
    repo="github.com/phorve/Nucleotide_Essentials.jl",
)