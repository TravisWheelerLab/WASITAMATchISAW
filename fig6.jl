using DelimitedFiles: writedlm
using SparseArrays: spzeros
using ProgressMeter: Progress, next!
using Base.Threads: @threads

include("source/chromosome.jl")
include("source/io.jl")
include("source/palindrome.jl")

pathchromosome = ARGS[1] # data/ncbi_chr22_sequence.fasta
pathannotation = ARGS[2] # data/ncbi_chr22_protein.txt, data/ncbi_chr22_rna.txt, data/ncbi_chr22_pseudo.txt
pathoutput = ARGS[3] # outputs/...

chromosome = loadchromosome(pathchromosome)
annotation = loadannotation(pathannotation)

frames = collect(annotate(chromosome, annotation))
framelengths = length.(frames)

lps = longestpalindromicsubstring.(frames)
lpslengths = (x::Tuple{Int, Int} -> x[2]-x[1]).(lps)

writedlm(pathoutput, hcat(framelengths, lpslengths))