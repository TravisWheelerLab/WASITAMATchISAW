using DelimitedFiles: writedlm
using Random: shuffle!
using ProgressMeter: Progress, next!
using Base.Threads: @threads

include("source/alignment.jl")
include("source/chromosome.jl")
include("source/io.jl")
include("source/palindrome.jl")
include("source/utils.jl")

pathchromosome = ARGS[1] # data/ncbi_chr22_sequence.fasta
pathannotation = ARGS[2] # data/ncbi_chr22_protein.txt, data/ncbi_chr22_rna.txt, data/ncbi_chr22_pseudo.txt
pathoutput = ARGS[3] # outputs/...
mode = ARGS[4]
doshuffle = "s" in ARGS

chromosome = loadchromosome(pathchromosome)
annotation = loadannotation(pathannotation)

frames = annotate(chromosome, annotation)
framelengths = length.(frames)
println(reduce(min, framelengths))
println(reduce(max, framelengths))
if doshuffle
    frames = shufflefast.(frames)
end

if mode == "lps"
    lps = longestpalindromicsubstring.(frames)
    lpslengths = (x::Tuple{Int, Int} -> x[2]-x[1]).(lps)
    writedlm(pathoutput, hcat(framelengths, lpslengths))
elseif mode == "alnrev"
    reversedframes = reverse.(frames)
    formatter = x::PairwiseAlignmentResult -> score(x)
    scores = align(Pairwise(), frames, reversedframes; formatter = formatter)
    writedlm(pathoutput, hcat(framelengths, scores))
end