using Base.Threads: @threads
using FASTX: sequence
using ProgressMeter: Progress, next!
using DelimitedFiles: writedlm
using Random: shuffle!
include("source/utils.jl")
include("source/io.jl");
include("source/trypsin.jl");
include("source/palindrome.jl");
# parse arguments
pathinput = ARGS[1]
pathoutput = ARGS[2]
minlength = parse(Int, ARGS[3])
maxlength = parse(Int, ARGS[4])
dopermute = "permute" in ARGS
# load the dataset
records = readfasta(pathinput)
# take sequences between lengths 50 and 2000
sequences = [seq for seq=sequence.(records) if 50 < length(seq) < 2000]
# splice out null characters
sequences = replace.(sequences, "X"=>"")
# reorder sequences - this helps the progress bar make more accurate estimates
shuffle!(sequences)
# collect the distribution of LPS lengths by sequence lengths
distribution = zeros(Int, (maxlength, maxlength))
if !dopermute
    trypticpalindromedistribution!(distribution, sequences, minlength, maxlength)
else
    permutedtrypticpalindromedistribution!(distribution, sequences, minlength, maxlength)
end
# 
writedlm(pathoutput, distribution)
