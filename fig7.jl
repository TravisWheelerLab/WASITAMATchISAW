using Base.Threads: @threads
using FASTX: sequence
using ProgressMeter: Progress, next!
using DelimitedFiles: writedlm
using Random: shuffle!

include("source/utils.jl")
include("source/io.jl");
include("source/trypsin.jl");
include("source/palindrome.jl");

pathinput = ARGS[1]
pathoutput = ARGS[2]
minlength = parse(Int, ARGS[3])
maxlength = parse(Int, ARGS[4])
doshuffle = "s" in ARGS

records = readfasta(pathinput)

sequences = [seq for seq=sequence.(records) if 50 < length(seq) < 2000]
shuffle!(sequences)
longestsequence = reduce(max, length.(sequences))

distribution = zeros(Int, (maxlength, maxlength))
if !doshuffle
    trypticpalindromedistribution!(distribution, sequences, minlength, maxlength)
else
    shuffledtrypticpalindromedistribution!(distribution, sequences, minlength, maxlength)
end

writedlm(pathoutput, distribution)
