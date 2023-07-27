# figure 7
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
longestpalindrome = parse(Int, ARGS[5])
if length(ARGS) > 5
    mode = ARGS[6]
else
    mode = "-"
end

records = readfasta(pathinput)

sequences = [seq for seq=sequence.(records) if 50 < length(seq) < 2000]
shuffle!(sequences)
longestsequence = reduce(max, length.(sequences))

distribution = zeros(Int, (longestsequence, longestpalindrome))
if mode == "-"
    trypticpalindromedistribution!(distribution, sequences, minlength, maxlength)
elseif mode == "s"
    shuffledtrypticpalindromedistribution!(distribution, sequences, minlength, maxlength)
end

writedlm(pathoutput, distribution)
