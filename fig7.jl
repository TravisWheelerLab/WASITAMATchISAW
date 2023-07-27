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
if length(ARGS) > 4
    mode = ARGS[5]
else
    mode = "-"
end

records = readfasta(pathinput)

sequences = [seq for seq=sequence.(records) if 50 < length(seq) < 2000]
shuffle!(sequences)
longestsequence = reduce(max, length.(sequences))

distribution = zeros(Int, (maxlength, maxlength))
if mode == "-"
    trypticpalindromedistribution!(distribution, sequences, minlength, maxlength)
elseif mode == "s"
    shuffledtrypticpalindromedistribution!(distribution, sequences, minlength, maxlength)
else
    @warn "unrecognized mode argument. the output matrix is blank."
end

writedlm(pathoutput, distribution)
