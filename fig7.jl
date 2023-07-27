# figure 7
using Base.Threads: @threads
using FASTX: sequence
using ProgressMeter: Progress, next!
using DelimitedFiles: writedlm

include("source/io.jl");
include("source/trypsin.jl");
include("source/palindrome.jl");
# method 1 - broadcast view
println("method 1 - broadcast view")
sprotall = [record for record=readfasta("data/sprot.fa") if 50 < length(sequence(record)) <= 2000]
sprotall = sprotall[1:12000]
nsprotall = length(sprotall)
palindromedistribution = zeros(Int, (2000, 2000))
minlength = 5
maxlength = 100
p = Progress(length(sprotall), 1, "Digesting...")
@threads for i=1:nsprotall
    seq = sequence(sprotall[i])
    trypticintervals = trypticpeptides(seq, minlength, maxlength)
    trypticsequences = view.(seq, (start:stop for (start, stop)=trypticintervals))
    lps = longestpalindromicsubstring.(trypticsequences)
    for (start, stop)=lps
        palindromelength = stop - start
        sequencelength = length(seq)
        palindromedistribution[palindromelength, sequencelength] += 1
    end
    next!(p)
end
writedlm("outputs/sprot_all-tryptic-palindrome-distribution.dlm", palindromedistribution)
# method 2 - lazily generate substrings
println("method 2 - lazily generate substrings")
sprotall = [record for record=readfasta("data/sprot.fa") if 50 < length(sequence(record)) <= 2000]
sprotall = sprotall[1:12000]
nsprotall = length(sprotall)
palindromedistribution = zeros(Int, (2000, 2000))
minlength = 5
maxlength = 100
p = Progress(length(sprotall), 1, "Digesting...")
@threads for i=1:nsprotall
    seq = sequence(sprotall[i])
    trypticintervals = trypticpeptides(seq, minlength, maxlength)
    trypticsequences = (seq[start:stop] for (start,stop)=trypticintervals)
    lps = longestpalindromicsubstring.(trypticsequences)
    for (start, stop)=lps
        palindromelength = stop - start
        sequencelength = length(seq)
        palindromedistribution[palindromelength, sequencelength] += 1
    end
    next!(p)
end
writedlm("outputs/sprot_all-tryptic-palindrome-distribution.dlm", palindromedistribution)
