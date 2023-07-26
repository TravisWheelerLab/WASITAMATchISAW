# figure 7
### distribution of LPS and LCS in tryptic peptides of Swiss-Prot
using Base.Threads: @threads
using FASTX: sequence
using ProgressMeter: Progress, next!
using DelimitedFiles: writedlm
### LPS and LCS
include("source/io.jl");
include("source/trypsin.jl");
include("source/palindrome.jl");
sprotall = [record for record=readfasta("data/sprot.fa") if 50 < length(sequence(record)) <= 2000]
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
