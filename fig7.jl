using FASTX: sequence, FASTARecord, FASTAWriter
using ProgressMeter: @showprogress, Progress, next!
using Base.Threads: @threads, nthreads
using Random: shuffle!
include("source/io.jl")
include("source/trypsin.jl")
include("source/palindrome.jl")
### partition Swiss-Prot into `nthreads` batches
sprot = [record for record=readfasta("data/sprot.fa") if 50 < length(sequence(record)) < 2000]
nsprot = length(sprot)
batch_size = ceil(Int, nsprot / nthreads())
sprot_batches = [sprot[1+batch_size*(i-1):min(nsprot, batch_size*i)] for i=1:nthreads()]
length.(sprot_batches)
### LPS and LCS
minlength = 5
maxlength = 100
p = Progress(length(sprot), 1, "Digesting...")
tempfiles = ["outputs/.temp_trypticpalindromes_thread$(t).txt" for t=1:nthreads()]
@threads for t=1:nthreads()
    open(tempfiles[t], "w") do writer
        for i=1:length(sprot_batches[t])
            seq = sequence(sprot_batches[t][i])
            trypticintervals = trypticpeptides(seq, minlength, maxlength)
            trypticsequences = (seq[start:stop] for (start,stop)=trypticintervals)
            lps = longestpalindromicsubstring.(trypticsequences)
            lps = (x->"\t$(x[1]),$(x[2])").(lps)
            write(writer, "$((t-1)*batch_size+i)")
            write.(writer, lps)
            write(writer, '\n')
            flush(writer)
            next!(p)
        end
    end
end
