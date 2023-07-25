using FASTX: sequence, FASTARecord, FASTAWriter
using ProgressMeter: @showprogress, Progress, next!
using Base.Threads: @threads, nthreads
println("threads: ", nthreads())
include("source/io.jl")
include("source/trypsin.jl")
include("source/palindrome.jl")
sprot = readfasta("data/sprot.fa");
try
    mkdir("data/tryptic")
catch IOError
    @warn "data/tryptic already exists"
end
minlength = 5
maxlength = 100
peptidename(i::Int, j::Int, minlength::Int, maxlength::Int) = "> sprot $(i) tryptic peptide $(j) minlength $(minlength) maxlength $(maxlength)\n"
p = Progress(length(sprot), 1, "Digesting...")
open("data/tryptic/sprot.fa", "w") do writer
   @threads for i=1:length(sprot)
        seq = sequence(sprot[i])
        tryptics = trypticpeptides(seq, minlength, maxlength)
        headers = peptidename.([i], 1:length(tryptics), [minlength], [maxlength])
        write.(writer, headers .* tryptics .* '\n')
        #for j=1:length(tryptic)
        #    header = peptidename(i, j, minlength, maxlength)
        #    write(writer, header * tryptic[j] * '\n')
        #end
        next!(p)
    end
end
