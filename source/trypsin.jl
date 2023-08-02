using Random: shuffle, shuffle!
using Base.Threads: @threads

trypsin_cleave_locations(q::AbstractString) = [i for i=1:length(q) if Char(q[i]) == 'K' || Char(q[i]) == 'R']

isproline(c::Char) = c == 'P'

# TODO - lazy implementation. maybe there's a functional expression?
"""
# sample all right bounded half open intervals [START, stop) for all 
# start, stop ∈ `locii` such that start < stop.
"""
function rightbounded_halfopen_subintervals(
    locii::Vector{Int}, 
    minlength::Int, 
    maxlength::Int
)
    n = length(locii)
    inbounds(a::Int, b::Int) = minlength <= interval_length((a, b)) <= maxlength
    [(locii[i], locii[j]-1) for i=1:n for j=i+1:n if inbounds(locii[i], locii[j]-1)]
end

function trypticpeptides(
    q::AbstractString,
    minlength::Int,
    maxlength::Int
)
    locii = trypsin_cleave_locations(q)
    proline = isproline.(collect(q))
    cleavable_locii = filter(i -> i==1 ? true : !proline[i-1], locii)
    rightbounded_halfopen_subintervals(cleavable_locii, minlength, maxlength)
end

function trypticpalindromedistribution!(
    distribution::Matrix{Int}, 
    seq::AbstractString, 
    minlength::Int, 
    maxlength::Int
)
    trypticintervals = trypticpeptides(seq, minlength, maxlength)
    trypticlengths = interval_length.(trypticintervals)
    n = length(trypticintervals)
    lazy_lps = (longestpalindromicsubstring(seq[start:stop]) 
                for (start,stop)=trypticintervals)
    palindromelengths = interval_length.(lazy_lps)
    for i=1:n
        distribution[trypticlengths[i], palindromelengths[i]] += 1
    end
end

function trypticpalindromedistribution!(
    distribution::Matrix{Int}, 
    seqs::Vector, 
    minlength::Int, 
    maxlength::Int
)
    nseqs = length(seqs)
    p = Progress(nseqs, 1, "Digesting...")
    @threads for i=1:nseqs
        trypticpalindromedistribution!(distribution, seqs[i], minlength, maxlength)
        next!(p)
    end
end

function permutedtrypticpalindromedistribution!(
    distribution::Matrix{Int}, 
    seq::AbstractString, 
    minlength::Int, 
    maxlength::Int
)
    trypticintervals = trypticpeptides(seq, minlength, maxlength)
    trypticlengths = interval_length.(trypticintervals)
    n = length(trypticintervals)
    lazy_lps = (longestpalindromicsubstring(shufflefast(seq[start:stop])) 
                for (start,stop)=trypticintervals)
    palindromelengths = interval_length.(lazy_lps)
    for i=1:n
        distribution[trypticlengths[i], palindromelengths[i]] += 1
    end
end

function permutedtrypticpalindromedistribution!(
    distribution::Matrix{Int}, 
    seqs::Vector, 
    minlength::Int, 
    maxlength::Int
)
    nseqs = length(seqs)
    p = Progress(nseqs, 1, "Digesting...")
    @threads for i=1:nseqs
        permutedtrypticpalindromedistribution!(distribution, seqs[i], minlength, maxlength)
        next!(p)
    end
end