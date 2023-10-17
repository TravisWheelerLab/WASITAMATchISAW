using Random: shuffle, shuffle!
using Base.Threads: @threads
using StringAlgorithms

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

function trypticpeptide_interval(
    seq::AbstractString,
    minlength::Int,
    maxlength::Int
)
    locii = trypsin_cleave_locations(seq)
    proline = isproline.(collect(seq))
    cleavable_locii = filter(i -> i==1 ? true : !proline[i-1], locii)
    rightbounded_halfopen_subintervals(cleavable_locii, minlength, maxlength)
end

function trypticpeptide(
    seq::AbstractString,
    minlength::Int,
    maxlength::Int,
)
    intervals = trypticpeptide_interval(seq, minlength, maxlength)
    (seq[start:stop] for (start,stop)=intervals)
end
