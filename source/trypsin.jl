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
    intervals = Tuple{Int, Int}[]
    for x=1:length(locii)-1
        start = -1
        stop = -1
        for i=x:length(locii)-1
            # the left bound is closed 
            start = locii[i]
            for j=i+1:length(locii)-1
                # the right bound is open
                stop = locii[j] - 1
                irvlength = interval_length((start, stop))
                if irvlength > maxlength
                    break
                elseif minlength <= irvlength # <= maxlength 
                    push!(intervals, (start, stop))
                end
            end
        end
    end
    intervals
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
    palindromelengths = interval_length.(
                            longestpalindromicsubstring.(
                                view.(seq, (x->x[1]:x[2]).(trypticintervals))))
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
    palindromelengths = interval_length.(
                            longestpalindromicsubstring.(
                                shufflefast.(
                                    view.(seq, (x->x[1]:x[2]).(trypticintervals)))))
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