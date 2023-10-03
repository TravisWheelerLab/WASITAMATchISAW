using Random: shuffle, shuffle!
using Base.Threads: @threads
using StringAlgorithms: longestcommonsubstring

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

function trypticLCSdistribution!(
    distribution::Matrix{Int},
    seq::AbstractString,
    minlength::Int,
    maxlength::Int
)
    trypticintervals = trypticpeptides(seq, minlength, maxlength)
    trypticlengths = interval_length.(trypticintervals)
    n = length(trypticintervals)
    lazy_peptides = (seq[start:stop] for (start,stop)=trypticintervals)
    lazy_shuffles = (shufflefast(seq[start:stop]) for (start,stop)=trypticintervals)
    lcslengths = longestcommonsubstring.(lazypeptides, lazyshuffles)#align(Pairwise(), 
        #lazy_peptides, 
        #lazy_shuffles; 
        #formatter=x::PairwiseAlignmentResult->Int(score(x)),
        #schema=SUBSTRING_ALIGNMENT_SCHEMA)
    for i=1:n
        # the distribution matrices are square with dim maxlength + 1
        # because we have to record zero-length LCS, we add one to the lcs length.
        # this has to accounted for later in the plotting function. 
        distribution[trypticlengths[i], lcslengths[i] + 1] += 1
    end
end

function trypticLCSdistribution!(
    distribution::Matrix{Int}, 
    seqs::Vector, 
    minlength::Int, 
    maxlength::Int;
    verbose=false
)
    nseqs = length(seqs)
    p = Progress(nseqs, 1, "Digesting (LCS)...")
    @threads for i=1:nseqs
        trypticLCSdistribution!(distribution, seqs[i], minlength, maxlength)
        if verbose 
            next!(p)
        end
    end
end

function trypticLCS2distribution!(
    distribution::Matrix{Int},
    seq::AbstractString,
    minlength::Int,
    maxlength::Int
)
    trypticintervals = trypticpeptides(seq, minlength, maxlength)
    trypticintervals2 = shuffle(trypticintervals)
    trypticlengths = interval_length.(trypticintervals)
    n = length(trypticintervals)
    lazy_peptides = (seq[start:stop] for (start,stop)=trypticintervals)
    lazy_peptides2 = (shufflefast(seq[start:stop]) for (start,stop)=trypticintervals2)
    lcslengths = longestcommonsubstring.(lazypeptides, lazy_peptides2)#align(Pairwise(), 
        #lazy_peptides, 
        #lazy_peptides2; 
        #formatter=x::PairwiseAlignmentResult->Int(score(x)),
        #schema=SUBSTRING_ALIGNMENT_SCHEMA)
    for i=1:n
        distribution[trypticlengths[i], lcslengths[i] + 1] += 1
    end
end

function trypticLCS2distribution!(
    distribution::Matrix{Int}, 
    seqs::Vector, 
    minlength::Int, 
    maxlength::Int;
    verbose=false
)
    nseqs = length(seqs)
    p = Progress(nseqs, 1, "Digesting (LCS2)...")
    @threads for i=1:nseqs
        trypticLCS2distribution!(distribution, seqs[i], minlength, maxlength)
        if verbose 
            next!(p)
        end
    end
end

function trypticLPSdistribution!(
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
    lpslengths = length.(lazy_lps)
    for i=1:n
        distribution[trypticlengths[i], lpslengths[i]] += 1
    end
end

function trypticLPSdistribution!(
    distribution::Matrix{Int}, 
    seqs::Vector, 
    minlength::Int, 
    maxlength::Int;
    verbose=false
)
    nseqs = length(seqs)
    p = Progress(nseqs, 1, "Digesting (LPS)...")
    @threads for i=1:nseqs
        trypticLPSdistribution!(distribution, seqs[i], minlength, maxlength)
        if verbose 
            next!(p)
        end
    end
end
