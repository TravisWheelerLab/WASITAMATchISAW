using Random: shuffle, shuffle!

trypsin_cleave_locations(q::AbstractString) = [i for i=1:length(q) if Char(q[i]) == 'K' || Char(q[i]) == 'R']

# TODO - lazy implementation. maybe there's a functional expression?
function nested_subintervals(
    locii::Vector{Int}, 
    minlength::Int, 
    maxlength::Int
)
    intervals = Tuple{Int, Int}[]
    for x=1:length(locii)-1
        start = -1
        stop = -1
        for i=x:length(locii)-1
            start = locii[i]
            for j=i+1:length(locii)-1
                stop = locii[j]
                interval_length = stop - start
                if minlength < interval_length < maxlength
                    push!(intervals, (start, stop))
                elseif interval_length > maxlength
                    break
                end
            end
        end
    end
    intervals
end

trypticpeptides(q::AbstractString, minlength::Int, maxlength::Int) = nested_subintervals(trypsin_cleave_locations(q), minlength, maxlength)

function trypticpalindromedistribution!(
    distribution::Matrix{Int}, 
    seq::AbstractString, 
    minlength::Int, 
    maxlength::Int
)
    trypticintervals = trypticpeptides(seq, minlength, maxlength)
    trypticlengths = interval_length.(trypticintervals)
    ntryptics = length(trypticintervals)
    trypticsequences = view.(seq, (start:stop for (start, stop)=trypticintervals))
    lps = longestpalindromicsubstring.(trypticsequences)
    palindromelengths = interval_length.(lps)
    for i=1:ntryptics
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

function shuffledtrypticpalindromedistribution!(
    distribution::Matrix{Int}, 
    seq::AbstractString, 
    minlength::Int, 
    maxlength::Int
)
    trypticintervals = trypticpeptides(seq, minlength, maxlength)
    trypticlengths = interval_length.(trypticintervals)
    ntryptics = length(trypticintervals)
    trypticsequences = (seq[start:stop] for (start, stop)=trypticintervals)
    lps = longestpalindromicsubstring.(shufflefast.(trypticsequences))
    palindromelengths = interval_length.(lps)
    for i=1:ntryptics
        distribution[trypticlengths[i], palindromelengths[i]] += 1
    end
end

function shuffledtrypticpalindromedistribution!(
    distribution::Matrix{Int}, 
    seqs::Vector, 
    minlength::Int, 
    maxlength::Int
)
    nseqs = length(seqs)
    p = Progress(nseqs, 1, "Digesting...")
    @threads for i=1:nseqs
        shuffledtrypticpalindromedistribution!(distribution, seqs[i], minlength, maxlength)
        next!(p)
    end
end