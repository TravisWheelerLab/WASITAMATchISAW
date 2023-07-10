using .Threads
using BioAlignments
using ProgressMeter

function pairalign_multithread(
    queries::Vector,
    references::Vector;
    verbose=true,
    mode=LocalAlignment(),
    cost=AffineGapScoreModel(BLOSUM62, -11, -1),
    formatter::Function=x::PairwiseAlignmentResult->x,
)
    n = length(queries)
    # sanity check
    @assert length(references) == n
    # what type does the formatter return? 
    ## `Base.return_types` returns a list of types. we unravel `...` this vector and recollect its elements as a `Union` of types.
    result_type = Union{Base.return_types(formatter)...}
    # create a place for each result
    results = Vector{result_type}(undef, n)
    # create a progress bar
    if verbose
        p = Progress(n, 0.1, "Aligning...")
    end
    # do the thing
    Threads.@threads for i=1:n
        res = BioAlignments.pairalign(mode, queries[i], references[i], cost)
        results[i] = formatter(res)
        if verbose
            next!(p)
        end
    end
    results
end
