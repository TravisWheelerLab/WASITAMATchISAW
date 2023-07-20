using Base.Iterators: Stateful, popfirst!
using Base.Threads: @threads
using BioAlignments: LocalAlignment, AffineGapScoreModel, BLOSUM62, pairalign
using ProgressMeter: Progress, next!

function paired_alignment(
    pairs::Base.Iterators.Zip,
    n::Int,
    model,
    schema,
    formatter,
    verbose,
)
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
    pairgenerator = Stateful(pairs)
    Threads.@threads for i=1:n
        (query, reference) = popfirst!(pairgenerator)
        res = pairalign(model, query, reference, schema)
        results[i] = formatter(res)
        if verbose
            next!(p)
        end
    end
    results
end