using Base.Iterators: Stateful, popfirst!, repeated
using Base.Threads: @threads
using BioAlignments: score, PairwiseAlignmentResult, LocalAlignment, GlobalAlignment, AffineGapScoreModel, BLOSUM62, pairalign
using ProgressMeter: Progress, next!, @showprogress

function _align(
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
        p = Progress(n, 1, "Aligning...")
    end
    # do the thing
    pairgenerator = Stateful(pairs)
    Threads.@threads for i=1:n
        (query, reference) = popfirst!(pairgenerator)
        res = pairalign(model, query, reference, schema)
        results[i] = formatter(res)
        if verbose# && (i % (floor(Int, n / 100)) == 0)
            next!(p)
        end
    end
    results
end

abstract type AlignmentProductSpace end
struct Pairwise <: AlignmentProductSpace end
struct OneToMany <: AlignmentProductSpace end
struct ManyToMany <: AlignmentProductSpace end

function align(
    ::Pairwise,
    queries,
    references;
    model=LocalAlignment(),
    schema=AffineGapScoreModel(BLOSUM62, -11, -1),
    formatter=x::PairwiseAlignmentResult->x,
    verbose=false,
)
    n = length(queries)
    if length(references) != n
        throw(ArgumentError("Pairwise alignment requires there be an equal number of query and reference sequences."))
    end
    _align(zip(queries, references), n, model, schema, formatter, verbose)
end

function align(
    ::OneToMany,
    query,
    references;
    model=LocalAlignment(),
    schema=AffineGapScoreModel(BLOSUM62, -11, -1),
    formatter=x::PairwiseAlignmentResult->x,
    verbose=false,
)
    n = length(references)
    _align(zip(repeated(query), references), n, model, schema, formatter, verbose)
end

function align(
    ::ManyToMany,
    queries,
    references;
    model=LocalAlignment(),
    schema=AffineGapScoreModel(BLOSUM62, -11, -1),
    formatter=x::PairwiseAlignmentResult->x,
    verbose=false,
)
    n = length(queries)
    m = length(references)
    result_type = Union{Base.return_types(formatter)...}
    results = Vector{result_type}(undef, n*m)
    if verbose
        p = Progress(n, 1, "Aligning...")
    end
    for i=1:n
        results[1+(i-1)*m:i*m] = align(
            OneToMany(),
            queries[i],
            references,
            model=model, schema=schema, formatter=formatter,
            verbose=false)
        if verbose
            next!(p)
        end
    end
    results
end