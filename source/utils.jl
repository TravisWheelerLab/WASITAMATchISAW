interval_length(interval::Tuple{Int, Int}) = 1 + interval[2] - interval[1] 

# https://stackoverflow.com/questions/36367482/unzip-an-array-of-tuples-in-julia
unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

function polysplit(
    str::AbstractString,
    delimiters::Union{AbstractVector, Tuple}
)
    placeholder = '_'
    for dlm=delimiters
        str = replace(str, dlm=>dlm*placeholder)
    end
    split(str, placeholder)[begin:end-1]
end

# https://stackoverflow.com/a/46998871
function shufflefast(s::AbstractString)
    ss = sizeof(s)
    l = length(s)

    ss == l && return String(shuffle!(copy(Vector{UInt8}(s))))

    v = Vector{Int}(l)
    i = start(s)
    for j in 1:l
        v[j] = i
        i = nextind(s, i)
    end

    p = pointer(s)
    u = Vector{UInt8}(ss)
    k = 1
    for i in randperm(l)
        for j in v[i]:(i == l ? ss : v[i+1]-1)
            u[k] = unsafe_load(p, j)
            k += 1
        end
    end
    String(u)
end

    
function namerecord_dir(
    pathtorecord::AbstractString;
    pathdlm='/',
)
    name, path... = reverse(split(pathtorecord, pathdlm))
    name = split(name, '.')[begin]
    join(path, pathdlm) * pathdlm * '.' * name * pathdlm
end

namerecord_file(pathtodir::AbstractString, i::Int) = pathtodir * "$(i).fa"

function individuate(
    pathtorecord::AbstractString;
    pathdlm='/',
)
    # load the records
    records = readfasta(pathtorecord)
    # create a hidden directory for the record files
    pathtorecord_dir = namerecord_dir(pathtorecord, pathdlm=pathdlm)
    mkdir(pathtorecord_dir)
    # iterate records by their index and write each to the directory
    for i=1:length(records)
        pathtorecord_file = namerecord(pathtorecord_dir, i)
        writefasta(pathtorecord_file, records[i])
    end
    pathtorecord_dir
end

function cleanup_individuate(
    pathtorecord::AbstractString;
    pathdlm='/',
)
    pathtorecord_dir = namerecord_dir(pathtorecord, pathdlm=pathdlm)
    rm(pathtorecord_dir, recursive=true)
end