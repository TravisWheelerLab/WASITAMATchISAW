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
