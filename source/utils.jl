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