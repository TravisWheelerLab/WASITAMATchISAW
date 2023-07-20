function readfasta(
    path::AbstractString;
    singleton=false,
)   
    records = collect(FASTAReader(open(path_to_record)))
    if singleton
        return records[1]
    else
        return records
    end
end