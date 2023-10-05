using FASTX: FASTAReader, FASTARecord, FASTAWriter, sequence
using DataFrames: DataFrame, dropmissing
using CSV
using DelimitedFiles: readdlm, writedlm

function parse_score_matrix_distribution(
    path_to_out::AbstractString;
    pattern::Regex = r":\n([.\n\w\d\s\-\t\r\~\(\)\=]*)\n\n"
)
    out = read(path_to_out, String)
    patterns = eachmatch(pattern, out)
    loader = x -> readdlm(IOBuffer(x.match))
    matrices = map(loader, patterns)

    alphabet = only.(matrices[1][2, :])
    k = length(alphabet)
    alpha_enum = Dict(alphabet .=> 1:k)
    frequencies = convert.(AbstractFloat, matrices[1][3, :])
    scores = convert.(AbstractFloat, matrices[2][3:end, 2:end])
    joints = convert.(AbstractFloat, matrices[4][3:end, 2:end])
    marginals = convert.(AbstractFloat, matrices[5][3, :])
    conditionals = convert.(AbstractFloat, matrices[6][3:end, 2:end])

    Dict("path" => path_to_out,
        "alphabet"     => alphabet,
        "enumerator"    => alpha_enum,
        "frequencies"   => frequencies,
        "scores"        => scores,
        "joints"        => joints,
        "marginals"     => marginals,
        "conditionals"  => conditionals)
end

function readsequences(
    pathtofasta::AbstractString
)
    String.(sequence.(readfasta(pathtofasta)))
end

function writesequences(
    pathtofasta::AbstractString,
    sequences::Vector;
    descriptions::Vector{String}=Vector{String}(undef, 0),
)
    n = length(sequences)
    if length(descriptions) != n
        if length(descriptions) > 0
            @warn "description vector and sequence vectors were different length; overwriting descriptions."
        end
        descriptions = string.(1:n)
    end
    writefasta(pathtofasta, FASTARecord.(descriptions, sequences))
end

function readfasta(
    pathtorecord::AbstractString;
    singleton=false,
)   
    records = collect(FASTAReader(open(pathtorecord)))
    if singleton
        return records[1]
    else
        return records
    end
end

function writefasta(
    pathtorecord::AbstractString,
    records::Vector{FASTARecord};
    singleton=false,
)
    FASTAWriter(open(pathtorecord, "w")) do writer
        for record=records
            write(writer, record)
        end
    end
end

function writefasta(
    pathtorecord::AbstractString,
    record::FASTARecord,
)
    writefasta(pathtorecord, [record])
end

function writefasta(
    pathtorecord::AbstractString,
    sequence::AbstractString;
    description="",
)
    writefasta(pathtorecord, FASTARecord(description, sequence))
end

function readtable(
    pathtotable::AbstractString,
)
    readdlm(pathtotable, '\t', Any, '\n')
end

function readtable(
    pathtotable::AbstractString,
    columns::AbstractVector,
)
    return DataFrame(readtable(pathtotable), columns)
end

function writetable(
    pathtotable::AbstractString,
    table::Matrix
)
    writedlm(pathtotable, table)
end

function writeframe(
    pathtoframe::AbstractString,
    frame::DataFrame,
)
    CSV.write(pathtoframe, frame)
end

function readframe(
    pathtoframe::AbstractString
)
    DataFrame(CSV.File(pathtoframe))
end

function readscopclass(
    pathtoscopclass::AbstractString
)
    scoptxt = readlines(open(pathtoscopclass))
    columns = split(scoptxt[6], ' ')[2:end]
    scopmatrix = permutedims(reduce(hcat, split.(scoptxt[7:end])))
    scoptable = DataFrame(scopmatrix, columns)
    scopcla = Tuple.(split.(scopmatrix[:, end], ','))
    scopcla = unzip(scopcla)
    scopcla = (x -> (y -> y[4:end]).(x)).(scopcla)
    scoptable."TP" = scopcla[1]
    scoptable."CL" = scopcla[2]
    scoptable."CF" = scopcla[3]
    scoptable."SF" = scopcla[4]
    scoptable."FA" = scopcla[5]
    return scoptable
end

readgenome(pathtogenome) = FASTXGenome(readfasta(pathtogenome))

readgtf(pathtogtf) = GTFAnnotation(readtable(pathtogtf, GTFCOLUMNS))
writegtf(pathtogtf, annotation::GTFAnnotation) = writetable(pathtogtf, Matrix(annotation.table))
