using FASTX: FASTAReader, FASTARecord, FASTAWriter
using DataFrames: DataFrame, dropmissing
using CSV
using DelimitedFiles: readdlm, writedlm

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

function writeframe(
    pathtoframe::AbstractString,
    frame::DataFrame,
)
    CSV.write(open(pathtoframe, "w"), frame)
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


"load chromosome record from a .fasta."
function loadchromosome(pathtorecord)
    record = readfasta(pathtorecord, singleton=true)
    metadata, sequence = split(string(record), '\n')
    gnav = split(metadata, ' ')[1][2:end]
    NCBINucleotideChromosome(record, sequence, gnav, metadata)
 end

"load a chromosome annotation file in the NCBI tabular format."
loadannotation(pathtoannotation) = NCBIGeneAnnotation(dropmissing(DataFrame(File(pathtoannotation))[:, ["GeneID", STARTPOS, ENDPOS, GNAV]]))