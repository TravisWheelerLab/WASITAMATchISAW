using FASTX: FASTAReader, FASTARecord, FASTAWriter
using DataFrames: DataFrame, dropmissing
using CSV: File
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
    sequence::AbstractString;
    description="",
)
    writefasta(pathtorecord, [FASTARecord(description, sequence)])
end

function readtable(
    pathtotable::AbstractString,
)
    readdlm(pathtotable, '\t', Any, '\n')
end

function readtable(
    pathtotable::AbstractString,
    columns::AbstractVector
)
    return DataFrame(readtable(pathtotable), columns)
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