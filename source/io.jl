using FASTX: FASTAReader
using DataFrames: DataFrame, dropmissing
using CSV: File

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

function readtable(
    pathtotable::AbstractString,
)
    readdlm(pathtotable, '\t', Any, '\n')
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