using FASTX: FASTARecord, description
import Base.\

abstract type Genome end

struct FASTXGenome <: Genome
    sequences::Vector{FASTARecord}
end

chromosome(n::Int, genome::FASTXGenome) = filter(
    record -> contains(description(record), "chr$(n)"), 
    genome.sequences)

abstract type Annotation end

struct GTFAnnotation <: Annotation
    table::DataFrame
end

\(annotation::GTFAnnotation, rowmask::BitVector) = GTFAnnotation(annotation.table[rowmask, :])

chromosome(n::Int, annotation::GTFAnnotation) = annotation \ (annotation.table.chromosome .== "chr$(n)")

exon(annotation::GTFAnnotation) = annotation \ (annotation.table.type .== "exon")        

intron(annotation::GTFAnnotation) = annotation \ (annotation.table.type .== "3UTR" .|| annotation.table.type == "5UTR")

intragenic(annotation::GTFAnnotation) = annotation \ (annotation.table.type .== "transcript")

const GTFCOLUMNS = ["chromosome","version","type","start","stop","dot","parity","zero","metadata"]
