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

chromosome(n::Int, annotation::GTFAnnotation) = annotation \ contains.(annotation.table.chromosome, "chr$(n)")

exon(annotation::GTFAnnotation) = annotation \ annotation.table.type .== "exon"

function intron(annotation::GTFAnnotation)
    return -1
    # every gene has a sequence of indices bracketed by the gene's start and stop locations and filled in the middle by the start and stop of each exon. introns are from the stop one of exon to the start of another, ie the first intron is from gene-start to intron-1-stop. we'll have to extract each of these indices, put them in order, then iterate pairs via 2-steps.
    # note - subsample! some introns are a million bp.
end

function intergenic(annotation::GTFAnnotation)
    # from gene-n-stop to gene-{n+1}-start
    # going to have to subsample bc intergenic regions are even bigger than introns.
end

const GTFCOLUMNS = ["chromosome","version","type","start","stop","dot","parity","zero","metadata"]

