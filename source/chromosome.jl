using DataFrames: DataFrame, dropmissing
using FASTX: FASTARecord

abstract type Chromosome end

struct NCBINucleotideChromosome <: Chromosome
    record::FASTARecord
    sequence::AbstractString
    gnav::AbstractString
    _metadata::AbstractString
end

record(chromosome::NCBINucleotideChromosome) = chromosome.record
sequence(chromosome::NCBINucleotideChromosome) = chromosome.sequence
gnav(chromosome::NCBINucleotideChromosome) = chromosome.gnav

abstract type Annotation end

struct NCBIGeneAnnotation <: Annotation 
    table::DataFrame
end

STARTPOS = "start_position_on_the_genomic_accession"
ENDPOS = "end_position_on_the_genomic_accession"
GNAV = "genomic_nucleotide_accession.version"

Base.length(annotation::NCBIGeneAnnotation) = size(annotation.table)[1]

startpos(annotation::NCBIGeneAnnotation) = annotation.table[:, STARTPOS]
endpos(annotation::NCBIGeneAnnotation) = annotation.table[:, ENDPOS]
gnav(annotation::NCBIGeneAnnotation) = annotation.table[:, GNAV]

"lazily iterates the substrings of `chromosome` induced by the intervals of `annotation`"
function annotate(
	chromosome::NCBINucleotideChromosome,
	annotation::NCBIGeneAnnotation,
)
	errors = gnav(annotation) .!= gnav(chromosome)
	if any(errors)
	    rate = sum(errors) / length(errors)
		@warn "the chromosome accession does not match the accession version of $(100 * rate)% annotated regions."
	end
	seq = sequence(chromosome)
	start = startpos(annotation)
	stop = endpos(annotation)
	nframes = length(annotation)
	(view(seq, start[i]:stop[i]) for i=1:nframes)
end