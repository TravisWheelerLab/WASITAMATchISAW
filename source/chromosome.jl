using CSV: File, DataFrames: DataFrame, dropmissing

struct Chromosome
    record::FASTARecord
    sequence::AbstractString
    gnav::AbstractString
    _metadata::AbstractString
end

"load chromosome record from a .fasta."
function loadchromosome(pathtorecord)
    record = readfasta(pathtorecord, singleton=true)
    metadata, sequence = split(string(record), '\n')
    gnav = split(metadata, ' ')[1][2:end]
    Chromosome(record, sequence, gnav, metadata)
 end

record(chromosome::Chromosome) = chromosome.record

sequence(chromosome::Chromosome) = chromosome.sequence

gnav(chromosome::Chromosome) = chromosome.gnav

struct Annotation
    table::DataFrame
end

STARTPOS = "start_position_on_the_genomic_accession"
ENDPOS = "end_position_on_the_genomic_accession"
GNAV = "genomic_nucleotide_accession.version"

"load a chromosome annotation file in the NCBI tabular format."
loadannotation(pathtoannotation) = Annotation(dropmissing(DataFrame(File(pathtoannotation))[:, [STARTPOS, ENDPOS, GNAV]]))

gnav(annotation::Annotation) = annotation[:, GNAV]

"extract the frames from an NCBI chromosome annotation as a vector of (startpos, endpos) tuples."
frames(annotation::Annotation) = eachrow(annotation.table[:, [STARTPOS, ENDPOS]])

segment(sequence, intervals::Vector{NTuple{2, Int}}) = [sequence[start:stop] for (start, stop) in intervals]

"extract regions in chromosome sequence using the frames from an NCBI annotation."
regions(chromosome::Chromosome, annotation::Annotation) = segment(sequence(chromosome), frames(annotation))

""
function annotate(chromosome::Chromosome, annotation::Annotation)
	errors = gnav(annotation) .!= gnav(chromosome)
	if any(errors)
	    rate = sum(errors) / length(errors)
		@warn "the chromosome accession does not match the accession version of $(rate) annotated regions."
	end
	regions(chromosome, annotation)	
end

