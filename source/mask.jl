using FASTX: FASTARecord
hardmask(sequence::AbstractString, mask) = replace(sequence, r"[a-z]" => mask)
hardmask(record::FASTARecord; mask='X') = FASTARecord(description(record), hard_mask(sequence(record)), mask)