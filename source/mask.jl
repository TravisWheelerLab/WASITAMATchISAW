using FASTX: FASTARecord, description, sequence
hardmask(sequence; mask='X') = String(replace(sequence, r"[a-z]" => mask))
hardmask(record::FASTARecord; mask='X') = FASTARecord(description(record), hardmask(sequence(record)); mask=mask)