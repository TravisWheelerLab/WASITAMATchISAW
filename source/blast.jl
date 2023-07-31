module BLAST

    execute(command::String) = read(Cmd(String.(split(command, ' '))), String)

    function blastp(
        pathtoqueryseq::String,
        pathtotargetdb::String;
        pathtomaskdb::String="",
        numthreads::Int=1,
    )
        command = "blastp -query $(pathtoqueryseq) -db $(pathtotargetdb) -word_size 4 -outfmt \"6 evalue score length qstart qend sstart send pctid\" -num_threads $(numthreads)"
        if length(pathtomaskdb) > 0
            command *= "-commanddb_hard_mask 100"
        end
        execute(command)
    end

    # TODO
    function blastn()
        @assert false "not yet implemented"
    end

    function makeblastdb(
        pathtotargetseq::String,
        pathtoresult::String;
        pathtomaskdb::String="",
    )
        command = "makeblastdb -in $(pathtotargetseq) -dbtype prot -out $(pathtoresult) -title $(pathtoresult)"
        if length(pathtomaskdb) > 0
            command *= " -mask_data $(pathtomaskdb)"
        end
        execute(command)
    end

    function convert2blastmask(
        pathtotargetmask::AbstractString,
        pathtoresult::AbstractString,
    )
        pathtomaskdb = pathtoresult * ".msk"
        command = "convert2blastmask -in $(pathtotargetmask) -out $(pathtomaskdb) -parse_seqids -masking_algorithm tantan -masking_options \"-p\" "
        execute(command)
    end

    function cleanup_makeblastdb(
        pathtoresult::AbstractString,
    )
        commands = [
            "rm $(pathtoresult * ".pto")",
            "rm $(pathtoresult * ".ptf")",
            "rm $(pathtoresult * ".psq")",
            "rm $(pathtoresult * ".pog")",
            "rm $(pathtoresult * ".pot")",
            "rm $(pathtoresult * ".pos")",
            "rm $(pathtoresult * ".pin")",
            "rm $(pathtoresult * ".phr")",
            "rm $(pathtoresult * ".pdb")"]
        for command=commands
            try
                execute(command)
            catch
                ; # throws an exception when one of the files it tries to rm doesn't exist
            end
        end
    end
    
    function cleanup_convert2blastmask(
        pathtoresult::AbstractString,
    )
        commands = [
            "rm $(pathtoresult * ".msk")",
            "rm $(pathtoresult * ".paa")",
            "rm $(pathtoresult * ".pab")",
            "rm $(pathtoresult * ".pac")"]
        for command=commands
            try
                execute(command)
            catch
                ;
            end
        end
    end
end