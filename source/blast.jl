using Distributed

execute(command::String) = read(Cmd(String.(split(command, ' '))), String)

const BLASTP_OUTFMT = "6 score evalue qseqid qstart qend sseqid sstart send"
const BLASTP_COLUMNS = split(BLASTP_OUTFMT, ' ')[2:end]

function blastp(
    pathtoqueryseq::String,
    pathtotargetdb::String;
    pathtomaskdb::String="",
    numthreads::Int=1,
)   
    if length(pathtomaskdb) == 0
        command = `blastp -query $pathtoqueryseq -db $pathtotargetdb -word_size 4 -outfmt $BLASTP_OUTFMT -num_threads $numthreads`
    else
        command = `blastp -query $pathtoqueryseq -db $pathtotargetdb -word_size 4 -outfmt $BLASTP_OUTFMT -num_threads $numthreads -db_hard_mask 100`
    end
    read(command, String)
end

function parse_blastp(
    result::AbstractString;
)
    if result == ""
        return zeros((0, length(BLASTP_COLUMNS)))
    else
        try
            return readdlm(IOBuffer(result))
        catch e
            print(e)
            print("RESULT ", result)
        end
    end
    
end

function parse_blastp(
    results::AbstractVector
)
    reduce(vcat, parse_blastp.(results))
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
#        "rm $(pathtoresult * ".pog")",
        "rm $(pathtoresult * ".pot")",
#        "rm $(pathtoresult * ".pos")",
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

function search(
    ::Pairwise,
    pathtoquery, 
    pathtoreference;
    pathdlm='/',
    verbose=false,
    careful=false,
    ntasks=1,
)
    if careful
        cleanup_individuate(pathtoquery)
        cleanup_individuate(pathtoreference)
    end
    query_dir = individuate(pathtoquery, pathdlm=pathdlm)
    reference_dir = individuate(pathtoreference, pathdlm=pathdlm)
    n = length(readdir(query_dir))
    result = Vector{String}(undef, n)
    p = Progress(n, 1)
    if verbose
        println("pairwise BLAST on $n sequence pairs\nquery $pathtoquery\nreference $pathtoreference")
    end
    nchunk = ceil(Int, n / ntasks)
    # thanks, tkf.
    # https://discourse.julialang.org/t/how-to-launch-several-run-cmd-in-parallel/68862/3
    @sync for tasknum=1:ntasks
        @async try
            for k=1:nchunk
                i = ((k-1)*ntasks)+tasknum
                if i > n
                    if verbose
                        println("task $tasknum stalled")
                    end
                    continue
                end
                query_file = namerecord_file(query_dir, i)
                reference_file = namerecord_file(reference_dir, i)
                reference_db = reference_dir * string(i)
                makeblastdb(reference_file, reference_db)
                result[i] = blastp(query_file, reference_db; numthreads=1)
                cleanup_makeblastdb(reference_db)
                if verbose
                    next!(p)
                end
            end
        catch e
            println(e)
        end
    end
    cleanup_individuate(pathtoquery, pathdlm=pathdlm)
    cleanup_individuate(pathtoreference, pathdlm=pathdlm)
    result
end