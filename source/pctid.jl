using BioAlignments: PairwiseAlignmentResult, alignment, cigar

OPCHRS = ('M', '=', 'X', 'I', 'D', 'N')

"""
    parsecigar(cigar::String)

Split a cigar string into vectors of values and operators.
"""
function parsecigar(
    cigar::AbstractString;
)
    operatorstrings = polysplit(cigar, OPCHRS)
    n = length(operatorstrings)
    operatornums = Vector{Int}(undef, n)
    operatorchrs = Vector{Char}(undef, n)
    for i=1:n
        opstr = operatorstrings[i]
        opnum = parse(Int, opstr[begin:end-1])
        opchr = opstr[end]
        operatornums[i] = opnum
        operatorchrs[i] = opchr
    end
    operatornums, operatorchrs
end

# these methods interface to the operator (character) vector
ismatch(c::Char) = (c == '=' || c == 'M')
ismismatch(c::Char) = c == 'X'
isinsert(c::Char) = c == 'I'
isdeletion(c::Char) = (c == 'D' || c == 'N')

"""
    percentid(operator_nums::Vector{Int}, operator_chrs::Vector{Char})

Calculate the portion of match operators to match and mismatch operators. 

Expects BioAlignment.jl's CIGAR implementation with '=' or 'M' as match and 
'X' as mismatch.
"""
function percentid(
    operatornums::Vector{<:Int},
    operatorchrs::Vector{<:Char},
)
    match = ismatch.(operatorchrs)
    nmatch = sum(match .* operatornums)
    mismatch = ismismatch.(operatorchrs)
    nmismatch = sum(mismatch .* operatornums)
    nmatch / (nmismatch + nmatch)
end

"""
    percent_id(res::PairwiseAlignmentResult)

Calculate the percent identity of a pairwise alignment.
"""
function percentid(
    result::PairwiseAlignmentResult
)
    aln = alignment(alignment(result))
    cig = cigar(aln)
    opnums, opchars = parsecigar(cig)
    percentid(opnums, opchars)
end
