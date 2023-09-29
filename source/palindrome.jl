using StringAlgorithms: manacher

"""
find the length of the longest palindromic substring (LPS) in a 
sequence.
"""
function longestpalindromicsubstring(
    sequence::AbstractString;
    mode=:manacher,
)
    if length(sequence) < 2
        return sequence
    else
        if mode == :naive
            radii = _naïve(sequence)
        elseif mode == :manacher
            radii = _manacher(sequence)
        elseif mode == :shoupu
            radii = _shoupu(sequence)
        end
        return longestpalindromicsubstring(sequence, radii)
    end
end

function longestpalindromicsubstring(sequence::AbstractString, radii::Vector{Int})
    augmentedargmax = argmax(radii)
    maxradii = radii[augmentedargmax]
    leftbound = ceil(Int, augmentedargmax / 2) - floor(Int, maxradii / 2)
    rightbound = leftbound + maxradii - 1
    sequence[leftbound:rightbound]
end

"""
find the size of the maximal palindrome centered at each loci of the sequence. 
quadratic wrt sequence length. 
"""
function palindromeradii(
    sequence::AbstractString;
    nullcharacter='_',
)
    augmentedcharacters = [nullcharacter]
    for character in sequence
        push!(augmentedcharacters, character, nullcharacter)
    end
    augmentedsequence = join(augmentedcharacters)
    n = length(augmentedsequence)
    @assert length(augmentedsequence) == n
    radii = zeros(Int, n)
    radius = 0
    center = 1
    while center <= n
        left = center - (radius + 1)
        right = center + (radius + 1)
        while (left >= 1
            && right <= n 
            && augmentedsequence[left] == augmentedsequence[right])
            radius += 1
            left = center - (radius + 1)
            right = center + (radius + 1)
        end
        radii[center] = radius
        oldradius = radius
        oldcenter = center
        center += 1
        radius = 0
        while center <= oldcenter + oldradius
            mirroredcenter = oldcenter - (center - oldcenter)
            maxmirroredradius = oldcenter + oldradius - center
            if radii[mirroredcenter] < maxmirroredradius
                radii[center] = radii[mirroredcenter]
                center += 1
            elseif radii[mirroredcenter] > maxmirroredradius
                radii[center] = maxmirroredradius
                center += 1
            else
                radius = maxmirroredradius
                break
            end
        end
    end
    radii
end

"""
using a naiive algorithm to find the size of each maximal palindrome centered 
at each loci of the sequence, find the left and right bounds of the LPS.
"""
function _naïve(
    sequence::AbstractString,
)
    palindromeradii(sequence)
end

"""
dl.acm.org/doi/10.1145/321892.321896
https://github.com/lucifer1004/StringAlgorithms.jl
"""
function _manacher(
    sequence::AbstractString,
)
    manacher(sequence)
end

"""
TODO - implement shoupu's algorithm for finding the LPS without constructing the
~2n augmented sequence.
arxiv.org/abs/2003.08211
"""
function _shoupu(
    sequence::AbstractString,
)
    @assert false "not implemented"
end