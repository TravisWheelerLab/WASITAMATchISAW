"""
find the left and right bounds of the longest palindromic substring (LPS) in a 
sequence.
"""
function longestpalindromicsubstring(
    sequence::AbstractString;
    mode=:naive,
)
    if length(sequence) == 0
        return nothing, nothing
    elseif length(sequence) == 1
        return 1, 1
    else
        if mode == :naive
            return _naïve(sequence)
        elseif mode == :manacher
            return _manacher(sequence)
        elseif mode == :shoupu
            return _shoupu(sequence)
        end
    end
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
    radii = palindromeradii(sequence)
    augmentedargmax = argmax(radii)
    maxradii = radii[augmentedargmax]
    leftbound = ceil(Int, augmentedargmax / 2) - floor(Int, maxradii / 2)
    rightbound = leftbound + maxradii - 1
    leftbound, rightbound
end

"""
TODO - implement manacher's algorithm for finding the LPS in linear time wrt 
sequence length.
dl.acm.org/doi/10.1145/321892.321896
"""
function _manacher(
    sequence::AbstractString,
)
    @assert false "not implemented"
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