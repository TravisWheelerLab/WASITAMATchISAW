using Random: shuffle!
using ProgressMeter: Progress, next!

"""
"""
function samplecomplement(
    space, 
    excluded_element::Int,
)
    n = length(space)
    i = rand(1:n)
    sampled_element = space[i]
    while sampled_element == excluded_element
        i = rand(1:n)
        sampled_element = space[i]
    end
    sampled_element
end

"""
"""
function samplecomplement(
    space, 
    excluded_element::Int, 
    weights::Vector{Int},
)
    n = length(space)
    @assert length(weights) == n
    weight_sum = sum(weights)
    sampled_element = excluded_element
    while sampled_element == excluded_element
        i = rand(1:n)
        p = weights[i] / weight_sum
        if rand() > p
            continue
        else
            sampled_element = space[i]
        end
    end
    sampled_element
end

"""
verify that the permutation ρ is disjoint wrt elements_by_cluster so that for 
all i elements_by_cluster[i] != ρ(elements_by_cluster)[i]. if more than half of
the elements belong to one cluster, no such ρ can exist.
"""
function isdisjoint(
    ρ::Vector{Int},
    elements_by_cluster::Vector{Int}
)
    n = length(elements_by_cluster)
    if length(ρ) != n
        return false
    end
    setρ = unique(ρ)
    if length(setρ) != n
        return false
    end
    if setρ == collect(1:n)
        return false
    end
    for i=1:n
        j = ρ[i]
        if elements_by_cluster[i] == elements_by_cluster[j]
            return false
        end
    end
    true
end

"""
construct a disjoint permutation ρ of elements_by_cluster. if the elements are
not uniformly distributed among the clusters, pass the kwarg flag weighted=true
to use a weighted sampling algorithm that is less likely to cause depletion. if
more than half of the elements belong to one cluster, no such ρ can exist.
"""
function disjointpermutation(
    elements_by_cluster::Vector{Int};
    weighted=false
)
    # extract cluster names and enumerate them
    nclusters = length(unique(elements_by_cluster))
    # construct a map from clusters to their elements
    clusters = [Int[] for _=1:nclusters]
    nelements = length(elements_by_cluster)
    for i=1:nelements
        a = elements_by_cluster[i]
        push!(clusters[a], i)
    end
    clusterlengths = length.(clusters)
    # shuffle to minimize sampling bias
    for a=1:nclusters
        shuffle!(clusters[a])
    end
    # construct the pairing
    ρ = fill(-1, nelements)
    while -1 in ρ
        # keep track of how many elements have been popped from each cluster.
        # TODO: replace lines 101-108 and 117-120 with logic using npops
        npops = zeros(Int, nclusters)
        for i=1:nelements
            # identify nonempty clusters and their lengths
            nonemptyclusters = filter(i -> npops[i] < clusterlengths[i], 1:nclusters)
            # check for depletion
            if length(nonemptyclusters) == 1 && elements_by_cluster[i] == nonemptyclusters[1]
                # there is only one non-empty cluster remaining, so pairing cannot proceed.
                # if one cluster contains more than 50%, pairing will always fail!
                @warn "depletion occurred after $(i) iterations"
                break
            end
            # retrieve the cluster of the element
            a = elements_by_cluster[i]
            # and select a complementary cluster
            b = -1
            if !weighted
                # if the elements are uniformly distributed into clusters, this method is faster
                # and (usually) won't cause depletion
                b = samplecomplement(1:nclusters, a)
                while npops[b] >= clusterlengths[b]
                    # loop until the sampled cluster is not depleted.
                    b = samplecomplement(1:nclusters, a)
                end
            else
                # if elements are nonuniformly distributed, weighted sampling is much more likely
                # to complete without depletion, although it can be slower.
                b = samplecomplement(1:nclusters, a, clusterlengths .- npops)
            end
            # finally, pair to a random element of the complementary cluster.
            npops[b] += 1
            j = npops[b]
            ρ[i] = clusters[b][j]
        end
    end
    ρ
end