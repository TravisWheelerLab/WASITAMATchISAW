using DelimitedFiles: readdlm
using StatsBase: mean, sample, Weights
using LinearAlgebra: diag
using FASTX: FASTARecord, sequence, description
using ProgressMeter

# modifies `distribution` so that `distribution[fixed_index] = adjusted_value` 
# with all other values adjusted proportionately 
function adjust_distribution!(
    distribution::AbstractArray{T},
    fixed_index::Int,
    adjusted_value::S
) where {T, S <: Real}
    original_value = distribution[fixed_index]
    # adjusts the remaining values, maintaining relative quantity w.r.t. each other
    relative_factor = (1 - adjusted_value) / (1 - original_value)
    distribution .*= relative_factor
    distribution[fixed_index] = adjusted_value
end

function adjust_distribution(
    distribution::AbstractArray{T},
    fixed_index::Int,
    adjusted_value::S
) where {T, S <: Real}
    new_distribution = copy(distribution)
    adjust_distribution!(new_distribution, fixed_index, adjusted_value)
    return new_distribution
end

# interface to the substitution matrix as a table of distributions
function adjust_distributions!(
    distributions::AbstractMatrix{T},
    fixed_indices::AbstractVector{Int},
    adjusted_values::AbstractVector{S}
) where {T, S <: Real}
    n, m = size(distributions)
    @assert length(fixed_indices) == m
    @assert 0 <= min(fixed_indices...) <= max(fixed_indices...) <= n
    for i=1:m
        distribution = view(distributions, i, :)
        adjust_distribution!(distribution, fixed_indices[i], adjusted_values[i])
    end
end

function adjust_distributions(
    distributions::AbstractMatrix{T},
    fixed_indices::AbstractVector{Int},
    adjusted_values::AbstractVector{S}
) where {T, S <: Real}
    new_distributions = copy(distributions)
    adjust_distributions!(new_distributions, fixed_indices, adjusted_values)
    return new_distributions
end

function adjust_model!(model::Dict, pctid::AbstractFloat)
    conditionals = model["conditionals"]
    n = size(conditionals)[1]
    old_diag = diag(conditionals)
    old_pctid = mean(old_diag)
    adjusted_diag = (pctid / old_pctid) .* old_diag
    adjust_distributions!(conditionals, 1:n, adjusted_diag)
    model["conditionals"] = conditionals
end

function adjust_model(model::Dict, pctid::AbstractFloat)
    model = deepcopy(model)
    adjust_model!(model, pctid)
    model
end

# x -> x + k(f - x)
step_towards(target::Real, position::Real, rate::Real) = position + rate * (target - position)

function _halve_mutation_likelihood!(matrix::AbstractMatrix{T}) where T <: Real
    conservation_likelihood = diag(matrix)
    adjusted_conservation_likelihood = step_towards.(1, conservation_likelihood, 0.5)
    n = size(matrix)[1]
    adjust_distributions!(matrix, 1:n, adjusted_conservation_likelihood)
end

function halve_mutation_likelihood!(substitution_model::Dict{String, Any})
    conditionals = substitution_model["conditionals"]
    _halve_mutation_likelihood!(conditionals)
    substitution_model["conditionals"] = conditionals;
end

function halve_mutation_likelihood(substitution_model::Dict{String, Any})
    new_substituion_model = deepcopy(substitution_model)
    halve_mutation_likelihood!(new_substituion_model)
    new_substituion_model
end

function _mutate(
    seq::AbstractString,
    substitution_model::Dict{String, Any}
)
    seq_vec = collect(seq)
    n = length(seq_vec)
    for i=1:n
        char = seq[i]
        if char in substitution_model["enumerator"].keys
            # if the character is known, find its transition vector
            char_idx = substitution_model["enumerator"][char]
            transition_probabilities = substitution_model["conditionals"][char_idx, :]
        else
            # otherwise, use the base frequencies of the known alphabet
            transition_probabilities = substitution_model["frequencies"]
        end
        seq_vec[i] = sample(substitution_model["alphabet"], Weights(transition_probabilities))
    end
    join(seq_vec)
end

function mutate(
    seq::AbstractString,
    substitution_model::Dict{String, Any};
    pctid::AbstractFloat=-1.0,
)
    if pctid != -1.0
        substitution_model = adjust_model(substitution_model, pctid)
    end
    _mutate(seq, substitution_model)
end

function mutate(
    seqs::Vector,
    substitution_model::Dict{String, Any};
    pctid::AbstractFloat=-1.0,
    verbose=false,
    nthreads=Base.Threads.nthreads(),
)
    if pctid != -1.0
        substitution_model = adjust_model(substitution_model, pctid)
    end
    if nthreads > 1
        n = length(seqs)
        mutated_seqs = Vector{String}(undef, n)
        if verbose
            modelname = substitution_model["path"]
            p = Progress(n, 1, "Mutating with model $modelname")
        end
        Base.Threads.@threads for i=1:n
            mutated_seqs[i] = mutate(seqs[i], substitution_model)
            if verbose
                next!(p)
            end
        end
        return mutated_seqs
    else
        return _mutate.(seqs, [substitution_model])
    end
end