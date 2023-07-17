using DelimitedFiles: readdlm
using StatsBase: sample, Weights

function parse_score_matrix_distribution(
    path_to_out::AbstractString;
    pattern::Regex = r":\n([.\n\w\d\s\-\t\r\~\(\)\=]*)\n\n"
)
    out = read(path_to_out, String)
    patterns = eachmatch(pattern, out)
    loader = x -> readdlm(IOBuffer(x.match))
    matrices = map(loader, patterns)

    alphabet = only.(matrices[1][2, :])
    k = length(alphabet)
    alpha_enum = Dict(alphabet .=> 1:k)
    frequencies = convert.(AbstractFloat, matrices[1][3, :])
    scores = convert.(AbstractFloat, matrices[2][3:end, 2:end])
    joints = convert.(AbstractFloat, matrices[4][3:end, 2:end])
    marginals = convert.(AbstractFloat, matrices[5][3, :])
    conditionals = convert.(AbstractFloat, matrices[6][3:end, 2:end])

    Dict("alphabet"     => alphabet,
        "enumerator"    => alpha_enum,
        "frequencies"   => frequencies,
        "scores"        => scores,
        "joints"        => joints,
        "marginals"     => marginals,
        "conditionals"  => conditionals)
end

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

function bl90_mut(pctid::AbstractFloat)
    model = parse_score_matrix_distribution(BLOSUM90)
    conditionals = model["conditionals"]
    n = size(conditionals)[1]
    old_diag = diag(conditionals)
    old_pctid = mean(old_diag)
    adjusted_diag = (pctid / old_pctid) .* old_diag
    adjust_distributions!(conditionals, 1:n, adjusted_diag)
    model["conditionals"] = conditionals
    return model
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

function mutate(
    seq::AbstractString,
    substitution_model::Dict{String, Any};
    n_mutations::Int = -1,
#    alphabet::Vector{Char},
#    enumerator::Dict{Char, Int},
#    conditionals::Matrix{T}
) #where T <: AbstractFloat
    seq_vec = collect(seq)
    if n_mutations == -1 || n_mutations ≥ length(seq)
        mut_idxs = 1:length(seq)
    elseif n_mutations > 0
        mut_idxs = sample(1:length(seq), n_mutations, replace=false)
    else
        return seq
    end
    for i=mut_idxs
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