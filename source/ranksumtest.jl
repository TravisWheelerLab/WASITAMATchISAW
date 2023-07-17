using HypothesisTests: MannWhitneyUTest, pvalue

ranksumtest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) = pvalue(MannWhitneyUTest(x, y))
