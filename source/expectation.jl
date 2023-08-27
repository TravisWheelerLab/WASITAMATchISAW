using StatsBase: countmap
coincidence(x::AbstractString) = sum(frequency^2 for frequency=values(countmap(x))./length(x))
ELCS(n::Int, lambda2::Real) = 2*log(n, lambda2)
ELCS(q::String) = ELCS(length(q), coincidence(q))
ELPS(n::Int, lambda2::Real; v::Int=0) = ELCS(n,lambda2)+v
ELPS(q::String; v::Int=0) = ELPS(length(q), coincidence(q); v=v)
