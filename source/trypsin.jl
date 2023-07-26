trypsin_cleave_locations(q::AbstractString) = [i for i=1:length(q) if Char(q[i]) == 'K' || Char(q[i]) == 'R']
function nested_subintervals(locii::Vector{Int}, minlength::Int, maxlength::Int)
    intervals = Tuple{Int, Int}[]
    for x=1:length(locii)-1
        start = -1
        stop = -1
        for i=x:length(locii)-1
            start = locii[i]
            for j=i+1:length(locii)-1
                stop = locii[j]
                interval_length = stop - start
                if minlength < interval_length < maxlength
                    push!(intervals, (start, stop))
                elseif interval_length > maxlength
                    break
                end
            end
        end
    end
    intervals
end
trypticpeptides(q::AbstractString, minlength::Int, maxlength::Int) = nested_subintervals(trypsin_cleave_locations(q), minlength, maxlength)