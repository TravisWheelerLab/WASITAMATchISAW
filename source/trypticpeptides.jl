trypsin_cleave_locations(q) = [i for i=1:length(q) if Char(q[i]) == 'K' || Char(q[i]) == 'R']

function trypticpeptides(
    q,
    minlength::Int,
    maxlength::Int
)
    peptides = String[]
    loc = trypsin_cleave_locations(q)
    for x=1:length(loc)-1
        peptide_start = -1
        peptide_end = -1
        for i=x:length(loc)-1
            peptide_start = loc[i]
            for j=i+1:length(loc)-1
                peptide_end = loc[j]
                dif = peptide_end - peptide_start
                # if dif < minlen continue else
                if minlength < dif < maxlength
                    pep = String(q[peptide_start:peptide_end])
                    push!(peptides, pep)
                elseif dif > maxlength
                    break
                end
            end
        end
    end
    peptides
end