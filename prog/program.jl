

function grammSchmidt(mat)
    # Funkcje pomocnicze
    function vectorNorm(vec1, vec2)
        return transpose(vec1)*(vec2)
    end
    
    function vecProjection(u, vec)
        return vectorNorm(u, vec)/vectorNorm(u,u) * u
    end
    
    function vecNormalize(vec)
        return sqrt(1 / vectorNorm(vec, vec)) * vec
    end

    baseU = []
    baseE = []

    # wyliczanie bazy E oraz U
    for i in 1:size(mat, 1)
        u = mat[:, i]
        for j in 1:(i-1)
            u = u - vecProjection(baseU[j], mat[:, i])
        end
        push!(baseU, u)
        push!(baseE, vecNormalize(u))
    end

    # TODO dokończenie wyznaczania macierzy z bazy
    baseE, baseU
end

print(grammSchmidt([1.0 2.0; 3.0 4.0]))
