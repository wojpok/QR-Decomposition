using LinearAlgebra # Biblioteka potrzebna żeby zwracać R jako macierz górnotrójkątną

# ------------------------------------------------------
#               Algorytm Gramma Schmidta
# ------------------------------------------------------

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

    baseU = Matrix{Float64}(undef, size(mat, 1), 0)
    baseE = Matrix{Float64}(undef, size(mat, 1), 0)

    # wyliczanie bazy E oraz U
    for i in 1:size(mat, 1)
        u = mat[:, i]
        for j in 1:(i-1)
            u = u - vecProjection(baseU[:, j], mat[:, i])
        end

        baseU = hcat(baseU, u)
        baseE = hcat(baseE, vecNormalize(u))
    end
    
    # wyliczanie macierzy R
    R = Matrix{Float64}(undef, size(mat, 1), size(mat, 1))
    for i in 1:size(mat, 1)
        for j in 1:size(mat, 1)
            R[i, j] = vectorNorm(baseE[:, i], mat[:, j])
        end
    end

    baseE, UpperTriangular(R)
end

# Prosta funkcja porównująca macierze
# Zwraca parę liczb: maksymalny i średni błąd macierzy result względem macierzy dest
function CompareMatrix(dest, result)
    maxErr, sumErr = 0, 0
    matSize = size(result, 1)

    for i in 1:matSize
        for j in 1:matSize
            nErr = abs(dest[i, j] - result[i, j])
            sumErr = sumErr + nErr
            maxErr = maxErr > nErr ? maxErr : nErr
        end
    end

    maxErr, sumErr/(matSize*matSize)
end

# Funckja oceniająca jakość algorytmu rozkładu QRalgorithm macierzy mat
# Używa Compare Matrix by określić jak dobrze Q^T Q = Id i QR = mat
function AlgorithmQuality(mat, QRalgoritm)
    println("===============================\nQR comparison")
    println("A = ",mat)

    q, r =  QRalgoritm(mat)

    id = transpose(q)*q
    em, ea = CompareMatrix(I, id)

    println("Q = ",q)
    println("Q^T*Q = ", id)
    println("Average Error: ",ea,"\nMax Error: ",em)

    qr = q*r
    em, ea = CompareMatrix(mat, qr)
    println("R = ",r)
    println("QR = ",qr)
    println("Average Error: ",ea,"\nMax Error: ",em)
    print("\n")
end

AlgorithmQuality([1 0; 6 5], grammSchmidt);
AlgorithmQuality([1 0 3; 6 5 7; 5 7 8], grammSchmidt);





