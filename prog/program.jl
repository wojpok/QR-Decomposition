using LinearAlgebra # Biblioteka potrzebna żeby zwracać R jako macierz górnotrójkątną

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

# ------------------------------------------------------
#               Algorytm Grama Schmidta
# ------------------------------------------------------

function GramSchmidt(mat)

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

# ------------------------------------------------------
#         Zmodyfikowany Algorytm Grama Schmidta
# ------------------------------------------------------

function MGS(mat)
    n = size(mat, 2)
    baseU = [mat[:, i] for i in 1:n]
    
    for i in 1:n
        for j in (i+1):n
            baseU[j] = baseU[j] - vecProjection(baseU[i], baseU[j])
        end
    end
    
    baseE = Matrix{Float64}(undef, n, 0)
    
    for i in 1:n
        baseE = hcat(baseE, vecNormalize(baseU[i]))
    end
    
    return baseE, transpose(baseE)*mat
    
end

# ------------------------------------------------------
#               Algorytm Householdera
# ------------------------------------------------------
function Householder(A)
    # rozmiar macierzy
    n = size(A, 2)
    
    # tworzy macierz Householdera ze wzoru
    function Hn(vec)
        vn = vec
        vn[1] -= sqrt(vectorNorm(vec, vec))
        return I - (2*vn*transpose(vn))/(transpose(vn)*vn)
    end
    
    # nakładanie macierzy H na Id, otrzymujemy H'
    function Hpn(mat)
        m = size(mat, 1)
        d = n - m
        res = convert(Matrix{Float64}, I[1:n, 1:n])
        for i in 1:m
            for j in 1:m
                res[d + i, d + j] = mat[i, j]
            end
        end
        return res
    end
    
    Q = I[1:n, 1:n]
    R = A
    # Iterowanie po rzędzach, branie podmacierzy i obliczanie kolejnych Hn
    for i = 1:1:n # jeśli dobrze rozumiem można iterować do n-1, dostaje się prawie identyczny wynik
                  # jak na wikipedii, oprócz elementu [3, 3] który ma znak przeciwny w macierzy R, ale odpowiednio wychodzi macierz Q
        Sub = R[i:n, i:n]
        Hp = Hpn(Hn(Sub[:, 1]))
        R = Hp*R
        Q = Q*Hp
    end
    
    return Q, UpperTriangular(R)
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
function AlgorithmQuality(mat, QRalgoritm; extended = false)
    println("===============================\nQR comparison")
    if extended; println("A = ",mat); end

    q, r =  QRalgoritm(mat)

    id = transpose(q)*q
    em, ea = CompareMatrix(I, id)

    if extended; println("Q = ",q); end
    if extended; println("Q^T*Q = ", id); end
    println("Average Error: ",ea,"\nMax Error: ",em)

    qr = q*r
    em, ea = CompareMatrix(mat, qr)
    if extended; println("R = ",r); end
    if extended; println("QR = ",qr); end
    println("Average Error: ",ea,"\nMax Error: ",em)
    print("\n")
end

t = transpose([12.0 6.0 -4.0; -51.0 167.0 24.0; 4.0 -68.0 -41.0])

AlgorithmQuality(t, GramSchmidt)
AlgorithmQuality(t, MGS)
AlgorithmQuality(t, Householder)


