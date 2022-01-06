using LinearAlgebra # Biblioteka potrzebna żeby zwracać R jako macierz górnotrójkątną
using Plots

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
    for i = 1:1:n-1 # jeśli dobrze rozumiem można iterować do n-1, dostaje się prawie identyczny wynik
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

# Losowa macierz
function RandomMatrix(size, rad)
    R = Matrix{Float64}(undef, size, 0)
    for i in 1:size 
        R = hcat(R, rand(size) * rad .- (rad/2))
        #R = hcat(R, [(j == i) ? Float64(1.0) : (abs(i - j) > 0) ? (i - j)*lambda : Float64(0.0) for j in 1:size]) 
    end
    return R
end 
R = (RandomMatrix(30, 10^10))




#add random noise
function addRandomNoise(M, amount=0.01)
    noise = rand(Float64, size(M))
    noise = noise .* 2 .- 1 
    noise = noise .* amount
    #uniform distribution from -amount to amount
    return M+noise
end

function checkChange(M, decomp; verbose = true, noise = true)
    noisy = noise
    if(noisy === true); noisy = addRandomNoise(M); end
    
    if(verbose); println("Total change in arguments ", sum(noisy-M)); end
    q, r = decomp(M)
    qn, rn = decomp(noisy)
    siz = size(M)[1]*size(M)[2]
    
    
    averageChange = (sum(M-q*r)-sum(noisy-qn*rn))/siz
    
    if(verbose)
        println("Total change in Q ", sum(qn-q))
        println("Average change in Q ", sum(qn-q)/siz)
        println("Total change in R ", sum(rn-r))
        println("Average change in R ", sum(rn-r)/siz)
        println("|ΔQ|/|Δarg| ", sum(qn-q)/sum(noisy-M))
        println("|Δq|/|Δarg| ", sum(qn-q)/(siz*sum(noisy-M)))
        println("|ΔR|/|Δarg| ", sum(rn-r)/sum(noisy-M))
        println("|Δr|/|Δarg| ", sum(rn-r)/(siz*sum(noisy-M)))
        println("Total change in M-Q*R ", sum(M-q*r)-sum(noisy-qn*rn))
        println("Average change in M-Q*R ", averageChange)
    end
    
    return averageChange

end

# =================================================
#           Wykresy
# =================================================

function StabilityRace(len, startRange, noise, name)
    
    t = RandomMatrix(10, startRange)
    org = t
    perturb = addRandomNoise(t, noise)
    
    xs, ys, zs = [], [], []
    
    for i in 1:len
        push!(xs, checkChange(t, GramSchmidt; verbose = false, noise = org))
        push!(ys, checkChange(t, MGS;         verbose = false, noise = org))
        push!(zs, checkChange(t, Householder; verbose = false, noise = org))
        
        t = perturb
        t = addRandomNoise(t, noise)
    end
    
    range = collect(1:len)
    
    p = plot(range, abs.(xs), title="Stabilność numeryczne - średnie błędy", xlabel="numer iteracji", label ="Gram-Schmidt")
    plot!(range, abs.(ys), label ="MGS")
    plot!(range, abs.(zs), label ="Householder")

    savefig(p, name)

end

function RandomRadiusRace(rad, it, tests, name)
    xs, ys, zs = [], [], []
    
    for i in 1:it
        push!(xs, 0)
        push!(ys, 0)
        push!(zs, 0)
        
        for j in 1:tests
            t = RandomMatrix(10, rad)
            
            Q, R = GramSchmidt(t)
            xs[i] += abs(sum(t - Q*R))
            
            Q, R = MGS(t)
            ys[i] += abs(sum(t - Q*R))
            
            Q, R = Householder(t)
            zs[i] += abs(sum(t - Q*R))
        end
        
        rad *= 2;
    end
    
    range = collect(1:it)
    
    p = plot(range,  log10.(xs ./ tests), title="Dokładność rozkładu losowych macierzy \\n- średnie rzędy błędów", xlabel="promień losowości", label ="Gram-Schmidt")
    plot!(range, log10.(ys ./ tests), label ="MGS")
    plot!(range, log10.(zs ./ tests), label ="Householder")

    savefig(p, name)
end
    
for i in 1:10
    StabilityRace(20, 100, 0.1, "stability"*string(i)*".png")
end

for i in 1:3
    RandomRadiusRace(5, 20, 20, "random"*string(i)*".png")
end