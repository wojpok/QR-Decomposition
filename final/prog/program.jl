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
    for i = 1:1:n-1 
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


# Losowa macierz
function RandomMatrix(size, rad)
    R = Matrix{Float64}(undef, size, 0)
    for i in 1:size 
        R = hcat(R, rand(size) * rad .- (rad/2))
        #R = hcat(R, [(j == i) ? Float64(1.0) : (abs(i - j) > 0) ? (i - j)*lambda : Float64(0.0) for j in 1:size]) 
    end
    return R
end 

#add random noise
function addRandomNoise(M, amount=0.01)
    noise = rand(Float64, size(M))
    noise = noise .* 2 .- 1 
    noise = noise .* amount
    #uniform distribution from -amount to amount
    return M+noise
end

function matNorm(M)
    sum = 0
    for i in M; sum += abs(i); end
    return sum/(size(M, 1)*size(M, 2))
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
    
    return matNorm(M - q*r - noisy + qn*rn)

end

# =================================================
#              Test stabilności numerycznej
# =================================================

function StabilityRace(len, startRange, noise, name)
    
    t = RandomMatrix(10, startRange)
    org = t
    perturb = addRandomNoise(t, noise)
    
    xs, ys, zs = [], [], []
    
    for i in 1:len
        push!(xs, checkChange(t, GramSchmidt; verbose = false, noise = perturb))
        push!(ys, checkChange(t, MGS;         verbose = false, noise = perturb))
        push!(zs, checkChange(t, Householder; verbose = false, noise = perturb))
        
        t = perturb
        perturb = addRandomNoise(t, noise)
    end
    
    range = collect(1:len)
    
    x = plot(range, log10.(abs.(xs)), title="Stabilność numeryczna - rzędy błędów", xlabel="numer iteracji", label ="Gram-Schmidt")
    plot!(range, log10.(abs.(ys)), label ="MGS")
    plot!(range, log10.(abs.(zs)), label ="Householder")

    savefig(x, name)
end

# =================================================
#              Test wielkości macierzy
# =================================================

function SizeRace(rad, it, tests, name)
    xs, ys, zs = [], [], []
    
    for i in 3:it
        push!(xs, 0)
        push!(ys, 0)
        push!(zs, 0)
        
        for j in 1:tests
            t = RandomMatrix(i, rad)
            
            Q, R = GramSchmidt(t)
            xs[i-2] += matNorm(t - Q*R)
            
            Q, R = MGS(t)
            ys[i-2] += matNorm(t - Q*R)
            
            Q, R = Householder(t)
            zs[i-2] += matNorm(t - Q*R)
        end
        
        #rad *= 2;
    end
    
    range = collect(3:it)
    
    p = plot(range,  log10.(xs ./ tests), title="Dokładność rozkładu losowych macierzy \\n- średnie rzędy błędów", xlabel="wielkość macierzy", label ="Gram-Schmidt")
    plot!(range, log10.(ys ./ tests), label ="MGS")
    plot!(range, log10.(zs ./ tests), label ="Householder")

    savefig(p, name)
end

function solveLinearSystem(M, ys, decomp=GramSchmidt)
    q, r = decomp(M)
    
    ys = transpose(q)*transpose(ys)
    
    solutions = zeros(Float64, size(ys))
    n = size(ys)[1]

    for i in n:-1:1
        sub = 0
        for j in n:-1:i+1
            sub += solutions[j]*r[i,j]
        end
        solutions[i]=(ys[i]-sub)/r[i,i]
    end
    
    #x = r\ys for debug purposes only
    return solutions
end

function GenerateRandomLinearSystem(siz, rad=10^8)
    coeffs = RandomMatrix(siz, rad)
    xs = rand(Float64, siz) .*rad .-(rad/2)
    ys = zeros(Float64, siz)
    for i in 1:siz
        for j in 1:siz
            ys[i] = ys[i]+coeffs[i,j]*xs[j]
        end
    end
    return coeffs, xs, transpose(ys)
end
function CompareSystemSolution(decomp, M, xs, ys, verbose=false)
    solution = solveLinearSystem(M, ys, decomp)
    if verbose==true
        println(xs-solution)
        println(sum(xs-solution)/size(xs)[1])
    end
    return sum(xs-solution)/size(xs)[1]
end

# =================================================
#         Test rozwiązań układów liniowych
# =================================================

function RandomRadiusSolutionRace(rad, it ,tests, name)
    xs, ys, zs = [], [], []
    
    for i in 1:it
        push!(xs, 0)
        push!(ys, 0)
        push!(zs, 0)
        
        for j in 1:tests
            M, a, b = GenerateRandomLinearSystem(10, rad)
            
            
            xs[i] += abs(CompareSystemSolution(GramSchmidt, M, a, b))#abs(sum(t - Q*R))
            
            ys[i] += abs(CompareSystemSolution(MGS, M, a, b))#abs(sum(t - Q*R))
            

            zs[i] += abs(CompareSystemSolution(Householder, M, a, b))#abs(sum(t - Q*R))
        end
        
        rad *= 2;
    end
    
    range = collect(1:it)
    
    p = plot(range,  log10.(xs ./ tests), label="Gram-Schmidt", title="Poprawność rozwiązań układów liniowych \\n- średnie rzędy błędów", xlabel="wykładnik promienia losowości")
    plot!(range, log10.(ys ./ tests), label="MGS")
    plot!(range, log10.(zs ./ tests), label="Householder")
    
    savefig(p, name)
end

# =================================================
#              Polecenia generowania
# =================================================

    
for i in 1:10
    StabilityRace(20, 100, 0.01, "wykresy/stability"*string(i)*".png")
end

for i in 1:10
    SizeRace(5, 30, 20, "wykresy/random"*string(i)*".png")
end

for i in 1:10
    RandomRadiusSolutionRace(0.001, 10, 100, "wykresy/linear"*string(i)*".png")
end
