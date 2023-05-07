module DN1

export RazprsenaMatrika
export sor

import Base.getindex
import Base.setindex!
import Base.firstindex
import Base.lastindex
import Base.*

using LinearAlgebra

"""
    RazprsenaMatrika(V, I)
    
Podatkovni tip za razpršeno matriko. V matriki 'V' se nahajajo neničelni elementi originalne matrike. V matriki 'I' pa so 
shranjeni indeksi stolpcev teh neničelnih elementov. Vsaka vrstica matrike 'V' vsebuje neničelne elemente iz iste vrstice 
kot originalna matrika. Če je originalna matrika dimenzije n x n, potem sta dimenziji 'V' in 'I' enaka n x m.
"""

mutable struct RazprsenaMatrika
    V
    I::Matrix{Int64}
end

"""
    RazprsenaMatrika(matrika)

Funkcija, ki vhodno matriko transformira v tip RazprsenaMatrika. 
Vhod je diagonalno dominantna razpršena matrika, izhod pa matriki 'V' in 'I'.
Primer:
    vhod: A =|3 0 -2 1 0|
          |0 3 1 2 0|
          |0 0 0 0 0|
          |-1 0 0 2 0|
          |0 4 0 0 7|
    izhod:
        V: |3 -2 1|
           |3 1 2|
           |0 0 0|
           |-1 2 0|
           |4 7 0|
        
        I: |1 3 4|
           |2 3 4|
           |0 0 0|
           |1 4 0|
           |2 5 0|
"""

function RazprsenaMatrika(matrika)
    max = 1
    for i=1:size(matrika, 1)
        st = 0
        for j=1:length(matrika[i, :])
            if matrika[i, j] != 0
                st = st + 1
            end
        end
        if st > max
            max = st
        end
    end

    stolpci = max
  
    V = zeros(size(matrika, 1), stolpci)
    I = zeros(size(matrika, 1), stolpci)
    br = 1
    for i=1:size(matrika, 1)
        count = 1
        flag = 0
        for j=1:length(matrika[i, :])
            if matrika[i, j] != 0
                V[br, count] = matrika[i, j]
                I[br, count] = j
                count = count + 1
                flag = 1
            end
        end
        br = br + 1
    end       

    return RazprsenaMatrika(V, I)
end


"""
    getindex(A::RazprsenaMatrika, i::Int, j::Int)

Funkcija, ki vrne element matrike na določenom indeksu. 
Vhod so: matrika tipa RazpršenaMatrika, indeks vrstice in indeks stolpca. Izhod pa je element na določenoj poziciji.
Primer:
    vhod: getindex(RazprsenaMatrika(A), 2 , 4)
    izhod: 2
"""

function getindex(A::RazprsenaMatrika, i::Int64, j::Int64)
    if i > size(A.V, 1) || j > size(A.V, 1) || i < 1 || j < 1
        throw(BoundsError(A, (i, j)))
    end

    for k=1:size(A.I, 2)
        if A.I[i, k] == j
            return A.V[i, k]
        end
    end
    return 0
end
 
"""
    setindex!(A::RazprsenaMatrika, vrednost, i::Int, j::Int)

Funkcija, ki shrani 'vrednost' na določenoj poziciji v matriki. 
Vhod so: matrika tipa RazpršenaMatrika, vrednost (ki jo želimo shraniti), indeks vrstice in indeks stolpca. 

Primer: setindex!(RazprsenaMatrika(A), 3, 2, 1) - element nove matrike v drugoj vrstici in prvem stolpcu ima vrednost 3.
"""

function setindex!(A::RazprsenaMatrika, vrednost, i::Int, j::Int)
    if i > size(A.V, 1) || j > size(A.V, 1) || i < 1 || j < 1
        throw(BoundsError(A, (i, j)))
    end

    #element na tom indeksu ima neku vrijednost
    for k=1:size(A.I, 2)
        if A.I[i, k] == j
            A.V[i, k] = vrednost
            return A
        end
    end

    #indeks nije do sad bio v V i I, na prvom mjestu
    if j < A.I[i, 1]
        tempV = zeros(size(A.V, 1), size(A.V, 2) + 1)
        tempI = zeros(size(A.V, 1), size(A.V, 2) + 1)
        for k = 1:size(tempV, 1)
            for n = 1: size(tempV, 2)
                if k == i && n == 1
                    tempV[k, n] = vrednost
                    tempI[k, n] = j
                elseif k == i
                    tempV[k, n] = A.V[k, n-1]
                    tempI[k, n] = A.I[k, n-1]
                elseif k != i && n == size(tempV, 2)
                    tempV[k, n] = 0
                    tempI[k, n] = 0
                else
                    tempV[k, n] = A.V[k, n]
                    tempI[k, n] = A.I[k, n]
                end
            end
        end
        A.V = tempV
        A.I = tempI
        return A
    end

    #indeks nije do sad bio v V i I, na zadnjem mjestu
    if j > A.I[i, size(A.I, 2)] && A.I[i, size(A.I, 2)] != 0
        tempV = zeros(size(A.V, 1), size(A.V, 2) + 1)
        tempI = zeros(size(A.V, 1), size(A.V, 2) + 1)
        for k = 1:size(tempV, 1)
            for n = 1: size(tempV, 2)
                if k == i && n == size(tempV, 2)
                    tempV[k, n] = vrednost
                    tempI[k, n] = j
                elseif k != i && n == size(tempV, 2)
                    tempV[k, n] = 0
                    tempI[k, n] = 0
                else
                    tempV[k, n] = A.V[k, n]
                    tempI[k, n] = A.I[n, n]
                end
            end
        end
        A.V = tempV
        A.I = tempI
        return A
    end 

    #indeks nije do sad bio v V i I, ima nula u toj vrstici
    if A.I[i, 1] == 0
        A.I[i, 1] = j
        A.V[i, 1] = vrednost
    else
        for k=1:size(A.I, 2)
            if A.I[i, k] > j
                A.I[i, k+1] = A.I[i, k]
                A.I[i, k] = j
                A.V[i, k+1] = A.V[i, k]
                A.V[i, k] = vrednost
                break
            end
        end
    end
    return A

end


"""
    firstindex(A::RazprsenaMatrika)

Funkcija vrne prvi indeks matrike, ki je tipa RazpršenaMatrika.  

Primer: vhod: firstindex(RazprsenaMatrika(A)), izhod (1, 1).
"""

function firstindex(A::RazprsenaMatrika)
    return 1, 1
end


"""
    lastindex(A::RazprsenaMatrika)

Funkcija vrne zadnji indeks matrike, ki je tipa RazpršenaMatrika.  

Primer: vhod: lastindex(RazprsenaMatrika(A)), izhod (5, 5).
"""

function lastindex(A::RazprsenaMatrika)
    return size(A.V, 1), size(A.V, 1)
end


"""
    *(A::RazprsenaMatrika, x::Vector)

Funkcija, ki pomnoži matriko tipa RazprsenaMatrika z vektorjem.
Primer: |1|
    v = |2|
        |3|
        |4|
        |5|

    vhod: *(RazprsenaMatrika(A), v)
    izhod vekor: |1|
                 |17|
                 |0|
                 |7|
                 |43|
"""

function *(A::RazprsenaMatrika, x::Vector)
    if size(A.V, 1) != length(x) 
        throw(error("Velikost matrike in vektorja nista združljiva!"))
    end

    b = zeros(length(x))
    print(b)
    for i=1:size(A.V, 1)
        for j=1:size(A.V, 2)
            if A.I[i, j] != 0
                b[i] = b[i] + A.V[i, j]*x[A.I[i, j]]
            end
        end
    end
    return b
end



"""
    sor(A::RazprsenaMatrika, b::Vector, x0::Vector, omega, tol=1e-10)

Funkcija, ki reši sistema Ax=b z uporabo metode - SOR iteracija (Successive-over-relaxation).
Vhod v funkcijo je:
    A - matrika tipa RazprsenaMatrika
    b - vektor ki je rešitev sistema
    x0 - začetni približek
    omega - relaksacijski parameter
    tol - toleranca (pogoj za vstavitev iteracij)
"""

function sor(A::RazprsenaMatrika, b::Vector, x0::Vector, omega, tol=1e-10)
    max_it = 500
    for it=1:max_it
        for i=1:size(A.V, 1)
            suma = 0
            for k=1:size(A.I, 2)
                if i != A.I[i, k] && A.I[i, k] != 0
                    suma = suma + A.V[i, k]*x0[A.I[i, k]]
                end
            end
            x0[i] = (1 - omega) * x0[i] + omega / getindex(A, i, i) * (b[i] - suma)
        end
        if norm(A*x0 - b) < tol
            return x0, it
        end
    end
    return x0, max_it
    
end

end # module DN1
