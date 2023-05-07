include("DN1.jl")

using LightGraphs
using Plots
using GraphPlot
using Cairo
using GraphRecipes

"""
    grafRavnina(g, omega)

Funkcija, ki dodeli koordinate vozliščima grafa, kako bi dobljeni graf predstavljao graf v ravnini.
Vhod v funkcijo je:
    g - preprost graf
    omega - relaksacijski parameter (parameter SOR iteracije).
Funkcija ustvari diagonalno dominantno razpršeno matriko. Naredi to na naslednji način: 
    - če vrstica predstavlja vozlišče, ki ni fiksirano, potem na diagonali ima 
    minus stopnjo vozlišča, v stolpcu pa je 1, če ta stolpec  predstavlja vozlišče, ki je njegov sosed, 
    - če ta vrstica predstavlja vozlišče, ki je fiksno, potem je 1 na diagonali.
Vsako drugo vozlišče se postavi da je fiksirano.
Funkcija na začetku inicializira koordinate vozlišč na naključne vrednosti. 
(za primer, ki sem uporabljala za testiranje parametra omega, sem uporabila tiste 
določene koordinate kot začetne, zapisane so dol v nadeljavnju)
Potem uporabi funkcijo za sor iteracijo da reši sistem.
"""

function grafRavnina(g, omega)
    A2 = adjacency_matrix(g)
    A = zeros(size(A2, 1), size(A2, 2))
    for i=1:size(A, 1) 
        for j=1:size(A, 2)
            A[i, j] = A2[i, j]
        end
    end

    br = 1
    for i=1:size(A, 1)
        if mod(br, 2) != 0
            suma = 0
            for j=1:size(A, 2)
                suma = suma + A[i, j]
            end
            A[i, i] = -suma
        else
            for j in 1:size(A, 2)
                A[i, j] = 0
            end
            A[i, i] = 1
        end
        br = br + 1
    end
    
    n = nv(g)
    
    #inicializira začetne koordinate na naključne vrednosti
    x_zacetni = rand(nv(g)).- 0.5
    y_zacetni = rand(nv(g)).- 0.5

    #alternativno za primer, ki sem uporabljala za testiranje parametra omega, sem uporabila tiste 
    #koordinate kot začetne
    #x = [0.0, 2.0, 2.0, 0.0, 1.0]
    #y = [0.0, 0.0, 2.0, 2.0, 1.0]

    #tukaj izračuna desno stran b
    bx = zeros(size(A, 1))
    by = zeros(size(A, 1))
    br = 1
    for i in 1:size(A, 1)
        if mod(br, 2) != 0
            bx[i] = x_zacetni[i]
            by[i] = y_zacetni[i]
        end
        br = br + 1
    end
    
    x, iter = sor(RazprsenaMatrika(A), bx, x_zacetni, omega)
    y, iter = sor(RazprsenaMatrika(A), by, y_zacetni, omega)
    
    return x, y, iter
end

#primer preprostega grafa
g = SimpleGraph(5)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 1)
add_edge!(g, 1, 3)
add_edge!(g, 1, 5)
add_edge!(g, 4, 5)

x, y, iter = grafRavnina(g, 1.5)

#vrednosti omega, za katere sem poskusila funkcijo grafRavnina
#in število iteracij, ki sem jih dobila za ustrezan omega
omega = [0.5, 1, 1.5, 1.7, 3, 6]
iter = [69, 17, 42, 42, 500, 500]

#vizualizacija odvisnosti števila iteracij od omega
plot(omega, iter, xlabel="omega", ylabel="število iteracij", title="Odvisnost hitrosti konvergence od omega")

#vizualizacija grafa v ravnini
scatter(x, y, title="Graf v ravnini", markersize=1)
for e in edges(g)
    i, j = src(e), dst(e)
    plot!([x[i], x[j]], [y[i], y[j]], color=:black, linewidth=2)
end
display(plot!(legend=false, xlabel="X", ylabel="Y"))
