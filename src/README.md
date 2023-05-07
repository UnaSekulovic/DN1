# SOR iteracija za razpršene matrike
V tej nalogi smo implementirali reševanje linearnega sistema enačb z uporabo SOR iteracije za razpršene matrike.
Na začetku smo definirali tip RazprsenaMatrika, ki diagonalno dominantno razpršeno matriko zaradi prostorskih zahtev hrani v dveh matrikah.
Za ta tip smo nato definirali metode, ki jih bomo uporabljali v funkciji za izračun SOR iteracije, in sicer metode za indeksiranje (getindex, setindex, firstindex, lastindex) in množenje matrike z desne z vektorjem.
Potem smo definirali funkcijo, ki rešuje sistem Ax=b s SOR iteracijo, kjer je A diagonalno dominantna razpršena matrika. Ta funkcija vrne rešitev sistema in število iteracij, ki so bile potrebne za dosego rešitve.

Na koncu smo še to metodo uporabili za vožitev grafa v ravnino. Definirali smo funkcijo, ki za dani graf fiksira vsako drugo vozlišče in nato uporabi funkcijo za SOR iteracijo, ki konstruira vožitev grafa v ravnino, tako da ostala vozlišča zavzamejo ravnovesno lego med fiksiranimi vozlišči. 

Za naš primer smo vzeli preprost graf, sestavljen iz 5 vozlišč, kjer so povezave v grafu naslednje: (1, 2), (2, 3), (3, 4), (4, 1), (1, 3), (1, 5), (4, 5). Za ta graf smo preizkusili nekaj vrednosti za parameter omega in pogledali, koliko iteracij je potrebnih, da konvergira SOR. Na sliki lahko vidimo odvisnost hitrosti konvergence od izbire parametra omega: 

![odvisnost](slike/odvisnost.svg)

Vidimo, da je v našem primeru najhitrejša konvergenca dosežena pri vrednosti omega  enako 1. Prav tako je razvidno, da je na splošno v tem primeru najhitrejša konvergenca dosežena za vrednost omega od 1 do okoli 1.8.