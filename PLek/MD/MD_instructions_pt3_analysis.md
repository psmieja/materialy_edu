# Dynamika Molekularna 3: opracowanie i analiza wyników

## Przypomnienie: formaty plików
- `gro` - współrzędne układu
- `top` - topologia układu
- `ndx` - indeksy atomów należących do poszczególnych grup, takich jak *Protein*, czy *Water_and_ions*
- `mdp` - (parametry symulacji MD takie, jak liczba iteracji, krok czasowy, algorytm całkowania, etc.
- `tpr` - plik binarny z topologią, współrzędnymi i prędkościami atomów, parametrami symulacji itp. Zawiera wszsytkie informacje potrzebne do uruchomienia symulacji. Przygotowujemy go w z wykorzystaniem
- `xtc` - binarny plik z trajektorią (współrzędnymi atomów w poszczególnych krokach czasowych)
- `edr` - plik z wartościami poszczególnych członów energii w  wybranych krokach czasowych.
 
Lista wszystkich formatów plików w gromacs:
https://manual.gromacs.org/2023.3/reference-manual/file-formats.html

## Postprodukcja

Symulacja odbywa się w periodycznych warunkach brzegowych i białko się przemieszcza, także często po pewnym czasie symulacji znajdzie się na brzegu "pudełka" symulacji, co sprawia, że w trajektorii atomy są rozrzucone po dwóch brzegach. Na potrzeby analizy najlepiej jest ustawić białko na środku:

```bash
gmx trjconv -s [plik tpr] -f [trajektoria wejściowa] -o [trajektoria wyjściowa] -center -pbc mol -ur compact
```
Wybieramy `Protein` dla `centering` i `System` jako `output`

To jednak nie chroni białka przed obracaniem się w czasie, co nie jest wygodne w analizie. Realnie interesuje nas rzeczywista ewolucja struktury i interakcji białko-ligand, a nie obrót i przesunięcie układu, które nie świadczą o rzeczywistej zmianie układu. Także możemy przesunąć (rotacja + translacja) w ten sposób, żeby kolejne klatki były mozliwie podobne do siebie:

```bash 
gmx trjconv -s [plik tpr] -f [trajektoria wejściowa] -o [trajektoria wyjściowa] -fit rot+trans
```
Wybieramy `Backbone` do `fit`owania i `System` jako `output`


### Energia interakcji

Interesuje nas energia interakcji białko-ligand. Te człony nie są jednakowoż wypisywane domyślnie w symulacji. Możemy uruchomić obliczenia ponownie, zmieniając jednak tylko to, co jest wypisywane do pliku `edr`, bez właściwych obliczeń.

Jak na poprzednich etapach jest nam potrzebny plik z parametrami (`ie.mdp`).Powinien być on identyczny z plikiem `mdp` z oryginalnej symulacji (tej, z której trajektorię wykorzystujemy), z tą różnicą, że na końcu sekcji `Output control` linijkę:
```
energygrps = Protein A1A
```

Używamy komendy `grompp`, żeby przygotować plik `tpr`
```bash
gmx grompp -f ie.mdp -c npt.gro -t npt.cpt -p SYSTEM.top -n index.ndx -o ie.tpr
```

Uruchamiamy

```bash
gmx mdrun -deffnm ie -rerun [trajektoria xtc] -nb cpu &
```

### Analiza wiązań wodorowych

Gromacs posiada wbudowane narzędzia do wykrywania wiązań wodorowych. Chcemy...

```bash
gmx hbond -f [trajektoria] -s [plik tpr] -n index.ndx -hbn hbond.ndx -hbm hbond.xpm -g hbond.log
```
Wybieramy dwie grupy: białko i ligand

*W zależności od wersji gromacs, zamiast `hbond` może być potrzebne użycie komendy `hbond-legacy`*


## Analiza wyników w Colabie (Jupyter notebook)

...

## Źródła, linki

Tutorial MD w gromacs
http://www.mdtutorials.com/gmx/complex/index.html

Creating movie in PyMol
https://www.blopig.com/blog/2018/12/turning-md-trajectories-into-movies-using-pymol/