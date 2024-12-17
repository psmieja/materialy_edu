# Dynamika Molekularna 2: opracowanie i analiza wyników

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

Interesuje nas energia interakcji 

Edytujemy plik mdp dodając na końcu sekcji `Output control` linijkę:
```
energygrps = Protein JZ4
```

Przygotowujemy i uruchamiamy obliczenia.

```bash
gmx grompp -f ie.mdp -c npt.gro -t npt.cpt -p SYSTEM.top -n index.ndx -o ie.tpr

nohup gmx mdrun -deffnm ie -rerun ../prod/md_0_200ns.xtc -nb cpu &
```



### Analiza wiązań wodorowych

Gromacs posiada wbudowane narzędzia do wykrywania wiązań wodorowych. Chcemy...

```bash
gmx hbond -f [trajektoria] -s [plik tpr] -n index.ndx -hbn hbond.ndx -hbm hbond.xpm -g hbond.log
```
Wybieramy dwie grupy: białko i ligand

*W zależności od wersji gromacs, zamiast `hbond` może być potrzebne użycie komendy `hbond_legacy`*


## Przygotowanie filmiku w PyMol

Otwieramy w PyMolu plik `SYSTEM.gro`:
```bash
pymol SYSTEM.gro
```

i ładujemy przygotowaną wcześniej trajektorię
```
PyMol> load_traj [plik xtc]
```

Jeżeli w badanym układzie cząsteczki wody nie odgrywają istotnej roli, można je wszsytkie usunąć (SYSTEM -> Action -> remove waters). Jeżeli chcemy zachować konkretne cząsteczki wody odgrywające znaczenie w badanej interakcji, możemy je wykluczyć z zaznaczenia...
```bash
remove resn WAT and not resi [indeks]
```

To samo dotyczy jonów

```bash
remove resn Cl-
```

Można przybliżyć widok na interesującą nas część układu, wybierając ją i klikając (Action -> zoom). Jeżeli chcemy, żeby stała się ona również środkiem rotacji przy przeciąganiu lewym przyciskiem myszy, należy wybrać (Action -> orient).

Jeżeli w pliku z trajektorią mamy np. co setną klatkę, można "wygładzić" ruchy, żeby filmik lepiej wyglądał. **Należy jednak zaznaczyć, że zmienia to pozycje atomów i uzyskana trajektoria nie powinna być wykorzystana do analizy!**
```
smooth SYSTEM, 30, 3
```
Warto zauważyć, że wygładzenie lepiej działa do demonstrowania dużych ruchów / zmian konformacyjnych i np obroty grup CH3 potrafią wyglądać niedorzecznie (i.e. wygładzenie prowadzi do utraty informacji w krótkich skalach czasowych). Można więc usunąć niepolarne wodory (Hide -> hydrogens -> nonpolar).

Jeśli między ligandem, a białkiem tworzą się wiązania wodorowe, można je zaznaczyć (Wizard -> Measurement -> Klikamy na pierwszy atom -> na drugi -> Done)

Następnie można wykonać ray tracing dla całej trajektorii (Movie -> Ray Trace Frames). Daje to bardziej realistycznie wyglądające oświetlenie układu, **zajmuje jednak bardzo dużo czasu, dlatego podczas zajęć można pominąć ten krok.**

Następnie można wyeksportować klatki (File -> Export Movie as -> PNG Images). PyMol pozwala od razu utworzyć plik MP4, ale w praktyce nie zawsze to działa. Na tym etapie powinniśmy dostać prompt z wyborem opcji `Ray` (Ray Tracing), bądź `Draw` (bez Ray Tracingu) oraz rozdzielczość filmiku. **Na potrzeby zajęć lepiej zrezygnować z ray tracingu.** Należy też na tym etapie wybrać rozdzielczość (Jeśli generowanie klatek jest zbyt wolne, można obniżyć rozdzielczość) i nazwę bazową pliku. Tzn. wybierając nazwę "mov" dostaniemy pliki nazwane `mov0001.png`, `mov0002.png`, itd.

Mając już w folderze pliki png można je złączyć w filmik programem `ffmpeg`
```bash
ffmpeg -i mov%04d.png -c:v libx264 mov .mp4
```


## Analiza wyników w Jupyter Notebook

Możemy uruchomić Jupyter Notebook wykonując komendę
```bash
jupyter notebook
```

Notatniki *Jupyter notebook* można tworzyć, otwierać i wykonywać również w Google Colab i zapisywać na dysku, co jest wygodne zwłaszcza jeśli nie potrzebujemy lokalnej przestrzeni dyskowej i mocy obliczeniowej, a problematyczne jest przygotowanie środowiska

## Źródła, linki

Tutorial MD w gromacs
http://www.mdtutorials.com/gmx/complex/index.html

Creating movie in PyMol
https://www.blopig.com/blog/2018/12/turning-md-trajectories-into-movies-using-pymol/