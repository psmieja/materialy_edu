# Dynamika Molekularna 2: symulacja i wizualizacja


## Obliczenia MD

### Minimalizacja

Wykorzystamy gotowe pliki z parametrami symulacji, wprowadzając jedynie niezbędne zmiany. Opis wszystkich parametrów w pliku mdp można znaleźć w dokumentacji (https://manual.gromacs.org/current/user-guide/mdp-options.html)

```bash
wget http://www.mdtutorials.com/gmx/complex/Files/em.mdp
```

```bash
gmx grompp -f em.mdp -c SYSTEM.gro -p SYSTEM.top -o em.tpr -maxwarn 2
gmx mdrun -v -deffnm em
```

### Przygotowanie NVT

Pobieramy plik

```bash
wget http://www.mdtutorials.com/gmx/complex/Files/nvt.mdp
```

Przygotowujemy grupy termostatowania

```shell
gmx make_ndx -f em.gro -o index.ndx
```

Przygotowujemy jedną grupę złożoną z białka i ligandu

```shell
> 1 | 13
Copied index group 1 'Protein'
Copied index group 13 'A1A'
Merged two groups with OR: 2614 22 -> 2636
```

Pojawia się nowa grupa ` 18 Protein_A1A :  2636 atoms`. Wszystko pozostałe chcemy wrzucić w drugą grupę:

```shell
> ! 18 
Copied index group 18 'Protein_A1A'
Complemented group: 24576 atoms
```

Możemy ją nazwać `Water_and_ions`

```shell
> name 19 Water_and_ions 
```

I wychodzimy

```shell
> q
```

Pojawił się plik `index.ndx` zawierający definicje grup. Teraz chcemy wejść w plik `nvt.mdp`, sprawdźmy linijkę `tc-grps`. Powinno być `Protein_A1A Water_and_ions`. W razie potrzeby edytujemy mdp, żeby poprawić grupy termostatowania

### Równoważenie NVT

Przygotowujemy plik symulacji...

```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p SYSTEM.top -n index.ndx -o nvt.tpr -maxwarn 2
```

Potem, jak już chcemy uruchomić symulację, robimy to poleceniem `gmx mdrun`. 

```bash
gmx mdrun -deffnm nvt &
```

### Równoważenie NPT

Pobieramy plik z parametrami symulacji.

```shell
wget http://www.mdtutorials.com/gmx/complex/Files/npt.mdp
```

Można na tym etapie wprowadzić niezbędne zmiany do pliku `mdp`. Przede wszsytkim należy zmienić nazwy grup termostatowania, jak w pliku do NVT. Następnie, znowu analogicznie do NVT:

```bash
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p SYSTEM.top -n index.ndx -o npt.tpr -maxwarn 2
```
Z tą różnicą, że teraz zaczynamy z "checkpointu" (`npt.cpt`)

Uruchomienie symulacji

```bash
gmx mdrun -deffnm npt
```

### Produkcyjne MD

Edytujemy w pliku mdp (`vim md.mdp` bądź inny edytor) czas symulacji i częstość wypisywania rzeczy. Warto zobaczyć, że uruchamiając `grompp` dostajemy informację ile przestrzeni dyskowej powinien zająć outupt!

```markup
; Run parameters
integrator              = md           ; leap-frog integrator
nsteps                  = [ile kroków] ;
dt                      = 0.002        ; 2 fs
; Output control
nstenergy               = [co ile kroków zapisać energię] 
nstlog                  = [co ile kroków wypisać logi]
nstxout-compressed      = [co ile kroków zapisać współrzędne]
```

**UWAGA:** Na ćwiczeniach ograniczmy czas symulacji do co najwyżej 100ps. Prawdziwa produkcyjna symulacja jest bardzo czasochłonna i nie bylibyśmy w stanie jej uruchomić na zajęciach.

Następnie uruchamiamy symulację

```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p SYSTEM.top -n index.ndx -o md.tpr -maxwarn 2
gmx mdrun -deffnm md
```
W wyniku dostajemy kilka plików, w tym:

- `md.log` - log symulacji (który odczytywaliśmy, żeby sprawdzić postęp symulacji
- `md.xtc` - plik trajektorii z zapisanymi położeniami atomów w poszczególnych klatkach.
- `md.edr` - plik zawierający informacje o wartościach energii, łącznie z poszczególnymi jej członami (elektrostatyczna, Van der Waalsa, etc)
  Pliki `.edr` i `.xtc` są plikami binarnymi i nie można ich obejrzeć w zwykłym edytorze tekstu i należy do tego użyć albo dedykowanych programów pakietu gromacs, albo innych pakietów do analizy trajektorii MD, co omówione zostanie na 2. zajęciach.




## Przygotowanie filmiku w PyMol

Otwieramy w PyMolu plik `SYSTEM.gro`:
```bash
pymol SYSTEM.gro
```

i ładujemy przygotowaną wcześniej trajektorię
```
PyMol> load_traj [plik xtc]
```

**Uwaga: w niektórych wersjach PyMol może być problem z załadowaniem trajektorii w ten sposób.** Zamiast tego można przygotować plik pdb, w którym kolejne klatki figurują jako osobne MODELe
```bash
gmx trjconv -s [plik tpr] -f [trajektoria xtc] -n [plik ndx] -o [plik wyjściowy pdb]
```
Możemy na tym etapie wybrać, które atomy mają zostać zapisane w pliku pdb. Jeśli w badanym układzie nie są dla nas istotne cząsteczki wody i jony, można jako `Output` wybrać tylko białko i ligand (`Protein_A1A`). W rezultacie dostajemy plik `pdb`, który możemy otworzyć w PyMol (co może chwilę zająć, bo wszystkie klatki muszą się załadować) i nie musimy już osobno ładować trajektorii.

Jeżeli w PyMolu jeszcze mamy cząsteczki wody, ale jednak ich nie potrzebujemy, możemy je wszysktie usunąć (SYSTEM -> Action -> remove waters). Jeżeli chcemy zachować konkretne cząsteczki wody odgrywające znaczenie w badanej interakcji, możemy je wykluczyć z zaznaczenia...
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

## Źródła, linki

Tutorial MD w gromacs
http://www.mdtutorials.com/gmx/complex/index.html

Creating movie in PyMol
https://www.blopig.com/blog/2018/12/turning-md-trajectories-into-movies-using-pymol/