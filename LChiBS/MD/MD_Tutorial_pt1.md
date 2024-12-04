# Dynamika Molekularna 1: przygotowanie układu i uruchomienie symulacji


## Instalacja narzędzi

Instalacja Ubuntu w WSL:
```bash
wsl --install -d Ubuntu-24.04
```
Pod koniec instalacji Ubuntu w WSL trzeba skonfigurować konto użytkownika i hasło. Terminal Ubuntu powinien uruchomić się automatycznie, ale możemy go też uruchomić w aplikacji Terminal w Winowsie, wybierając `Ubuntu-24.04` z menu drop-down przy zakładkach kart, bądź wybierając program `Ubuntu-24.04` z menu start.

**Reszta komend już w Ubuntu!**

Do tworzenia i zarządzania środowiskami wirtualnymi i instalacji niektórych pakietów oprogramowania wykorzystamy **Micromambę**
```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```
Po zainstalowaniu i skonfigurowaniu Mamby możemy uruchomić ponownie terminal Ubuntu i stworzyć środowisko wirtualne na nasze potrzeby
```bash
micromamba create --name MD
```
Aktywujemy utworzone środowisko
```bash
micromamba activate MD
```
I instalujemy potrzebne programy

Najważniejszy jest pakiet AmberTools23 do przygotowania symulacji w pakiecie Amber. My będziemy używać pakietu gromacs do obliczeń, ale użyjemy pola siłowego Amber, dlatego pliki przygotujemy w AmberTools
https://ambermd.org/GetAmber.php#ambertools
```bash
micromamba install -c conda-forge ambertools=23
```
Poza tym konieczne są pakiet PDB2PQR do ustalenia sprotonowania reszt aminokwasowych i dodania wodorów do białka oraz pakiet parmed do konwersji plików
```bash
pip install pdb2pqr parmed
```

Program PyMol do wizualizacji można zainstalować z użyciem `condy` / `mamby`, ale też `apt`, albo używać instalacji w systemie Windows
```bash
micromamba install pymol-open-source
```
albo
```
apt install pymol
```

Same obliczenia uruchamiać będziemy w pakiecie **gromacs**, który możemy zainstalować `apt`-em
```bash
apt install gromacs
```

W przypadku posiadania nietypowego systemu/sprzętu, konieczności optymalizacji oprogramowania, bądź potrzeby użycia innej wersji pakietu, można pobrać kod źródłowy i skompilować go z odpowiednią konfiguracją.


## Przygotowanie plików

Pobieramy depozyt 3HTB z rcsb.org 

```bash
wget https://files.rcsb.org/download/3htb.pdb
```

### Przygotowanie ligandu

Interesuje nas ligand JZ4
```bash
awk '$1=="HETATM"' 3htb.pdb | awk '$4=="JZ4"' > jz4.pdb
```

W pliku pdb nie ma wodorów, które są istotne w MD. W pakiecie Amber jest program `reduce`, ale nie zawsze działa tak jak trzeba. Można zawsze zrobić to ręcznie w pyMolu.
W tym celu należy uruchomić pymola z plikiem pdb (`pymol jz4.pdb`), uruchomić **Builder**. Można, jeżeli to konieczne, poprawić typy wiązań, a następnie dodać automatycznie wodory **AddH** i ewentualnie poprawić ręcznie cokolwiek, co jest niepoprawne.

```bash
# przenumerowanie atomow i reszt
pdb4amber -i jz4_h.pdb -o jz4_h_renum.pdb

# przypisanie ładunków, wygenerowanie mol2
antechamber -fi pdb -fo mol2 -i jz4_h_renum.pdb -o jz4.mol2 -c bcc -pf y -nc 0

# wygenerowanie pliku frcmod (parametrów pola siłowego)
parmchk2 -i jz4.mol2 -o jz4.frcmod -f mol2
```

### Przygotowanie białka

```bash
awk '$1=="ATOM"' 3htb.pdb > protein.pdb
```

#### Sprotonowanie reszt i dodanie wodorów.

Są serwery które to potrafią (i mają dużo opcji), np.
- H++ (https://newbiophysics.cs.vt.edu/H++/)
- PDB2PQR (https://server.poissonboltzmann.org/pdb2pqr)

Z czym PDB2PQR jest też dostępne jako program standalone i tak możemy też go użyć. Musimy tu przede wszystkim ustalić pH układu, żeby poprawnie ustalić sprotonowanie reszt białka (możemy np wziąć pH eksperymentu z publikacji).

```bash
pdb2pqr --ffout AMBER --with-ph 6.5 --pdb-output protein_H.pdb protein.pdb protein_H.pdb

# od razu przenumerujmy...
pdb4amber -i protein_H.pdb -o protein_H_renum.pdb
```


### Przygotowanie białka

Trzeba teraz białko i ligand wrzucić w jeden plik

```bash 
cat protein_H_renum.pdb jz4_h_renum.pdb | awk '$1!="END"' > system.pdb
pdb4amber -i system.pdb -o system_renum.pdb
```

### Przygotowanie układu w `tleap`


Tworzymy plik wejściowy do programu tleap `tleap.in` zawierający:

```
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2

loadamberparams jz4.frcmod

JZ4 = loadmol2 jz4.mol2
mol = loadpdb system_renum.pdb

solvatebox mol TIP3PBOX 10.0
addions mol Na+ 0
addions mol Cl- 0
savepdb mol SYSTEM_solv.pdb
saveamberparm mol SYSTEM.prmtop SYSTEM.inpcrd
quit
```

I urachamiamy skrypt w `tleap`:

```bash
tleap -f tleap.in 
```

### Konwersja plików symulacji do formatu GROMACS

Same obliczenia MD będziemy odpalać w GROMACS, a nie Amber, dlatego zmienimy format plików z topologią i współrzędnymi. Można do tego użyć pakietu `parmed`.

W interpreterze Pythona:

```python
import parmed
parm = parmed.load_file('./SYSTEM.prmtop', './SYSTEM.inpcrd')
parm.save('SYSTEM.top', format='gromacs')
parm.save('SYSTEM.gro')
```

Teraz mamy pliki `.top` i `.gro` potrzebne do uruchomienia symulacji


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

```bash
cd ..
mkdir nvt
cp min/SYSTEM.top nvt
cp min/em.gro nvt
wget http://www.mdtutorials.com/gmx/complex/Files/nvt.mdp
```

Przygotowujemy grupy termostatowania
```
gmx make_ndx -f em.gro -o index.ndx
```

Przygotowujemy jedną grupę złożoną z białka i ligandu
```
> 1 | 13

Copied index group 1 'Protein'
Copied index group 13 'JZ4'
Merged two groups with OR: 2614 22 -> 2636
```

Pojawia się nowa grupa ` 18 Protein_JZ4 :  2636 atoms`. Wszystko pozostałe chcemy wrzucić w drugą grupę:
```
> ! 18 

Copied index group 18 'Protein_JZ4'
Complemented group: 24576 atoms
```

Możemy ją nazwać `Water_and_ions`

```
> name 19 Water_and_ions 
```

I wychodzimy
```
> q
```

Pojawił się plik `index.ndx` zawierający definicje grup. Teraz chcemy wejść w plik `nvt.mdp`, sprawdźmy linijkę `tc-grps`. Powinno być `Protein_JZ4 Water_and_ions`. W razie potrzeby edytujemy mdp, żeby poprawić grupki termostatowania

### Równoważenie NVT

Przygotowujemy plik symulacji...

```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p SYSTEM.top -n index.ndx -o nvt.tpr -maxwarn 2
```

Potem, jak już chcemy uruchomić symulację, robimy to poleceniem `gmx mdrun`. Żeby komenda działała w tle, trzeba dorzucić na końcu `&`. Niestety, jeśli rozłączymy się z serwerem, celowo lub przypadkowo, program zostanie przerwany. Żeby temu zapobiec odpalamy go z `nohup`:

```bash
nohup gmx mdrun -deffnm nvt &
```

### Równoważenie NPT

Najpierw przygotowujemy wszystkie potrzebne pliki
```
mkdir ../npt
cd ../npt
cp ../nvt/nvt.cpt .
cp ../nvt/nvt.gro .
cp ../nvt/index.ndx .
cp ../nvt/SYSTEM.top .
wget http://www.mdtutorials.com/gmx/complex/Files/npt.mdp
```

Można na tym etapie wprowadzić niezbędne zmiany do pliku `mdp`

```bash
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p SYSTEM.top -n index.ndx -o npt.tpr -maxwarn 2
```

Uruchomienie symulacji
```bash
nohup gmx mdrun -deffnm npt &
```


### Produkcyjne MD 

```bash
mkdir ../prod
cd ../prod
cp ../npt/npt.cpt .
cp ../npt/npt.gro .
cp ../npt/index.ndx .
cp ../npt/SYSTEM.top .
wget http://www.mdtutorials.com/gmx/complex/Files/md.mdp
```

Edytujemy w pliku mdp (`vim md.mdp` bądź inny edytor) czas symulacji i częstość wypisywania rzeczy tak, żeby dostać 200ns symulacji i jednocześnie nie mieć za dużo outputu. Warto zobaczyć, że uruchamiając `grompp` dostajemy informację ile przestrzeni dyskowej powinien zająć outupt!

```
; Run parameters
integrator              = md         ; leap-frog integrator
nsteps                  = 100000000  ; 2 * 100000000 = 200000 ps (200 ns)
dt                      = 0.002      ; 2 fs
; Output control
nstenergy               = 50000      ; save energies every 100.0 ps
nstlog                  = 50000      ; update log file every 100.0 ps
nstxout-compressed      = 50000      ; save coordinates every 100.0 ps
```

A następnie odpalamy symulację

```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p SYSTEM.top -n index.ndx -o md_0_200ns.tpr -maxwarn 2

nohup gmx mdrun -deffnm md_0_200ns &
```


Możemy sprawdzać postęp symulacji odczytując logi
```
tail -n 30 md_0_200ns.log
```

Możemy jednocześnie sprawdzać wykorzystanie zasobów poleceniem `htop` (dla karty NVIDIA `nvidia-smi`).


W wyniku dostajemy kilka plików, w tym:
- `md_0_200ns.log` - log symulacji (który odczytywaliśmy, żeby sprawdzić postęp symulacji
- `md_0_200ns.xtc` - plik trajektorii z zapisanymi położeniami atomów w poszczególnych klatkach.
- `md_0_200ns.edr` - plik zawierający informacje o wartościach energii, łącznie z poszczególnymi jej członami (elektrostatyczna, Van der Waalsa, etc)

Pliki `.edr` i `.xtc` są plikami binarnymi i nie można ich obejrzeć w zwykłym edytorze tekstu i należy do tego użyć albo dedykowanych programów pakietu gromacs, albo innych pakietów do analizy trajektorii MD, co omówione zostanie na 2. zajęciach.


## Źródła, linki

Depozyt 3HTB, lisozym + ligand JZ4
https://www.rcsb.org/structure/3HTB

Publikacja dot. depozytu
https://pmc.ncbi.nlm.nih.gov/articles/PMC2788029/pdf/nihms149280.pdf

Tutorial MD w gromacs (stąd układ i pliki mdp)
http://www.mdtutorials.com/gmx/complex/index.html

Tutorial Amberowy do modelowania P450 z hemem, stąd część informacji
https://ambermd.org/tutorials/advanced/tutorial20/mcpbpy_heme.php

Dokumentacja AmberTools23
https://ambermd.org/doc12/Amber23.pdf

PDB2PQR do protonacji i dodawania wodorów do białka
https://pdb2pqr.readthedocs.io/en/latest/using/index.html

#### Gromacs

Wszystkie opcje w pliku mdp
https://manual.gromacs.org/current/user-guide/mdp-options.html

Program `gmx grompp`
https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html#gmx-grompp