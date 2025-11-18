## Wstep

Znowu pracujemy z proteazą HIV-1

https://en.wikipedia.org/wiki/HIV-1_protease

Wykorzystamy depozyt wykrystalizowany już z lekiem

https://www.rcsb.org/structure/2BQV

**Uwaga:** Jeśli masz wątpliwości, co do parametrów danego programu, możesz użyć help, co zazwyczaj można uzyskać poprzez flagę `-h`, albo uruchamiając program bez żadnych argumentów, e.g.

```shell
pdb4amber -h
```

### Załadowanie modułów

Pamiętaj, że na początku musisz załadować potrzebne moduły

Będziemy wykorzystywać 

- `gromacs/2024.4`

- `ambertools/24`

### Pobranie depozytu

Pobierz **2BQV** depozyt w formacie `pdb`

### Przygotowanie leku

Zapisz współrzędne leku (odpowiednie rekordy `HETATM` z depozytu do nowego pliku `A1A_raw.pdb` 



Wykorzustując program `reduce` z pakietu `AmberTools` dodaj wodory do leku. Zapisz uzyskaną cząsteczkę do nowego pliku `A1A_H.pdb`  

```shell
reduce A1A_raw.pdb > A1A_H.pdb
```

Obejrzyj plik w PyMol, zobacz, czy wygląda poprawnie



Wykorzystaj program `pdb4amber`, żeby poprawić numerację w pliku. Zapisz plik wynikowy jako `A1A_renum.pdb`

```shell
pdb4amber -i A1A_H.pdb -o A1A_renum.pdb
```

W następnej kolejności trzeba przygotować plik `mol2` który zawiera informacje o ładunkach i wiązniach. 

```shell
antechamber -fi pdb -fo mol2 -i A1A_renum.pdb -o A1A.mol2 -c bcc -pf y -nc 0
```

Konieczne jest przygotowanie również pliku `frcmod`, który zawiera dodatkowe parametry pola siłowego pot

```shell
parmchk2 -i A1A.mol2 -o A1A.frcmod -f mol2
```

### Przygotowanie białka

Wyciągamy atomy białka z pliku depozytu (wszystkie rekordy `ATOM` ) i zapisujemy je do plik `protein_raw.pdb`

```shell
awk '$1=="ATOM"' 2BQV.pdb > protein_raw.pdb
```

Ustalenie protonacji reszt i dodanie wodoru do białka 

Istnieje kilka rozwiązań.

- Serwer H++ http://newbiophysics.cs.vt.edu/H++/

- Serwer PDB2PQR https://server.poissonboltzmann.org/pdb2pqr

- narzędzie CMD PDB2PQR

```shell
pdb2pqr --ffout AMBER --titration-state-method=propka --with-ph 5.5 --pdb-output protein_H.pdb ./protein_raw.pdb protein_H.pqr
```

Numeracja (analogicznie)

```shell
pdb4amber -i protein_H.pdb -o protein_renum.pdb
```

### System

Złączenie plików

```shell
cat protein_renum.pdb A1A_renum.pdb | awk '$1!="END"' > system.pdb
```

Numeracja (analogicznie)

```shell
cat protein_renum.pdb A1A_renum.pdb | awk '$1!="END"' > system.pdb
```

Utwórz plik `tleap.in` o treści:

```shell
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2

loadamberparams A1A.frcmod

A1A = loadmol2 A1A.mol2
mol = loadpdb system_renum.pdb

solvatebox mol TIP3PBOX 10.0
addions mol Na+ 0
addions mol Cl- 0
savepdb mol SYSTEM_solv.pdb
saveamberparm mol SYSTEM.prmtop SYSTEM.inpcrd
quit
```

Uruchomienie tleap

```shell
tleap -f tleap.in 
```



Uzyskane pliki `SYSTEM.prmtop` i `SYSTEM.inpcrd` są dostosowane do programu Amber. My obliczenia jednak uruchamiać będziemy w Gromacsie, zatem trzeba przekonwertować pliki.



W tym celu, w  interpreterze Python wykonaj poniższy kod. Jeśli nie jest dostępny `parmed`, można zainstalować go komendą `pip install parmed`. Jeśli nie ma takiej możliwości lokalnie, można to uczynić w Colabie

```python
import parmed
parm = parmed.load_file('./SYSTEM.prmtop', './SYSTEM.inpcrd')
parm.save('SYSTEM.top', format='gromacs')
parm.save('SYSTEM.gro')
```

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

Pojawia się nowa grupa ` 18 Protein_JZ4 :  2636 atoms`. Wszystko pozostałe chcemy wrzucić w drugą grupę:

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

Pojawił się plik `index.ndx` zawierający definicje grup. Teraz chcemy wejść w plik `nvt.mdp`, sprawdźmy linijkę `tc-grps`. Powinno być `Protein_A1A Water_and_ions`. W razie potrzeby edytujemy mdp, żeby poprawić grupki termostatowania

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

```shell
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

```markup
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

```shell
tail -n 30 md_0_200ns.log
```

Możemy jednocześnie sprawdzać wykorzystanie zasobów poleceniem `htop` (dla karty NVIDIA `nvidia-smi`).
W wyniku dostajemy kilka plików, w tym:

- `md_0_200ns.log` - log symulacji (który odczytywaliśmy, żeby sprawdzić postęp symulacji
- `md_0_200ns.xtc` - plik trajektorii z zapisanymi położeniami atomów w poszczególnych klatkach.
- `md_0_200ns.edr` - plik zawierający informacje o wartościach energii, łącznie z poszczególnymi jej członami (elektrostatyczna, Van der Waalsa, etc)
  Pliki `.edr` i `.xtc` są plikami binarnymi i nie można ich obejrzeć w zwykłym edytorze tekstu i należy do tego użyć albo dedykowanych programów pakietu gromacs, albo innych pakietów do analizy trajektorii MD, co omówione zostanie na 2. zajęciach.
