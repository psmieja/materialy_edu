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
- `pymol`

### Pobranie depozytu

Pobierz **2BQV** depozyt w formacie `pdb`

### Przygotowanie leku

Zapisz współrzędne leku (odpowiednie rekordy `HETATM` z depozytu do nowego pliku `A1A_raw.pdb` 



Wykorzustując program `reduce` z pakietu `AmberTools` dodaj wodory do leku. Zapisz uzyskaną cząsteczkę do nowego pliku `A1A_H.pdb`  

```shell
reduce [plik wejściowy] > [plik wyjściowy]
```

Obejrzyj plik w PyMol, zobacz, czy wygląda poprawnie



Wykorzystaj program `pdb4amber`, żeby poprawić numerację w pliku `A1A_H.pdb`. Zapisz plik wynikowy jako `A1A_renum.pdb`

```shell
pdb4amber -i [plik wejściowy] -o [plik wyjściowy]
```

W następnej kolejności trzeba przygotować plik w formacie `mol2` który zawiera informacje o ładunkach i wiązniach. W tym celu wykorzystaj program `antechamber` wykonując poniższą komendę

```shell
antechamber -fi pdb -fo mol2 -i A1A_renum.pdb -o A1A.mol2 -c bcc -pf y -nc 0
```

Konieczne jest przygotowanie również pliku `frcmod`, który zawiera dodatkowe parametry pola siłowego dla leku. Wykonaj poniższą komendę.

```shell
parmchk2 -i A1A.mol2 -o A1A.frcmod -f mol2
```


### Przygotowanie białka

W następnej kolejności należy wyciągnąć same atomy białka z pliku depozytu (wszystkie rekordy `ATOM`). Zapisz je do pliku
`protein_raw.pdb`

Kolejnym krokiem jest ustalenie protonacji reszt białka i dodanie do odpowiednio atomów wodoru. Na tym etapie trzeba ustalić w jakim pH odbywa się symulacja. Na potrzeby tego ćwiczenia proszę przyjąć pH 5.5

Istnieje kilka rozwiązań.
- Serwer H++ http://newbiophysics.cs.vt.edu/H++/
- Serwer PDB2PQR https://server.poissonboltzmann.org/pdb2pqr
- narzędzie w lini komend PDB2PQR

W ostatnim plik `pdb` z poprawnymi typami reszt w formacie Amber uzyskać możemy wykonując następującą komendę:
```shell
pdb2pqr --ffout AMBER --titration-state-method=propka --with-ph 5.5 --pdb-output protein_H.pdb ./protein_raw.pdb protein_H.pqr
```

Poprawa numeracji białka (analogicznie do ligandu)

```shell
pdb4amber -i protein_H.pdb -o protein_renum.pdb
```

### System

Utwórz nowy plik `system.pdb` zawierający zarówno atomy białka jak i ligandu zawarte w utworzonych wcześniej plikach `protein_renum.pdb` i `A1A_renum.pdb`. **Uwaga:** Nie umieszczaj w nowym pliku rekordu `END`

Ponownie uruchom program `pdb4amber`, tym razem dla całego systemu, tworząc plik `system_renum.pdb`

### Przygotowanie symulacji: tleap
Następnym krokiem jest przygotowanie plików wejściowych do symulacji. Wymaga to między innymi zdefiniowania pól siłowych, załadowania współrzędnych oraz dodania rozpuszczalnika. W przypadku pól siłowych Amber robimy to w programie leap

Utwórz nowy plik `tleap.in` o treści:

```
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

Pojawił się plik `index.ndx` zawierający definicje grup. Teraz chcemy wejść w plik `nvt.mdp`, sprawdźmy linijkę `tc-grps`. Powinno być `Protein_A1A Water_and_ions`. W razie potrzeby edytujemy mdp, żeby poprawić grupy termostatowania

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
nohup gmx mdrun -deffnm npt &
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

A następnie odpalamy symulację

```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p SYSTEM.top -n index.ndx -o md.tpr -maxwarn 2
nohup gmx mdrun -deffnm md &
```

Możemy sprawdzać postęp symulacji odczytując logi

```shell
tail -n 30 md.log
```

Możemy jednocześnie sprawdzać wykorzystanie zasobów poleceniem `htop` (dla karty NVIDIA `nvidia-smi`).
W wyniku dostajemy kilka plików, w tym:

- `md.log` - log symulacji (który odczytywaliśmy, żeby sprawdzić postęp symulacji
- `md.xtc` - plik trajektorii z zapisanymi położeniami atomów w poszczególnych klatkach.
- `md.edr` - plik zawierający informacje o wartościach energii, łącznie z poszczególnymi jej członami (elektrostatyczna, Van der Waalsa, etc)
  Pliki `.edr` i `.xtc` są plikami binarnymi i nie można ich obejrzeć w zwykłym edytorze tekstu i należy do tego użyć albo dedykowanych programów pakietu gromacs, albo innych pakietów do analizy trajektorii MD, co omówione zostanie na 2. zajęciach.
