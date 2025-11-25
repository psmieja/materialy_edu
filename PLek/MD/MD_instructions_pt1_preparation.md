# Dynamika Molekularna 1: przygotowanie układu

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

Pamiętaj, że musisz ładować potrzebne moduły (polecenie `module load`)
Będziemy wykorzystywać:
- `gromacs/2024.4`
- `ambertools/24`
- `pymol`
- opcjonalnie `python`

**Uwaga:** Poszczególne programy mogą nie działać jeśli są załadowane inne moduły (mogą wystąpić konflikty zależności na danym urządzeniu). Można wtedy np ustawić osobne okna terminala do poszczególnych narzędzi, albo przed użyciem konkretnego programu "odładować" wszystkie moduły (`module reset`) i załadować tylko te konieczne do uruchomienia danego programu.

**Uwaga:** W ćwiczeniu tym tworzymy wiele plików. Polecamy, w związku z tym, stworzyć nowy katalog na to zadanie.


### Pobranie depozytu

Pobierz **2BQV** depozyt w formacie `pdb`

### Przygotowanie leku

Zapisz współrzędne leku (odpowiednie rekordy `HETATM`) z depozytu do nowego pliku `A1A_raw.pdb`. Upewnij się, że w pliku są tylko porządane współrzędne. Możesz uruchomić uzyskany plik w PyMol i go obejrzeć.

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

Kolejnym krokiem jest ustalenie protonacji reszt białka i dodanie do odpowiednio atomów wodoru. Co istotne Amber ma własne kody reszt zależne od sprotonowania (np. zamiast `HIS` może być `HIE`).

Na tym etapie trzeba ustalić w jakim pH odbywa się symulacja. Na potrzeby tego ćwiczenia proszę przyjąć pH 5.5. 

Istnieje kilka rozwiązań.
- Serwer H++ http://newbiophysics.cs.vt.edu/H++/
- Serwer PDB2PQR https://server.poissonboltzmann.org/pdb2pqr
- narzędzie w lini komend PDB2PQR

Polecamy ostatnie, jako najprostsze; można tam uzyskać plik `pdb` z poprawnymi typami reszt w formacie Amber uzyskać możemy wykonując następującą komendę:
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
