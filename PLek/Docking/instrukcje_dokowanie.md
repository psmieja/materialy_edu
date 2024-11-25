## Oprogramowanie

Wykorzystujemy

- AutoDockTools: narzędzie GUI do przygotowania plików do dokowania AutoDock i Viną ( instalowane jako część MGLTools https://ccsb.scripps.edu/mgltools/ )

- Vina: Oprogramowanie CLI do dokowania ( https://vina.scripps.edu/downloads/ )

- PyMol do inspekcji i wizualizacji wyników

## Białko

Będziemy pracować z proteazą HIV-1 i jej inhibitorami.

https://en.wikipedia.org/wiki/HIV-1_protease

Ma bardzo dużo struktur z różnymi ligandami. Na nasze potrzeby użyjemy depozytu 1HPV. Póki co jednak przejdźmy do ligandu.

## Wybór i przygotowanie ligandu

Interesują nas leki stanowiące inhibitory tego białka

> There are ten HIV-1 PR inhibitors that are currently approved by the [Food and Drug Administration](https://en.wikipedia.org/wiki/Food_and_Drug_Administration "Food and Drug Administration"): [indinavir](https://en.wikipedia.org/wiki/Indinavir "Indinavir"), [saquinavir](https://en.wikipedia.org/wiki/Saquinavir "Saquinavir"), [ritonavir](https://en.wikipedia.org/wiki/Ritonavir "Ritonavir"), [nelfinavir](https://en.wikipedia.org/wiki/Nelfinavir "Nelfinavir"), [lopinavir](https://en.wikipedia.org/wiki/Lopinavir "Lopinavir"), [amprenavir](https://en.wikipedia.org/wiki/Amprenavir "Amprenavir"), [fosamprenevir](https://en.wikipedia.org/wiki/Fosamprenavir "Fosamprenavir"), [atazanavir](https://en.wikipedia.org/wiki/Atazanavir "Atazanavir"), [tipranavir](https://en.wikipedia.org/wiki/Tipranavir "Tipranavir"), and [darunavir](https://en.wikipedia.org/wiki/Darunavir "Darunavir").



Można wziąć jeden z nich i pobrać z

 https://pubchem.ncbi.nlm.nih.gov/

Najlepiej pobrać konfomer 3D. Jeśli w PubChemie nie ma konfomerów 3D można je utworzyć w **RDKit** (jak na zajęciach), albo np. za pomocą serwera

https://www.novoprolabs.com/tools/smiles2pdb

Mogą być też dostępne do pobrania na RCSB

https://www.rcsb.org/



Mając strukturę ligandu, możemy ją załadować do AutoDockTools. Jeśli plik wymaga konwersji do innego formatu można to zrobić w PyMolu (otwierając plik, potem File -> Export Molecule -> Save -> Wybieramy typ i nazwę -> Save)



#### W AutoDockTools

Mając ligand w ADT:

Brązowe menu -> Ligand -> Input -> Choose -> Zaznaczamy Ligand -> Select Molecule For AutoDock4

Brązowe menu -> Ligand -> Output -> Save as PDBQT -> Zapisujemy jako ligand.pdbqt



## Przygotowanie białka

Pobieramy depozyt 1HPV

```shell
wget https://files.rcsb.org/download/1HPV.pdb
```

Wyciągamy same atomy białka. Nie interesuje nas woda i obecny ligand

```bash
grep ATOM 1HPV.pdb > protein.pdb
```

#### W ADT

Otwieramy białko w ADT i:

Edit -> Hydrogens -> Add -> Zaznaczamy "Polars only""

Edit -> Charges -> Add Kollman Charges

Grid -> Macromolecule -> Choose

Grid -> Grid Box -> ustawiamy parametry wg. uznania  suwaczkami tak, żeby zaznaczyć całą kieszeń aktywną -> File -> Output Grid Dimensions File -> zapisujemy jako grid.txt

#### Grid Box Config

W `grid.txt` mamy format, którego nie przyjmuje `vina`, także należy utworzyć plik `config.txt` określający grid box w poprawnym formacie

```
center_x = 
center_y = 
center_z = 
size_x = 
size_y = 
size_z = 
```

Wartości są w Angstromach, także jeśli mamy "`npts 50 ... ... ...`" przy domyślnym spacing 0.375, oznacza to, że `size_x` wynosi $50 \cdot 0.375 \rm{\AA} = 18.75 \rm{\AA}$

## Dokowanie

Mając wszystkie pliki, uruchamiamy dokowanie w `vina`:

```bash
vina --receptor protein.pdbqt --ligand ligand.pdbqt \
     --config config.txt --exhaustiveness=32 \
     --num_modes 20 --out out_32_20.pdbqt
```

`out_[exhaustiveness]_[num_modes].pdbqt `

## Dodatkowe materiały

Zamiast AutoDockTools możemy do przygotowania układu użyć narzędzi CLI takich, jak Mekko ( https://meeko.readthedocs.io/en/release ). Tutorial:

https://autodock-vina.readthedocs.io/en/latest/docking_basic.html

Można dokowanie wykorzystać do HTS (High-Throughput screening), uruchamiając je dla całej bazy cząsteczek, np. podzbioru bazy ZINC
