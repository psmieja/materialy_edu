# Dynamika Molekularna 2: opracowanie i analiza wyników

### Postprodukcja

Symulacja odbywa się w periodycznych warunkach brzegowych i białko się przemieszcza, także często po pewnym czasie symulacji znajdzie się na brzegu "pudełka" symulacji, co sprawia, że w trajektorii atomy są rozrzucone po dwóch brzegach. Na potrzeby analizy najlepiej jest ustawić białko na środku:

```bash
gmx trjconv -s md_0_200ns.tpr -f md_0_200ns.xtc -o md_0_200ns_center.xtc -center -pbc mol -ur compact
```
Wybieramy `Protein` dla `centering` i `System` jako `output`

To jednak nie chroni białka przed obracaniem się w czasie, co nie jest wygodne w analizie. Realnie interesuje nas rzeczywista ewolucja struktury i interakcji białko-ligand, a nie obrót i przesunięcie układu, które nie świadczą o rzeczywistej zmianie układu. Także możemy przesunąć (rotacja + translacja) w ten sposób, żeby kolejne klatki były mozliwie podobne do siebie:

```bash 
gmx trjconv -s md_0_200ns.tpr -f md_0_200ns_center.xtc -o md_0_200ns_fit.xtc -fit rot+trans
```
Wybieramy `Backbone` do `fit`owania i `System` jako `output`


#### Energia interakcji


```bash
mkdir ie
cp npt/npt.gro ie
cp npt/npt.cpt ie
cp npt/SYSTEM.top ie
cp npt/index.ndx ie
cp prod/md.mdp ie
cd ie
mv md.mdp ie.mdp
```

Edytujemy plik mdp dodając na końcu sekcji `Output control` linijkę:
```
energygrps = Protein JZ4
```

Przygotowujemy i uruchamiamy obliczenia.

```bash
gmx grompp -f ie.mdp -c npt.gro -t npt.cpt -p SYSTEM.top -n index.ndx -o ie.tpr

nohup gmx mdrun -deffnm ie -rerun ../prod/md_0_200ns.xtc -nb cpu &
```



#### Analiza wiązań wodorowych

Gromacs posiada wbudowane narzędzia do wykrywania wiązań wodorowych. Chcemy...

```bash
gmx hbond-legacy -f md_0_200ns_fit.xtc -s md_0_200ns.tpr -n index.ndx -hbn hbond.ndx -hbm hbond.xpm -g hbond.log
```
Wybieramy dwie grupy: białko i ligand




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