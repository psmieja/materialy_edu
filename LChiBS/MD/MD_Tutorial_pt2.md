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
gmx hbond -f md_0_200ns_fit.xtc -s md_0_200ns.tpr -n index.ndx -hbn hbond.ndx -hbm hbond.xpm -g hbond.log
```
Wybieramy dwie grupy: białko i ligand

*W nowszych wersjach gromacs zamiast komendy hbond należy użyć `hbond_legacy`*


## Przygotowanie filmiku w PyMol

Uruchamiamy w PyMolu plik `SYSTEM.gro` i ładujemy trajektorię

```
PyMol> load_traj [plik z trajektorią]
```

Na tym etapie można odpowiednio ustawić widok, zaznaczyć wiązania, etc.

Jeżeli w pliku z trajektorią mamy np. co setną klatkę, można "wygładzić" ruchy
```
smooth SYSTEM, 30, 3
```

Należy jednak zaznaczyć, że zmienia to pozycje atomów i uzyskana trajektoria nie powinna być wykorzystana do analizy!

Następnie można wykonać ray tracing dla całej trajektorii (Movie -> Ray Trace Frames) i przygotować filmik (File -> Export Movie as -> PNG Images). PyMol pozwala oficjalnie od razu utworzyć plik MP4, ale w praktyce niekoniecznie to działa. Mając natomiast pliki png można je złączyć w filmik programem `ffmpeg`

```
ffmpeg -i *.png -c:v libx264 movie.mp4 
```



## Źródła, linki

Tutorial MD w gromacs
http://www.mdtutorials.com/gmx/complex/index.html

Creating movie in PyMol
https://www.blopig.com/blog/2018/12/turning-md-trajectories-into-movies-using-pymol/