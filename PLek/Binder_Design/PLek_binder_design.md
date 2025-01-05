
https://pymolwiki.org/index.php/Create
https://pymolwiki.org/index.php/Cealign

Nałożenie strukturalne w Bioshell?

## Nasz target

https://en.wikipedia.org/wiki/PD-1_and_PD-L1_inhibitors
https://www.youtube.com/watch?v=2v9POD652Ts

https://www.nature.com/articles/s41586-023-06415-8

Interesuje nas binder albo do PD-1, albo do PD-L1, gdzieś na ich interfejsie. Załóżmy, że szukamy bindera do PD-L1. Bierzemy tylko to 

To jest całkiem modny układ. Znaleziono małe cząsteczki, peptydy i białka stanowiące inhibitory tej interakcji. Projektowanie białek wiążących się do PD-L1 jest już dobrze opracowanym tematem; m.in. papery z Baker Labu i Kuhlman Labu...

Paper z Kuhlman Lab o Binder Design Pipeline wprowadzający EvoPro; na podstawie PD-L1 
https://www.biorxiv.org/content/10.1101/2023.05.03.539278v1.full

Można zobaczyć kompleks PD-1/PD-L1 w depozycie 3BIK (Łańcuch A to PD-L1, łańcuch B to PD-1)

Chcemy zobaczyć depozyt 5O45; mamy tam białko PD-L1 (łańcuch A) z peptydowym inhibitorem (łańcuch B)
Depozyt
https://www.rcsb.org/structure/5O45

Dobrze jest spojrzeć na znany binder i wykorzystać...

```bash
wget https://files.rcsb.org/download/5O45.pdb
```

```
sele chain A and resi 56+115+123
set_name sele, hotspots
show sticks, hotspots
```

Tworzymy plik `target_atoms.pdb`, gdzie są tylko rekordy ATOM z łańcucha A.

## RFDiffusion

RFDiffusion Colab
https://colab.research.google.com/github/sokrypton/ColabDesign/blob/main/rf/examples/diffusion.ipynb

Github RFDiffusion, z całkiem dobrym Readme i instrukcjami:
https://github.com/RosettaCommons/RFdiffusion?tab=readme-ov-file#binder-design

#### Contigs

Cytując z dokumentacji RFDiffusion
https://github.com/RosettaCommons/RFdiffusion?tab=readme-ov-file#binder-design
>Anything prefixed by a letter indicates that this is a motif, with the letter corresponding to the chain letter in the input pdb files. E.g. A10-25 pertains to residues ('A',10),('A',11)...('A',25) in the corresponding input pdb
Anything not prefixed by a letter indicates protein to be built. This can be input as a length range. These length ranges are randomly sampled each iteration of RFdiffusion inference.
To specify chain breaks, we use /0 .
In more detail, if we want to scaffold a motif, the input is just like RFjoint Inpainting, except needing to navigate the hydra config input. If we want to scaffold residues 10-25 on chain A a pdb, this would be done with 'contigmap.contigs=[5-15/A10-25/30-40]'. This asks RFdiffusion to build 5-15 residues (randomly sampled at each inference cycle) N-terminally of A10-25 from the input pdb, followed by 30-40 residues (again, randomly sampled) to its C-terminus.

'contigmap.contigs=[B1-100/0 100-100]'
łańcuch A, od 1 do 100 reszty ... długość od 100 do 100

17-134, ale możemy też wziąć wszystko (1-145)
długość bindera 70-100

W colabie inaczej się to zapisuje...

czyli contigs = A1-145:70-100
hotspots A5,A115,A123

'ppi.hotspot_res=[A30,A33,A34]'
#### ... 

Można wybrać tylko do 32 projektów w colabie. W praktyce chcemy dziesiątki tysięcy. Autorzy RFDiff sugerują >=10k, inni 50k
Tak czy inaczej trwa to trochę czasu, także lepiej ograniczyć się do kilku designów.

Pole **pdb** zostawiamy puste, powinien się pojawić prompt do wgrania pliku z dysku. Tam możemy wgrać przygotowane wcześniej `target_atoms.pdb`

RFDiffusion ma silną skłonność do generowania struktur z wieloma helisami. Jeśli chcemy coś dłuższego i więcej wstęg, można zaznaczyć "use_beta_model", czyli model, który powinien preferować $\beta$-kartki. 

### Dodatkowe opcje (niekoniecznie w Colabie)

Auxillary potentials
https://github.com/RosettaCommons/RFdiffusion/blob/main/rfdiffusion/potentials/potentials.py

`monomer_ROG`: radius of gyration 
Dodatkowy potencja ograniczający rozrzucenie atomów daleko od środka masy. 


#### Trajektorie RFDiffusion



## ProteinMPNN

W notebooku z RFDiffusion jest wbudowany ProteinMPNN, ale są też dedykowane notebooki. Ogólnie lepiej jest zrobić to samemu...

ProteinMPNN są najmniej kosztowną częścią obliczeń co do zasady; można je policzyć na CPU i nie trwa to długo.
ProteinMPBB daje nam ileś wyników i metryk wg. których możemy je sortować. Natomiast lepiej niż na samych metrykach proteinmpnn, popatrzeć na wyniki refoldingu w AlphaFold

https://colab.research.google.com/github/dauparas/ProteinMPNN/blob/main/colab_notebooks/quickdemo.ipynb

Należy odznaczyć opcję `homomer`, ustalić designed_chain i fixed_chain

Dostaniemy wyniki w stylu:

```
Generating sequences...
>tmp, score=2.8057, fixed_chains=['A'], designed_chains=['B'], model_name=v_48_020
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
>T=0.1, sample=0, score=1.1183, seq_recovery=0.0366
SAPVLMIKIEPTENREFDLLVSGGKVVAAVQEVSPNNKVKLDPAAIEPLLARAREFLSTVTIDAVLRATETGVEAFPWNPEP
>T=0.1, sample=0, score=1.0830, seq_recovery=0.0610
MKPVLGVIETPTKKLEFRLYVSGGKVVAAFLGYSKNNKVKIDPEALEPLLKKAKEFLTNTKIDEVLVATESGVVSVPFDPNP
>T=0.1, sample=0, score=1.1398, seq_recovery=0.0366
SEPVLAIAELPTKNIETRFLVSGGKVVAVVQAVSKDNKVKLDPEALEPLRARAEELLSGVRIDEVLVATPEKVVRLPFDPAP
>T=0.1, sample=0, score=1.0478, seq_recovery=0.0610
AAPVLLIVEIPTEKLEFRLYVSGGKVVAAAIGYSKKNKVKLDEEAIEPLLAKGKEFLSTVTIDAVLVATETGFVAVPFDPNP
```

Możemy to sobie póki co wrzucić do pliku...
Powtarzamy to dla wszystkich designów
...
Bierzemy num_designs: 4 
Wybierzmy design o najlepszym score (to niekoniecznie najlepszy, ale whatever)


Structures can be aligned using Clustal via the UniProt website
https://www.uniprot.org/align



## Refolding

Podstawowym kolejnym krokiem jest *refolding*, czyli predykcja struktury dla danej sekwencji.

Można teoretycznie przewidzieć cały kompleks, ale zważywszy na to, że strukturę targetu już znamy, chcemy się skupić na strukturze bindera. Interesuje nas wtedy 

To również możemy zrobić w notebooku, tym razem ColabFold.
Mamy kilka do wyboru. Generalnie, ze względu na szybkość polecamy połączenie OpenFold+mmseq...

https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb

RELAXATION!!!
num_relax
ALE to nie ma sensu, bo relaxowanie bez targetu jest mało użyteczne. Dosyć istotne jest bowiem ustawienie reszt na interfejsie target-binder. 

Metryka plDDT

> AlphaFold produces a per-residue estimate of its confidence on a scale from 0 - 100. This confidence measure is called pLDDT and corresponds to the model’s predicted score on the lDDT-Cα metric. It is stored in the B-factor fields of the mmCIF and PDB files available for download (although unlike a B-factor, higher pLDDT is better). pLDDT is also used to colour-code the residues of the model in the 3D structure viewer. The following rules of thumb provide guidance on the expected reliability of a given region:
> - Regions with pLDDT > 90 are expected to be modelled to high accuracy. These should be suitable for any application that benefits from high accuracy (e.g. characterising binding sites).
> - Regions with pLDDT between 70 and 90 are expected to be modelled well (a generally good backbone prediction).
> - Regions with pLDDT between 50 and 70 are low confidence and should be treated with caution.
> - The 3D coordinates of regions with pLDDT < 50 often have a ribbon-like appearance and should not be interpreted. We show in our paper that pLDDT < 50 is a reasonably strong predictor of disorder, i.e. it suggests such a region is either unstructured in physiological conditions or only structured as part of a complex. (Note: this relationship has typically been tested in the context of well-studied proteins, which may have more evolutionarily-related sequences available than a randomly chosen UniProt entry.)

### Nałożenie, RMSD

Mając nowe struktury chcemy je nałożyć na oryginalne i obliczyć RMSD. 

Komendy `align` i `super` wymagają informacji o sekwencji.
Bazują na uliniowieniu sekwencyjnym, które nie jest możliwe, bo dla struktury z RFDiffusion mamy de fakto poliglicynę, czyli nie uliniowimy jej w ten sposób. 

Są też algorytmy nakładające str

W naszym przypadku mamy korespondencję 1-1 dla węgli alfa, czyli optymalnie powinniśmy wykorzystać rototranslację, która da minimalne RMSD.

Algorytmy nałożenia struktur
https://pymolwiki.org/index.php/Kabsch
https://pymolwiki.org/index.php/Pair_fit

RMSD daje nam informacje o zgodności przewi

Ale warto też popatrzeć na kompleks -- czy nie ma jakichś clashów; jak wygląda interfejs. Jeśli oczywiście robimy to dla kilkudziesięciu tysięcy, to nie robimy tego na tym etapie, tylko filtrujemy to, np. na podstawie energii interakcji (np z klasycznej Rosetty, która uwzględni człon energii VdW.

 + pLDDT z AlphaFolda dają informację o stabilności foldu, także generalnie preferujemy sekwencje, które dają struktury z wysokim pLDDT


## Filtrowanie, metryki pipeline

Z RFDiffusion

Z ProteinMPNN

Z https://www.biorxiv.org/content/10.1101/2023.05.03.539278v1.full:
> We orthogonally verified the interface quality of our designs using Rosetta-based scoring metrics. Following FastRelax(27) rotamer optimization and backbone minimization, the InterfaceAnalyzer mover(28) was applied to calculate interface energy between the binder and target protein. The interface quality metric dG_separated/dSASA × 100 (dG/dSASA) was used to score the 6,375 designed AiDs. The score distribution shifts towards lower (i.e., better) energy interfaces after EvoPro optimization, as compared to the interfaces of the starting miniprotein scaffolds (Fig. 1F). This metric was then used to sort and select the 100 best designs. For these, the resulting models from AF2/EvoPro were then compared to predictions made with another DL-based method, OmegaFold(29). These predictions were made with a 28-amino acid flexible linker connecting the C-terminus of HA-PD1 to the N-terminus of the miniprotein to resemble the final autoinhibited complexes we aimed to build. Sequences with similar AF2 and OmegaFold models (as determined by RMSD < 3Å) were further filtered to reduce redundancy of their starting scaffolds to yield a final set of 23 designs for experimental characterization (Table S2). The final TEV protease-cleavable linker for each construct was optimized using OmegaFold to generate the 23 masked antagonist (MA) constructs (Table S3, see Methods).



#### Rosetta Interface Energy

https://docs.rosettacommons.org/docs/latest/application_documentation/analysis/interface-analyzer

---

AF2 pLDDT > 90
RMSD tro diff backbone < 2A
AF2 pae_interface < 4
RF2 pLDDT > 80


https://www.biorxiv.org/content/10.1101/2023.05.03.539278v1.full


## Usprawnienie...

ThermoMPNN


ColabFold all notebooks
https://github.com/sokrypton/ColabFold

AF2
https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb

AF3
https://alphafoldserver.com/


Extra notebooks




ThermoMPNN
https://colab.research.google.com/drive/1OcT4eYwzxUFNlHNPk9_5uvxGNMVg3CFA

DiffDock
https://colab.research.google.com/github/hgbrian/biocolabs/blob/master/DiffDock.ipynb