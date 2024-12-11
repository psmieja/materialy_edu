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