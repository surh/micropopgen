# Software installed in fraser-server

# AntiSMASH
Istalled as conda environment. Conflicted with something on sur evironment.
Created dedicated module which can be loaded with modules. Current version 4.1

```bash
module load antismash
```

or via anaconda

```bash
module load anaconda/3.6
source activate antismash
```

## CheckM

Python2 program, installed in conda environment

```bash
#module load anaconda
module load fraserconda
source activate python2
module load hmmer
module load prodigal
module load pplacer
```

Data directory is '/godot/users/sur/data/checkm/'

Should load it. Test with:

```bash
checkm test ~/test.out
```

# Clustal Omega

Multiple sequence alignment. Compiled on fraserv

```bash
module load clustalO
```

# EggNOG mapper

Tool for functional annotation. Installed with bacteria, eukaryotic,
human, mammal, and archaeal databases. Data is at
/opt/modules/pkgs/eggnog/1.0.3/data. Can be loaded with:

```bash
module load eggnog
```

# MAFFT
Multiple sequence alignment. Compiled on fraserv

```bash
module load mafft
```

## MIDAS


# Prokka
Installed as conda package. Seemed to conflict with roary but it works
after updating.

```bash
module load anaconda
source activate sur
```

## Roary

```bash
module load anaconda
source activate sur
```
