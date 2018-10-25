# Software installed in fraser-server

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
## MIDAS

## Roary

```bash
module load anaconda
source activate sur
```


# MAFFT
Multiple sequence alignment. Compiled on fraserv

```bash
module load mafft
```

# Clustal Omega

Multiple sequence alignment. Compiled on fraserv

```bash
module load clustalO
```

# Prokka
Installed as conda package. Seemed to conflict with roary but it works
after updating.

```bash
module load anaconda
source activate sur
```
