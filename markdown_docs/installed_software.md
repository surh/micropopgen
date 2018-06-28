# Software installed in fraser-server

## CheckM

Python2 program, installed in conda environment

```bash
module load anaconda
source activate python2
module load hhmer
module load prodigal
module load pplacer 
```

Should load it. Test with:

```bash
checkm test ~/test.out
```

Error
ImportError: /usr/local/anaconda3/envs/python2/lib/python2.7/lib-dynload/../../libz.so.1: version `ZLIB_1.2.9' not found (required by /usr/lib/libpng16.so.16)

## MIDAS
