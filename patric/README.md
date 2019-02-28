# Interacting with the PATRIC database

Scripts to interact with the [PATRIC database](https://www.patricbrc.org)
and the files downloaded from it.

## Downloading data from PATRIC.

The following coede require a module named fraserconda which contains the
fraserconda conda environment.

* The script `patric_download_genomes.py` connects to the PATRIC ftp site
and downloads genomes whose IDs where given. See full documentation
with:

`patric_download_genomes.py --help`

* The pipeline `patric_download_genomes.nf` splits a table of genomes into
many separate processes and sends them independently to
`patric_download_genomes.py`. Besides the parameters from the python script,
the nextflow pipeline allows you to define the number of genomes per batch
(default 300), the number of simultaneous download processes (default 2,
not recommended to increase), and the waiting time before establishing
a new connection (default 30 seconds). It depends on
[sutilspy](https://wwww.github.com/surh/sutilspy) executable scripts
`split_tables.py` and `cat_tables.py` which are managed via subrepo at
relative path `../sutilspy/bin/` from the pipeline projectDir.

* The script `patric_check_downloads.py` checks that a directory contains an
.fna file and that all the `*.gff` files and `*.features.tab` files contain
only features that exist within the sequences that exist in the `.fna` file.
It can also delete directories that have no .fna file
