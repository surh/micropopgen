#!/bin/bash

# Usage:
#   $ backup_dir.sh indir <outfile>

# Check directory to compress
if [ -d $1 ]; then
  indir=$1
else
  echo -e "ERROR: directory ($1) does not exist"
  exit 1
fi

# Prepare outfile name
if [ "$2" == "" ]; then
  prefix=`date +%Y-%m-%d`
    suffix=`basename $indir`
  outfile="$prefix.$suffix.tar.bz2"
else
  outfile=$2
fi

# Avoid overwriting files
if [ -f $outfile ]; then
  echo -e "ERROR: outfile ($outfile) already exists"
  exit 2
fi

# echo -e $indir
# echo -e $outfile

# compress
tar -cvvjf $outfile $indir
