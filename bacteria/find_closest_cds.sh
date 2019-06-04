#!/usr/bin/env bash

# $ find_closest_cds.sh <features.tab.gz>


if [ ! -f feats.bed ] && [ ! -f closest.txt ];
then
  zcat $1 | awk '($6 == "CDS")' | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t.\t" $5}' > feats.bed
  closestBed -D a -N -s -a feats.bed -b feats.bed > closest.txt
  rm feats.bed
else
  echo "File feats.bed or closest.txt already exists" && exit 1
fi
