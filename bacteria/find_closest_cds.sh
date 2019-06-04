#!/usr/bin/env bash
// Copyright (C) 2019 Sur Herrera Paredes

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


# $ find_closest_cds.sh <features.tab.gz>


if [ ! -f feats.bed ] && [ ! -f closest.txt ];
then
  zcat $1 | awk '($6 == "CDS")' | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t.\t" $5}' > feats.bed
  closestBed -D a -N -s -a feats.bed -b feats.bed > closest.txt
  rm feats.bed
else
  echo "File feats.bed or closest.txt already exists" && exit 1
fi
