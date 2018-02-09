# Pipeline for MKtest comparisons

1. Generate list of samples from output directory. 1 per line. For example:

```bash
ll -d /godot/hmp/midas/out/SRS* | awk '{print $9}' | awk 'BEGIN{FS ="/"}{print $6}' > samples.txt
```

2. Get sample body site metadata from database.
Use data_management/get_sample_map.py script. For example:

```bash
 ~/micropopgen/src/micropopgen/data_management/get_sample_map.py --samples samples.txt --db ~/micropopgen/data/database/metagenomes.db --outfile map.txt
 ```
3. Edit the map file from previous step to include headers. There should be
a column named 'ID' which contains the sample ID matching the sample
names in the MIDAS output, and a column named 'Group' indicating the
group (e.g. body site) to which a sample belongs.

4. Run find_comparisons.r. For example:

```bash
~/micropopgen/src/rsutiles/stitch_file.r ~/micropopgen/src/micropopgen/midas/analyze_midas_snps/find_comparisons.r coverage.txt map.txt 3
```
