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
