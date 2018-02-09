# Pipeline for MKtest comparisons

1. Generate list of samples from output directory. 1 per line. For example

```bash
ll -d /godot/hmp/midas/out/SRS* | awk '{print $9}' | awk 'BEGIN{FS ="/"}{print $6}' > samples.txt
```

data_management/get_sample_map.py
