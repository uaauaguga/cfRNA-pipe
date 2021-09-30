


### Splicing Analysis

```bash
ls output/test_pe/bam/*/genome.sorted.bam | awk 'BEGIN{ORS=","}{print}' | sed 's/,$/\n/'


rmats.py --b1 SRR5687235.txt --gtf reference/gtf/gencode.v38.annotation.gtf -t paired --readLength 50 --nthread 4 --od SRR5687235 --tmp SRR5687235-tmp --task prep --readLength 150 --variable-read-length 
rmats.py --b1 SRR5687235.txt --b2 SRR5687236.txt --gtf reference/gtf/gencode.v38.annotation.gtf -t paired --readLength 150 --nthread 1 --tmp tmp --od test-output --task post
```
