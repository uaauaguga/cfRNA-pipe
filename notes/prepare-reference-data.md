## Reference data preparation

### Reference genome and annotation

- reference genome `reference/fasta/hg38.fa`
```bash
wget -P reference/source http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
zcat reference/source/GRCh38.primary_assembly.genome.fa.gz >  reference/fasta/hg38.fa
samtools faidx reference/fasta/hg38.fa
```

- gencode annotations

```bash
wget  -P reference/source  http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
mkdir reference/gtf
zcat reference/source/gencode.v38.annotation.gtf.gz > reference/gtf/gencode.v38.annotation.gtf 
mkdir reference/bed
gffread reference/gtf/gencode.v38.annotation.gtf --bed -o reference/bed/gencode.v38.annotation.bed
run/extract-3putr.py -b reference/bed/gencode.v38.annotation.bed -u reference/bed/gencode.v38.annotation.3putr.bed
bedtools subtract -a reference/bed/gencode.v38.annotation.3putr.bed -b reference/bed/gencode.v38.annotation.3putr.bed -S > reference/bed/gencode.v38.annotation.3putr.subtracted.bed
run/select-3putr.py -i reference/bed/gencode.v38.annotation.3putr.subtracted.bed -o reference/bed/gencode.v38.annotation.3putr.subtracted.selected.bed
```

- Prepare 3p UTR annotation
```bash
bedtools sort -faidx reference/fasta/hg38.fa.fai -i test.3putr.subtracted.selected.bed > test.3putr.subtracted.selected.sorted.bed
```


- Build STAR index

```bash
STAR --runMode genomeGenerate --genomeDir reference/star-index/hg38 --genomeFastaFiles reference/fasta/hg38.fa --sjdbGTFfile reference/gtf/gencode.v38.annotation.gtf --runThreadN 8
```


- rRNA sequences `reference/fasta/rRNA.fa`

```bash
## Download ncbi annotations (to obtain rRNA sequence ids)
wget -P reference/source https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
## Chromsome name mapping between Refseq and gencode
wget -P reference/source https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_RefSeq2UCSC.txt
wget -P reference/source https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_UCSC2gencode.txt
## Get rRNA sequences
run/get-rRNA-sequences.py
```

- circRNA sequences `reference/fasta/circRNA.fa`

```bash
# Download putative spliced sequence from circBase
wget -P reference/source http://www.circbase.org/download/human_hg19_circRNAs_putative_spliced_sequence.fa.gz
# Extract back-spliced junction sequences
run/get-circRNA-junctions.py
samtools faidx reference/fasta/circRNA.fa
```

- UniVec sequence

```bash
wget -P reference/source https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec 
cp reference/source/UniVec  reference/fasta/UniVec.fa
samtools faidx  reference/fasta/UniVec.fa
```

- ERCC splike in sequence

```bash
wget -P reference/source https://assets.thermofisher.cn/TFS-Assets/LSG/manuals/cms_095047.txt 
cat reference/source/cms_095047.txt | awk 'BEGIN{FS="\t"}{print ">"$1"|"$2"\n"$5}' > reference/fasta/spikein.fa
```


- Known RNA editing sites

```bash
wget -P reference/source http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg38.txt.gz
```


- kraken2's database

```bash
wget  -P reference/source ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz
tar xvzf reference/source/minikraken2_v2_8GB_201904.tgz -C reference/kraken2db
mv reference/kraken2db/minikraken2_v2_8GB_201904_UPDATE reference/kraken2db/minikraken2_v2
```
