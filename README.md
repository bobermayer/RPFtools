#RPFtools

RPFtools is a collection of scripts for the analysis of ribosome profiling data.

## Prerequisites

These tools run on Python 2.7.11 with numpy, scipy, pandas, twobitreader, and pysam.

## Usage

### convert GTF to bed file
```
python gtf_to_bed.py -i annotation.gtf(.gz) -a annotation.tsv > annotation.bed
```
produces a 12-column bed file from the gtf (e.g., from Gencode or Ensembl) with all transcripts and their ORFs (if present), as well as an annotation table containing geneIDs and gene/transcript biotypes

### get longest ORF per gene
if only one ORF (the longest) per gene is required, this tool creates a reduced bed file with the geneID in the 3rd column
```
python get_longest_ORF_per_gene.py -b annotation.bed -a annotation.tsv -o annotation.collapsed.bed
```

### find all ORFs in all transcripts
```
python get_all_ORFs.py -b annotation.bed -G genome.2bit -s orf_stats.csv -o annotation.all_ORFs.bed
```
given a bed file with transcript definitions and a 2bit file with the genome, finds all ORFs of minimum length ``--minlength`` [12] and writes a (large) bed file with their coordinates. the 3rd column now contains an ID combined from transcript ID, coordinates, and a hash value computed from the ORF sequence.
``orf_stats.csv`` is a table with statistics (ORF and UTR lengths for each ORF)

### analyze phasing of ribosome profiling data
```
python get_RPF_phasing.py -b annotation.collapsed.bed -B RPF.bam -o RPF_phasing.pdf 
```
given a bam file with mapped reads, collects read counts in all ORFs of the given bed file and checks for frame bias, assuming P-sites 12nt downstream of 5'end (can be changed using ``--offset``)

### analyze coverage profiles around start and stop codons
```
python get_RPF_profiles.py -b annotation.collapsed.bed -B RPF.bam -L "29,30" -o "12,12" -N 150 > RPF_profiles.out
```
takes all reads of length 29 and 30nt with 12nt offset over the ORFs of the bed file and outputs for each gene the P-site density (relative to all reads in that gene) in 2x150nt windows around start and stop codon. These can then be averaged over genes for meta-gene plots.

### calculate ORFscores
```
python get_ORFscores.py -b annotation.all_ORFs.bed -B RPF.bam -L "29,30" -o "12,12" > ORFscores.out
```
calculates ORFscores [Bazzini et al.](http://emboj.embopress.org/content/33/9/981.long) for each ORF in the input bed file, using reads of length 29 and 30nt with 12nt offset. Output contains reads per million mapped over transcript and ORF, ORFscore, % of p0 positions covered and fraction of multimappers contributing to these scores. If multiple bam files are given using ``-B RPF_1.bam,RPF_2.bam``, scores are also calculated using pooled reads (use ``-L "29,30|29,30"`` and ``-o "12,12|12,12" to specify lengths and offsets).

### analyze codon coverage
```
python get_RPF_codon_counts.py -b annotation.collapsed.bed -B RPF.bam -L "29,30" -o "12,12" -G genome.2bit > RPF_codon_counts.out
```
aggregates for each ORF in the input bed file all P-site counts (for A-site counts, use ``-o 15``) over all codons. First and last codons can be excluded using ``-e`` (default: 25,1). 
