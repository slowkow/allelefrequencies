---
title: "Allele Frequencies Net Database"
author: "Kamil Slowikowski"
date: "2023-03-14"
output: github_document
---

# Allele Frequencies Net Database

Here, we share a script `allelefrequencies.py` that downloads allele
frequencies for HLA, KIR, MIC, and cytokine genes from this website:

- http://allelefrequencies.net

Then, the script writes .tsv files like this:

```
ls allelefrequencies.net/hla/
A.tsv
B.tsv
C.tsv
DPA1.tsv
DPB1.tsv
DQA1.tsv
DQB1.tsv
DRB1.tsv
```

Each file looks like this:

```
head allelefrequencies.net/A.tsv

allele   population                               allele_freq  sample_size
A*01:01  Argentina Rosario Toba                   0.0760       86
A*01:01  Armenia combined Regions                 0.1250       100
A*01:01  Australia Cape York Peninsula Aborigine  0.0530       103
A*01:01  Australia Groote Eylandt Aborigine       0.0270       75
A*01:01  Australia New South Wales Caucasian      0.1870       134
A*01:01  Australia Yuendumu Aborigine             0.0080       191
A*01:01  Austria                                  0.1460       200
A*01:01  Azores Central Islands                   0.0800       59
A*01:01  Azores Oriental Islands                  0.1150       43
```

Please cite the latest manuscript about Allele Frequency Net Database:

Gonzalez-Galarza FF, McCabe A, Santos EJMD, Jones J, Takeshita L, Ortega-Rivera
ND, et al. [Allele frequency net database (AFND) 2020 update: gold-standard
data classification, open access genotype data and new query
tools.](https://pubmed.ncbi.nlm.nih.gov/31722398) Nucleic Acids Res. 2020;48:
D783–D788. doi:10.1093/nar/gkz1029

And thanks to David A. Wells for sharing
[scrapeAF](https://github.com/DAWells/scrapeAF), which I forked to create this
script.
