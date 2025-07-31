# HBG1/2 - ONT-analysis Pipeline

This GitHub repository contains Python scripts for the classification and summarization of long-read sequencing data (e.g., Oxford Nanopore) targeting the **HBG1/2 locus** after CRISPR-Cas9 editing.

The pipeline was developed as part of a functional assessment of deletion and indel outcomes following genome editing at the β-globin locus.

---

## 📂 Repository structure
```
├── barcode01/
│ └── barcode01.sorted.bam
├── barcode02/
│ └── barcode02.sorted.bam
├── global_HBG1-2_analyse_30.07.25.py
├── merge_crispr_results.30.07.25.py
└── Global_resumed_outcomes_CRISPR.csv ← output summary table
...
```

## 🔁 Preprocessing & Alignment

Run the following Bash script to merge and align your `.fastq.gz` files:

```bash
bash Script_EPI2ME_like_31.07.25.sh
```
Requirements: 'minimap2', 'samtools'


## ⚙️ Analyze

- Python ≥ 3.7  
- [pysam](https://pysam.readthedocs.io/en/latest/)  
- [pandas](https://pandas.pydata.org/)
  
Tested with:  
pysam 0.22.0   
pandas 2.3.1  
Python 3.13.2  

1. **global_HBG1-2_analyse_30.07.25.py**
  
This script processes .sorted.bam files aligned to chr11, classifies reads based on:

- 5 kb deletion size and breakpoints
- presence of small indels near cut sites
- truncated or unmapped reads
- overall alignment length

It produces one .csv summary per barcode/sample (e.g., barcode01_resumed_outcomes_CRISPR.csv)

Run from project root:
```bash
python3 global_HBG1-2_analyse_30.07.25.py
```
Note: Each subfolder (e.g., barcode01/) should contain a file named *.sorted.bam.

2. **merge_crispr_results.30.07.25.py**  
  
This script searches all subfolders (starting with barcode) for *_resumed_outcomes_CRISPR.csv files and merges selected columns into a global summary file: Global_resumed_outcomes_CRISPR.csv

Run with:

```bash
python3 merge_crispr_results.30.07.25.py
```

Output:
Each per-sample CSV and the final merged table contain the following columns:
- sample
- total_reads
- WT
- small_indels (≤10bp near cut sites)
- del_5kb (expected large deletions)
- truncated_reads (short reads)
- Artifact (very short / unmapped)
- Unclassified_reads (did not match any category)





## 👩‍🔬 Author

Camille Bergès  
PharmD, PhD candidate  
camille.berges@u-bordeaux.fr  

Date of last update: 30 July 2025

📄 License
MIT License – free to use, adapt, and cite. Please acknowledge the original author in any derivative work.

---
