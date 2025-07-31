#!/bin/bash

"""
Script: Script_EPI2ME_like_31.07.25.sh
Description: Merge and align fastq files from ONT sequencing to obtain BAM ; script based on the Epi2Me (v5.1.9) workflow and the wf-alignment pipeline (v1.1.2).
Author: Camille Bergès
Date: 31.07.2025
Dependencies: Bash (v5.1.16), , minimap2, samtools
Usage: Run Usage: Run from the root directory containing subfolders (e.g., barcode01/)
       with raw sequencing files in fastq.gz format.
"""


REF="combined_refs.mmi" # pre-indexed reference with minimap2
THREADS=6  # ← adaptable

for barcode_dir in barcode*/; do
    barcode_name=$(basename "$barcode_dir")
    merged_fastq="${barcode_name}_merged.fastq.gz"
    sorted_bam="${barcode_name}.sorted.bam"

    echo "▶️ Traitement de $barcode_name..."

    # Merging of fastq.gz files
    cat "$barcode_dir"/*.fastq.gz > "$merged_fastq"

    # Alignement with minimap2, classificiation and indexation
    minimap2 -t 16 -ax map-ont --cap-kalloc 100m --cap-sw-mem 50m "$REF" "$merged_fastq" \
    | samtools sort -@ 16 --write-index -o "$sorted_bam"

    echo "✅ Done for $barcode_name : $sorted_bam"
done

echo "🏁 All barcodes processed"
