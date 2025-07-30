"""
Script: classify_crispr_reads.py
Description: This script processes long-read sequencing data from CRISPR-edited samples 
             targeting the HBG1/2 locus (chr11). It classifies each aligned read into 
             categories (WT, 5kb deletion, small indels, truncated, etc.) based on CIGAR 
             string interpretation and positional criteria. Results are exported as 
             individual summary tables (one per sample).
Author: Camille Bergès
Date: 30.07.2025
Dependencies: pysam 0.23.3, pandas 2.3.1, Python  3.13.2
Usage: Run in a directory containing subfolders with sorted BAM files (e.g., barcode01/*.sorted.bam).
"""
*

import pysam
import pandas as pd
from collections import defaultdict
import os
import re

# === CLASSIFICATION PARAMETERS ===

# PCR Long Range HBG1/2 (9.5kb) Locus 
LOCUS_START = 5_244_000
LOCUS_END   = 5_259_000

# --- Deletion boundaries  ---

# Expected CRISPR product (∆5kb)
DEL_CRISPR_START = 5_250_097
DEL_CRISPR_END   = 5_255_025
DEL_MARGIN       = 60  # tolerance

# Alternative 5kb product (as visualized in IGV)
DEL_CRISPR_alt_START = 5_249_028
DEL_CRISPR_alt_END   = 5_253_928
DEL_alt_MARGIN       = 60 # tolerance

# Threshold length
WT_MIN_LEN   = 6000
WT_MAX_LEN   = 10500
DEL_MIN_LEN  = 4600
DEL_MAX_LEN  = 5400
SHORT_READ_LEN = 1800

# Window around cut sites to detect small indels
INDEL_WINDOW = 10  # ±10 nt
CUT_SITES = [5_250_100, 5_255_040]



# === READ CLASSIFICATION FUNCTION ===
def classify_read(read):
    if read.is_unmapped or read.is_secondary:
        return "ignored"

    if read.mapping_quality < 10:
        return "ignored"

    if read.reference_name != "chr11":
        return "ignored"

    ref_start = read.reference_start
    ref_end   = read.reference_end
    align_len = ref_end - ref_start
    qlen      = read.query_length

    # 0. The read must be included in the PCR locus
    if ref_start < LOCUS_START or ref_end > LOCUS_END:
        return "ignored"

    # 1. Very small reads = Artifact
    if qlen < SHORT_READ_LEN:
        return "Artifact"

    # 2. 5kb deletion via CIGAR
    if read.cigartuples:
        ref_pos = ref_start
        for op, length in read.cigartuples:
            if op in (0, 7, 8):  # match / equal / mismatch
                ref_pos += length

            elif op in (2, 3):  # deletion (D) or skipped region (N)
                del_start = ref_pos
                del_end   = ref_pos + length

                # 🔹 fine detection specific for the expected loci
                if 4600 <= length <= 5400:
                    if abs(del_start - DEL_CRISPR_START) <= DEL_MARGIN and \
                       abs(del_end   - DEL_CRISPR_END)   <= DEL_MARGIN:
                        return "del_5kb"
                    if abs(del_start - DEL_CRISPR_alt_START) <= DEL_alt_MARGIN and \
                       abs(del_end   - DEL_CRISPR_alt_END)   <= DEL_alt_MARGIN:
                        return "del_5kb"

                    # 🔸 If loci does not match perfectly but product size does (4928bp +/-5)
                    if 4923 <= length <= 4933:
                        return "del_5kb"

                ref_pos += length

    # 3. Small indels at the cut site
    if read.cigartuples:
        ref_pos = read.reference_start
        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                ref_pos += length

            elif op == 2:  # small deletion
                if length <= 10:
                    for cut in CUT_SITES:
                        if abs(ref_pos - cut) <= INDEL_WINDOW:
                            return "small_indels"
                ref_pos += length

            elif op == 1:  # small insertion
                if length <= 10:
                    for cut in CUT_SITES:
                        if abs(ref_pos - cut) <= INDEL_WINDOW:
                            return "small_indels"
                # # ref_pos remains unchanged here. 
                        # Note: ONT sequencing does not allow reliable detection of small indels due to high background noise.


    # 4. WT: sorted thanks to length
    if WT_MIN_LEN <= align_len <= WT_MAX_LEN:
        return "WT"

    # 5. Truncated reads 
    if 2000 <= align_len < WT_MIN_LEN:
        return "truncated_reads"

    # 6. Unclassified reads
    return "Unclassified_reads"




# === BAM Processing ===
for root, dirs, files in os.walk("."):
    for file in files:
        if file.endswith(".sorted.bam"):
            bam_path = os.path.join(root, file)
            sample = os.path.basename(root)  # barcodeXX
            print(f"📦 Treatment for {sample}...")

            try:
                bamfile = pysam.AlignmentFile(bam_path, "rb")
            except ValueError as e:
                print(f"⚠️  Error, could not open {bam_path} : {e}")
                continue

            counts = defaultdict(int)

            try:
                for read in bamfile.fetch("chr11"):
                    cat = classify_read(read)
                    if cat != "ignored":
                        counts[cat] += 1
            except ValueError:
                print(f"⚠️  'chr11' not found in {file}.")
                continue

            counts["total_reads"] = sum(counts.values())
            counts["sample"] = sample

            #  DataFrame columns
            columns_order = [
                "sample", "total_reads",
                "WT", "truncated_reads",
                "small_indels", "del_5kb",
                "Artifact", "Unclassified_reads"
            ]
            for col in columns_order:
                if col not in counts:
                    counts[col] = 0

            # Individual DataFrame building
            df_indiv = pd.DataFrame([counts])[columns_order]

            # Export CSV in the current directory 
            output_csv = os.path.join(root, f"{sample}_resumed_outcomes_CRISPR.csv")
            df_indiv.to_csv(output_csv, index=False)
            print(f"✅ Resumed outcomes exported in {output_csv}")



