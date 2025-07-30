"""
Script: merge_crispr_results.py
Description: Collects per-sample classification results from subfolders (barcodeXX),
             keeps selected columns, and concatenates them into a global summary CSV.
Author: Camille Bergès
Date: 30.07.2025
Dependencies: pandas 2.3.1, Python  3.13.2, os (standard library)
Usage: Run Usage: Run from the root directory containing subfolders (e.g., barcode01/)
       with files named "_resumed_outcomes_CRISPR.csv".
"""

import os
import pandas as pd

# Columns of interest
columns_to_keep = ["sample", "total_reads", "WT", "small_indels", "del_5kb"]

# Dataframes
all_data = []

# Run files
for folder in os.listdir():
    folder_path = os.path.join(".", folder)
    if os.path.isdir(folder_path) and folder.startswith("barcode"):
        for file in os.listdir(folder_path):
            if file.endswith("_resumed_outcomes_CRISPR.csv"):
                csv_path = os.path.join(folder_path, file)
                try:
                    df = pd.read_csv(csv_path)
                    df = df[columns_to_keep]
                    all_data.append(df)
                    print(f"📥 Imported file for {df['sample'].iloc[0]} : {file}")
                except Exception as e:
                    print(f"⚠️ Erreur avec {csv_path} : {e}")

# Concatenation
if all_data:
    df_final = pd.concat(all_data, ignore_index=True)
    df_final.to_csv("_resumed_outcomes_CRISPR.csv", index=False)
    print("\n✅ Global file exported: Global_resumed_outcomes_CRISPR.csv")
else:
    print("❌ No .csv file found to concatenate.")

