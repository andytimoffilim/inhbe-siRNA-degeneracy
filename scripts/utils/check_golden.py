# check_golden.py
import pandas as pd
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]

df = pd.read_csv(r"E:\datasets\mRNA\pcsk9_sirna_golden.csv")
print(f"Shape: {df.shape}")
print(f"Columns: {df.columns.tolist()}")

print("\n--- IC50 stats ---")
print(df['ic50_pm'].describe())
print(f"log10(IC50): mean={np.log10(df['ic50_pm']).mean():.2f}, "
      f"std={np.log10(df['ic50_pm']).std():.2f}")

print("\n--- Site coverage ---")
print(f"site_start range: {df['site_start'].min()}–{df['site_start'].max()}")
print(f"Missing sites: {df['site_start'].isna().sum()}")

print("\n--- Clusters ---")
cluster_A = df[(df['site_start'] >= 3520) & (df['site_start'] <= 3570)]
cluster_B = df[(df['site_start'] >= 3472) & (df['site_start'] <= 3490)]
outside = df[~df.index.isin(cluster_A.index) & ~df.index.isin(cluster_B.index)]
print(f"A (3520-3570): n={len(cluster_A)}, median IC50={cluster_A['ic50_pm'].median():.1f} pM")
print(f"B (3472-3490): n={len(cluster_B)}, median IC50={cluster_B['ic50_pm'].median():.1f} pM")
print(f"Outside:       n={len(outside)}, median IC50={outside['ic50_pm'].median():.1f} pM")

print("\n--- Sequences ---")
print(f"Missing sense_seq: {df['sense_seq'].isna().sum()}")
print(f"Missing antisense_seq: {df['antisense_seq'].isna().sum()}")
print(df[['duplex_name', 'site_start', 'site_end',
          'sense_seq', 'antisense_seq', 'ic50_pm']].head(10))