#!/usr/bin/env python3
"""
Off-target анализ гидов siRNA против транскриптома человека.
BLAST-база: e:\datasets\blast_tmp\human_rna
"""

import os
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REPO_ROOT = Path(__file__).resolve().parents[2]

# ---------- 1. НАСТРОЙКИ ----------
BLAST_DB = r"e:\datasets\blast_tmp\human_rna"
OUTDIR = "blast_work"
os.makedirs(OUTDIR, exist_ok=True)

# ---------- 2. ГИДЫ ----------
GUIDES = {
    "DSR-0104292_SSR-0108323": "UAAGACGGCAGAAUGGAAAGAUU",
    "DSR-0104108_SSR-0108095": "UAAGACGCGAGAUGGAAGGAUU",
    "DSR-0104068_SSR-0108000": "UAAGACGCAGAAUGGAAGGAUU",
    "DSR-0104083_SSR-0108024": "UCUCAGACAAAGAUGGCCCAUU",
    "DSR-0104075_SSR-0108016": "UAAGACGCAGAAUGGAAAGAUU",
    "DSR-0104071_SSR-0108012": "UUAGCAAAGCUGUAGACCUCCUU",
    "DSR-0104072_SSR-0108013": "UUAGCAGUAGCAAAGCUGUAGUU",
    "DSR-0104091_SSR-0108032": "UUAAUCGUUGUCUCCCAGUGUU",
    "DSR-0104099_SSR-0108040": "UUGAAGUGGAGUCUGUACAGUU",
}

# ---------- 3. РУЧНОЙ REVCOMP В ДНК ----------
def revcomp_dna(seq_rna: str) -> str:
    """РНК → ДНК → обратный комплемент. Без Bio.Seq, чтобы избежать сюрпризов."""
    seq_dna = seq_rna.upper().replace("U", "T")
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(comp[b] for b in reversed(seq_dna))

# ---------- 4. ЗАПИСЬ FASTA (чистый заголовок!) ----------
def write_guides_fasta(guides: dict, outpath: str):
    records = []
    for name, seq_rna in guides.items():
        sense_dna = revcomp_dna(seq_rna)
        # В description НЕТ нуклеотидных символов!
        records.append(SeqRecord(
            Seq(sense_dna),
            id=name,
            description="siRNA guide sense strand (DNA)"
        ))
    SeqIO.write(records, outpath, "fasta")
    return outpath

guides_fasta = write_guides_fasta(GUIDES, os.path.join(OUTDIR, "guides.fasta"))

# Диагностика: печатаем содержимое FASTA
print("=== FASTA CONTENT ===")
with open(guides_fasta) as f:
    print(f.read())
print("=====================\n")

# ---------- 5. BLAST ----------
def run_blast(query_fasta, db_path, out_tsv):
    result = subprocess.run([
        "blastn",
        "-query", query_fasta,
        "-db", db_path,
        "-out", out_tsv,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                   "qstart qend sstart send evalue bitscore",
        "-word_size", "7",
        "-evalue", "100",          # <-- permissive, чтобы ничего не потерять
        "-max_target_seqs", "500",
        "-num_threads", "4",
    ], check=True, capture_output=True, text=True)
    if result.stderr.strip():
        print(f"[BLAST stderr]\n{result.stderr.strip()[:1000]}\n")
    return out_tsv

offtarget_tsv = run_blast(
    guides_fasta, BLAST_DB,
    os.path.join(OUTDIR, "offtarget.tsv")
)

# ---------- 6. АНАЛИЗ ----------
cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

df = pd.read_csv(offtarget_tsv, sep="\t", names=cols)

df["is_long"] = df["length"] >= 15
df_long = df[df["is_long"]].copy()
df_long["site"] = df_long.apply(
    lambda r: f"{min(r.sstart, r.send)}-{max(r.sstart, r.send)}", axis=1
)

# ---------- 7. СВОДКА ----------
print("\n=== СВОДКА: значимые совпадения (length >= 15 нт) ===\n")
summary = []
for guide_id in GUIDES:
    sub = df_long[df_long["qseqid"] == guide_id]
    inhbe = sub[sub["sseqid"].str.contains("NM_031479", na=False)]
    offtargets = sub[~sub["sseqid"].str.contains("NM_031479", na=False)]
    summary.append({
        "guide": guide_id,
        "total_hits": len(sub),
        "inhbe_hits": len(inhbe),
        "offtarget_hits": len(offtargets),
    })
summary_df = pd.DataFrame(summary)
print(summary_df.to_string(index=False))

# ---------- 8. ДЕТАЛИ OFF-TARGET ----------
print("\n=== ДЕТАЛИ: off-target (не INHBE) ===\n")
for guide_id in GUIDES:
    sub = df_long[
        (df_long["qseqid"] == guide_id) &
        (~df_long["sseqid"].str.contains("NM_031479", na=False))
    ]
    if len(sub) == 0:
        print(f"{guide_id}: off-target не обнаружено")
        continue
    print(f"\n{guide_id}: {len(sub)} off-target")
    for _, row in sub.head(10).iterrows():
        print(f"  {row.sseqid:25s} site={row.site:12s} "
              f"len={row.length} mism={row.mismatch} "
              f"pident={row.pident:.1f}% evalue={row.evalue:.2e}")

# ---------- 9. ЭКСПОРТ ----------
out_xlsx = os.path.join(OUTDIR, "offtarget_analysis.xlsx")
with pd.ExcelWriter(out_xlsx) as writer:
    summary_df.to_excel(writer, sheet_name="Summary", index=False)
    df_long.to_excel(writer, sheet_name="All_hits", index=False)

print(f"\n[+] Результаты сохранены в {out_xlsx}")