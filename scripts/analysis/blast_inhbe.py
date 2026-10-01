#!/usr/bin/env python3
"""
BLAST-анализ гидов siRNA против мРНК INHBE (человек и мышь).

Требования:
    pip install biopython pandas
    BLAST+ (makeblastdb, blastn) в PATH

Использование:
    python blast_inhbe.py
"""

import os
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REPO_ROOT = Path(__file__).resolve().parents[2]

# ---------- 1. ВХОДНЫЕ ДАННЫЕ ----------

# Последовательности гидов (antisense strands), 5'->3'
GUIDES = {
    "DSR-0104292 (SSR-0108323)": "UAAGACGGCAGAAUGGAAAGAUU",
    "DSR-0104108 (SSR-0108095)": "UAAGACGCGAGAUGGAAGGAUU",
    "DSR-0104068 (SSR-0108000)": "UAAGACGCAGAAUGGAAGGAUU",
    "DSR-0104083 (SSR-0108024)": "UCUCAGACAAAGAUGGCCCAUU",
    "DSR-0104075 (SSR-0108016)": "UAAGACGCAGAAUGGAAAGAUU",
    "DSR-0104071 (SSR-0108012)": "UUAGCAAAGCUGUAGACCUCCUU",
    "DSR-0104072 (SSR-0108013)": "UUAGCAGUAGCAAAGCUGUAGUU",
    "DSR-0104091 (SSR-0108032)": "UUAAUCGUUGUCUCCCAGUGUU",
    "DSR-0104099 (SSR-0108040)": "UUGAAGUGGAGUCUGUACAGUU",
}

# Референсные мРНК (скачать заранее или через Biopython Entrez)
# Человек: NM_031479.4 (INHBE)
# Мышь:   NM_008382.3 (Inhbe)
#REF_HUMAN = "data/NM_031479.4.fasta"
REF_HUMAN = "data/NM_031479.5.fasta"


REF_MOUSE = "data/NM_008382.3.fasta"

WORKDIR = "blast_work"
os.makedirs(WORKDIR, exist_ok=True)

# ---------- 2. СКАЧИВАНИЕ РЕФЕРЕНСОВ (если нет локально) ----------

def download_refseq(accession: str, outpath: str):
    """Скачивает FASTA из NCBI, если файла нет."""
    if os.path.exists(outpath):
        print(f"[+] {outpath} уже существует")
        return
    from Bio import Entrez
    Entrez.email = "your.email@example.com"  # <-- ЗАМЕНИТЕ на свой email
    print(f"[*] Скачиваю {accession}...")
    with Entrez.efetch(db="nucleotide", id=accession,
                       rettype="fasta", retmode="text") as handle:
        data = handle.read()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, "w") as f:
        f.write(data)
    print(f"[+] Сохранено в {outpath}")

download_refseq("NM_031479.5", REF_HUMAN) #download_refseq("NM_031479.4", REF_HUMAN)
download_refseq("NM_008382.3", REF_MOUSE)

# ---------- 3. ПОДГОТОВКА BLAST-БАЗ ----------

def make_blast_db(fasta_path: str, db_name: str):
    """Создаёт BLAST-базу из FASTA."""
    db_path = os.path.join(WORKDIR, db_name)
    subprocess.run(
        ["makeblastdb", "-in", fasta_path, "-dbtype", "nucl",
         "-out", db_path, "-parse_seqids"],
        check=True, capture_output=True
    )
    return db_path

db_human = make_blast_db(REF_HUMAN, "inhbe_human")
db_mouse = make_blast_db(REF_MOUSE, "inhbe_mouse")

# ---------- 4. ЗАПИСЬ ГИДОВ В FASTA ----------

def write_guides_fasta(guides: dict, outpath: str):
    records = []
    for name, seq in guides.items():
        # Гид — антисмысловая цепь. Для поиска сайта в мРНК
        # нужен обратный комплемент гида (то есть смысловая цепь).
        sense = str(Seq(seq).reverse_complement())
        records.append(SeqRecord(Seq(sense), id=name.replace(" ", "_"),
                                 description=f"guide={seq}; sense={sense}"))
    SeqIO.write(records, outpath, "fasta")
    return outpath

guides_fasta = write_guides_fasta(GUIDES, os.path.join(WORKDIR, "guides.fasta"))

# ---------- 5. ЗАПУСК BLAST ----------

def run_blast(query_fasta: str, db_path: str, out_tsv: str):
    """Запускает blastn с выводом в TSV."""
    subprocess.run([
        "blastn",
        "-query", query_fasta,
        "-db", db_path,
        "-out", out_tsv,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                   "qstart qend sstart send evalue bitscore",
        "-word_size", "7",
        "-evalue", "10",
        "-strand", "plus",
        "-max_target_seqs", "10"
    ], check=True, capture_output=True)
    return out_tsv

cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

human_tsv = run_blast(guides_fasta, db_human,
                      os.path.join(WORKDIR, "human.tsv"))
mouse_tsv = run_blast(guides_fasta, db_mouse,
                      os.path.join(WORKDIR, "mouse.tsv"))

# ---------- 6. АНАЛИЗ РЕЗУЛЬТАТОВ ----------

def parse_blast(tsv_path: str) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t", names=cols)
    # Извлекаем позицию связывания (1-based) для плюс-цепи мРНК
    df["site"] = df.apply(
        lambda r: f"{min(r.sstart, r.send)}-{max(r.sstart, r.send)}", axis=1
    )
    return df

human_df = parse_blast(human_tsv)
mouse_df = parse_blast(mouse_tsv)

print("\n=== BLAST против INHBE человека (NM_031479.5) ===")
print(human_df[["qseqid", "site", "pident", "length", "mismatch",
                "evalue"]].to_string(index=False))

print("\n=== BLAST против INHBE мыши (NM_008382.3) ===")
print(mouse_df[["qseqid", "site", "pident", "length", "mismatch",
                "evalue"]].to_string(index=False))

# ---------- 7. OFF-TARGET АНАЛИЗ (полный транскриптом) ----------
# Для полного off-target анализа нужен BLAST против всего RefSeq.
# Это делается отдельно: скачать human.rna.fna (~100 МБ),
# сделать makeblastdb и запустить blastn с -max_target_seqs 1000.

print("\n[i] Для off-target анализа скачайте RefSeq RNA:")
print("    https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RNA/")
print("    и запустите blastn против всей базы.")