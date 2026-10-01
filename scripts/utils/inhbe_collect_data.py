#!/usr/bin/env python3
"""
inhbe_collect_data.py
Сбор golden-подмножества INHBE-siRNA для переноса модели с PCSK9.

Источники:
  - WO 2026/055502 A3 (PDF, если есть)
  - Fallback: known guides из отчёта Wave Life Sciences

Выход:
  - inhbe_golden.csv    (dsr_name, sense, antisense, site, ic50_pm)
  - inhbe_helm_raw.json (сырые HELM-строки)
  - inhbe_blast.tsv     (BLAST-результаты)

Требования:
  pip install biopython pandas pdfplumber
  BLAST+ (makeblastdb, blastn) в PATH
"""

import os
import re
import csv
import json
import subprocess
from pathlib import Path
from collections import defaultdict

REPO_ROOT = Path(__file__).resolve().parents[2]

# ============================================================
# КОНФИГ
# ============================================================

BASE       = REPO_ROOT
OUT_DIR    = BASE / "data"
WORK_DIR   = OUT_DIR / "blast_work"
DATA_DIR   = BASE / "data"

OUT_DIR.mkdir(parents=True, exist_ok=True)
WORK_DIR.mkdir(parents=True, exist_ok=True)
# DATA_DIR уже существует — не создаём, чтобы не перезаписать

PATENT_PDF   = BASE / "WO2026055502A3.pdf"
PATENT_TXT   = BASE / "WO2026055502A3.txt"    # если есть OCR-текст
REF_HUMAN    = DATA_DIR / "NM_031479.5.fasta"
REF_MOUSE    = DATA_DIR / "NM_008382.3.fasta"

# Страницы патента, где лежат HELM и IC50 (из отчёта)
HELM_PAGES   = [23, 24, 25, 65]
IC50_PAGES   = list(range(230, 245)) + list(range(279, 290))

# ============================================================
# FALLBACK: known data из отчёта
# ============================================================

KNOWN_GUIDES = {
    "DSR-0104292": "UAAGACGGCAGAAUGGAAAGAUU",  # SSR-0108323
    "DSR-0104108": "UAAGACGCGAGAUGGAAGGAUU",   # SSR-0108095
    "DSR-0104068": "UAAGACGCAGAAUGGAAGGAUU",   # SSR-0108000
    "DSR-0104083": "UCUCAGACAAAGAUGGCCCAUU",   # SSR-0108024
    "DSR-0104075": "UAAGACGCAGAAUGGAAAGAUU",   # SSR-0108016
    "DSR-0104071": "UUAGCAAAGCUGUAGACCUCCUU",  # SSR-0108012
    "DSR-0104072": "UUAGCAGUAGCAAAGCUGUAGUU",  # SSR-0108013
    "DSR-0104091": "UUAAUCGUUGUCUCCCAGUGUU",   # SSR-0108032
    "DSR-0104099": "UUGAAGUGGAGUCUGUACAGUU",   # SSR-0108040
}

KNOWN_HELM = {
    "SSR-0108051": (
        "RNA1{p.m(U)[Ssp].m(C)p.m(U)p.m(U)p.m(U)p.m(C)p."
        "[fl2r](C)p.m(A)p.[fl2r](U)p.[fl2r](U)p.[fl2r](C)p."
        "m(U)p.m(G)p.m(C)p.m(C)p.m(G)p.m(U)p.m(C)p.m(U)p.m(U)"
        "[Ssp].m(A)}"
    ),
    "SSR-0108323": (
        "RNA2{[d5m](U)[Ssp].[fl2r](A)[Rsp].m(A)[n001S]."
        "[fl2r](G)p.m(A)p.[fl2r](C)p.m(G)p.m(G)p.m(C)p.m(A)"
        "[n001R].[fl2r](G)p.m(A)p.m(A)p.[fl2r](U)p.m(G)p."
        "[fl2r](G)p.m(A)p.m(A)p.m(A)p.m(G)p.m(A)[Ssp].m(U)"
        "[Ssp].m(U)}"
    ),
}

KNOWN_SITES = {
    "DSR-0104292": (1107, 1126),
    "DSR-0104075": (1107, 1126),
    "DSR-0104068": (1107, 1126),
    "DSR-0104072": (509, 522),
    "DSR-0104091": (1857, 1869),
    "DSR-0104083": (1432, 1454),
    "DSR-0104071": (499, 519),
    "DSR-0104099": (519, 539),
}

# ============================================================
# УТИЛИТЫ
# ============================================================

def read_fasta(path):
    if not Path(path).exists():
        return None
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def helm_to_rna_seq(helm_str):
    """Извлекает последовательность нуклеотидов из HELM-строки."""
    m = re.search(r"RNA\d*\{(.+?)\}", helm_str)
    if not m:
        return None
    body = m.group(1)

    seq = []
    i = 0
    while i < len(body):
        c = body[i]
        if c == "[":
            j = body.find("]", i)
            i = (j + 1) if j != -1 else i + 1
            continue
        if c == "(":
            j = body.find(")", i)
            i = (j + 1) if j != -1 else i + 1
            continue
        if c in "AUGC":
            seq.append(c)
        i += 1
    return "".join(seq) if seq else None


def rna_to_dna(seq_rna):
    return (seq_rna or "").upper().replace("U", "T")


def revcomp_dna(seq_dna):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq_dna.upper()))


def find_site(mrna, seq_dna, min_len=17):
    """Ищет seq_dna в мРНК, при неудаче укорачивает с 3'-конца."""
    seq_dna = re.sub(r"[^ACGT]", "", seq_dna.upper())
    if len(seq_dna) < min_len:
        return None
    for trim in range(0, len(seq_dna) - min_len + 1):
        sub = seq_dna[: len(seq_dna) - trim] if trim else seq_dna
        p = mrna.find(sub)
        if p != -1:
            return (p + 1, p + len(sub), len(sub), trim)
    return None


# ============================================================
# ПАРСИНГ PDF (HELM + IC50)
# ============================================================

def extract_helm_from_pdf(pdf_path, pages):
    try:
        import pdfplumber
    except ImportError:
        print("[!] pdfplumber не установлен — PDF-парсинг пропущен")
        return {}

    if not Path(pdf_path).exists():
        print(f"[!] PDF не найден: {pdf_path}")
        return {}

    helm_blocks = {}
    with pdfplumber.open(pdf_path) as pdf:
        for pnum in pages:
            if pnum > len(pdf.pages):
                continue
            text = pdf.pages[pnum - 1].extract_text() or ""
            for m in re.finditer(
                r"(DSR-\d+|SSR-\d+)[^\n]{0,80}?(RNA\d?\{[^\n]+\})",
                text,
            ):
                name = m.group(1)
                helm = m.group(2)
                helm_blocks.setdefault(name, []).append(helm)
    return helm_blocks


def extract_ic50_from_pdf(pdf_path, pages):
    try:
        import pdfplumber
    except ImportError:
        return {}

    if not Path(pdf_path).exists():
        return {}

    ic50 = defaultdict(list)
    with pdfplumber.open(pdf_path) as pdf:
        for pnum in pages:
            if pnum > len(pdf.pages):
                continue
            text = pdf.pages[pnum - 1].extract_text() or ""
            for m in re.finditer(
                r"(DSR-\d+)[^\n]{0,200}?"
                r"(?:IC50|ED50|IC\s*50)[^\d\-]{0,20}"
                r"([\d.]+)\s*(pM|nM|µM|uM|M|mg/kg)",
                text,
                re.IGNORECASE,
            ):
                ic50[m.group(1)].append({
                    "value": float(m.group(2)),
                    "unit": m.group(3),
                    "page": pnum,
                })
    return dict(ic50)


def ic50_to_pm(records):
    """Среднее по всем записям, приведённое к пМ."""
    if not records:
        return None
    pm_vals = []
    for r in records:
        v, u = r["value"], r["unit"].lower()
        if u == "pm":
            pm_vals.append(v)
        elif u == "nm":
            pm_vals.append(v * 1000)
        elif u in ("µm", "um", "mkm"):
            pm_vals.append(v * 1e6)
    return sum(pm_vals) / len(pm_vals) if pm_vals else None


# ============================================================
# BLAST
# ============================================================

def ensure_blast_db(fasta_path, db_name):
    db_path = WORK_DIR / db_name
    if not (Path(str(db_path) + ".nin").exists() or
            Path(str(db_path) + ".nsq").exists()):
        subprocess.run(
            ["makeblastdb", "-in", str(fasta_path), "-dbtype", "nucl",
             "-out", str(db_path), "-parse_seqids"],
            check=True, capture_output=True,
        )
    return str(db_path)


def write_guides_fasta(guides, outpath):
    """guides: {name: rna_seq}. Записывает sense (revcomp гида без overhang)."""
    with open(outpath, "w") as f:
        for name, rna in guides.items():
            guide_dna = rna_to_dna(rna)
            core = guide_dna[:-2] if guide_dna.endswith("TT") else guide_dna
            sense_dna = revcomp_dna(core)
            safe_name = re.sub(r"[^\w\-]", "_", name)
            f.write(f">{safe_name}\n{sense_dna}\n")
    return outpath


def run_blast(query_fasta, db_path, out_tsv, evalue=1.0):
    subprocess.run([
        "blastn", "-query", str(query_fasta), "-db", str(db_path),
        "-out", str(out_tsv),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                   "qstart qend sstart send evalue bitscore",
        "-word_size", "7", "-evalue", str(evalue),
        "-max_target_seqs", "50",
    ], check=True, capture_output=True)
    return out_tsv


def parse_blast_tsv(tsv_path):
    rows = []
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch",
            "gapopen", "qstart", "qend", "sstart", "send",
            "evalue", "bitscore"]
    with open(tsv_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == len(cols):
                rows.append(dict(zip(cols, parts)))
    return rows


# ============================================================
# СБОРКА
# ============================================================

def build_golden():
    mrna = read_fasta(REF_HUMAN)
    if mrna:
        print(f"[+] NM_031479.5: {len(mrna)} nt")
    else:
        print(f"[!] Нет {REF_HUMAN}")
        print(f"    Скачайте: https://www.ncbi.nlm.nih.gov/nuccore/NM_031479.5")
        print(f"    и положите в {DATA_DIR}")
        mrna = ""

    mrna_mouse = read_fasta(REF_MOUSE)
    if mrna_mouse:
        print(f"[+] NM_008382.3: {len(mrna_mouse)} nt")

    # --- 1. Пытаемся извлечь из PDF ---
    helm_from_pdf = extract_helm_from_pdf(PATENT_PDF, HELM_PAGES)
    ic50_from_pdf = extract_ic50_from_pdf(PATENT_PDF, IC50_PAGES)
    print(f"[+] HELM из PDF:  {len(helm_from_pdf)} блоков")
    print(f"[+] IC50 из PDF:  {len(ic50_from_pdf)} записей")

    # Сохраняем сырые HELM
    raw_helm_path = OUT_DIR / "inhbe_helm_raw.json"
    with open(raw_helm_path, "w", encoding="utf-8") as f:
        json.dump({
            "from_pdf": helm_from_pdf,
            "fallback": KNOWN_HELM,
        }, f, ensure_ascii=False, indent=2)
    print(f"[+] Сырые HELM -> {raw_helm_path}")

    # --- 2. Собираем golden ---
    golden = []
    for dsr_name, guide_rna in KNOWN_GUIDES.items():
        guide_dna = rna_to_dna(guide_rna)
        core = guide_dna[:-2] if guide_dna.endswith("TT") else guide_dna
        sense_dna = revcomp_dna(core)

        # Site: сначала из мРНК, потом fallback
        site = None
        if mrna:
            site = find_site(mrna, sense_dna)
        if not site and dsr_name in KNOWN_SITES:
            s, e = KNOWN_SITES[dsr_name]
            site = (s, e, e - s + 1, 0)

        # IC50
        ic50_pm = None
        if dsr_name in ic50_from_pdf:
            ic50_pm = ic50_to_pm(ic50_from_pdf[dsr_name])

        golden.append({
            "dsr_name":       dsr_name,
            "sense_seq":      sense_dna,
            "antisense_seq":  guide_dna,
            "sense_len":      len(sense_dna),
            "antisense_len":  len(guide_dna),
            "site_start":     site[0] if site else "",
            "site_end":       site[1] if site else "",
            "site_length":    site[2] if site else "",
            "ic50_pm":        round(ic50_pm, 3) if ic50_pm else "",
            "source":         "WO 2026/055502 A3",
            "notes":          "pdf" if helm_from_pdf else "fallback_from_report",
        })

    # --- 3. BLAST-верификация ---
    if (REF_HUMAN.exists() or REF_MOUSE.exists()):
        print("\n[BLAST] Подготовка...")
        guides_fasta = WORK_DIR / "guides.fasta"
        write_guides_fasta(KNOWN_GUIDES, guides_fasta)

        blast_results = []
        if REF_HUMAN.exists():
            db_h = ensure_blast_db(REF_HUMAN, "inhbe_human")
            tsv_h = run_blast(guides_fasta, db_h, WORK_DIR / "human.tsv")
            for r in parse_blast_tsv(tsv_h):
                r["reference"] = "human"
                blast_results.append(r)
        if REF_MOUSE.exists():
            db_m = ensure_blast_db(REF_MOUSE, "inhbe_mouse")
            tsv_m = run_blast(guides_fasta, db_m, WORK_DIR / "mouse.tsv")
            for r in parse_blast_tsv(tsv_m):
                r["reference"] = "mouse"
                blast_results.append(r)

        if blast_results:
            blast_csv = OUT_DIR / "inhbe_blast.tsv"
            with open(blast_csv, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(blast_results[0].keys()))
                w.writeheader()
                w.writerows(blast_results)
            print(f"[+] BLAST -> {blast_csv} ({len(blast_results)} hits)")

            # Обновляем site из BLAST для человека
            for r in blast_results:
                if r["reference"] != "human":
                    continue
                if not r["sseqid"].startswith("NM_031479"):
                    continue
                for g in golden:
                    if g["dsr_name"] == r["qseqid"]:
                        s = min(int(r["sstart"]), int(r["send"]))
                        e = max(int(r["sstart"]), int(r["send"]))
                        g["site_start"] = s
                        g["site_end"] = e
                        g["site_length"] = int(r["length"])
                        if "blast_verified" not in g["notes"]:
                            g["notes"] = (g["notes"] + ";blast_verified").strip(";")

    # --- 4. Запись golden ---
    out_csv = OUT_DIR / "inhbe_golden.csv"
    fields = ["dsr_name", "sense_seq", "antisense_seq",
              "sense_len", "antisense_len",
              "site_start", "site_end", "site_length",
              "ic50_pm", "source", "notes"]
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(golden)

    # --- 5. Сводка ---
    n_site = sum(1 for r in golden if r["site_start"])
    n_ic50 = sum(1 for r in golden if r["ic50_pm"])
    print(f"\n=== INHBE GOLDEN ===")
    print(f"Всего дуплексов:        {len(golden)}")
    print(f"С сайтом в NM_031479.5: {n_site}")
    print(f"С IC50 (пМ):            {n_ic50}")
    print(f"\nCSV: {out_csv}")

    if n_ic50 == 0:
        print(f"\n[!] IC50 не извлечены. Варианты:")
        print(f"    1. Скачать PDF патента WO 2026/055502 A3 в {BASE}")
        print(f"    2. Если PDF — скан, прогнать OCR (см. IC50_probe_scans.py)")
        print(f"    3. Использовать клинические PD-данные Phase 1 как прокси-метки:")
        print(f"       75 мг -> 61.3%, 240 мг -> 80.1%, 400 мг -> 85.8% (Activin E)")


if __name__ == "__main__":
    build_golden()