"""
Собирает golden INHBE: 6 дуплексов с IC50, последовательностями и сайтами.
Сайт ищется в окне ±30 nt вокруг известного hot-spot 1107.
"""
import csv
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

BASE = REPO_ROOT
GOLDEN = BASE / "data" / "inhbe_golden.csv"

DSR_PAIRS = {
    "DSR-0104292": ("SSR-0108323", "SSR-0108051"),
    "DSR-0104290": ("SSR-0108321", "SSR-0108051"),
    "DSR-0104455": ("SSR-0108301", "SSR-0108001"),
    "DSR-0104212": ("SSR-0108302", "SSR-0108001"),
    "DSR-0104522": ("SSR-0108682", "SSR-0108051"),
    "DSR-0104523": ("SSR-0108697", "SSR-0108051"),
}

SEQUENCES = {
    "SSR-0108051": "UCUUUCCAUUCUGCCGUCUUA",
    "SSR-0108001": "UCCUUCCAUUCUGCCGUCUUA",
    "SSR-0108323": "UAAGACGGCAGAAUGGAAAGAUU",
    "SSR-0108321": "UAAGACGGCAGAAUGGAAAGAUU",
    "SSR-0108682": "UAAGACGGCAGAAUGGAAAGAUU",
    "SSR-0108697": "UAAGACGGCAGAAUGGAAAGAUU",
    "SSR-0108301": "UAAGACGGCAGAAUGGAAGGAUU",
    "SSR-0108302": "UAAGACGGCAGAAUGGAAGGAUU",
}

# В parse_and_build_golden.py заменить IC50:
IC50 = {
    "DSR-0104290": {"human": 0.74},
    "DSR-0104292": {"human": 0.94},
    "DSR-0104455": {"human": 0.55},
    "DSR-0104212": {"human": 0.90},
    "DSR-0104522": {"human": 0.54},
    "DSR-0104523": {"human": 0.55},
}

HOTSPOT = 1107        # известная позиция из BLAST
WINDOW = 30           # ±30 nt


def read_fasta(path):
    if not path.exists():
        return ""
    seq = []
    for line in path.read_text().splitlines():
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq).upper()


def find_site_window(mrna, seq, center=HOTSPOT, radius=WINDOW):
    """Ищет сайт в окне ±radius вокруг center. Возвращает (start, end, len, mm)."""
    seq = re.sub(r"[^ACGT]", "", seq.upper().replace("U", "T"))
    if not mrna or len(seq) < 17:
        return None
    L = len(seq)
    lo = max(0, center - radius - 1)
    hi = min(len(mrna), center + radius + L)

    best = None
    for i in range(lo, hi - L + 1):
        window = mrna[i:i + L]
        mm = sum(1 for a, b in zip(seq, window) if a != b)
        if mm <= 2:
            if best is None or mm < best[3]:
                best = (i + 1, i + L, L, mm)
                if mm == 0:
                    break
    return best


def main():
    mrna = read_fasta(BASE / "data" / "NM_031479.5.fasta")
    print(f"NM_031479.5: {len(mrna)} nt")
    print(f"Hot-spot: {HOTSPOT}, window ±{WINDOW}\n")

    rows = []
    for dsr, (anti_id, sense_id) in DSR_PAIRS.items():
        sense = SEQUENCES.get(sense_id, "")
        anti = SEQUENCES.get(anti_id, "")

        site = find_site_window(mrna, sense) if sense and mrna else None
        site_start = site[0] if site else ""
        site_end = site[1] if site else ""
        site_len = site[2] if site else ""
        site_mm = site[3] if site else ""

        ic = IC50.get(dsr, {})
        golden_ic50 = ic.get("nhp") or ic.get("human") or ic.get("rat") or ic.get("mouse")
        golden_source = ("nhp" if ic.get("nhp") else
                        "human" if ic.get("human") else
                        "rat" if ic.get("rat") else
                        "mouse" if ic.get("mouse") else "")

        rows.append({
            "dsr_name": dsr,
            "sense_id": sense_id,
            "antisense_id": anti_id,
            "sense_seq": sense,
            "antisense_seq": anti,
            "sense_len": len(sense),
            "antisense_len": len(anti),
            "site_start": site_start,
            "site_end": site_end,
            "site_length": site_len,
            "site_mismatch": site_mm,
            "ic50_nm": f"{golden_ic50:.3f}" if golden_ic50 else "",
            "ic50_pm": f"{golden_ic50 * 1000:.1f}" if golden_ic50 else "",
            "ic50_source": golden_source,
            "ic50_mouse_nm": ic.get("mouse", ""),
            "ic50_human_nm": ic.get("human", ""),
            "ic50_rat_nm": ic.get("rat", ""),
            "ic50_nhp_nm": ic.get("nhp", ""),
            "source": "WO 2026/055502 A3",
        })

    rows.sort(key=lambda r: float(r["ic50_nm"]) if r["ic50_nm"] else 9999)

    fields = list(rows[0].keys())
    with open(GOLDEN, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print(f"{'DSR':<15} {'site':>12} {'mm':>4} {'IC50 nM':>9} {'src':>7} "
          f"{'sense':<12} {'anti':<13}")
    print("-" * 90)
    for r in rows:
        site = f"{r['site_start']}-{r['site_end']}" if r["site_start"] else "—"
        print(f"{r['dsr_name']:<15} {site:>12} {str(r['site_mismatch']):>4} "
              f"{r['ic50_nm']:>9} {r['ic50_source']:>7} "
              f"{r['sense_id']:<12} {r['antisense_id']:<13}")

    print(f"\nCSV: {GOLDEN}")
    print(f"Всего: {len(rows)}, с IC50: {sum(1 for r in rows if r['ic50_nm'])}")
    print(f"На hot-spot 1107: {sum(1 for r in rows if r['site_start'] and abs(int(r['site_start']) - HOTSPOT) <= WINDOW)}")


if __name__ == "__main__":
    main()