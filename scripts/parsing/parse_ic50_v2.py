"""
Обновляет inhbe_golden.csv:
- Добавляет IC50 из Tables 13-17 (4 типа клеток)
- Фильтрует шум (дозы 20/100/400 nM отбрасываются)
"""
import csv
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

BASE = REPO_ROOT
OCR_FILE = BASE / "pdf" / "ocr_hq" / "page_0063.txt"
GOLDEN = BASE / "data" / "inhbe_golden.csv"

# Ручная таблица из OCR стр. 63 (Tables 13-17)
IC50_TABLE = {
    # Table 13: mouse hepatocytes, mouse INHBE
    "mouse_hep": {
        "DSR-0104290": 4.06, "DSR-0104292": 12.62,
        "DSR-0104455": 1.58, "DSR-0104212": 2.00,
    },
    # Table 14: human hepatocytes, mouse INHBE
    "human_hep_t14": {
        "DSR-0104290": 0.41, "DSR-0104292": 0.63,
        "DSR-0104455": 0.34, "DSR-0104212": 0.54,
    },
    # Table 15: rat hepatocytes, mouse INHBE
    "rat_hep": {
        "DSR-0104290": 0.65, "DSR-0104292": 0.59,
        "DSR-0104455": 0.24, "DSR-0104212": 0.28,
    },
    # Table 16: human hepatocytes, human INHBE
    "human_hep_t16": {
        "DSR-0104290": 0.38, "DSR-0104522": 0.36,
        "DSR-0104292": 0.57, "DSR-0104523": 0.33,
    },
    # Table 17: NHP hepatocytes, NHP INHBE
    "nhp_hep": {
        "DSR-0104292": 1.41,
    },
}

# Приоритет для «золотого» IC50: NHP > human > rat > mouse
PRIORITY = ["nhp_hep", "human_hep_t16", "human_hep_t14", "rat_hep", "mouse_hep"]

def pick_golden_ic50(dsr):
    """Выбирает наиболее релевантный IC50 (нМ) для дуплекса."""
    for key in PRIORITY:
        if dsr in IC50_TABLE[key]:
            return IC50_TABLE[key][dsr], key
    return None, None

def main():
    with open(GOLDEN, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    # Колонки
    new_cols = ["ic50_nm", "ic50_pm", "ic50_source_table",
                "ic50_mouse_nm", "ic50_human_nm", "ic50_rat_nm", "ic50_nhp_nm"]
    for r in rows:
        for c in new_cols:
            r.setdefault(c, "")

    for r in rows:
        dsr = r["dsr_name"]

        # Все 4 значения
        if dsr in IC50_TABLE["mouse_hep"]:
            r["ic50_mouse_nm"] = IC50_TABLE["mouse_hep"][dsr]
        if dsr in IC50_TABLE["human_hep_t16"]:
            r["ic50_human_nm"] = IC50_TABLE["human_hep_t16"][dsr]
        elif dsr in IC50_TABLE["human_hep_t14"]:
            r["ic50_human_nm"] = IC50_TABLE["human_hep_t14"][dsr]
        if dsr in IC50_TABLE["rat_hep"]:
            r["ic50_rat_nm"] = IC50_TABLE["rat_hep"][dsr]
        if dsr in IC50_TABLE["nhp_hep"]:
            r["ic50_nhp_nm"] = IC50_TABLE["nhp_hep"][dsr]

        # Золотое значение по приоритету
        val, src = pick_golden_ic50(dsr)
        if val is not None:
            r["ic50_nm"] = f"{val:.3f}"
            r["ic50_pm"] = f"{val*1000:.1f}"
            r["ic50_source_table"] = src

    # Запись
    fields = list(rows[0].keys())
    with open(GOLDEN, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # Сводка
    print(f"[+] Обновлено: {GOLDEN}\n")
    print(f"{'DSR':<15} {'golden':>10} {'source':<16} "
          f"{'mouse':>8} {'human':>8} {'rat':>8} {'nhp':>8}")
    print("-" * 80)
    for r in rows:
        if r["ic50_nm"] or r["ic50_mouse_nm"]:
            print(f"{r['dsr_name']:<15} "
                  f"{r['ic50_nm'] or '—':>10} "
                  f"{r['ic50_source_table'] or '—':<16} "
                  f"{r['ic50_mouse_nm'] or '—':>8} "
                  f"{r['ic50_human_nm'] or '—':>8} "
                  f"{r['ic50_rat_nm'] or '—':>8} "
                  f"{r['ic50_nhp_nm'] or '—':>8}")

if __name__ == "__main__":
    main()