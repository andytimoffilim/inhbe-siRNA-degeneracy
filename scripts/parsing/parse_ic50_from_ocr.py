"""
Парсит IC50 из OCR-файлов патента и обновляет inhbe_golden.csv.
"""
import csv
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

BASE = REPO_ROOT
OCR_DIR = BASE / "pdf" / "ocr_hq"
GOLDEN = BASE / "data" / "inhbe_golden.csv"

# Известные IC50 из OCR (нМ)
IC50_DATA = {
    "DSR-0104292": [
        {"value_nm": 12.62, "ci_low": 7.19, "ci_high": 22.38,
         "cell": "mouse hepatocytes", "table": "Table 13", "page": 63},
        {"value_nm": 1.41, "ci_low": 0.65, "ci_high": 2.97,
         "cell": "NHP hepatocytes", "table": "Table 17", "page": 63},
    ],
}

# ED50 (мг/кг) для мышиного суррогата
ED50_DATA = {
    "DSR-0104212": {"value_mg_kg": "< 1", "source": "Summary, p. 86"},
}

def update_golden():
    with open(GOLDEN, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    # Добавляем новые колонки
    new_cols = ["ic50_nm_mean", "ic50_nm_ci", "ic50_cell_type",
                "ic50_table", "ed50_mg_kg", "ed50_source"]
    for r in rows:
        for c in new_cols:
            r.setdefault(c, "")

    for r in rows:
        name = r["dsr_name"]

        # IC50
        if name in IC50_DATA:
            records = IC50_DATA[name]
            # Среднее по всем записям (в нМ)
            vals = [rec["value_nm"] for rec in records]
            mean_nm = sum(vals) / len(vals)

            # pM
            mean_pm = mean_nm * 1000

            r["ic50_nm_mean"] = f"{mean_nm:.3f}"
            r["ic50_pm"] = f"{mean_pm:.1f}"
            r["ic50_nm_ci"] = "; ".join(
                f"{rec['ci_low']}–{rec['ci_high']}" for rec in records
            )
            r["ic50_cell_type"] = "; ".join(
                rec["cell"] for rec in records
            )
            r["ic50_table"] = "; ".join(
                rec["table"] for rec in records
            )

        # ED50
        if name in ED50_DATA:
            r["ed50_mg_kg"] = ED50_DATA[name]["value_mg_kg"]
            r["ed50_source"] = ED50_DATA[name]["source"]

    # Запись
    fields = list(rows[0].keys())
    with open(GOLDEN, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # Сводка
    n_ic50 = sum(1 for r in rows if r["ic50_pm"])
    print(f"[+] Обновлено: {GOLDEN}")
    print(f"    С IC50: {n_ic50} дуплексов")
    for r in rows:
        if r["ic50_pm"]:
            print(f"    {r['dsr_name']}: {r['ic50_pm']} pM "
                  f"({r['ic50_nm_mean']} nM, {r['ic50_cell_type']})")


def search_more_ic50():
    """Ищет дополнительные IC50 в OCR-файлах."""
    print("\n=== Поиск IC50 в OCR-файлах ===")
    pattern = re.compile(
        r"(DSR-\d+)[^\n]{0,60}?(\d+\.\d+)\s*(?:nM|нМ)"
    )
    for txt in sorted(OCR_DIR.glob("*.txt")):
        text = txt.read_text(encoding="utf-8")
        for m in pattern.finditer(text):
            print(f"  {txt.name}: {m.group(1)} = {m.group(2)} nM")


if __name__ == "__main__":
    update_golden()
    search_more_ic50()