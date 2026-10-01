"""
Ищет ВСЕ IC50 в OCR-файлах: таблицы, текст, разные форматы.
"""
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

OCR_DIRS = [
    REPO_ROOT / "ocr_text",
    REPO_ROOT / "ocr_text",
    REPO_ROOT / "ocr_text",
]

PATTERNS = [
    # Табличный: DSR-0104290 4.06 to 17.37
    re.compile(r"(DSR-\d+)\s+([\d.]+)\s+(?:to\s+)?([\d.]+)?"),
    # Текстовый: DSR-0104292 ... IC50 ... 12.62 nM
    re.compile(r"(DSR-\d+)[^\n]{0,200}?IC50[^\d]{0,20}([\d.]+)\s*(?:nM|нМ)",
               re.IGNORECASE),
    # Табличный с "IC50(nM)": DSR-0104292 12.62
    re.compile(r"(DSR-\d+)\s+([\d.]+)\s*\|", re.IGNORECASE),
]

def scan():
    all_hits = {}
    for d in OCR_DIRS:
        if not d.exists():
            continue
        for txt in sorted(d.glob("*.txt")):
            text = txt.read_text(encoding="utf-8", errors="ignore")
            for pat in PATTERNS:
                for m in pat.finditer(text):
                    dsr = m.group(1)
                    val = m.group(2)
                    try:
                        f = float(val)
                    except ValueError:
                        continue
                    if 0.01 <= f <= 1000:
                        all_hits.setdefault(dsr, []).append({
                            "value": f,
                            "file": txt.name,
                            "match": m.group(0)[:120],
                        })

    print(f"=== Найдено IC50 для {len(all_hits)} дуплексов ===\n")
    for dsr in sorted(all_hits):
        print(f"{dsr}:")
        for h in all_hits[dsr]:
            print(f"  {h['value']:>8.3f} nM  ({h['file']})")
            print(f"    {h['match']!r}")
        print()

if __name__ == "__main__":
    scan()