"""
Ищет в OCR-файлах контекст вокруг имён новых DSR-дуплексов.
Показывает строки с последовательностями (ACGTU) рядом с именами.
"""
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

OCR_DIRS = [
    REPO_ROOT / "ocr_text",
    REPO_ROOT / "ocr_text",
    REPO_ROOT / "ocr_text",
]

TARGETS = [
    "DSR-0104290", "DSR-0104455", "DSR-0104212",
    "DSR-0104522", "DSR-0104523",
]

# Последовательность: 15+ букв ACGTU
SEQ_RE = re.compile(r"[ACGTUacgtu]{15,}")

def main():
    for d in OCR_DIRS:
        if not d.exists():
            continue
        print(f"\n=== {d.name} ===\n")
        for txt in sorted(d.glob("*.txt")):
            text = txt.read_text(encoding="utf-8", errors="ignore")
            for target in TARGETS:
                if target not in text:
                    continue
                print(f"--- {txt.name}: {target} ---")
                # Показать все строки с target
                for line in text.split("\n"):
                    if target in line:
                        print(f"  [{line[:200]}]")
                # Показать 10 строк после target
                idx = text.find(target)
                if idx != -1:
                    context = text[idx:idx+1500]
                    print(f"  Контекст:")
                    for line in context.split("\n")[:15]:
                        line = line.strip()
                        if line:
                            # Помечаем строки с последовательностями
                            seqs = SEQ_RE.findall(line)
                            mark = "  <<< SEQ" if seqs else ""
                            print(f"    {line[:150]}{mark}")
                print()

if __name__ == "__main__":
    main()