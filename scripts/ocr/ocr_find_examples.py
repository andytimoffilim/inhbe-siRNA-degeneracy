"""
Быстрый OCR-скан APBDY для поиска страниц с Examples и IC50.
DPI=150 — быстро, для разведки.
"""
import re
import time
from pathlib import Path
import pytesseract
from pdf2image import convert_from_path

REPO_ROOT = Path(__file__).resolve().parents[2]

POPPLER = r"C:\Release-26.09.0-0\Library\bin"
TESSERACT = r"C:\Program Files\Tesseract-OCR\tesseract.exe"
pytesseract.pytesseract.tesseract_cmd = TESSERACT

PDF = REPO_ROOT / "pdf" / "WO2026055502-APBDY-20260312-0011.pdf"
OUT = REPO_ROOT / "ocr_text"
OUT.mkdir(parents=True, exist_ok=True)

# Сканируем все страницы в низком DPI
PAGES = range(1, 148)
DPI = 150

MARKERS = [
    "Example 7", "Example 8", "Example 10", "Example 11",
    "Example 12", "Example 17",
    "Table 13", "Table 14", "Table 15", "Table 16",
    "Table 17", "Table 18", "Table 19",
    "IC50", "IC 50", "ED50",
    "DSR-0104292",
]

def ocr_page(pdf, pnum, dpi=150):
    imgs = convert_from_path(
        str(pdf), dpi=dpi,
        first_page=pnum, last_page=pnum,
        poppler_path=POPPLER,
    )
    if not imgs:
        return ""
    return pytesseract.image_to_string(imgs[0], lang="eng")

def main():
    print(f"Диагностика {PDF.name}: {len(PAGES)} страниц, DPI={DPI}")
    print("Ищем: Examples 7, 8, 10, 17, Tables 13-19, IC50\n")

    hits = {}
    t0 = time.time()

    for pnum in PAGES:
        try:
            text = ocr_page(PDF, pnum, DPI)
        except Exception as e:
            print(f"[err] p{pnum}: {e}")
            continue

        (OUT / f"p{pnum:04d}.txt").write_text(text, encoding="utf-8")

        found = [m for m in MARKERS if m.lower() in text.lower()]

        if found:
            hits[pnum] = found
            print(f"p{pnum:>4}: {len(text):>5} chars  >>> {', '.join(found)}")
        elif pnum % 20 == 0:
            print(f"p{pnum:>4}: {len(text):>5} chars")

    dt = time.time() - t0
    print(f"\nГотово за {dt/60:.1f} мин")
    print(f"Страниц с маркерами: {len(hits)}")
    print(f"Результаты: {OUT}")

    # Группировка: где какие маркеры
    print("\n=== КАРТА МАРКЕРОВ ===")
    for pnum, marks in sorted(hits.items()):
        print(f"p{pnum:>4}: {marks}")

if __name__ == "__main__":
    main()