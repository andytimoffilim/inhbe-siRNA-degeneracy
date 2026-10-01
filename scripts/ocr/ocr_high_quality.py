"""
Высококачественный OCR для страниц с IC50/ED50.
DPI=400 для лучшего качества на мелком шрифте таблиц.
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

# Все страницы с DSR-0104292 или Examples/Tables
TARGET_PAGES = [
    9, 10, 13,
    17,
    22,
    41,
    60, 61, 63,          # Example 7, 8, Tables 13-17, IC50
    65, 66, 67, 68,      # IUPAC, Table 18
    69, 70, 71,          # Example 10, 11, 12, Table 19
    76, 77, 78, 79, 80, 81, 82,
    84, 85, 86, 87, 88,  # Example 17, ED50
]

DPI = 400   # выше для таблиц
LANG = "eng"

def ocr_page(pdf, pnum, dpi=400):
    imgs = convert_from_path(
        str(pdf), dpi=dpi,
        first_page=pnum, last_page=pnum,
        poppler_path=POPPLER,
    )
    if not imgs:
        return ""
    # --psm 6 = assume uniform block of text (хорошо для таблиц)
    return pytesseract.image_to_string(imgs[0], lang=LANG, config="--psm 6")

def main():
    print(f"HQ OCR {len(TARGET_PAGES)} страниц, DPI={DPI}")
    print(f"Выход: {OUT}\n")

    t0 = time.time()
    for i, pnum in enumerate(TARGET_PAGES, 1):
        txt_path = OUT / f"page_{pnum:04d}.txt"
        if txt_path.exists() and txt_path.stat().st_size > 200:
            print(f"[skip] p{pnum} — уже есть")
            continue

        try:
            t1 = time.time()
            text = ocr_page(PDF, pnum, DPI)
            dt = time.time() - t1
        except Exception as e:
            print(f"[err]  p{pnum}: {e}")
            continue

        txt_path.write_text(text, encoding="utf-8")

        # маркеры
        low = text.lower()
        marks = []
        for kw in ["IC50", "IC 50", "ED50", "Example", "Table 1"]:
            if kw.lower() in low:
                marks.append(kw)

        flags = f"  >>> {', '.join(marks)}" if marks else ""
        print(f"[{i:>2}/{len(TARGET_PAGES)}] p{pnum:>4}: "
              f"{len(text):>6} chars ({dt:.1f}s){flags}")

    print(f"\nГотово за {(time.time()-t0)/60:.1f} мин")
    print(f"Результаты: {OUT}")

if __name__ == "__main__":
    main()