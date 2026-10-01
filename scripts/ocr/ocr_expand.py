"""
OCR всех страниц с DSR- и Table-маркерами (60-88, 17, 22, 41, 76-82).
"""
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

# Все страницы с DSR-/Table-/Example-маркерами из диагностики
TARGET_PAGES = [
    17, 22, 41,
    *range(60, 89),   # 60..88 — Examples 7-17, Tables 11-19
    100, 120, 140,    # для полноты
]

DPI = 400
LANG = "eng"

def ocr_page(pdf, pnum, dpi=400):
    imgs = convert_from_path(
        str(pdf), dpi=dpi,
        first_page=pnum, last_page=pnum,
        poppler_path=POPPLER,
    )
    if not imgs:
        return ""
    return pytesseract.image_to_string(imgs[0], lang=LANG, config="--psm 6")

def main():
    print(f"OCR {len(TARGET_PAGES)} страниц, DPI={DPI}")
    t0 = time.time()

    for i, pnum in enumerate(TARGET_PAGES, 1):
        txt = OUT / f"page_{pnum:04d}.txt"
        if txt.exists() and txt.stat().st_size > 200:
            print(f"[skip] p{pnum}")
            continue
        try:
            text = ocr_page(PDF, pnum, DPI)
            txt.write_text(text, encoding="utf-8")
            print(f"[{i:>2}/{len(TARGET_PAGES)}] p{pnum}: {len(text)} chars")
        except Exception as e:
            print(f"[err] p{pnum}: {e}")

    print(f"\nГотово за {(time.time()-t0)/60:.1f} мин -> {OUT}")

if __name__ == "__main__":
    main()