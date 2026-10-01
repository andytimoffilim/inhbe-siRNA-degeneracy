"""
OCR только целевых страниц патента WO 2026/055502 A3.
PDF — скан, поэтому используем pdf2image + Tesseract.
"""
import re
import time
from pathlib import Path
import pytesseract
from pdf2image import convert_from_path

REPO_ROOT = Path(__file__).resolve().parents[2]

# ---------- пути ----------
POPPLER = r"C:\Release-26.09.0-0\Library\bin"
TESSERACT = r"C:\Program Files\Tesseract-OCR\tesseract.exe"
pytesseract.pytesseract.tesseract_cmd = TESSERACT

PDF = REPO_ROOT / "pdf" / "WO2026055502-APBDY-20260312-0011.pdf"
OUT = REPO_ROOT / "ocr_text"
OUT.mkdir(parents=True, exist_ok=True)

# ---------- целевые страницы ----------
TARGET_PAGES = (
    list(range(10, 15))       # III-a/III-b
    + [23, 24, 25]            # Table 1A (HELM)
    + [65]                    # IUPAC
    + list(range(230, 246))   # Examples 7, 8, 10 + Tables 13-17
    + list(range(279, 290))   # Example 17
)

DPI = 300          # 300 — оптимально для сканов патентов
LANG = "eng"       # английский; если есть русский — "eng+rus"

# ---------- маркеры для поиска ----------
MARKERS = [
    "DSR-", "SSR-", "IC50", "ED50", "IC 50", "EC50",
    "Example 7", "Example 8", "Example 10", "Example 17",
    "Table 13", "Table 14", "Table 15", "Table 16", "Table 17", "Table 19",
    "GalNAc", "HELM", "RNA1", "RNA2", "INHBE", "Activin",
]

# ---------- OCR одной страницы ----------
def ocr_page(pdf_path, pnum, dpi=300, lang="eng"):
    imgs = convert_from_path(
        str(pdf_path), dpi=dpi,
        first_page=pnum, last_page=pnum,
        poppler_path=POPPLER,
    )
    if not imgs:
        return ""
    return pytesseract.image_to_string(imgs[0], lang=lang)

# ---------- main ----------
def main():
    print(f"OCR {len(TARGET_PAGES)} страниц из {PDF.name}")
    print(f"DPI={DPI}, lang={LANG}")
    print(f"Выход: {OUT}\n")

    t_start = time.time()

    for i, pnum in enumerate(TARGET_PAGES, 1):
        txt_path = OUT / f"page_{pnum:04d}.txt"

        # пропускаем уже сделанные
        if txt_path.exists() and txt_path.stat().st_size > 50:
            print(f"[skip] p{pnum} — уже есть")
            continue

        try:
            t0 = time.time()
            text = ocr_page(PDF, pnum, DPI, LANG)
            dt = time.time() - t0
        except Exception as e:
            print(f"[err]  p{pnum}: {e}")
            continue

        txt_path.write_text(text, encoding="utf-8")

        # ищем маркеры
        low = text.lower()
        found = [m for m in MARKERS if m.lower() in low]

        flags = f"  >>> {', '.join(found[:5])}" if found else ""
        print(f"[{i:>2}/{len(TARGET_PAGES)}] "
              f"p{pnum:>4}: {len(text):>6} chars ({dt:.1f}s){flags}")

    total = time.time() - t_start
    print(f"\nГотово за {total/60:.1f} мин")

    # ---------- сводка: где IC50 и DSR ----------
    print("\n=== ГДЕ ЧТО НАЙДЕНО ===\n")

    for pnum in TARGET_PAGES:
        txt_path = OUT / f"page_{pnum:04d}.txt"
        if not txt_path.exists():
            continue
        text = txt_path.read_text(encoding="utf-8")

        # IC50-паттерны
        ic50_hits = re.findall(
            r"IC\s*50[^\d]{0,15}([\d.]+)\s*(pM|nM|µM|uM|M)",
            text, re.IGNORECASE,
        )
        # DSR/SSR-номера
        dsr_hits = set(re.findall(r"DSR-\d+", text))
        ssr_hits = set(re.findall(r"SSR-\d+", text))
        # HELM-блоки
        helm_hits = re.findall(r"RNA[12]\{", text)
        # Примеры
        examples = re.findall(r"Example\s+\d+", text)

        if not (ic50_hits or dsr_hits or ssr_hits or helm_hits or examples):
            continue

        print(f"p{pnum}:")
        if examples:
            print(f"  Примеры:  {sorted(set(examples))}")
        if dsr_hits:
            print(f"  DSR:      {sorted(dsr_hits)[:10]}")
        if ssr_hits:
            print(f"  SSR:      {sorted(ssr_hits)[:10]}")
        if helm_hits:
            print(f"  HELM:     {len(helm_hits)} блоков")
        if ic50_hits:
            print(f"  IC50:     {ic50_hits[:10]}")
        print()

if __name__ == "__main__":
    main()