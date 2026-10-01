# probe_scans.py
# Разведка сканов: метаданные + OCR первых 2 страниц.
# Windows 11, Python 3.
# Требуется: pdfplumber, pdf2image, pillow, pytesseract
# Внешние: poppler (C:\Release-26.09.0-0), tesseract (C:\Program Files\Tesseract-OCR)

import os
import sys
from pathlib import Path
import pdfplumber

REPO_ROOT = Path(__file__).resolve().parents[2]

# --- ваши пути ---
POPPLER_PATH  = r"C:\Release-26.09.0-0\Library\bin"
TESSERACT_CMD = r"C:\Program Files\Tesseract-OCR\tesseract.exe"
OCR_LANG      = "eng"          # или "eng+rus"
DPI           = 200
MAX_PAGES_OCR = 2

BASE = Path(r"E:\datasets\mRNA\new_patents")
OUT  = Path(r"E:\datasets\mRNA\diagnose\scan_probe")
OUT.mkdir(parents=True, exist_ok=True)

SCAN_FILES = [
    "WO2008147824A2.pdf",
    "WO2010002626A2.pdf",
    "WO2010048605A1.pdf",
]


def check_binaries():
    """Проверяет, что внешние утилиты на месте."""
    problems = []

    if not Path(TESSERACT_CMD).exists():
        problems.append(f"Tesseract не найден: {TESSERACT_CMD}")

    if POPPLER_PATH:
        p = Path(POPPLER_PATH)
        if not p.exists():
            problems.append(f"POPPLER_PATH не существует: {p}")
        else:
            for exe in ("pdftoppm.exe", "pdfinfo.exe", "pdftocairo.exe"):
                if not (p / exe).exists():
                    problems.append(f"В {p} нет {exe}")

    # Пробуем импорт python-зависимостей
    for mod in ("pdfplumber", "pdf2image", "PIL", "pytesseract"):
        try:
            __import__(mod)
        except Exception as e:
            problems.append(f"Не установлен пакет {mod}: {e}")

    return problems


def read_metadata(pdf_path):
    try:
        with pdfplumber.open(pdf_path) as pdf:
            return pdf.metadata or {}
    except Exception as e:
        return {"__error__": str(e)}


def ocr_pages(pdf_path, pages, dpi=DPI, lang=OCR_LANG):
    """Возвращает {page_number: ocr_text}."""
    results = {}

    from pdf2image import convert_from_path
    import pytesseract
    pytesseract.pytesseract.tesseract_cmd = TESSERACT_CMD

    for pn in pages:
        try:
            kwargs = dict(dpi=dpi, first_page=pn, last_page=pn)
            if POPPLER_PATH:
                kwargs["poppler_path"] = POPPLER_PATH
            imgs = convert_from_path(str(pdf_path), **kwargs)
            for img in imgs:
                png_path = OUT / f"{pdf_path.stem}__page_{pn:04d}.png"
                img.save(png_path, "PNG")
                try:
                    txt = pytesseract.image_to_string(img, lang=lang)
                except Exception as e:
                    txt = f"[OCR error: {e}]"
                results[pn] = txt
                (OUT / f"{pdf_path.stem}__page_{pn:04d}.txt").write_text(
                    txt, encoding="utf-8"
                )
                print(f"  [OCR] {png_path.name} -> "
                      f"{pdf_path.stem}__page_{pn:04d}.txt ({len(txt)} симв.)")
        except Exception as e:
            print(f"  [skip] page {pn}: {e}")

    return results


def main():
    print("Проверка окружения:")
    problems = check_binaries()
    if problems:
        print("  !!! Найдены проблемы:")
        for p in problems:
            print(f"      - {p}")
        print("  Скрипт продолжит работу, но OCR может не сработать.")
    else:
        print("  всё на месте: poppler, tesseract, пакеты Python\n")

    if not BASE.exists():
        print(f"Нет папки: {BASE}")
        sys.exit(1)

    for fname in SCAN_FILES:
        pdf_path = BASE / fname
        print(f"\n=== {fname} ===")
        if not pdf_path.exists():
            print("  файл не найден")
            continue

        md = read_metadata(pdf_path)
        if md:
            printed = False
            for k in ("Title", "Author", "Subject", "Creator", "Producer"):
                if md.get(k):
                    print(f"  meta.{k}: {md[k]}")
                    printed = True
            if not printed:
                print("  (метаданные пусты)")

        pages = list(range(1, MAX_PAGES_OCR + 1))
        ocr = ocr_pages(pdf_path, pages)

        for pn in pages:
            txt = (ocr.get(pn) or "").strip()
            if not txt:
                continue
            low = txt.lower()
            marks = [kw for kw in
                     ("pcsk9", "proprotein convertase", "sirna", "dsrna",
                      "alnylam", "ic50", "hypercholesterol", "antisense")
                     if kw in low]
            print(f"  --- page {pn} OCR (первые 500 символов) ---")
            print("  " + " ".join(txt.split())[:500])
            if marks:
                print(f"  >>> маркеры: {', '.join(marks)}")

    print(f"\nВсё в: {OUT}")


if __name__ == "__main__":
    main()