"""
Извлекает HELM-нотации для 6 дуплексов из Table 1A (стр. 23-25).
Парсит HELM → ACGTU-последовательности.
"""
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

OCR_DIR = REPO_ROOT / "ocr_text"
OUT_DIR = Path(r"C:\Users\user\PycharmProjects\DSR-0104292_LOOP\inhbe")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Пары из стр. 22 (DSR → antisense-гид, sense-пассажир)
DSR_PAIRS = {
    "DSR-0104292": ("SSR-0108323", "SSR-0108051"),
    "DSR-0104290": ("SSR-0108321", "SSR-0108051"),
    "DSR-0104455": ("SSR-0108301", "SSR-0108001"),
    "DSR-0104212": ("SSR-0108302", "SSR-0108001"),
    "DSR-0104522": ("SSR-0108682", "SSR-0108051"),
    "DSR-0104523": ("SSR-0108697", "SSR-0108051"),
}

def helm_to_seq(helm):
    """Парсит HELM → последовательность оснований."""
    m = re.search(r"RNA\d*\{(.+?)\}", helm)
    if not m:
        return None
    body = m.group(1)
    seq = []
    i = 0
    while i < len(body):
        c = body[i]
        if c == "[":
            j = body.find("]", i)
            i = (j + 1) if j != -1 else i + 1
            continue
        if c == "(":
            j = body.find(")", i)
            i = (j + 1) if j != -1 else i + 1
            continue
        if c in "AUGC":
            seq.append(c)
        i += 1
    return "".join(seq) if seq else None

def extract_helm_for_ssr(text, ssr_id):
    """
    Ищет HELM-блок для SSR-ID.
    Формат в Table 1A: <SSR-ID> RNA1{...}RNA2{...}
    или: DSR-XXX RNA1{...}|RNA2{...}
    """
    # Паттерн: SSR-ID за которым следует HELM с RNA1/RNA2 в пределах 500 символов
    pattern = re.compile(
        rf"{re.escape(ssr_id)}\s*[\.\s]*(RNA\d?\{{[^}}]+\}})",
        re.DOTALL,
    )
    m = pattern.search(text)
    if m:
        return m.group(1)
    return None

def main():
    # Читаем все OCR-файлы из обеих папок
    ocr_files = list(OCR_DIR.glob("page_00*.txt")) + \
                list(OCR_DIR.glob("page_00*.txt")) + \
                list(Path(OCR_DIR.parent / "ocr_hq").glob("page_00*.txt")) + \
                list(Path(OCR_DIR.parent / "ocr_full").glob("page_00*.txt"))

    # Собираем весь текст
    full_text = ""
    for f in sorted(set(ocr_files)):
        full_text += f"\n\n=== {f.name} ===\n" + f.read_text(encoding="utf-8", errors="ignore")

    print(f"Всего символов: {len(full_text)}")
    n_helm = len(re.findall(r"RNA1\s*\{", full_text))
    print(f"HELM-блоков всего: {n_helm}\n")

    # Для каждого дуплекса
    results = {}
    for dsr, (anti, sense) in DSR_PAIRS.items():
        print(f"=== {dsr} (antisense={anti}, sense={sense}) ===")

        helm_anti = extract_helm_for_ssr(full_text, anti)
        helm_sense = extract_helm_for_ssr(full_text, sense)

        seq_anti = helm_to_seq(helm_anti) if helm_anti else None
        seq_sense = helm_to_seq(helm_sense) if helm_sense else None

        if seq_anti:
            print(f"  {anti}: {seq_anti}")
        else:
            print(f"  {anti}: HELM НЕ НАЙДЕН")

        if seq_sense:
            print(f"  {sense}: {seq_sense}")
        else:
            print(f"  {sense}: HELM НЕ НАЙДЕН")

        results[dsr] = {
            "antisense_id": anti,
            "sense_id": sense,
            "antisense_seq": seq_anti,
            "sense_seq": seq_sense,
            "helm_anti": helm_anti,
            "helm_sense": helm_sense,
        }
        print()

    # Сохраняем
    import json
    out = OUT_DIR / "inhbe_helm_parsed.json"
    out.write_text(json.dumps(results, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    print(f"Сохранено: {out}")

if __name__ == "__main__":
    main()