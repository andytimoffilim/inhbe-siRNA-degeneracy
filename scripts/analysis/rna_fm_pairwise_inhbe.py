"""
RNA-FM pairwise distance для 6 INHBE-дуплексов.
Аналог статьи по PCSK9: d_E(sense, mRNA), d_E(antisense, mRNA) vs log10(IC50).
"""
import csv
import re
import warnings
from pathlib import Path
from collections import defaultdict

import numpy as np
import torch
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

REPO_ROOT = Path(__file__).resolve().parents[2]

warnings.filterwarnings("ignore")

# ============================================================
# КОНФИГ
# ============================================================

BASE = REPO_ROOT
GOLDEN = BASE / "data" / "inhbe_golden.csv"
FASTA = BASE / "data" / "NM_031479.5.fasta"
RNA_FM_WEIGHTS = BASE / "weights" / "RNA-FM_pretrained.pth"

WINDOW_NT = 50         # размер мРНК-окна
PCA_K = 10             # число главных компонент
SEED = 42

torch.manual_seed(SEED)
np.random.seed(SEED)


# ============================================================
# 1. RNA-FM — загрузка
# ============================================================

def load_rna_fm():
    """Загружает RNA-FM. Пробует пакет fm, потом локальный чекпоинт."""
    # Попытка 1: rna-fm пакет
    try:
        import fm
        model, alphabet = fm.pretrained.rna_fm_t12()
        model.eval()
        print("[+] RNA-FM загружен через пакет fm")
        return model, alphabet
    except ImportError:
        pass

    # Попытка 2: локальный чекпоинт
    if not RNA_FM_WEIGHTS.exists():
        print(f"[!] Нет {RNA_FM_WEIGHTS}")
        print("    Установите: pip install rna-fm")
        print("    Или скачайте RNA-FM_pretrained.pth")
        print("    Fallback: k-mer эмбеддинги (менее точно)")
        return None, None

    try:
        import fm
        model, alphabet = fm.pretrained.load_model_and_alphabet_local(
            str(RNA_FM_WEIGHTS)
        )
        model.eval()
        print(f"[+] RNA-FM загружен из {RNA_FM_WEIGHTS}")
        return model, alphabet
    except Exception as e:
        print(f"[!] Не удалось загрузить RNA-FM: {e}")
        return None, None


# ============================================================
# 2. Fallback: k-mer эмбеддинги (если RNA-FM нет)
# ============================================================

def kmer_embed(seq, k=3, dim=64):
    """Простейший fallback: частоты k-mer'ов."""
    seq = seq.upper().replace("T", "U")
    vocab = {}
    idx = 0
    for i in range(k):
        for a in "AUGC":
            for b in "AUGC":
                for c in "AUGC" if k == 3 else [" "]:
                    kmer = (a + b + c)[:k]
                    if kmer not in vocab:
                        vocab[kmer] = idx
                        idx += 1
    vec = np.zeros(len(vocab), dtype=np.float32)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in vocab:
            vec[vocab[kmer]] += 1
    if vec.sum() > 0:
        vec /= vec.sum()
    # Проекция на dim
    if len(vec) >= dim:
        vec = vec[:dim]
    else:
        vec = np.pad(vec, (0, dim - len(vec)))
    return vec


# ============================================================
# 3. Эмбеддинги
# ============================================================

def rna_fm_embed(model, alphabet, seq, max_len=1024):
    """Mean-pooled эмбеддинг из RNA-FM."""
    seq = seq.upper().replace("T", "U")
    seq = re.sub(r"[^AUGC]", "N", seq)[:max_len]
    data = [("rna", seq)]        # ← ИСПРАВЛЕНО: (label, sequence)
    batch_converter = alphabet.get_batch_converter()
    _, _, tokens = batch_converter(data)
    with torch.no_grad():
        out = model(tokens, repr_layers=[12], return_contacts=False)
        emb = out["representations"][12]
        # Mean pool по длине (без [CLS], [EOS])
        emb = emb[0, 1:len(seq)+1, :].mean(dim=0).numpy()
    return emb


def get_embedding(model, alphabet, seq):
    if model is not None:
        return rna_fm_embed(model, alphabet, seq)
    return kmer_embed(seq)


# ============================================================
# 4. Утилиты
# ============================================================

def read_fasta(path):
    seq = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq).upper()


def get_mrna_window(mrna, site_start, site_end, window=WINDOW_NT):
    """Окно window nt, центрированное на сайте."""
    center = (site_start + site_end) // 2
    half = window // 2
    lo = max(0, center - half - 1)
    hi = min(len(mrna), center + half)
    return mrna[lo:hi]


def euclidean(a, b):
    return float(np.linalg.norm(a - b))


def cosine(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


# ============================================================
# 5. Bootstrap CI для Spearman
# ============================================================

def bootstrap_spearman(x, y, n_iter=2000, seed=SEED):
    rng = np.random.RandomState(seed)
    n = len(x)
    rhos = []
    for _ in range(n_iter):
        idx = rng.choice(n, n, replace=True)
        if len(set(idx)) < 3:
            continue
        rho, _ = spearmanr(x[idx], y[idx])
        if not np.isnan(rho):
            rhos.append(rho)
    rhos = np.array(rhos)
    return np.percentile(rhos, [2.5, 97.5])


# ============================================================
# 6. MAIN
# ============================================================

def main():
    # Читаем golden
    with open(GOLDEN, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    print(f"[+] Golden: {len(rows)} дуплексов")

    mrna = read_fasta(FASTA)
    print(f"[+] NM_031479.5: {len(mrna)} nt\n")

    # Загрузка RNA-FM
    model, alphabet = load_rna_fm()
    model_name = "RNA-FM" if model else "k-mer (fallback)"

    # --- эмбеддинги ---
    data = []
    for r in rows:
        dsr = r["dsr_name"]
        sense = r["sense_seq"].replace("T", "U")
        anti = r["antisense_seq"].replace("T", "U")
        s_start = int(r["site_start"]) if r["site_start"] else 0
        s_end = int(r["site_end"]) if r["site_end"] else 0

        if not s_start:
            continue

        window_seq = get_mrna_window(mrna, s_start, s_end, WINDOW_NT)

        v_s = get_embedding(model, alphabet, sense)
        v_a = get_embedding(model, alphabet, anti)
        v_m = get_embedding(model, alphabet, window_seq)

        ic50 = float(r["ic50_nm"])
        data.append({
            "dsr": dsr,
            "sense": sense,
            "anti": anti,
            "window": window_seq,
            "site": f"{s_start}-{s_end}",
            "ic50_nm": ic50,
            "log_ic50": np.log10(ic50),
            "v_s": v_s,
            "v_a": v_a,
            "v_m": v_m,
        })

    print(f"[+] Эмбеддинги готовы ({model_name}): {len(data)} дуплексов\n")

    # --- PCA на все эмбеддинги вместе ---
    all_vecs = np.vstack([d["v_s"] for d in data] +
                         [d["v_a"] for d in data] +
                         [d["v_m"] for d in data])
    pca = PCA(n_components=min(PCA_K, all_vecs.shape[0], all_vecs.shape[1]))
    pca.fit(all_vecs)
    print(f"[+] PCA: {PCA_K} компонент, "
          f"explained var = {pca.explained_variance_ratio_.sum():.3f}\n")

    for d in data:
        d["v_s_pca"] = pca.transform(d["v_s"].reshape(1, -1))[0]
        d["v_a_pca"] = pca.transform(d["v_a"].reshape(1, -1))[0]
        d["v_m_pca"] = pca.transform(d["v_m"].reshape(1, -1))[0]

    # --- расстояния ---
    for d in data:
        d["d_E_sense"] = euclidean(d["v_s_pca"], d["v_m_pca"])
        d["d_E_anti"] = euclidean(d["v_a_pca"], d["v_m_pca"])
        d["cos_sense"] = cosine(d["v_s_pca"], d["v_m_pca"])
        d["cos_anti"] = cosine(d["v_a_pca"], d["v_m_pca"])

    # --- вывод таблицы ---
    print(f"{'DSR':<15} {'site':>11} {'IC50':>6} {'dE_s':>7} "
          f"{'dE_a':>7} {'cos_s':>7} {'cos_a':>7}")
    print("-" * 75)
    for d in sorted(data, key=lambda x: x["ic50_nm"]):
        print(f"{d['dsr']:<15} {d['site']:>11} {d['ic50_nm']:>6.3f} "
              f"{d['d_E_sense']:>7.3f} {d['d_E_anti']:>7.3f} "
              f"{d['cos_sense']:>7.3f} {d['cos_anti']:>7.3f}")

    # --- корреляции ---
    print("\n=== КОРРЕЛЯЦИИ (Spearman) ===")
    for key, name in [
        ("d_E_sense", "d_E(sense, mRNA)"),
        ("d_E_anti",  "d_E(antisense, mRNA)"),
        ("cos_sense", "cos(sense, mRNA)"),
        ("cos_anti",  "cos(antisense, mRNA)"),
    ]:
        x = np.array([d[key] for d in data])
        y = np.array([d["log_ic50"] for d in data])
        rho, p = spearmanr(x, y)
        ci = bootstrap_spearman(x, y, n_iter=2000)
        print(f"  {name:<30} ρ = {rho:+.3f}  p = {p:.4f}  "
              f"95% CI [{ci[0]:+.3f}, {ci[1]:+.3f}]")

    # --- сравнение hot-spot vs other ---
    print("\n=== ГРУППЫ ===")
    hot = [d for d in data if 1100 <= int(d["site"].split("-")[0]) <= 1140]
    print(f"  Hot-spot 1107: {len(hot)}/{len(data)}")

    # --- сохранить ---
    out = BASE / "data" / "rna_fm_pairwise_results.csv"
    with open(out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["dsr", "site", "ic50_nm", "log_ic50",
                    "d_E_sense", "d_E_anti", "cos_sense", "cos_anti"])
        for d in sorted(data, key=lambda x: x["ic50_nm"]):
            w.writerow([d["dsr"], d["site"], d["ic50_nm"], d["log_ic50"],
                        d["d_E_sense"], d["d_E_anti"],
                        d["cos_sense"], d["cos_anti"]])
    print(f"\n[+] CSV: {out}")


if __name__ == "__main__":
    main()