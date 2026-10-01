"""
Проверка: различает ли RNA-FM разные последовательности?
"""
import re
import numpy as np
import torch

import fm

REPO_ROOT = Path(__file__).resolve().parents[2]
model, alphabet = fm.pretrained.rna_fm_t12()
model.eval()

def embed(seq):
    seq = seq.upper().replace("T", "U")
    seq = re.sub(r"[^AUGC]", "N", seq)
    _, _, tokens = alphabet.get_batch_converter()([("rna", seq)])
    with torch.no_grad():
        out = model(tokens, repr_layers=[12], return_contacts=False)
    emb = out["representations"][12][0, 1:len(seq)+1, :].mean(dim=0).numpy()
    return emb

# 6 уникальных последовательностей из датасета
seqs = {
    "sense_A":   "UCUUUCCAUUCUGCCGUCUUA",
    "sense_B":   "UCCUUCCAUUCUGCCGUCUUA",
    "anti_A":    "UAAGACGGCAGAAUGGAAAGAUU",
    "anti_B":    "UAAGACGGCAGAAUGGAAGGAUU",
    "mRNA_1107": "UAAA",  # placeholder, заменить на реальное окно
}

print("Попарные евклидовы расстояния между эмбеддингами:")
keys = list(seqs.keys())
embs = {k: embed(v) for k, v in seqs.items()}
for i, k1 in enumerate(keys):
    for k2 in keys[i+1:]:
        d = np.linalg.norm(embs[k1] - embs[k2])
        print(f"  {k1:<12} vs {k2:<12}: {d:.4f}")