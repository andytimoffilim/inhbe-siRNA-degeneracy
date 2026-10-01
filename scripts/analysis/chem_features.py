"""
Гибрид: RNA-FM + химия. RF на малом датасете (LOO).
"""
import csv
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import LeaveOneOut, cross_val_predict

REPO_ROOT = Path(__file__).resolve().parents[2]

GOLDEN = Path(r"C:\Users\user\PycharmProjects\DSR-0104292_LOOP"
              r"\inhbe\inhbe_golden.csv")

# Химические признаки из патента
CHEM = {
    # 5'-модификация антисмысла: ptz = 1 если есть [cp][tz], иначе 0
    "DSR-0104292": {"5ptz": 1, "pn_count": 2, "n_ps": 6, "mismatch": 1},
    "DSR-0104290": {"5ptz": 0, "pn_count": 2, "n_ps": 4, "mismatch": 1},
    "DSR-0104455": {"5ptz": 0, "pn_count": 2, "n_ps": 4, "mismatch": 2},
    "DSR-0104212": {"5ptz": 0, "pn_count": 2, "n_ps": 4, "mismatch": 2},
    "DSR-0104522": {"5ptz": 0, "pn_count": 2, "n_ps": 4, "mismatch": 1},
    "DSR-0104523": {"5ptz": 0, "pn_count": 2, "n_ps": 4, "mismatch": 1},
}

def main():
    with open(GOLDEN, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    X, y, names = [], [], []
    for r in rows:
        d = r["dsr_name"]
        c = CHEM.get(d, {})
        X.append([
            c.get("5ptz", 0),
            c.get("pn_count", 0),
            c.get("n_ps", 0),
            c.get("mismatch", 0),
            int(r["sense_len"]),
            int(r["antisense_len"]),
        ])
        y.append(np.log10(float(r["ic50_nm"])))
        names.append(d)

    X = np.array(X, dtype=float)
    y = np.array(y)

    print(f"{'DSR':<15} {'5ptz':>5} {'pn':>4} {'ps':>4} {'mm':>4} "
          f"{'len_s':>6} {'len_a':>6} {'logIC50':>8}")
    print("-" * 75)
    for n, x, yy in zip(names, X, y):
        print(f"{n:<15} {int(x[0]):>5} {int(x[1]):>4} {int(x[2]):>4} "
              f"{int(x[3]):>4} {int(x[4]):>6} {int(x[5]):>6} {yy:>8.3f}")

    # LOO
    model = RandomForestRegressor(n_estimators=200, max_depth=3,
                                   random_state=42)
    loo = LeaveOneOut()
    y_pred = cross_val_predict(model, X, y, cv=loo)
    rho, p = spearmanr(y, y_pred)

    print(f"\n=== LOO Spearman ===")
    print(f"  ρ = {rho:+.3f}, p = {p:.4f}")
    print(f"\n{'DSR':<15} {'actual':>8} {'pred':>8} {'err':>8}")
    for n, yt, yp in zip(names, y, y_pred):
        print(f"  {n:<15} {yt:>8.3f} {yp:>8.3f} {yp-yt:>+8.3f}")

    # Feature importance
    model.fit(X, y)
    print(f"\n=== Feature importance ===")
    feat_names = ["5ptz", "pn_count", "n_ps", "mismatch", "len_s", "len_a"]
    for n, imp in sorted(zip(feat_names, model.feature_importances_),
                          key=lambda x: -x[1]):
        print(f"  {n:<12} {imp:.3f}")


if __name__ == "__main__":
    main()