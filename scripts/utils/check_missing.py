import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

CSV = REPO_ROOT / "data"

with open(CSV, encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f))

print(f"{'DSR':<15} {'site_start':>10} {'site_end':>9} {'len':>5}  notes")
print("-" * 70)
for r in rows:
    print(f"{r['dsr_name']:<15} {r['site_start']:>10} {r['site_end']:>9} "
          f"{r['site_length']:>5}  {r['notes']}")

missing = [r["dsr_name"] for r in rows if not r["site_start"]]
print(f"\nБез сайта: {missing}")