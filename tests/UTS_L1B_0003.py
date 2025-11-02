# tests/UTS_L1B_0003.py
import os, sys
import numpy as np

# --- root del repo (este script está en tests/)
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from config.globalConfig import globalConfig
from common.io.writeToa import readToa

# RUTAS (ajusta la de OUT si quieres apuntar a myoutputs_eq/noeq)
OUT_DIR = "C:/Users/Kolarov/Desktop/EarthObservation/Codes/S3/EODP-TS-L1B/myoutputs_noeq".replace("\\","/")
REF_DIR = "C:/Users/Kolarov/Desktop/EarthObservation/Codes/S3/EODP-TS-L1B/output".replace("\\","/")

cfg = globalConfig()
bands = list(cfg.bands)

tol_abs = 1.0e-4               # como en ISM
three_sigma_frac = 1 - 0.997   # 3σ

print("\n" + "="*60)
print("=== RESUMEN TEST L1B (ISM-0003) ===")
print(f"REF: {REF_DIR}")
print(f"OUT: {OUT_DIR}")
print("-"*60)

def compare(name):
    out_p = os.path.join(OUT_DIR, name)
    ref_p = os.path.join(REF_DIR, name)
    if not os.path.isfile(out_p) or not os.path.isfile(ref_p):
        return "SKIP", 0, 0
    A = readToa(OUT_DIR, name)
    B = readToa(REF_DIR, name)
    diff = np.abs(A - B)
    n_out = int(np.sum(diff > tol_abs))
    total = diff.size
    status = "OK" if n_out < total * three_sigma_frac else "NOK"
    return status, n_out, total

ok = nok = skip = 0
for b in bands:
    # Si activaste equalización, compara también el intermedio igualado:
    name_eq = f"l1b_toa_eq_{b}.nc"
    st, n, tot = compare(name_eq)
    if st != "SKIP":
        pct = 100*n/max(1,tot)
        print(f"{b:<8} {'l1b_toa_eq':<14} -> {st:>3s} | fuera tol: {n}/{tot} ({pct:.3f}%)")
        ok += (st=="OK"); nok += (st=="NOK")
    else:
        print(f"{b:<8} {'l1b_toa_eq':<14} -> SKIP (no encontrado)")
        skip += 1

    # Producto final en radiancias:
    name = f"l1b_toa_{b}.nc"
    st, n, tot = compare(name)
    if st != "SKIP":
        pct = 100*n/max(1,tot)
        print(f"{b:<8} {'l1b_toa':<14}    -> {st:>3s} | fuera tol: {n}/{tot} ({pct:.3f}%)")
        ok += (st=="OK"); nok += (st=="NOK")
    else:
        print(f"{b:<8} {'l1b_toa':<14}    -> SKIP (no encontrado)")
        skip += 1

print("-"*60)
print(f"Totales → OK: {ok} | NOK: {nok} | SKIP: {skip}")
print("="*60)
