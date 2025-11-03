# tests/UTS_E2E_0004.py
import os, sys
import numpy as np

# --- root del repo (este script está en tests/)
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from config.globalConfig import globalConfig
from common.io.writeToa import readToa

# ====== RUTAS EXACTAS ======
OUT_ISM_DIR = "C:/Users/Kolarov/Desktop/EarthObservation/Codes/S3/EODP-TS-ISM/myoutput".replace("\\","/")
OUT_L1B_DIR = "C:/Users/Kolarov/Desktop/EarthObservation/Codes/S3/EODP-TS-L1B/myoutputs_noeq".replace("\\","/")

REF_E2E_ROOT = "C:/Users/Kolarov/Desktop/EarthObservation/Codes/S3/EODP-TS-E2E".replace("\\","/")
REF_ISM_DIR  = f"{REF_E2E_ROOT}/ism_out"
REF_L1B_DIR  = f"{REF_E2E_ROOT}/l1b_out"

cfg = globalConfig()
bands = list(cfg.bands)

# Tolerancias
tol_abs = 1.0e-4
three_sigma_frac = 1 - 0.997

print("\n" + "="*64)
print("=== RESUMEN TEST E2E (final ISM DN + L1B radiancias) ===")
print(f"REF ISM : {REF_ISM_DIR}")
print(f"REF L1B : {REF_L1B_DIR}")
print(f"OUT ISM : {OUT_ISM_DIR}")
print(f"OUT L1B : {OUT_L1B_DIR}")
print("-"*64)

def compare(dir_out, dir_ref, name):
    out_p = os.path.join(dir_out, name)
    ref_p = os.path.join(dir_ref, name)
    if not os.path.isfile(out_p) or not os.path.isfile(ref_p):
        return "SKIP", 0, 0
    A = readToa(dir_out, name)
    B = readToa(dir_ref, name)
    if A.shape != B.shape:
        print(f"  [SKIP] shape mismatch: OUT {A.shape} vs REF {B.shape} -> {name}")
        return "SKIP", 0, 0
    diff = np.abs(A - B)
    n_out = int(np.sum(diff > tol_abs))
    total = int(diff.size)
    status = "OK" if n_out < total * three_sigma_frac else "NOK"
    return status, n_out, total

ok = nok = skip = 0

# ---- 1) ISM final (DN) ----
print("— ISM final (DN)")
for b in bands:
    name = f"ism_toa_{b}.nc"
    st, n, tot = compare(OUT_ISM_DIR, REF_ISM_DIR, name)
    if st == "SKIP":
        print(f"{b:<8} ism_toa         -> SKIP (no encontrado)")
        skip += 1
    else:
        pct = 100*n/max(1,tot)
        print(f"{b:<8} ism_toa         -> {st:>3s} | fuera tol: {n}/{tot} ({pct:.3f}%)")
        ok += (st=="OK"); nok += (st=="NOK")
print("-"*64)

# ---- 2) L1B (radiancias) ----
print("— L1B (radiancias)")
for b in bands:
    name = f"l1b_toa_{b}.nc"
    st, n, tot = compare(OUT_L1B_DIR, REF_L1B_DIR, name)
    if st == "SKIP":
        print(f"{b:<8} l1b_toa         -> SKIP (no encontrado)")
        skip += 1
    else:
        pct = 100*n/max(1,tot)
        print(f"{b:<8} l1b_toa         -> {st:>3s} | fuera tol: {n}/{tot} ({pct:.3f}%)")
        ok += (st=="OK"); nok += (st=="NOK")

print("-"*64)
print(f"Totales → OK: {ok} | NOK: {nok} | SKIP: {skip}")
print("="*64)
