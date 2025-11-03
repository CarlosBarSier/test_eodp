import os, sys
import numpy as np

# --- root del repo (este script está en tests/)
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from config.globalConfig import globalConfig
from common.io.writeToa import readToa

# ====== RUTAS (ajústalas si hace falta; usa slashes) ======
OUT_ISM_DIR = "C:/Users/Kolarov/Desktop/EarthObservation/Codes/S3/EODP-TS-ISM/myoutput".replace("\\","/")
OUT_L1B_DIR = "C:/Users/Kolarov/Desktop/EarthObservation/Codes/S3/EODP-TS-L1B/myoutputs_noeq".replace("\\","/")

REF_E2E_ROOT = "C:/Users/Kolarov/Desktop/EarthObservation/Codes/S3/EODP-TS-E2E".replace("\\","/")
REF_ISM_DIR  = f"{REF_E2E_ROOT}/ism_out"
REF_L1B_DIR  = f"{REF_E2E_ROOT}/l1b_out"

# ====== Parámetros ======
bands = list(globalConfig().bands)
tol_abs = 1.0e-4                 # tolerancia absoluta
three_sigma_frac = 1 - 0.997     # 3σ (~0.3% fuera de tol)

def center_crop(A, target_shape):
    r, c = A.shape
    tr, tc = target_shape
    sr = max((r - tr)//2, 0)
    sc = max((c - tc)//2, 0)
    return A[sr:sr+tr, sc:sc+tc]

def compare(dir_out, dir_ref, name, allow_crop=True, label=""):
    out_p = os.path.join(dir_out, name)
    ref_p = os.path.join(dir_ref, name)
    if not os.path.isfile(out_p) or not os.path.isfile(ref_p):
        return "SKIP", 0, 0, f"{label}: SKIP (no encontrado OUT/REF)"
    try:
        A = readToa(dir_out, name)
        B = readToa(dir_ref, name)
    except Exception as e:
        return "SKIP", 0, 0, f"{label}: SKIP ({e})"

    if A.shape != B.shape:
        if not allow_crop:
            return "SKIP", 0, 0, f"{label}: SKIP (shape OUT {A.shape} vs REF {B.shape})"
        tr, tc = min(A.shape[0], B.shape[0]), min(A.shape[1], B.shape[1])
        A = center_crop(A, (tr, tc))
        B = center_crop(B, (tr, tc))
        note = f" (recorte centrado {A.shape})"
    else:
        note = ""

    diff = np.abs(A - B)
    n_out = int(np.sum(diff > tol_abs))
    total = int(diff.size)
    status = "OK" if n_out < total * three_sigma_frac else "NOK"
    pct = 100.0 * n_out / max(1, total)
    return status, n_out, total, f"{label}: {status} | fuera tol: {n_out}/{total} ({pct:.3f}%){note}"

# ================== ISM (DN) ==================
print("\n" + "="*68)
print("=== E2E · ISM final (DN) — con recorte centrado si hace falta ===")
print(f"OUT ISM : {OUT_ISM_DIR}")
print(f"REF ISM : {REF_ISM_DIR}")
print("-"*68)

ok = nok = skip = 0
for b in bands:
    name = f"ism_toa_{b}.nc"
    st, n, tot, msg = compare(OUT_ISM_DIR, REF_ISM_DIR, name, allow_crop=True, label=f"{b:<8} ism_toa")
    print(msg)
    if st == "OK": ok += 1
    elif st == "NOK": nok += 1
    else: skip += 1

print("-"*68)
print(f"Totales ISM → OK: {ok} | NOK: {nok} | SKIP: {skip}")

# ================== L1B (radiancias) ==================
print("\n" + "="*68)
print("=== E2E · L1B (radiancias) — con recorte centrado si hace falta ===")
print(f"OUT L1B : {OUT_L1B_DIR}")
print(f"REF L1B : {REF_L1B_DIR}")
print("-"*68)

ok = nok = skip = 0
for b in bands:
    name = f"l1b_toa_{b}.nc"
    st, n, tot, msg = compare(OUT_L1B_DIR, REF_L1B_DIR, name, allow_crop=True, label=f"{b:<8} l1b_toa")
    print(msg)
    if st == "OK": ok += 1
    elif st == "NOK": nok += 1
    else: skip += 1

print("-"*68)
print(f"Totales L1B → OK: {ok} | NOK: {nok} | SKIP: {skip}")
print("="*68)
