"""
EODP-TS-ISM-0001 — OPTICAL STAGE
Comparación entre resultados generados (myoutput) y dataset de referencia (output)
"""

import os
import numpy as np
from netCDF4 import Dataset

# === CONFIGURACIÓN ===
# Rutas exactas de tu entorno
ref_dir = r"C:\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-ISM\output"
out_dir = r"C:\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-ISM\myoutput"
bands = ["VNIR-0", "VNIR-1", "VNIR-2", "VNIR-3"]

# === FUNCIONES AUXILIARES ===

def read_nc(filepath):
    """Lee un archivo .nc y devuelve el array principal"""
    with Dataset(filepath, mode="r") as ds:
        # asume que el primer campo 2D o 3D del dataset es la imagen
        for var in ds.variables.values():
            if var.ndim >= 2:
                return np.array(var[:])
    raise ValueError(f"No se encontró variable válida en {filepath}")

def pct_rel_err(a, b, eps=1e-12):
    a = np.asarray(a)
    b = np.asarray(b)
    denom = np.maximum(np.abs(b), eps)
    return 100.0 * np.abs(a - b) / denom

def assert_within_threshold(actual, reference, pct_threshold=0.01, coverage=0.997):
    per = pct_rel_err(actual, reference)
    finite = np.isfinite(per)
    per = per[finite]
    ok = (per <= pct_threshold)
    frac = ok.mean()
    passed = frac >= coverage
    return passed, frac * 100.0

# === TEST ===

def run_test_optical_stage():
    print("=== EODP-TS-ISM-0001 — OPTICAL STAGE ===")
    all_passed = True
    for band in bands:
        print(f"\n--- {band} ---")
        ref_isrf = os.path.join(ref_dir, f"ism_toa_isrf_{band}.nc")
        out_isrf = os.path.join(out_dir, f"ism_toa_isrf_{band}.nc")
        ref_opt = os.path.join(ref_dir, f"ism_toa_optical_{band}.nc")
        out_opt = os.path.join(out_dir, f"ism_toa_optical_{band}.nc")

        # --- ISRF ---
        if not os.path.exists(ref_isrf) or not os.path.exists(out_isrf):
            print(f"[ISRF] ❌ Faltan archivos: {ref_isrf} o {out_isrf}")
            all_passed = False
        else:
            ref = read_nc(ref_isrf)
            got = read_nc(out_isrf)
            passed, frac = assert_within_threshold(got, ref)
            print(f"[ISRF] {'✅' if passed else '❌'} {frac:.3f}% de los puntos dentro del 0.01% permitido")
            all_passed &= passed

        # --- OPTICAL ---
        if not os.path.exists(ref_opt) or not os.path.exists(out_opt):
            print(f"[OPTICAL] ❌ Faltan archivos: {ref_opt} o {out_opt}")
            all_passed = False
        else:
            ref = read_nc(ref_opt)
            got = read_nc(out_opt)
            passed, frac = assert_within_threshold(got, ref)
            print(f"[OPTICAL] {'✅' if passed else '❌'} {frac:.3f}% de los puntos dentro del 0.01% permitido")
            all_passed &= passed

    if all_passed:
        print("\n✅ TEST EODP-TS-ISM-0001 PASADO")
    else:
        print("\n❌ TEST EODP-TS-ISM-0001 FALLÓ")

# === EJECUCIÓN DIRECTA ===
if __name__ == "__main__":
    run_test_optical_stage()
