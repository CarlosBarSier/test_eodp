"""
UTS_ISM_0001 - Optical Path (ISRF & MTF)
Checks EODP-ALG-ISM-1010/1020/1030
"""

import os
from netCDF4 import Dataset
import numpy as np

# Rutas base (ajústalas a tu PC)
proj = r"C:\Users\Kolarov\Desktop\EarthObservation\Codes\S3\test_eodp"
prod_dir = os.path.join(proj, "outputs", "ism")
ref_dir  = os.path.join(proj, "auxiliary", "EODP-TS-ISM", "reference")

# Lista de bandas a verificar
bands = ["VNIR-0", "VNIR-1", "VNIR-2", "VNIR-3"]

def read_matrix(fn):
    ds = Dataset(fn)
    for k in ds.variables:
        if ds.variables[k].ndim == 2:
            return ds.variables[k][:].astype("float64")
    raise RuntimeError(f"No 2D var found in {fn}")

def compare(tag):
    f_prod = os.path.join(prod_dir, f"{tag}.nc")
    f_ref  = os.path.join(ref_dir,  f"{tag}.nc")
    A = read_matrix(f_prod)
    B = read_matrix(f_ref)
    diff = A - B
    mae = np.nanmean(np.abs(diff))
    rms = np.sqrt(np.nanmean(diff**2))
    rel3sigma = (3*np.nanstd(diff)) / (np.nanmean(np.abs(B))+1e-12) * 100
    passed = rel3sigma < 0.01
    print(f"{tag}: PASS={passed} | MAE={mae:.3e}, RMSE={rms:.3e}, 3σ_rel={rel3sigma:.4e}%")
    return passed

if __name__ == "__main__":
    all_pass = True
    for b in bands:
        for stage in ["ism_toa_isrf_", "ism_toa_optical_"]:
            tag = stage + b
            ok = compare(tag)
            if not ok:
                all_pass = False
    print("=== RESULT: UTS_ISM_0001", "PASS" if all_pass else "FAIL", "===")