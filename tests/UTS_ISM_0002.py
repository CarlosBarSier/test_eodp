import os, sys, glob
import numpy as np

# --- Añadimos la raíz del proyecto al path (sube desde tests/)
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# Import correctos con tu estructura
from common.io.writeToa import readToa
from config.globalConfig import globalConfig


# ====== RUTAS CORRECTAS, DERIVADAS DEL PROYECTO ======
# --- Rutas absolutas correctas (usa slashes) ---
OUT_DIR = "C:/Users/Kolarov/Desktop/EarthObservation/Codes/S3/EODP-TS-ISM/myoutput"
REFERENCE_DIR = "C:/Users/Kolarov/Desktop/EarthObservation/Codes/S3/EODP-TS-ISM/output"


# Permite sobreescribir con variables de entorno si quieres
OUT_DIR = os.environ.get("ISM_OUT_DIR", OUT_DIR)
REFERENCE_DIR = os.environ.get("ISM_REF_DIR", REFERENCE_DIR)

# Normaliza a forward slashes (Windows friendly)

print("=== ISM-0002 · Detection & Video Conversion ===")
print(f"REF: {REFERENCE_DIR}")
print(f"OUT: {OUT_DIR}")
print("------------------------------------------------")

# ==============================================================
# Parámetros del test (spec ISM-0002)
# ==============================================================
cfg = globalConfig()
bands = cfg.bands                     # p.ej. ["VNIR-0","VNIR-1","VNIR-2","VNIR-3"]

tol_abs = 1.0e-4                      # 0.01% típico del dataset (en unidades del fichero)
three_sigma_frac = 1 - 0.997          # 3σ ~ 99.7% -> umbral puntos fuera

def check_pair(outdir, refdir, filename, tol):
    """
    Lee OUT y REF, devuelve (#puntos fuera, total, boolean OK).
    """
    # Comprobar existencia
    f_out = os.path.join(outdir, filename)
    f_ref = os.path.join(refdir, filename)
    if not os.path.isfile(f_out):
        raise FileNotFoundError(f"No existe salida: {f_out}")
    if not os.path.isfile(f_ref):
        raise FileNotFoundError(f"No existe referencia: {f_ref}")

    a = readToa(outdir, filename)          # numpy 2D
    b = readToa(refdir, filename)          # numpy 2D

    # Diferencias absolutas
    diff = np.abs(a - b)
    n_out = int(np.sum(diff > tol))
    total = int(diff.size)
    ok = (n_out < total * three_sigma_frac)
    return n_out, total, ok

def pretty(filename, ok, n_out, total):
    status = "OK" if ok else "NOK"
    pct = 100.0 * n_out / max(1, total)
    return f"{filename}: {status}  |  fuera tol: {n_out}/{total} ({pct:.3f}%)"

def main():
    print("=== ISM-0002 · Detection & Video Conversion ===")
    print(f"REF: {REFERENCE_DIR}")
    print(f"OUT: {OUT_DIR}")
    print("------------------------------------------------")

    # —— Ficheros de DETECCIÓN (según tu pipeline):
    #     - ism_toa_e_<BAND>.nc           (electrones tras PRNU/DS/bad/dead)
    #     - ism_toa_detection_<BAND>.nc   (solo conversión + efectos agregados)
    #     - ism_toa_ds_<BAND>.nc          (dark signal)
    #     - ism_toa_prnu_<BAND>.nc        (prnu)
    det_lists = {
        "ism_toa_e":         [f"ism_toa_e_{b}.nc" for b in bands],
        "ism_toa_detection": [f"ism_toa_detection_{b}.nc" for b in bands],
        "ism_toa_ds":        [f"ism_toa_ds_{b}.nc" for b in bands],
        "ism_toa_prnu":      [f"ism_toa_prnu_{b}.nc" for b in bands],
    }

    print("— Detection stage")
    for label, namelist in det_lists.items():
        for fname in namelist:
            try:
                n_out, total, ok = check_pair(OUT_DIR, REFERENCE_DIR, fname, tol_abs)
                print(pretty(fname, ok, n_out, total))
            except FileNotFoundError as e:
                print(f"{fname}: SKIP ({e})")
        print("------------------------------------------------")

    # —— Fichero final de VCU (DN):
    #     tu cadena guarda el producto final con el prefijo global "ism_toa_"
    #     (escrito en ism.processModule -> writeToa(globalConfig.ism_toa + band, ...))
    print("— Video Chain (VCU) – producto final en DN")
    for b in bands:
        fname = f"ism_toa_{b}.nc"   # salida final (DN)
        try:
            n_out, total, ok = check_pair(OUT_DIR, REFERENCE_DIR, fname, tol_abs)
            print(pretty(fname, ok, n_out, total))
        except FileNotFoundError as e:
            print(f"{fname}: SKIP ({e})")
    print("------------------------------------------------")

    # —— Métrica adicional pedida por la spec: % saturados por banda (en DN)
    print("— Saturación (DN)")
    levels = (2**12) - 1  # si quieres, lee bit depth de tu config de ismConfig
    for b in bands:
        fname = f"ism_toa_{b}.nc"
        try:
            dn = readToa(OUT_DIR, fname)
            sat_pct = float(np.mean(dn >= levels) * 100.0)
            print(f"{fname}: pixels saturados = {sat_pct:.3f}% (levels_max={levels})")
        except FileNotFoundError as e:
            print(f"{fname}: SKIP ({e})")

if __name__ == "__main__":
    main()

# === RESUMEN SIMPLE EN TERMINAL ===
print("\n" + "="*60)
print("=== RESUMEN TEST ISM-0002 ===")
ok = nok = skip = 0

# recorrer todos los archivos procesados y contar OK/NOK/SKIP
for b in bands:
    for prefix in ["ism_toa_e", "ism_toa_detection", "ism_toa_ds", "ism_toa_prnu", "ism_toa"]:
        fname = f"{prefix}_{b}.nc"
        fpath = os.path.join(OUT_DIR, fname)
        if not os.path.isfile(fpath):
            print(f"{b:<8} {prefix:<20} -> SKIP (no encontrado)")
            skip += 1
            continue
        try:
            out_arr = readToa(OUT_DIR, fname)
            ref_arr = readToa(REFERENCE_DIR, fname)
            diff = np.abs(out_arr - ref_arr)
            n_out = int(np.sum(diff > tol_abs))
            total = diff.size
            status = "OK" if n_out < total * three_sigma_frac else "NOK"
            pct = 100 * n_out / total
            print(f"{b:<8} {prefix:<20} -> {status:>3s} | fuera tol: {n_out}/{total} ({pct:.3f}%)")
            if status == "OK": ok += 1
            else: nok += 1
        except Exception as e:
            print(f"{b:<8} {prefix:<20} -> ERROR ({e})")
            skip += 1

print("-"*60)
print(f"Totales → OK: {ok} | NOK: {nok} | SKIP/ERROR: {skip}")
print("="*60)