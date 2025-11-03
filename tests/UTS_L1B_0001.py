# tests/UTS_L1B_0001.py  — L1B Equalization & Restoration (con plots y múltiples carpetas)
import os, sys
import numpy as np
import matplotlib.pyplot as plt

# --- raíz del repo (este script está en tests/) ---
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from config.globalConfig import globalConfig
from common.io.writeToa import readToa

# ====== RUTAS: ajusta estas cuatro según tu máquina ======
IN_DIR        = r"C:\\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-L1B\input"           # contiene ism_toa_isrf_*.nc
OUT_EQ_DIR    = r"C:\\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-L1B\myoutputs_eq"    # salida L1B con equalization=True
OUT_NOEQ_DIR  = r"C:\\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-L1B\myoutputs_noeq"  # salida L1B con equalization=False (opcional)
REF_DIR       = r"C:\\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-L1B\output"          # referencias oficiales

# ====== Parámetros de validación (según especificación) ======
tol_abs = 1.0e-4               # 0.01% en magnitud de TOA en los datasets de referencia
three_sigma_frac = 1 - 0.997   # 3σ

cfg = globalConfig()
bands = list(cfg.bands)

print("\n" + "="*70)
print("=== TEST L1B-0001 (Equalization & Restoration) ===")
print(f"IN:  {IN_DIR}")
print(f"OUT(eq):    {OUT_EQ_DIR}")
print(f"OUT(noeq):  {OUT_NOEQ_DIR}  (usado si existe)")
print(f"REF: {REF_DIR}")
print("-"*70)

def compare(prod_name, out_dir, band):
    """Compara un producto (p.ej. l1b_toa_) contra la referencia con criterio 3σ."""
    out_name = f"{prod_name}_{band}.nc"
    out_p = os.path.join(out_dir, out_name)
    ref_p = os.path.join(REF_DIR, out_name)
    if not (os.path.isfile(out_p) and os.path.isfile(ref_p)):
        return "SKIP", 0, 0
    A = readToa(out_dir, out_name)
    B = readToa(REF_DIR, out_name)
    diff = np.abs(A - B)
    n_out = int(np.sum(diff > tol_abs))
    total = diff.size
    status = "OK" if n_out < total * three_sigma_frac else "NOK"
    return status, n_out, total

def plot_midrow_curve(y, label, out_png, title, xlab="Pixel ACT [-]", ylab="Radiance [mW/m²/sr]"):
    """Guarda plot simple de una línea en el ACT de la fila central."""
    plt.figure()
    plt.plot(y, label=label)
    plt.title(title)
    plt.xlabel(xlab); plt.ylabel(ylab); plt.grid(True)
    if label: plt.legend()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()

# Totales
ok = nok = skip = 0

for b in bands:
    print(f"\n>>> Banda {b}")
    # 1) Comparación cuantitativa con referencia (ambos productos)
    for prod in ["l1b_toa_eq", "l1b_toa"]:
        st, n, tot = compare(prod, OUT_EQ_DIR, b)
        if st != "SKIP":
            pct = 100*n/max(1,tot)
            print(f"  {prod:<12} -> {st:>3s} | fuera tol: {n}/{tot} ({pct:.3f}%)")
            ok += (st=="OK"); nok += (st=="NOK")
        else:
            print(f"  {prod:<12} -> SKIP (no encontrado)")
            skip += 1

    # 2) Plots exigidas por la guía: l1b_toa (EQ) vs ism_toa_isrf (input del test) en la fila ALT central
    try:
        toa_l1b = readToa(OUT_EQ_DIR, f"l1b_toa_{b}.nc")         # radiancias restauradas
        toa_isrf = readToa(IN_DIR,     f"ism_toa_isrf_{b}.nc")    # tras ISRF (irradiancias/energía en detector vs radiancia final — ver guía)
        mid = toa_l1b.shape[0] // 2
        y_l1b = toa_l1b[mid, :].astype(float)
        # Para comparar formas, si las unidades difieren, puedes normalizar:
        # y_isrf = (toa_isrf[mid, :] * scale) si quisieras escalar. Para la guía, la comparación es cualitativa del perfil.
        y_isrf = toa_isrf[mid, :].astype(float)

        out_png = os.path.join(OUT_EQ_DIR, f"plot_L1B_vs_ISRF_{b}.png")
        plt.figure()
        plt.plot(y_l1b, label="l1b_toa (EQ)")
        plt.plot(y_isrf, label="ism_toa_isrf")
        plt.title(f"Banda {b} — Perfil ALT central: L1B restaurado vs ISRF")
        plt.xlabel("Pixel ACT [-]"); plt.ylabel("Arbitrary / Comparable units"); plt.grid(True); plt.legend()
        os.makedirs(os.path.dirname(out_png), exist_ok=True)
        plt.savefig(out_png, dpi=150, bbox_inches="tight"); plt.close()
        print(f"  Plot guardado: {out_png}")
    except Exception as e:
        print(f"  [WARN] No se pudo graficar L1B vs ISRF para {b}: {e}")

    # 3) Comparación EQ vs NOEQ (si existe OUT_NOEQ_DIR)
    try:
        if os.path.isdir(OUT_NOEQ_DIR):
            toa_eq  = readToa(OUT_EQ_DIR,   f"l1b_toa_{b}.nc")
            toa_noq = readToa(OUT_NOEQ_DIR, f"l1b_toa_{b}.nc")
            mid = toa_eq.shape[0] // 2
            y_eq  = toa_eq[mid, :].astype(float)
            y_noq = toa_noq[mid, :].astype(float)
            out_png = os.path.join(OUT_EQ_DIR, f"plot_EQ_vs_NOEQ_{b}.png")
            plt.figure()
            plt.plot(y_eq,  label="l1b_toa (EQ=True)")
            plt.plot(y_noq, label="l1b_toa (EQ=False)")
            plt.title(f"Banda {b} — Perfil ALT central: EQ vs NOEQ")
            plt.xlabel("Pixel ACT [-]"); plt.ylabel("Radiance [mW/m²/sr]"); plt.grid(True); plt.legend()
            os.makedirs(os.path.dirname(out_png), exist_ok=True)
            plt.savefig(out_png, dpi=150, bbox_inches="tight"); plt.close()
            print(f"  Plot guardado: {out_png}")
        else:
            print("  (Sin OUT_NOEQ_DIR: omito plot EQ vs NOEQ)")
    except Exception as e:
        print(f"  [WARN] No se pudo graficar EQ vs NOEQ para {b}: {e}")

print("\n" + "-"*70)
print(f"Totales → OK: {ok} | NOK: {nok} | SKIP: {skip}")
print("="*70)
