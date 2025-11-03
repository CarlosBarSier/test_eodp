import os, sys
import numpy as np
import matplotlib.pyplot as plt

# --- root del repo (este script está en tests/) ---
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from config.globalConfig import globalConfig
from common.io.writeToa import readToa

# ====== RUTAS (alineadas con tus mains E2E) ======
# ISM E2E: coincide con mainIsm.py
OUT_ISM_DIR = "C:/Users/Kolarov/Desktop/EarthObservation/Codes/S3/EODP-TS-ISM/myoutput_E2E".replace("\\","/")
# L1B E2E: coincide con mainL1b.py
OUT_L1B_DIR = "C:/Users/Kolarov/Desktop/EarthObservation/Codes/S3/EODP-TS-L1B/myoutputs_E2E".replace("\\","/")

# Referencias E2E (carpetas estándar del paquete E2E)
REF_E2E_ROOT = "C:/Users/Kolarov/Desktop/EarthObservation/Codes/S3/EODP-TS-E2E".replace("\\","/")
REF_ISM_DIR  = f"{REF_E2E_ROOT}/ism_out"
REF_L1B_DIR  = f"{REF_E2E_ROOT}/l1b_out"

# (Opcional) Si has corrido el GM por línea de comandos y quieres graficar DEM/latitudes:
GM_DIR = None  # por ejemplo: "C:/.../EODP-TS-E2E/gm_out". Déjalo en None si no lo tienes.

# ====== Parámetros ======
bands = list(globalConfig().bands)
tol_abs = 1.0e-4               # tolerancia absoluta (~0.01% en las referencias)
three_sigma_frac = 1 - 0.997   # 3σ (~0.3% de puntos pueden irse)
PLOT_DIR = os.path.join(OUT_L1B_DIR, "plots_e2e")  # donde guardamos PNGs

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
        return "SKIP", 0, 0, f"{label}: SKIP (no encontrado OUT/REF: {out_p} | {ref_p})"
    try:
        A = readToa(dir_out, name)
        B = readToa(dir_ref, name)
    except Exception as e:
        return "SKIP", 0, 0, f"{label}: SKIP ({e})"

    note = ""
    if A.shape != B.shape:
        if not allow_crop:
            return "SKIP", 0, 0, f"{label}: SKIP (shape OUT {A.shape} vs REF {B.shape})"
        tr, tc = min(A.shape[0], B.shape[0]), min(A.shape[1], B.shape[1])
        A = center_crop(A, (tr, tc))
        B = center_crop(B, (tr, tc))
        note = f" (recorte centrado {A.shape})"

    diff = np.abs(A - B)
    n_out = int(np.sum(diff > tol_abs))
    total = int(diff.size)
    status = "OK" if n_out < total * three_sigma_frac else "NOK"
    pct = 100.0 * n_out / max(1, total)
    return status, n_out, total, f"{label}: {status} | fuera tol: {n_out}/{total} ({pct:.3f} %){note}"

def save_profile_plot(ylist, labels, title, out_png, xlab="Pixel ACT [-]", ylab="Radiance / arbitrary", legend=True):
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.figure()
    for y, lab in zip(ylist, labels):
        plt.plot(y, label=lab)
    plt.title(title)
    plt.xlabel(xlab); plt.ylabel(ylab); plt.grid(True)
    if legend: plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()

# ================== ISM (comparaciones y plots intermedios) ==================
print("\n" + "="*76)
print("=== E2E · ISM (incluye intermedios si existen) — con recorte centrado si hace falta ===")
print(f"OUT ISM : {OUT_ISM_DIR}")
print(f"REF ISM : {REF_ISM_DIR}")
print("-"*76)

ok = nok = skip = 0
for b in bands:
    # Final del ISM (DN)
    name_final = f"ism_toa_{b}.nc"
    st, n, tot, msg = compare(OUT_ISM_DIR, REF_ISM_DIR, name_final, allow_crop=True, label=f"{b:<8} ism_toa")
    print(msg)
    ok += (st=="OK"); nok += (st=="NOK"); skip += (st=="SKIP")

    # Intermedios ISM: ISRF y óptico (si existen) — solo informe (no hay REF oficial para ambos en E2E)
    for inter in ["ism_toa_isrf", "ism_toa_optical"]:
        f = os.path.join(OUT_ISM_DIR, f"{inter}_{b}.nc")
        if os.path.isfile(f):
            try:
                A = readToa(OUT_ISM_DIR, f"{inter}_{b}.nc")
                mid = A.shape[0]//2
                y = A[mid,:].astype(float)
                out_png = os.path.join(PLOT_DIR, f"{inter}_{b}.png")
                save_profile_plot([y], [inter], f"{b} — {inter} (perfil ALT central)", out_png,
                                  ylab="Irradiance / DN / comparable (según etapa)")
                print(f"  Plot intermedio guardado: {out_png}")
            except Exception as e:
                print(f"  [WARN] No se pudo plotear {inter} {b}: {e}")

print("-"*76)
print(f"Totales ISM → OK: {ok} | NOK: {nok} | SKIP: {skip}")

# ================== L1B (radiancias) ==================
print("\n" + "="*76)
print("=== E2E · L1B (radiancias) — con recorte centrado si hace falta ===")
print(f"OUT L1B : {OUT_L1B_DIR}")
print(f"REF L1B : {REF_L1B_DIR}")
print("-"*76)

ok = nok = skip = 0
for b in bands:
    # Validación cuantitativa contra referencia
    name_l1b = f"l1b_toa_{b}.nc"
    st, n, tot, msg = compare(OUT_L1B_DIR, REF_L1B_DIR, name_l1b, allow_crop=True, label=f"{b:<8} l1b_toa")
    print(msg)
    ok += (st=="OK"); nok += (st=="NOK"); skip += (st=="SKIP")

    # Plot requerido: L1B (radiancias) vs ISRF (del ISM) — perfil ALT central (comparación cualitativa)
    try:
        l1b = readToa(OUT_L1B_DIR, f"l1b_toa_{b}.nc")
        isrf = readToa(OUT_ISM_DIR, f"ism_toa_isrf_{b}.nc")
        mid = min(l1b.shape[0], isrf.shape[0])//2
        # recorte centrado si columnas difieren
        tc = min(l1b.shape[1], isrf.shape[1])
        y_l1b  = center_crop(l1b,  (2*mid+1, tc))[mid,:].astype(float)
        y_isrf = center_crop(isrf, (2*mid+1, tc))[mid,:].astype(float)
        out_png = os.path.join(PLOT_DIR, f"L1B_vs_ISRF_{b}.png")
        save_profile_plot([y_l1b, y_isrf], ["l1b_toa (radiances)", "ism_toa_isrf"],
                          f"{b} — Perfil ALT central: L1B vs ISRF", out_png,
                          ylab="Radiance / comparable")
        print(f"  Plot L1B vs ISRF guardado: {out_png}")
    except Exception as e:
        print(f"  [WARN] No se pudo plotear L1B vs ISRF para {b}: {e}")

print("-"*76)
print(f"Totales L1B → OK: {ok} | NOK: {nok} | SKIP: {skip}")
print("="*76)

# ================== (Opcional) GM: latitudes / DEM ==================
if GM_DIR:
    print("\n" + "="*76)
    print("=== E2E · GM (opcional) — plots de latitudes/altitudes DEM ===")
    try:
        # Ejemplos de nombres típicos; cámbialos si tus ficheros difieren
        lat = readToa(GM_DIR, "lat.nc")      # matriz 2D de latitudes
        alt = readToa(GM_DIR, "altitude.nc") # DEM interpolado
        os.makedirs(PLOT_DIR, exist_ok=True)
        # Latitudes (perfil central o imagen rápida)
        plt.figure(); plt.imshow(lat, aspect='auto'); plt.title("GM — Latitudes"); plt.colorbar();
        plt.savefig(os.path.join(PLOT_DIR, "GM_latitudes.png"), dpi=150); plt.close()
        # DEM altitudes
        plt.figure(); plt.imshow(alt, aspect='auto'); plt.title("GM — Altitudes DEM"); plt.colorbar();
        plt.savefig(os.path.join(PLOT_DIR, "GM_altitudes.png"), dpi=150); plt.close()
        print(f"  Plots GM guardados en: {PLOT_DIR}")
    except Exception as e:
        print(f"  [WARN] GM no ploteado: {e}")
