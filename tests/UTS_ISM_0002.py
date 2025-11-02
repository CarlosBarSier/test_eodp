import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat

# --- IMPORTA TUS CLASES EXISTENTES ---
from ism.src.detectionPhase import detectionPhase
from ism.src.videoChainPhase import videoChainPhase

# raíz = carpeta test_eodp (sube desde tests/)
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# SALIDAS del Test 2 (puedes dejarlas aquí, o donde quieras)
OUTDIR = os.path.join(ROOT, "myoutput")
os.makedirs(OUTDIR, exist_ok=True)

# **ENTRADA**: salida de la ÓPTICA del Test 1
OPT_OUTDIR = r"C:\Users\Kolarov\Desktop\EarthObservation\Codes\S3\EODP-TS-ISM\myoutput"
BAND = "VNIR-0"
IRRAD_FILE = os.path.join(OUTDIR, "ism_toa_optical_" + BAND + ".mat")
d = loadmat(IRRAD_FILE)['toa']

os.makedirs(OUTDIR, exist_ok=True)

# --------- LOADER GENÉRICO SIN readToa ---------
def load_toa_irradiance(fname):
    """
    Carga la matriz de irradiancia post-óptica.
    Admite:
      - .mat: busca una variable 2D ('toa' si existe; si no, la primera 2D)
      - .npy/.npz: carga array 2D (o primera 2D si .npz)
    Devuelve ndarray 2D en mW/m^2 (tal y como suele guardarse tras óptica).
    """
    if fname.lower().endswith(".mat"):
        d = loadmat(fname)
        # intenta 'toa' primero
        cand = None
        if "toa" in d and isinstance(d["toa"], np.ndarray) and d["toa"].ndim == 2:
            cand = d["toa"]
        else:
            # coge la primera 2D "real"
            for k, v in d.items():
                if k.startswith("__"):
                    continue
                if isinstance(v, np.ndarray) and v.ndim == 2 and np.isrealobj(v):
                    cand = v
                    break
        if cand is None:
            raise ValueError(f"No encontré matriz 2D válida en {fname}")
        return np.array(cand, dtype=np.float64)

    if fname.lower().endswith(".npy"):
        arr = np.load(fname)
        if arr.ndim != 2:
            raise ValueError("El .npy no contiene un array 2D")
        return arr.astype(np.float64)

    if fname.lower().endswith(".npz"):
        npz = np.load(fname)
        for k in npz.files:
            v = npz[k]
            if isinstance(v, np.ndarray) and v.ndim == 2:
                return v.astype(np.float64)
        raise ValueError("El .npz no tiene ninguna matriz 2D")

    raise ValueError(f"Extensión no soportada: {fname}")

# --------- PLOTS RÁPIDOS (si no tienes utilidades) ---------
def save_img(mat, title, xlabel, ylabel, outdir, basename):
    plt.figure(figsize=(7,5))
    im = plt.imshow(mat, origin="upper")
    plt.title(title)
    plt.xlabel(xlabel); plt.ylabel(ylabel)
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, basename + ".png"), dpi=150)
    plt.close()

def save_cut_alt(mat, title, outdir, basename):
    ialt = mat.shape[0] // 2
    plt.figure(figsize=(7,4))
    plt.plot(mat[ialt, :])
    plt.title(f"{title} (ALT cut @ row {ialt})")
    plt.xlabel("act_columns")
    plt.ylabel("value")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, basename + f"_alt{ialt}.png"), dpi=150)
    plt.close()
    return ialt

# ================== TEST 2: DETECTION + VCU ==================
# 1) Carga irradiancia post-óptica (normalmente guardada en mW/m^2)
toa_irrad_mWm2 = load_toa_irradiance(IRRAD_FILE)

# 2) Detection
det = detectionPhase(AUXDIR, INDIR, OUTDIR)
# OJO: tu detectionPhase debe convertir internamente mW/m^2 -> W/m^2 (*1e-3)
toa_e = det.compute(toa_irrad_mWm2, BAND)   # devuelve e- tras PRNU/DS/bad-dead

# Plots y guardado (electrones)
save_img(toa_e, "TOA after detection [e-]", "act_columns", "alt_lines", OUTDIR, f"ism2_toa_e_{BAND}")
ialt = save_cut_alt(toa_e, "TOA after detection [e-]", OUTDIR, f"ism2_toa_e_{BAND}")
savemat(os.path.join(OUTDIR, f"ism2_toa_e_{BAND}.mat"), {"Ne": toa_e})

# 3) Video Chain
vcu = videoChainPhase(AUXDIR, INDIR, OUTDIR)
toa_dn = vcu.compute(toa_e, BAND)           # devuelve DN (uint)

# Plots y guardado (DN)
save_img(toa_dn, "TOA after VCU [DN]", "act_columns", "alt_lines", OUTDIR, f"ism2_toa_dn_{BAND}")
save_cut_alt(toa_dn, "TOA after VCU [DN]", OUTDIR, f"ism2_toa_dn_{BAND}")
savemat(os.path.join(OUTDIR, f"ism2_toa_dn_{BAND}.mat"), {"DN": toa_dn})

# 4) % saturados (para el informe del test)
bit_depth = det.ismConfig.bit_depth if hasattr(det, "ismConfig") else 12
levels = (2**bit_depth) - 1
sat_pct = float(np.mean(toa_dn >= levels) * 100.0)
print(f"[ISM-0002] {BAND}  Saturated pixels = {sat_pct:.3f}%")
savemat(os.path.join(OUTDIR, f"ism2_metrics_{BAND}.mat"),
        {"sat_pct": np.array([[sat_pct]], dtype=np.float64)})
