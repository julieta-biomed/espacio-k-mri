"""
mri.py — Fantasma cerebral sintético con verdad conocida, ecuación de señal
spin-echo, operaciones de espacio-k y métricas de calidad / detectabilidad.

Todo el módulo trabaja en unidades físicas explícitas:
  - T1, T2, TR, TE en milisegundos
  - FOV en milímetros
  - El fantasma es una etiqueta entera por píxel; los mapas PD/T1/T2 se
    construyen a partir de la tabla TEJIDOS.

Autor: Julieta Enríquez Fernández
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Parámetros de tejido a 1.5 T (valores nominales de literatura, usados como
# ENTRADA del modelo; todas las afirmaciones del artículo se miden sobre las
# imágenes que este modelo produce, no se toman de la literatura).
# PD relativa al agua libre.
# ---------------------------------------------------------------------------
TEJIDOS = {
    #  etiqueta: (nombre,              PD,    T1 ms,  T2 ms)
    0: ("aire",                        0.00,     1.0,    1.0),
    1: ("grasa subcutanea",            1.00,   250.0,   70.0),
    2: ("hueso cortical",              0.10,   500.0,   20.0),
    3: ("liquido cefalorraquideo",     1.00,  4000.0, 2000.0),
    4: ("materia gris",                0.85,  1000.0,  100.0),
    5: ("materia blanca",              0.75,   650.0,   80.0),
    6: ("lesion (edema)",              0.90,  1200.0,  150.0),
}

NOMBRES = {k: v[0] for k, v in TEJIDOS.items()}

# Protocolos clínicos típicos (TR, TE) en ms
PROTOCOLOS = {
    "T1w":  (500.0, 14.0),
    "DPw":  (2500.0, 14.0),   # densidad protonica
    "T2w":  (4000.0, 100.0),
}


# ---------------------------------------------------------------------------
# Fantasma
# ---------------------------------------------------------------------------
def _elipse(N, cy, cx, ry, rx, theta=0.0):
    """Máscara booleana de una elipse. Coordenadas en píxeles."""
    y, x = np.ogrid[:N, :N]
    yy = (y - cy).astype(float)
    xx = (x - cx).astype(float)
    ct, st = np.cos(theta), np.sin(theta)
    yr = yy * ct + xx * st
    xr = -yy * st + xx * ct
    return (yr / ry) ** 2 + (xr / rx) ** 2 <= 1.0


def fantasma_etiquetas(N=256, lesion_chica=True, lesion_grande=True,
                       fov_mm=240.0, d_chica_mm=3.0, d_grande_mm=12.0):
    """
    Corte axial esquemático de cabeza, etiquetado por tejido.

    Devuelve un array (N, N) de enteros con las etiquetas de TEJIDOS.
    Las lesiones se pueden encender/apagar de forma independiente: eso permite
    construir pares presente/ausente para la tarea de detección.
    """
    mm_px = fov_mm / N
    lab = np.zeros((N, N), dtype=np.int16)
    c = N / 2.0

    # Cuero cabelludo (grasa) -> cráneo -> espacio subaracnoideo -> parénquima
    lab[_elipse(N, c, c, 0.42 * N, 0.36 * N)] = 1          # grasa
    lab[_elipse(N, c, c, 0.400 * N, 0.340 * N)] = 2        # hueso
    lab[_elipse(N, c, c, 0.378 * N, 0.318 * N)] = 3        # LCR
    lab[_elipse(N, c, c, 0.368 * N, 0.308 * N)] = 4        # materia gris

    # Materia blanca central
    lab[_elipse(N, c, c, 0.300 * N, 0.240 * N)] = 5

    # Ventrículos laterales (LCR)
    lab[_elipse(N, c - 0.02 * N, c - 0.085 * N, 0.115 * N, 0.040 * N, 0.22)] = 3
    lab[_elipse(N, c - 0.02 * N, c + 0.085 * N, 0.115 * N, 0.040 * N, -0.22)] = 3

    # Núcleos grises profundos: dan estructura interna al parénquima
    lab[_elipse(N, c + 0.02 * N, c - 0.150 * N, 0.070 * N, 0.030 * N, 0.15)] = 4
    lab[_elipse(N, c + 0.02 * N, c + 0.150 * N, 0.070 * N, 0.030 * N, -0.15)] = 4

    # Lesiones, ambas rodeadas de materia blanca homogénea
    if lesion_grande:
        r = (d_grande_mm / 2.0) / mm_px
        lab[_elipse(N, c + 0.135 * N, c - 0.060 * N, r, r)] = 6
    if lesion_chica:
        r = (d_chica_mm / 2.0) / mm_px
        lab[_elipse(N, c - 0.150 * N, c + 0.045 * N, r, r)] = 6

    return lab


def centros_lesiones(N=256, fov_mm=240.0, d_chica_mm=3.0, d_grande_mm=12.0):
    """Centros (fila, columna) y radios en píxeles de ambas lesiones."""
    mm_px = fov_mm / N
    c = N / 2.0
    return {
        "grande": (c + 0.135 * N, c - 0.060 * N, (d_grande_mm / 2.0) / mm_px),
        "chica":  (c - 0.150 * N, c + 0.045 * N, (d_chica_mm / 2.0) / mm_px),
    }


def mapas_tejido(lab, t2_lesion=None, t1_lesion=None, pd_lesion=None):
    """
    Convierte etiquetas en mapas continuos de PD, T1 y T2.

    Los parámetros de la lesión (etiqueta 6) se pueden sobrescribir: así se
    controla la DIFICULTAD de la tarea por contraste, que es lo que varía
    clínicamente, en vez de subir el ruido a niveles irreales.
    """
    pd = np.zeros(lab.shape, float)
    t1 = np.ones(lab.shape, float)
    t2 = np.ones(lab.shape, float)
    for k, (_, p, a, b) in TEJIDOS.items():
        m = lab == k
        pd[m], t1[m], t2[m] = p, a, b
    m = lab == 6
    if pd_lesion is not None:
        pd[m] = pd_lesion
    if t1_lesion is not None:
        t1[m] = t1_lesion
    if t2_lesion is not None:
        t2[m] = t2_lesion
    return pd, t1, t2


def senal_spin_echo(pd, t1, t2, TR, TE):
    """
    Magnitud de la señal de un spin-echo estándar:
        S = PD * (1 - exp(-TR/T1)) * exp(-TE/T2)
    """
    return pd * (1.0 - np.exp(-TR / t1)) * np.exp(-TE / t2)


def imagen_protocolo(lab, protocolo="T2w", **kw_lesion):
    TR, TE = PROTOCOLOS[protocolo]
    pd, t1, t2 = mapas_tejido(lab, **kw_lesion)
    return senal_spin_echo(pd, t1, t2, TR, TE)


def sigma_para_snr(snr_objetivo, img, lab, etiqueta=5, N=None):
    """
    Sigma de espacio-k que da el SNR por píxel pedido en el tejido indicado,
    con muestreo completo. Para una IFFT2 normalizada en 1/N^2 el ruido de
    imagen tiene desviación sigma/N por componente.
    """
    N = N or img.shape[0]
    return float(img[lab == etiqueta].mean() * N / snr_objetivo)


# ---------------------------------------------------------------------------
# Espacio-k
# ---------------------------------------------------------------------------
def a_kspace(img):
    """Imagen -> espacio-k, con la frecuencia cero al centro."""
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(img)))


def a_imagen(k):
    """Espacio-k (centrado) -> imagen compleja."""
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(k)))


def ruido_kspace(shape, sigma, rng):
    """Ruido gaussiano complejo, i.i.d. por muestra de espacio-k."""
    return (rng.normal(0, sigma, shape) + 1j * rng.normal(0, sigma, shape))


def mascara_lineas_centrales(N, n_lineas, eje=0):
    """
    Máscara que conserva las n_lineas centrales del eje indicado
    (eje=0 -> codificación de fase por convención de este repo).
    """
    m = np.zeros((N, N), bool)
    c = N // 2
    lo = c - n_lineas // 2
    hi = lo + n_lineas
    if eje == 0:
        m[lo:hi, :] = True
    else:
        m[:, lo:hi] = True
    return m


def mascara_anillo(N, r_int, r_ext):
    """Máscara anular en el espacio-k, radios en píxeles."""
    c = N // 2
    y, x = np.ogrid[:N, :N]
    r = np.hypot(y - c, x - c)
    return (r >= r_int) & (r < r_ext)


def radio_k(N):
    c = N // 2
    y, x = np.ogrid[:N, :N]
    return np.hypot(y - c, x - c)


def energia_acumulada_radial(k, n_bins=128):
    """Fracción acumulada de |k|^2 dentro de un radio dado."""
    N = k.shape[0]
    r = radio_k(N)
    e = np.abs(k) ** 2
    r_max = N / 2.0
    bordes = np.linspace(0, r_max, n_bins + 1)
    idx = np.clip(np.digitize(r.ravel(), bordes) - 1, 0, n_bins - 1)
    acum = np.bincount(idx, weights=e.ravel(), minlength=n_bins).cumsum()
    return bordes[1:], acum / e.sum()


# ---------------------------------------------------------------------------
# Métricas
# ---------------------------------------------------------------------------
def nrmse(ref, test):
    return float(np.sqrt(np.mean((ref - test) ** 2)) / (ref.max() - ref.min()))


def psnr(ref, test):
    mse = np.mean((ref - test) ** 2)
    if mse == 0:
        return float("inf")   # reconstrucción idéntica a la referencia
    return float(10 * np.log10((ref.max() - ref.min()) ** 2 / mse))


def contraste(img, lab, a, b):
    """Contraste relativo de Michelson entre dos tejidos."""
    ma, mb = img[lab == a].mean(), img[lab == b].mean()
    return float((ma - mb) / (ma + mb))


def cnr(img, lab, a, b, ruido_std):
    return float(abs(img[lab == a].mean() - img[lab == b].mean()) / ruido_std)


def disco(N, cy, cx, r):
    y, x = np.ogrid[:N, :N]
    return np.hypot(y - cy, x - cx) <= r


def observador_roi(img, cy, cx, r, factor_fondo=(2.0, 3.2)):
    """
    Observador tipo radiólogo, no ideal: diferencia entre la media en un disco
    centrado en la lesión y la media en un anillo de fondo inmediato.
    No prewhitening, no conoce la correlación del ruido.
    """
    N = img.shape[0]
    d = disco(N, cy, cx, r)
    fondo = disco(N, cy, cx, r * factor_fondo[1]) & ~disco(N, cy, cx, r * factor_fondo[0])
    return float(img[d].mean() - img[fondo].mean())


def auc(scores_pos, scores_neg):
    """AUC por el estadístico de Mann-Whitney (exacto, sin binning)."""
    pos = np.asarray(scores_pos)
    neg = np.asarray(scores_neg)
    todos = np.concatenate([pos, neg])
    rangos = todos.argsort().argsort() + 1.0
    # promediar rangos empatados
    orden = np.sort(todos)
    for v in np.unique(orden):
        m = todos == v
        if m.sum() > 1:
            rangos[m] = rangos[m].mean()
    n1, n2 = len(pos), len(neg)
    r1 = rangos[:n1].sum()
    return float((r1 - n1 * (n1 + 1) / 2) / (n1 * n2))


def ic_auc(scores_pos, scores_neg, n_boot=2000, semilla=0):
    """Intervalo de confianza al 95 % del AUC por bootstrap no paramétrico."""
    rng = np.random.default_rng(semilla)
    pos = np.asarray(scores_pos)
    neg = np.asarray(scores_neg)
    vals = np.empty(n_boot)
    for i in range(n_boot):
        a = rng.choice(pos, len(pos), replace=True)
        b = rng.choice(neg, len(neg), replace=True)
        vals[i] = auc(a, b)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def ssim(ref, test, L=None, k1=0.01, k2=0.03, ventana=7):
    """
    SSIM global con ventana uniforme (implementación propia para que el repo
    no dependa de scikit-image en la ruta crítica).
    """
    from scipy.ndimage import uniform_filter
    if L is None:
        L = ref.max() - ref.min()
    c1, c2 = (k1 * L) ** 2, (k2 * L) ** 2
    mu_x = uniform_filter(ref, ventana)
    mu_y = uniform_filter(test, ventana)
    sxx = uniform_filter(ref * ref, ventana) - mu_x ** 2
    syy = uniform_filter(test * test, ventana) - mu_y ** 2
    sxy = uniform_filter(ref * test, ventana) - mu_x * mu_y
    num = (2 * mu_x * mu_y + c1) * (2 * sxy + c2)
    den = (mu_x ** 2 + mu_y ** 2 + c1) * (sxx + syy + c2)
    return float(np.mean(num / den))


# ---------------------------------------------------------------------------
# Tiempo de adquisición
# ---------------------------------------------------------------------------
def tiempo_adquisicion_s(TR_ms, n_lineas_fase, NEX=1, n_cortes=1,
                         cortes_entrelazados=True):
    """
    Tiempo de una secuencia spin-echo convencional.
    Si los cortes se entrelazan dentro del mismo TR (práctica estándar),
    n_cortes no multiplica el tiempo mientras quepan en el TR.
    """
    base = TR_ms * n_lineas_fase * NEX / 1000.0
    return base if cortes_entrelazados else base * n_cortes
