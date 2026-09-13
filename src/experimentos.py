"""
experimentos.py — Rutinas de medición usadas por el notebook y por los scripts
de figuras. Todas devuelven diccionarios o arrays de NÚMEROS; no dibujan nada.
"""

from __future__ import annotations

import numpy as np

import mri


# ---------------------------------------------------------------------------
def energia_radial(protocolo="T2w", N=256):
    lab = mri.fantasma_etiquetas(N)
    img = mri.imagen_protocolo(lab, protocolo)
    k = mri.a_kspace(img)
    r, acum = mri.energia_acumulada_radial(k)
    out = {"r_px": r, "acum": acum, "N": N}
    for frac in (0.01, 0.02, 0.05, 0.10, 0.25, 0.50):
        r_obj = frac * (N / 2)
        i = np.searchsorted(r, r_obj)
        out[f"energia_r{frac:.2f}"] = float(acum[min(i, len(acum) - 1)])
    # fracción de MUESTRAS dentro de cada radio, para comparar
    rr = mri.radio_k(N)
    out["muestras"] = {frac: float((rr <= frac * N / 2).mean())
                       for frac in (0.01, 0.02, 0.05, 0.10, 0.25, 0.50)}
    return out


# ---------------------------------------------------------------------------
def reconstruir(img_limpia, n_lineas, sigma, rng, N=None):
    """
    Adquisición simulada: pasa a espacio-k, añade ruido gaussiano complejo
    i.i.d. por muestra, conserva n_lineas centrales de codificación de fase,
    rellena con ceros y reconstruye la magnitud.
    """
    N = N or img_limpia.shape[0]
    k = mri.a_kspace(img_limpia)
    if sigma > 0:
        k = k + mri.ruido_kspace(k.shape, sigma, rng)
    m = mri.mascara_lineas_centrales(N, n_lineas, eje=0)
    return np.abs(mri.a_imagen(k * m))


def barrido_truncamiento(protocolo="T2w", N=256,
                         lineas=(256, 192, 160, 128, 96, 64, 48, 32, 24, 16)):
    """Métricas de imagen completa (sin ruido) frente a la referencia total."""
    lab = mri.fantasma_etiquetas(N)
    img = mri.imagen_protocolo(lab, protocolo)
    ref = np.abs(mri.a_imagen(mri.a_kspace(img)))
    rng = np.random.default_rng(0)
    filas = []
    for n in lineas:
        rec = reconstruir(img, n, 0.0, rng, N)
        filas.append({
            "lineas": n,
            "frac": n / N,
            "ssim": mri.ssim(ref, rec),
            "psnr": mri.psnr(ref, rec),
            "nrmse": mri.nrmse(ref, rec),
            "res_mm": 240.0 / n,
            "t_s": mri.tiempo_adquisicion_s(mri.PROTOCOLOS[protocolo][0], n),
        })
    return filas


# ---------------------------------------------------------------------------
def tarea_deteccion(cual="grande", protocolo="T2w", N=256, sigma=1.0,
                    n_real=300, lineas=(256, 192, 128, 96, 64, 48, 32, 24, 16),
                    semilla=1234, t2_lesion=None, pd_lesion=None):
    """
    Tarea sí/no de lesión conocida en posición conocida.
    Se generan n_real realizaciones de ruido con lesión y n_real sin lesión;
    el observador de ROI puntúa cada imagen y se calcula el AUC.
    """
    kw = dict(lesion_chica=(cual == "chica"), lesion_grande=(cual == "grande"))
    lab_pres = mri.fantasma_etiquetas(N, **kw)
    lab_aus = mri.fantasma_etiquetas(N, lesion_chica=False, lesion_grande=False)
    img_pres = mri.imagen_protocolo(lab_pres, protocolo, t2_lesion=t2_lesion,
                                    pd_lesion=pd_lesion)
    img_aus = mri.imagen_protocolo(lab_aus, protocolo)
    cy, cx, r = mri.centros_lesiones(N)[cual]

    rng = np.random.default_rng(semilla)
    filas = []
    for n in lineas:
        sp, sn = [], []
        for _ in range(n_real):
            sp.append(mri.observador_roi(reconstruir(img_pres, n, sigma, rng, N), cy, cx, r))
            sn.append(mri.observador_roi(reconstruir(img_aus, n, sigma, rng, N), cy, cx, r))
        sp, sn = np.array(sp), np.array(sn)
        # d' equivalente, útil para ver la separación en unidades de ruido
        dprime = (sp.mean() - sn.mean()) / np.sqrt(0.5 * (sp.var() + sn.var()))
        lo, hi = mri.ic_auc(sp, sn)
        filas.append({"lineas": n, "auc": mri.auc(sp, sn), "ic95": (lo, hi),
                      "dprime": float(dprime), "sd_fondo": float(sn.std())})
    return filas


def calibrar_t2_lesion(cual, sigma, auc_objetivo=0.90, protocolo="T2w", N=256,
                       n_real=250, lo=80.5, hi=400.0, iteraciones=14,
                       pd_lesion=0.75):
    """
    Busca por bisección el T2 de la lesión que deja la tarea en el AUC objetivo
    con muestreo completo. Se ajusta el CONTRASTE, no el ruido: una tarea
    saturada en AUC 1.0 no puede discriminar entre esquemas de muestreo, y
    subir el ruido hasta saturar la lleva a un régimen de magnitud Rician que
    no representa ninguna adquisición clínica real.
    (T2 de la materia blanca de fondo = 80 ms.)
    """
    for _ in range(iteraciones):
        mid = 0.5 * (lo + hi)
        a = tarea_deteccion(cual, protocolo, N, sigma, n_real, lineas=(N,),
                            semilla=7, t2_lesion=mid, pd_lesion=pd_lesion)[0]["auc"]
        if a > auc_objetivo:
            hi = mid
        else:
            lo = mid
    return float(0.5 * (lo + hi))


# ---------------------------------------------------------------------------
def gibbs_escalon(N=256, alturas=(1.0,), lineas=(256, 128, 64, 32, 16)):
    """
    Sobreimpulso de Gibbs medido en un escalón ideal 1D truncado en frecuencia.
    El valor teórico asintótico es 8.949 % de la altura del escalón.
    """
    x = np.zeros(N)
    x[N // 2:] = 1.0
    filas = []
    for n in lineas:
        X = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(x)))
        m = np.zeros(N, bool)
        c = N // 2
        m[c - n // 2:c - n // 2 + n] = True
        y = np.real(np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(X * m))))
        pico = y[N // 2:].max()
        filas.append({"lineas": n, "sobreimpulso_pct": float((pico - 1.0) * 100)})
    return filas


def gibbs_en_fantasma(protocolo="T2w", N=256, n_lineas=64):
    """Perfil a través de la interfaz LCR/parénquima, con y sin truncamiento."""
    lab = mri.fantasma_etiquetas(N)
    img = mri.imagen_protocolo(lab, protocolo)
    rng = np.random.default_rng(0)
    ref = np.abs(mri.a_imagen(mri.a_kspace(img)))
    rec = reconstruir(img, n_lineas, 0.0, rng, N)
    col = N // 2
    return {"ref": ref[:, col], "trunc": rec[:, col], "col": col}
