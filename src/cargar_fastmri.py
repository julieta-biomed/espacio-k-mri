"""
cargar_fastmri.py — Puente entre un archivo crudo de fastMRI y el pipeline de
este repositorio.

ESTE CÓDIGO NO SE HA EJECUTADO. El entorno donde se hicieron las mediciones del
artículo no tiene acceso de red a fastmri.med.nyu.edu, así que todos los
resultados publicados vienen del fantasma sintético. Esto queda escrito para
que cualquiera con los datos pueda repetir el experimento sobre anatomía real.

Cómo conseguir los datos
------------------------
1. Solicitar acceso en https://fastmri.med.nyu.edu/ (registro gratuito, se
   acepta un acuerdo de uso de datos).
2. Descargar el conjunto `knee_singlecoil_val` o `brain_multicoil_val`.
3. Dejar los `.h5` en la carpeta `datos/` de este repositorio.

Uso
---
    python src/cargar_fastmri.py datos/file1000000.h5 --corte 20

Estructura del archivo fastMRI
------------------------------
    kspace : (n_cortes, n_bobinas, n_fase, n_lectura) complejo, o
             (n_cortes, n_fase, n_lectura) en las adquisiciones de una bobina
    reconstruction_rss / reconstruction_esc : reconstrucción de referencia

Nota importante sobre la comparación
------------------------------------
Con datos reales NO hay verdad conocida: la referencia es la reconstrucción a
muestreo completo, que ya trae ruido. Por eso las métricas dejan de medir error
absoluto y pasan a medir diferencia entre dos reconstrucciones. Es exactamente
la limitación que el fantasma sintético evita, y la razón de haber empezado por
ahí. La tarea de detección con AUC directamente no se puede reproducir sobre
datos reales sin insertar lesiones sintéticas en el espacio-k.
"""

from __future__ import annotations

import argparse

import numpy as np


def leer_kspace(ruta, corte=None):
    """
    Devuelve (kspace, referencia) para un corte.

    kspace: array complejo (n_fase, n_lectura) ya combinado entre bobinas por
            raíz de la suma de cuadrados en el dominio de la imagen.
    """
    import h5py  # importado aquí para que el resto del repo no lo requiera

    with h5py.File(ruta, "r") as f:
        k = np.asarray(f["kspace"])
        ref = np.asarray(f["reconstruction_rss"]) if "reconstruction_rss" in f \
            else (np.asarray(f["reconstruction_esc"]) if "reconstruction_esc" in f
                  else None)

    corte = k.shape[0] // 2 if corte is None else corte
    ks = k[corte]

    if ks.ndim == 3:                      # multibobina: (bobina, fase, lectura)
        imgs = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(ks, axes=(-2, -1))),
                               axes=(-2, -1))
        img = np.sqrt((np.abs(imgs) ** 2).sum(axis=0))   # RSS
        # se vuelve al espacio-k de la imagen combinada para poder truncar
        kcomb = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(img)))
    else:                                  # una bobina
        kcomb = np.fft.fftshift(ks)

    return kcomb, (ref[corte] if ref is not None else None)


def recortar_central(k, n=320):
    """Recorta la región central n x n del espacio-k (fastMRI trae sobremuestreo
    en la dirección de lectura)."""
    cy, cx = k.shape[0] // 2, k.shape[1] // 2
    h = n // 2
    return k[cy - h:cy + h, cx - h:cx + h]


def barrido_truncamiento_real(ruta, corte=None, n_lineas=(320, 240, 160, 120,
                                                          80, 60, 40, 30, 20)):
    """
    Repite la tabla de truncamiento del artículo sobre datos reales.
    La referencia es la reconstrucción a muestreo completo, no una verdad.
    """
    import sys, os
    sys.path.insert(0, os.path.dirname(__file__))
    import mri

    k, _ = leer_kspace(ruta, corte)
    k = recortar_central(k, 320)
    N = k.shape[0]
    ref = np.abs(mri.a_imagen(k))

    filas = []
    for n in n_lineas:
        m = mri.mascara_lineas_centrales(N, n, eje=0)
        rec = np.abs(mri.a_imagen(k * m))
        filas.append({"lineas": n, "frac": n / N,
                      "ssim": mri.ssim(ref, rec),
                      "psnr": mri.psnr(ref, rec),
                      "nrmse": mri.nrmse(ref, rec)})
    return filas


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("ruta", help="archivo .h5 de fastMRI")
    ap.add_argument("--corte", type=int, default=None)
    a = ap.parse_args()

    print(f'{"lineas":>7} {"frac":>6} {"SSIM":>8} {"PSNR":>8} {"NRMSE":>9}')
    for r in barrido_truncamiento_real(a.ruta, a.corte):
        print(f'{r["lineas"]:7d} {r["frac"]:6.2f} {r["ssim"]:8.4f} '
              f'{r["psnr"]:8.2f} {r["nrmse"]:9.5f}')
    print("\nPredicción del artículo: el SSIM debe caer MÁS RÁPIDO que en el "
          "fantasma sintético, porque la anatomía real tiene más energía en "
          "alta frecuencia.")
