"""
figuras.py — Genera las figuras del artículo. Cada función guarda un PNG en
figuras/ y devuelve la ruta.
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import mri
import experimentos as E

DIR = os.path.join(os.path.dirname(__file__), "..", "figuras")
os.makedirs(DIR, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 130, "font.size": 9,
    "axes.titlesize": 9.5, "axes.labelsize": 9, "savefig.bbox": "tight",
})

N = 256
SIGMA_SNR15 = None  # se calcula al vuelo


def _sin_ejes(ax):
    ax.set_xticks([]); ax.set_yticks([])


# ---------------------------------------------------------------------------
def fig1_contraste():
    lab = mri.fantasma_etiquetas(N)
    pd, t1, t2 = mri.mapas_tejido(lab)

    fig = plt.figure(figsize=(10.4, 5.9))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.2], hspace=0.30, wspace=0.42)

    paneles = [
        ("Densidad protónica\n(lo análogo a un CT)", pd, None),
        ("T1w  TR=500 TE=14", mri.senal_spin_echo(pd, t1, t2, 500, 14), None),
        ("DPw  TR=2500 TE=14", mri.senal_spin_echo(pd, t1, t2, 2500, 14), None),
        ("T2w  TR=4000 TE=100", mri.senal_spin_echo(pd, t1, t2, 4000, 100), None),
    ]
    for j, (tit, im, _) in enumerate(paneles):
        ax = fig.add_subplot(gs[0, j])
        ax.imshow(im, cmap="gray")
        c = mri.contraste(im, lab, 6, 5)
        ax.set_title(f"{tit}\nC(lesión/MB) = {c:+.3f}")
        _sin_ejes(ax)
        # flecha a la lesión grande
        cy, cx, r = mri.centros_lesiones(N)["grande"]
        ax.add_patch(plt.Circle((cx, cy), r * 2.1, fill=False, ec="#e04a2f", lw=1.1))

    # Mapa de contraste
    TRs = np.linspace(100, 5000, 160)
    TEs = np.linspace(2, 200, 160)
    C = np.empty((len(TRs), len(TEs)))
    for i, TR in enumerate(TRs):
        for j, TE in enumerate(TEs):
            s = mri.senal_spin_echo(pd, t1, t2, TR, TE)
            a, b = s[lab == 6].mean(), s[lab == 5].mean()
            C[i, j] = (a - b) / (a + b)

    ax = fig.add_subplot(gs[1, :2])
    vm = np.abs(C).max()
    im = ax.pcolormesh(TEs, TRs, C, cmap="RdBu_r", vmin=-vm, vmax=vm, shading="auto")
    ax.contour(TEs, TRs, C, levels=[0.0], colors="k", linewidths=2.0)
    ax.annotate("contraste = 0\nla lesión es invisible", xy=(30, 1900),
                xytext=(72, 1150), fontsize=8.2,
                arrowprops=dict(arrowstyle="->", lw=1.1, color="k"))
    for nom, (TR, TE) in mri.PROTOCOLOS.items():
        ax.plot(TE, TR, "o", ms=6, mfc="w", mec="k", mew=1.2)
        ax.annotate(nom, (TE, TR), textcoords="offset points", xytext=(7, 5),
                    fontsize=8.5, fontweight="bold")
    ax.set_xlabel("TE (ms)"); ax.set_ylabel("TR (ms)")
    ax.set_title("Contraste lesión / materia blanca en función de TR y TE\n"
                 "el mismo tejido, el mismo escáner", pad=6)
    cb = fig.colorbar(im, ax=ax, pad=0.03)
    cb.set_label("contraste de Michelson", fontsize=8.5)

    # Curvas de señal
    ax = fig.add_subplot(gs[1, 2:])
    TEv = np.linspace(0, 200, 200)
    for et, col in [(5, "#2b6cb0"), (4, "#38a169"), (6, "#e04a2f"), (3, "#805ad5")]:
        _, p, a, b = mri.TEJIDOS[et]
        ax.plot(TEv, p * (1 - np.exp(-4000 / a)) * np.exp(-TEv / b),
                color=col, lw=1.8, label=mri.NOMBRES[et])
    ax.axvline(100, color="k", ls=":", lw=1)
    ax.annotate("TE = 100 ms (T2w)", (100, 0.03), textcoords="offset points",
                xytext=(6, 2), fontsize=8)
    ax.set_xlabel("TE (ms)"); ax.set_ylabel("señal relativa")
    ax.set_title("Decaimiento T2 a TR = 4000 ms:\nel contraste se construye esperando", pad=6)
    ax.legend(fontsize=7.8, frameon=False, loc="upper right", ncol=1,
              handlelength=1.4, borderaxespad=0.2)
    ax.set_ylim(0, 1.02)

    ruta = os.path.join(DIR, "fig1_contraste.png")
    fig.savefig(ruta); plt.close(fig)
    return ruta


# ---------------------------------------------------------------------------
def fig2_kspace():
    lab = mri.fantasma_etiquetas(N)
    img = mri.imagen_protocolo(lab, "T2w")
    k = mri.a_kspace(img)

    r_centro = 0.05 * (N / 2)
    m_c = mri.mascara_anillo(N, 0, r_centro)
    m_p = ~m_c
    ic = np.abs(mri.a_imagen(k * m_c))
    ip = np.abs(mri.a_imagen(k * m_p))

    fig, axs = plt.subplots(1, 5, figsize=(15.2, 3.5),
                            gridspec_kw={"wspace": 0.32,
                                         "width_ratios": [1, 1, 1, 1, 1.45]})

    axs[0].imshow(img, cmap="gray"); _sin_ejes(axs[0])
    axs[0].set_title("Imagen T2w\n(verdad conocida)", fontsize=9)

    axs[1].imshow(np.log10(np.abs(k) + 1e-3), cmap="magma"); _sin_ejes(axs[1])
    axs[1].add_patch(plt.Circle((N / 2, N / 2), r_centro * 4, fill=False,
                                ec="#4ade80", lw=1.4))
    axs[1].set_title("Espacio-k\nlog₁₀|S(kx,ky)|", fontsize=9)

    axs[2].imshow(ic, cmap="gray"); _sin_ejes(axs[2])
    axs[2].set_title("Solo el centro\n0.20 % muestras · 84.5 % energía", fontsize=9)

    axs[3].imshow(ip, cmap="gray"); _sin_ejes(axs[3])
    axs[3].set_title("Solo la periferia\n99.8 % muestras · 15.5 % energía", fontsize=9)

    r, acum = mri.energia_acumulada_radial(k)
    ax = axs[4]
    ax.plot(r / (N / 2) * 100, acum * 100, lw=2, color="#2b6cb0")
    ax.axvline(5, color="#e04a2f", ls="--", lw=1.2)
    ax.annotate("r = 5 %  →  84.5 %\n(0.20 % de las muestras)", (5, 84.5),
                xytext=(40, -42), textcoords="offset points", fontsize=8.2,
                arrowprops=dict(arrowstyle="->", lw=1.1, color="#e04a2f"))
    ax.set_xlabel("radio en el espacio-k (% del máx.)", fontsize=8.5)
    
    ax.set_title("Energía acumulada (%)\nel contraste vive en el centro", fontsize=9)
    ax.set_xlim(0, 100); ax.set_ylim(0, 101)
    ax.grid(alpha=0.25)

    ruta = os.path.join(DIR, "fig2_kspace.png")
    fig.savefig(ruta); plt.close(fig)
    return ruta


# ---------------------------------------------------------------------------
def fig3_truncamiento(lineas=(256, 96, 48, 24, 16)):
    lab = mri.fantasma_etiquetas(N)
    img = mri.imagen_protocolo(lab, "T2w")
    sigma = mri.sigma_para_snr(15, img, lab)
    rng = np.random.default_rng(3)
    ref = np.abs(mri.a_imagen(mri.a_kspace(img)))

    cy, cx, rr = mri.centros_lesiones(N)["grande"]
    cy2, cx2, rr2 = mri.centros_lesiones(N)["chica"]

    fig, axs = plt.subplots(3, len(lineas), figsize=(11.0, 7.0),
                            gridspec_kw={"hspace": 0.07, "wspace": 0.05,
                                         "height_ratios": [1, 1, 1]})
    for j, n in enumerate(lineas):
        rec = E.reconstruir(img, n, sigma, rng, N)
        limpia = E.reconstruir(img, n, 0.0, rng, N)
        axs[0, j].imshow(rec, cmap="gray", vmin=0, vmax=ref.max())
        _sin_ejes(axs[0, j])
        t = mri.tiempo_adquisicion_s(4000, n)
        axs[0, j].set_title(f"{n} líneas · {t/60:.1f} min\n"
                            f"Δy {240/n:.2f} mm · SSIM {mri.ssim(ref, limpia):.2f}",
                            fontsize=8.8)
        for i, (a, b, c) in enumerate([(cy, cx, 30), (cy2, cx2, 22)]):
            sub = rec[int(a - c):int(a + c), int(b - c):int(b + c)]
            axs[i + 1, j].imshow(sub, cmap="gray", vmin=0, vmax=ref.max())
            _sin_ejes(axs[i + 1, j])
    axs[0, 0].set_ylabel("corte completo", fontsize=9)
    axs[1, 0].set_ylabel("lesión 12 mm", fontsize=9)
    axs[2, 0].set_ylabel("lesión 3 mm", fontsize=9)
    fig.suptitle("Truncar el espacio-k: menos líneas, menos tiempo, menos resolución.\n"
                 "SNR por píxel = 15 con muestreo completo; SSIM medido sobre la "
                 "reconstrucción sin ruido", y=0.98, fontsize=10)
    ruta = os.path.join(DIR, "fig3_truncamiento.png")
    fig.savefig(ruta); plt.close(fig)
    return ruta


# ---------------------------------------------------------------------------
def fig4_ssim_vs_auc(datos):
    trunc = {d["lineas"]: d for d in E.barrido_truncamiento()}
    lin = [d["lineas"] for d in datos["grande"]]

    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    x = np.arange(len(lin))

    ax.plot(x, [trunc[n]["ssim"] for n in lin], "s--", color="#718096",
            lw=1.8, ms=5, label="SSIM contra la imagen completa")

    for cual, col, mk in [("grande", "#2b6cb0", "o"), ("chica", "#e04a2f", "^")]:
        a = np.array([d["auc"] for d in datos[cual]])
        lo = np.array([d["ic95"][0] for d in datos[cual]])
        hi = np.array([d["ic95"][1] for d in datos[cual]])
        ax.errorbar(x, a, yerr=[a - lo, hi - a], fmt=mk + "-", color=col, lw=2,
                    ms=5.5, capsize=3,
                    label=f"AUC detección lesión {'12 mm' if cual=='grande' else '3 mm'}")

    ax.axhline(0.5, color="k", ls=":", lw=1)
    ax.annotate("azar (AUC = 0.5)", (0.05, 0.512), fontsize=8, color="#444")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{n}\n{4000*n/1000/60:.1f} min" for n in lin], fontsize=8)
    ax.set_xlabel("líneas de codificación de fase adquiridas  /  tiempo de adquisición")
    ax.set_ylabel("SSIM   |   AUC")
    ax.set_ylim(0.43, 1.05)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8.4, frameon=True, framealpha=0.95, edgecolor="none",
              loc="lower right", borderaxespad=0.6)

    # banda que marca la zona de "misma detectabilidad" para la lesión grande
    i48 = lin.index(48)
    ax.axvspan(-0.4, i48 + 0.4, color="#2b6cb0", alpha=0.06)
    ax.annotate("AUC de la lesión de 12 mm sin cambio significativo\n"
                "SSIM cae de 1.00 a 0.73",
                xy=(i48 / 2, 0.545), fontsize=8.6, ha="center", color="#2b6cb0")
    ax.set_xlim(-0.4, len(lin) - 0.6)
    ax.set_title("El SSIM se desploma y la detectabilidad no: dependen del tamaño de la lesión")

    ruta = os.path.join(DIR, "fig4_ssim_vs_auc.png")
    fig.savefig(ruta); plt.close(fig)
    return ruta


# ---------------------------------------------------------------------------
def fig5_gibbs():
    g = E.gibbs_en_fantasma(n_lineas=64)
    fig, axs = plt.subplots(1, 2, figsize=(10.2, 3.8), gridspec_kw={"wspace": 0.26})

    y = np.arange(len(g["ref"]))
    axs[0].plot(y, g["ref"], color="#718096", lw=1.4, label="256 líneas")
    axs[0].plot(y, g["trunc"], color="#e04a2f", lw=1.6, label="64 líneas")
    axs[0].set_xlim(40, 110)
    axs[0].set_ylim(0.16, 0.40)
    # oscilación dentro de una meseta que la verdad dice constante
    a, b = 40, 51
    nivel = g["ref"][45]
    seg = g["trunc"][a:b]
    lo_p = (seg.min() - nivel) / nivel * 100
    hi_p = (seg.max() - nivel) / nivel * 100
    axs[0].axhspan(seg.min(), seg.max(), xmin=0.0, xmax=(b - 40) / 70,
                   color="#e04a2f", alpha=0.10)
    axs[0].annotate(f"tejido homogéneo:\nla verdad es plana,\nla imagen oscila "
                    f"{lo_p:+.0f} % / {hi_p:+.0f} %",
                    xy=(float(np.argmax(seg)) + a, seg.max()),
                    xytext=(20, 8), textcoords="offset points", fontsize=8.2,
                    arrowprops=dict(arrowstyle="->", lw=1))
    axs[0].set_xlabel("posición en y (píxeles)")
    axs[0].set_ylabel("señal")
    axs[0].set_title("Perfil por la interfaz LCR / parénquima\nel timbre de Gibbs son bandas reales en la imagen", fontsize=9.5)
    axs[0].legend(fontsize=8, frameon=False)
    axs[0].grid(alpha=0.25)

    # convergencia del sobreimpulso
    def over(Nn, n):
        x = np.zeros(Nn); x[Nn // 2:] = 1.0
        X = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(x)))
        m = np.zeros(Nn, bool); c = Nn // 2
        m[c - n // 2:c - n // 2 + n] = True
        yv = np.real(np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(X * m))))
        return (yv[Nn // 2:].max() - 1) * 100

    fracs = np.array([1 / 4, 1 / 8, 1 / 16, 1 / 32, 1 / 64, 1 / 128, 1 / 256, 1 / 512])
    for Nn, col in [(256, "#e04a2f"), (2048, "#38a169"), (32768, "#2b6cb0")]:
        ok = [(f, over(Nn, int(Nn * f))) for f in fracs if int(Nn * f) >= 4]
        axs[1].semilogx([1 / f for f, _ in ok], [v for _, v in ok], "o-", ms=4,
                        color=col, lw=1.6, label=f"malla de {Nn} muestras")
    axs[1].axhline(8.949, color="k", ls="--", lw=1.3)
    axs[1].annotate("8.949 %  —  Gibbs, límite continuo", (10, 8.949),
                    textcoords="offset points", xytext=(6, -22), fontsize=8.4)
    axs[1].set_xlabel("factor de truncamiento  N / n")
    axs[1].set_ylabel("sobreimpulso (%)")
    axs[1].set_title("Verificación numérica del sobreimpulso\ncon 256 muestras la malla se queda corta", fontsize=9.5)
    axs[1].legend(fontsize=8, frameon=False, loc="lower right")
    axs[1].grid(alpha=0.25, which="both")
    axs[1].set_ylim(6, 11)

    ruta = os.path.join(DIR, "fig5_gibbs.png")
    fig.savefig(ruta); plt.close(fig)
    return ruta


if __name__ == "__main__":
    import json
    datos = json.load(open(os.path.join(os.path.dirname(__file__), "..", "_auc.json")))
    for f in (fig1_contraste, fig2_kspace, fig3_truncamiento, fig5_gibbs):
        print(f())
    print(fig4_ssim_vs_auc(datos))
