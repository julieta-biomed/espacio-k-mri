# El espacio-k: dónde vive el contraste y dónde vive el borde

Medición de cuánto espacio-k hace falta realmente en una resonancia magnética, sobre un fantasma cerebral sintético con **verdad conocida**.

**Resultado principal:** recortar de 256 a 48 líneas de codificación de fase (17.1 → 3.2 min de adquisición) hace caer el SSIM de 1.00 a 0.73, **no cambia** el AUC de detección de una lesión de 12 mm (0.933 → 0.924, intervalos superpuestos) y hunde el de una de 3 mm (0.907 → 0.831).

Artículo: [`el-espacio-k-donde-vive-el-contraste.md`](el-espacio-k-donde-vive-el-contraste.md)

---

## Resultados

### El contraste no es densidad

Mismo tejido, mismo escáner, distinto TR y TE. Contraste de Michelson.

| Protocolo | TR (ms) | TE (ms) | lesión/blanca | gris/blanca | LCR/blanca |
|---|---|---|---|---|---|
| solo densidad protónica | — | — | +0.091 | — | — |
| T1w | 500 | 14 | −0.095 | −0.075 | −0.487 |
| DPw | 2500 | 14 | +0.076 | +0.048 | −0.144 |
| T2w | 4000 | 100 | +0.350 | +0.178 | +0.474 |

Barrido completo TR ∈ [100, 5000], TE ∈ [2, 200]: el contraste de la lesión va de **−0.190 a +0.583**, y vale **cero** en TR = 1454 ms, TE = 10 ms. Con ese protocolo la lesión es invisible en una imagen que se ve perfecta.

### Energía radial del espacio-k

| Radio (% del máx.) | % de muestras | % de energía |
|---|---|---|
| 1 % | 0.01 % | 70.13 % |
| 2 % | 0.03 % | 73.20 % |
| **5 %** | **0.20 %** | **84.45 %** |
| 10 % | 0.78 % | 85.98 % |
| 25 % | 4.90 % | 94.65 % |
| 50 % | 19.61 % | 97.50 % |

### Truncamiento: resolución, tiempo y métricas de imagen

FOV 240 mm, matriz 256, T2w con TR = 4000 ms. Métricas contra la reconstrucción completa **sin ruido**.

| Líneas | Frac. | Δy (mm) | t (s) | t (min) | SSIM | PSNR (dB) | NRMSE |
|---|---|---|---|---|---|---|---|
| 256 | 1.00 | 0.94 | 1024 | 17.1 | 1.0000 | ∞ | 0.00000 |
| 192 | 0.75 | 1.25 | 768 | 12.8 | 0.9406 | 32.47 | 0.02379 |
| 160 | 0.62 | 1.50 | 640 | 10.7 | 0.9385 | 30.63 | 0.02941 |
| 128 | 0.50 | 1.88 | 512 | 8.5 | 0.8987 | 28.83 | 0.03620 |
| 96 | 0.38 | 2.50 | 384 | 6.4 | 0.8632 | 26.79 | 0.04574 |
| 64 | 0.25 | 3.75 | 256 | 4.3 | 0.8112 | 24.74 | 0.05785 |
| **48** | 0.19 | 5.00 | 192 | **3.2** | **0.7341** | 22.91 | 0.07148 |
| 32 | 0.12 | 7.50 | 128 | 2.1 | 0.6513 | 20.60 | 0.09327 |
| 24 | 0.09 | 10.00 | 96 | 1.6 | 0.6303 | 19.69 | 0.10358 |
| 16 | 0.06 | 15.00 | 64 | 1.1 | 0.6242 | 18.92 | 0.11318 |

### Detectabilidad (600 realizaciones por condición, IC95 por bootstrap)

SNR por píxel = 15 en muestreo completo. Dificultad calibrada por contraste, no por ruido: lesión de 12 mm con T2 = 83.14 ms (contraste +0.0065), lesión de 3 mm con T2 = 85.37 ms (contraste +0.0222), fondo de materia blanca T2 = 80 ms.

| Líneas | t (min) | AUC 12 mm | AUC 3 mm | SSIM |
|---|---|---|---|---|
| 256 | 17.1 | 0.933 [0.919–0.946] | 0.907 [0.891–0.922] | 1.000 |
| 192 | 12.8 | 0.927 [0.913–0.940] | 0.907 [0.890–0.923] | 0.941 |
| 128 | 8.5 | 0.927 [0.912–0.941] | 0.888 [0.870–0.906] | 0.899 |
| 96 | 6.4 | 0.930 [0.916–0.944] | 0.881 [0.861–0.900] | 0.863 |
| 64 | 4.3 | 0.918 [0.903–0.932] | 0.847 [0.825–0.868] | 0.811 |
| **48** | **3.2** | **0.924 [0.909–0.938]** | **0.831 [0.808–0.852]** | **0.734** |
| 32 | 2.1 | 0.904 [0.887–0.920] | 0.746 [0.719–0.773] | 0.651 |
| 24 | 1.6 | 0.913 [0.897–0.927] | 0.725 [0.696–0.753] | 0.630 |
| 16 | 1.1 | 0.857 [0.836–0.876] | 0.659 [0.629–0.689] | 0.624 |

La explicación física: a 48 líneas la lesión de 12 mm conserva el **89 %** de su amplitud (la PSF no la alcanza) mientras el ruido baja; la de 3 mm conserva el **41 %** y pierde señal más rápido de lo que gana en ruido.

| Líneas | Señal, 12 mm | Señal, 3 mm |
|---|---|---|
| 256 | 0.2312 (100 %) | 0.2312 (100 %) |
| 128 | 0.2236 (97 %) | 0.1990 (86 %) |
| 48 | 0.2054 (89 %) | 0.0946 (41 %) |
| 24 | 0.1805 (78 %) | 0.0429 (19 %) |
| 16 | 0.1335 (58 %) | 0.0277 (12 %) |

### Timbre de Gibbs: verificación numérica

Valor teórico en el límite continuo: **8.949 %**.

| N / n | malla 256 | malla 32768 |
|---|---|---|
| 2 | 6.83 % | — |
| 4 | 8.44 % | — |
| 8 | 8.88 % | 8.819 % |
| 16 | 9.18 % | 8.916 % |
| 64 | — | **8.947 %** |
| 256 | — | **8.953 %** |
| 1024 | — | 9.014 % |

Una malla de 256 muestras **no basta** para reproducir el valor teórico: el máximo del sobreimpulso continuo cae entre muestras. Con truncamientos moderados sobre malla fina, el valor converge a 8.949 % por ambos lados.

En anatomía la cifra relevante es otra: en un borde aislado en 2D con 64 de 256 líneas el sobreimpulso es **5.04 %**, pero sobre una meseta que la verdad dice constante la reconstrucción oscila **−6 % / +10 %** de su propio valor, por superposición del timbre de bordes vecinos.

### Tiempos de adquisición (un corte, NEX = 1)

| Protocolo | TR (ms) | 256 líneas | 128 líneas | 48 líneas |
|---|---|---|---|---|
| T1w | 500 | 2.1 min | 1.1 min | 0.4 min |
| DPw | 2500 | 10.7 min | 5.3 min | 2.0 min |
| T2w | 4000 | 17.1 min | 8.5 min | 3.2 min |

Un CT helicoidal de cráneo: del orden de 5–10 s de adquisición.

---

## Figuras

| | |
|---|---|
| `figuras/fig1_contraste.png` | Ponderaciones T1/DP/T2, mapa de contraste en el plano TR–TE con el contorno de contraste nulo |
| `figuras/fig2_kspace.png` | Espacio-k, reconstrucciones solo-centro y solo-periferia, energía acumulada radial |
| `figuras/fig3_truncamiento.png` | Reconstrucciones a 256/96/48/24/16 líneas con recortes sobre cada lesión |
| `figuras/fig4_ssim_vs_auc.png` | **La figura principal**: SSIM y AUC sobre las mismas imágenes |
| `figuras/fig5_gibbs.png` | Perfil por la interfaz LCR/parénquima y convergencia del sobreimpulso |

---

## Cómo correrlo

```bash
git clone https://github.com/USUARIO/espacio-k-mri.git
cd espacio-k-mri
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
jupyter notebook notebooks/espacio_k.ipynb
```

El notebook está en el repositorio **ya ejecutado, con todas las salidas guardadas**. Correrlo completo toma unos 5 minutos; la celda cara es el barrido de AUC (600 realizaciones × 9 condiciones × 2 lesiones, con bootstrap de 2000 remuestreos por punto).

Para regenerar solo las figuras:

```bash
cd src && python figuras.py
```

## Estructura

```
notebooks/espacio_k.ipynb    notebook ejecutado con salidas
src/mri.py                   fantasma, señal spin-echo, espacio-k, métricas, AUC
src/experimentos.py          rutinas de medición (devuelven números, no dibujan)
src/figuras.py               generación de las cinco figuras
src/cargar_fastmri.py        puente a datos reales — ESCRITO, NO EJECUTADO
figuras/                     PNG generados
resultados/resultados.json   todas las cifras en crudo
```

## Datos reales

Los resultados publicados salen **todos** del fantasma sintético, porque solo con verdad conocida se puede medir error absoluto y montar una tarea de detección con lesión conocida. `src/cargar_fastmri.py` lee un `.h5` de [fastMRI](https://fastmri.med.nyu.edu/) (NYU Langone), combina bobinas por raíz de la suma de cuadrados y repite la tabla de truncamiento sobre anatomía real. **No lo he ejecutado**: el entorno donde se hicieron las mediciones no tiene acceso de red a ese dominio.

Predicción concreta para quien lo corra: el SSIM debería caer **más rápido** que en el fantasma, porque la anatomía real tiene bastante más energía en alta frecuencia que un dibujo de elipses.

## Limitaciones

Están desarrolladas en la sección "Dónde falla esto" del artículo. Las tres que más pesan:

1. **El observador no es ideal.** Es una diferencia disco-menos-anillo con posición conocida. Un observador con prewhitening no mostraría la planicie: como el truncamiento es lineal, quien conoce la correlación del ruido no gana nada descartando datos. La planicie es una propiedad de observadores no ideales — la clase a la que pertenece un humano.
2. **El fantasma es un dibujo.** Sin textura, sin inhomogeneidad de B0/B1, sin perfiles de bobina. Las cifras absolutas de AUC dependen de eso; la forma de las curvas debería transferir.
3. **Relleno con ceros, no reconstrucción.** Sin compressed sensing, sin regularización, sin redes. Es el piso, no el estado del arte.

## Licencia

MIT. Ver [LICENSE](LICENSE).
