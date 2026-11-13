# Asteroid Impact Probability

Simulación Monte Carlo en C para estimar la **probabilidad de impacto de un asteroide con la Tierra**. Las órbitas de los asteroides nunca se conocen con precisión perfecta — las observaciones telescópicas tienen error de medición. En lugar de calcular una única trayectoria determinista, el proyecto muestrea miles de trayectorias plausibles consistentes con esa incertidumbre y cuenta cuántas de ellas intersectan la Tierra.

---

## Tabla de contenidos

- [Estructura del proyecto](#estructura-del-proyecto)
- [Modelo físico](#modelo-físico)
- [Método Monte Carlo](#método-monte-carlo)
- [Cómo compilar y ejecutar](#cómo-compilar-y-ejecutar)
- [Salida del programa](#salida-del-programa)
- [Visualizador interactivo](#visualizador-interactivo)
- [Parámetros configurables](#parámetros-configurables)
- [Fuentes de datos de referencia](#fuentes-de-datos-de-referencia)

---

## Estructura del proyecto

```
Asteroid-Impact-Probability/
├── src/
│   ├── main.c          # Motor Monte Carlo: muestreo, loop de simulación, salida
│   └── physics.c       # Integrador Velocity Verlet, vectores, energías
├── include/
│   └── physics.h       # Structs (Vector3, Body), constantes astronómicas, prototipos
├── data/
│   └── trajectories.csv  # Trayectorias generadas (1 de cada 10 simulaciones)
├── visualizer.html       # Visualizador 3D interactivo (N-cuerpos en JavaScript)
├── visualize_results.py  # Visualizador alternativo 3D con matplotlib (lee el CSV)
├── Makefile              # Compilación: gcc -O3 -Wall
└── nbody-analisis.md     # Notas de referencia sobre el modelo físico
```

---

## Modelo físico

### Gravitación N-cuerpos

El sistema está compuesto por **Sol, Tierra, Luna y el asteroide**. La aceleración gravitacional sobre cualquier cuerpo `i` debida a los demás es:

$$
\mathbf{a}_i = G \sum_{j \neq i} m_j \, \frac{\mathbf{r}_j - \mathbf{r}_i}{\left|\mathbf{r}_j - \mathbf{r}_i\right|^3}
$$

Sol, Tierra y Luna se integran **entre sí** en cada paso de tiempo usando N-cuerpos completo (`step_verlet_bodies`). El asteroide se integra bajo la gravedad de los tres cuerpos actualizados (`step_verlet`).

### Integrador Velocity Verlet

Cada paso de tiempo `Δt = 1 día (86 400 s)` aplica el esquema simétrico de Verlet:


$$
\begin{aligned}
\mathbf{r}(t + \Delta t) &= \mathbf{r}(t) + \mathbf{v}(t)\,\Delta t + \tfrac{1}{2}\,\mathbf{a}(t)\,\Delta t^2 && \text{(nueva posición)} \\
\mathbf{a}(t + \Delta t) &= \mathbf{a}\big(\mathbf{r}(t + \Delta t)\big) && \text{(nueva aceleración)} \\
\mathbf{v}(t + \Delta t) &= \mathbf{v}(t) + \tfrac{1}{2}\,\big[\mathbf{a}(t) + \mathbf{a}(t + \Delta t)\big]\,\Delta t && \text{(nueva velocidad)}
\end{aligned}
$$


El promedio de aceleraciones en el paso de velocidad es lo que distingue a Velocity Verlet de métodos más simples: **conserva la energía** a largo plazo sin disipación artificial.

### Condiciones iniciales de los cuerpos principales

Las velocidades orbitales se calculan analíticamente a partir de la condición de órbita circular `v = √(G·M/r)`:

| Cuerpo | Posición inicial | Velocidad inicial |
|--------|-----------------|-------------------|
| Sol    | (0, 0, 0) | (0, 0, 0) |
| Tierra | (1.496×10¹¹ m, 0, 0) | (0, ~29 784 m/s, 0) |
| Luna   | (Tierra + 3.844×10⁸ m, 0, 0) | (0, v_tierra + ~1 022 m/s, 0) |

Cada simulación Monte Carlo **resetea los cuerpos a estas condiciones iniciales**, garantizando que el sistema físico de fondo es idéntico para todas las muestras.

## Método Monte Carlo

La incertidumbre de medición del telescopio se modela como una distribución gaussiana sobre la posición y velocidad base del asteroide. Para cada simulación `i` de `N = 1000`:

$$
\begin{aligned}
\mathbf{r}_i &= \mathbf{r}_{\text{base}} + \mathcal{N}\!\left(\mathbf{0},\, \sigma_{\text{pos}}^2 \mathbf{I}\right), & \sigma_{\text{pos}} &= 1\,000\ \text{km} \\
\mathbf{v}_i &= \mathbf{v}_{\text{base}} + \mathcal{N}\!\left(\mathbf{0},\, \sigma_{\text{vel}}^2 \mathbf{I}\right), & \sigma_{\text{vel}} &= 100\ \text{m/s}
\end{aligned}
$$

Los números aleatorios se generan con el método **Box-Muller**:

$$
z = \sqrt{-2 \ln u_1}\,\cos\!\left(2\pi u_2\right), \qquad u_1, u_2 \sim \mathcal{U}(0, 1)
$$

Cada muestra se integra durante **365 pasos (1 año)**. Si la distancia entre el asteroide y la Tierra cae por debajo de `RADIUS_EARTH = 6.371×10⁶ m`, se registra un impacto.

La probabilidad final es:

```
P(impacto) = N_impactos / N_simulaciones
```

---

## Salida del programa

```
Semilla fija: 42 (resultados reproducibles)
Iniciando Montecarlo en C: 1000 simulaciones, 365 dias cada una...
Progreso: 100/1000 simulaciones completadas.
...

╔══════════════════════════════════════╗
║        RESULTADOS FINALES            ║
╠══════════════════════════════════════╣
║  Semilla usada       : 42            ║
║  Simulaciones        : 1000          ║
║  Duracion por sim.   : 365      dias ║
╠══════════════════════════════════════╣
║  Impactos detectados : 0             ║
║  Probabilidad        : 0.0000 %      ║
╠══════════════════════════════════════╣
║  CONSERVACION DE ENERGIA (sim. 0)    ║
║  E inicial : -4.984186e+18 J         ║
║  E final   : -4.986417e+18 J         ║
║  Error max : 0.046787     %          ║
║  Integrador: ACEPTABLE (<1.00%)      ║
╚══════════════════════════════════════╝
```

Además genera `data/trajectories.csv` con el estado completo (posición, velocidad, energías) de 1 de cada 10 simulaciones, paso a paso.

---

## Visualizador interactivo

El archivo [`visualizer.html`](visualizer.html) es un visualizador autocontenido que puede usar los calculos contenidos en el CSV **O** implementar el mismo integrador Velocity Verlet N-cuerpos directamente en JavaScript y correr la simulación en el navegador.

```bash
open visualizer.html
```

## Parámetros configurables

Todos los parámetros de la simulación C se definen en [`src/main.c`](src/main.c):

| Parámetro | Valor actual | Descripción |
|---|---|---|
| `N_SIMULATIONS` | 1 000 | Número de muestras Monte Carlo |
| `DT` | 86 400 s | Paso de tiempo (1 día) |
| `MAX_STEPS` | 365 | Duración de cada simulación (1 año) |
| `UNCERTAINTY_POS` | 1×10⁶ m | Incertidumbre de posición (1 000 km, σ gaussiano) |
| `UNCERTAINTY_VEL` | 100 m/s | Incertidumbre de velocidad (σ gaussiano) |

---
