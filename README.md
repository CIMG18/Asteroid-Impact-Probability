# Asteroid Impact Probability

Simulación del sistema **Sol–Tierra–Luna** usando gravitación N-cuerpos en C puro. El programa integra las órbitas de los tres cuerpos durante 10 años con el método Velocity Verlet y exporta las trayectorias a un CSV que se visualiza como animación 2D en Python.

---

## Estructura del proyecto

```
Asteroid-Impact-Probability/
├── src/
│   ├── main.c          # Integración N-cuerpos Sol-Tierra-Luna, exporta CSV
│   └── physics.c       # Integrador Velocity Verlet, vectores, energías
├── include/
│   └── physics.h       # Structs (Vector3, Body), constantes astronómicas, prototipos
├── data/
│   └── trajectories.csv  # Posiciones de los 3 cuerpos por cada día simulado
├── visualize_results.py  # Animación 2D con matplotlib.animation.FuncAnimation
└── Makefile              # Compilación: gcc -O3 -Wall
```

---

## Modelo físico

### Gravitación N-cuerpos

El sistema está compuesto por **Sol, Tierra y Luna**. La aceleración gravitacional sobre cualquier cuerpo `i` debida a los demás es:

$$
\mathbf{a}_i = G \sum_{j \neq i} m_j \, \frac{\mathbf{r}_j - \mathbf{r}_i}{\left|\mathbf{r}_j - \mathbf{r}_i\right|^3}
$$

Los tres cuerpos se integran **entre sí** en cada paso de tiempo usando N-cuerpos completo (`step_verlet_bodies`): el Sol perturba a la Luna, la Tierra perturba al Sol, etc. No hay cuerpos fijos.

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

### Condiciones iniciales

Las velocidades orbitales se calculan analíticamente a partir de la condición de órbita circular `v = √(G·M/r)`:

| Cuerpo | Posición inicial | Velocidad inicial |
|--------|-----------------|-------------------|
| Sol    | (0, 0, 0) | (0, 0, 0) |
| Tierra | (1.496×10¹¹ m, 0, 0) | (0, ~29 784 m/s, 0) |
| Luna   | (Tierra + 3.844×10⁸ m, 0, 0) | (0, v_tierra + ~1 022 m/s, 0) |

---

## Cómo compilar y ejecutar

### Requisitos

- `gcc` (incluido en macOS con Xcode Command Line Tools)
- `make`
- Python: `pip install matplotlib pandas numpy`

### Compilar y simular

```bash
make && ./main
```

Genera `data/trajectories.csv` con las posiciones de los 3 cuerpos en cada uno de los 3 650 pasos (10 años × 365 días).

### Visualizar la animación

```bash
python3 visualize_results.py
```

### Limpiar binarios

```bash
make clean
```

---

## Salida del programa

```
Simulando sistema Sol-Tierra-Luna: 10 anios (3650 pasos de 1 dia)...

╔══════════════════════════════════════╗
║   SISTEMA SOL-TIERRA-LUNA            ║
╠══════════════════════════════════════╣
║  Duracion        : 10         anios  ║
║  Pasos           : 3650              ║
║  DT              : 1 dia (86400 s)   ║
╠══════════════════════════════════════╣
║  CONSERVACION DE ENERGIA             ║
║  E inicial : -2.679957e+33 J         ║
║  Error max : 0.000028     %          ║
║  Integrador: EXCELENTE (<0.01%)      ║
╚══════════════════════════════════════╝

Datos guardados en data/trajectories.csv
```

El diagnóstico de conservación de energía verifica que el integrador Velocity Verlet funciona correctamente — un error < 0.01% en 10 años de simulación confirma que las órbitas son físicamente fiables.

---

## Visualización

La animación `visualize_results.py` muestra dos paneles simultáneos:

| Panel | Contenido |
|---|---|
| **Izquierdo** — Vista heliocéntrica | Sol en el centro, Tierra trazando su órbita progresivamente |
| **Derecho** — Referencial geocéntrico | Tierra fija en el origen, Luna orbitando a su alrededor |

La estela de la Tierra se acumula desde el día 0, mostrando cuántas vueltas completas ha dado. La Luna usa una ventana deslizante de 60 días para mantener la imagen legible.

---

## Parámetros configurables

**En [`src/main.c`](src/main.c):**

| Parámetro | Valor actual | Descripción |
|---|---|---|
| `N_YEARS` | 10 | Años de simulación |
| `DT` | 86 400 s | Paso de tiempo (1 día) |

**En [`visualize_results.py`](visualize_results.py):**

| Parámetro | Valor actual | Descripción |
|---|---|---|
| `TRAIL_LEN` | 60 | Días de estela visible (Sol y Luna) |
| `INTERVAL` | 20 ms | Tiempo entre frames (~50 fps) |
| `STEP_SKIP` | 2 | Pasos del CSV por frame (velocidad de reproducción) |

---

## Fuentes de datos de referencia

- **JPL Small-Body Database (SBDB)** — elementos orbitales con incertidumbres 1-sigma:
  `https://ssd-api.jpl.nasa.gov/sbdb.api`

- **JPL Sentry System** — objetos near-Earth monitoreados para posible impacto:
  `https://ssd-api.jpl.nasa.gov/sentry.api`

- **JPL Horizons DE421** — efemérides para las condiciones iniciales del sistema Sol-Tierra-Luna.
