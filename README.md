# Asteroid Impact Probability

Simulación Monte Carlo en C puro para estimar la **probabilidad de impacto de un asteroide con la Tierra**. Las órbitas de los asteroides nunca se conocen con precisión perfecta — las observaciones telescópicas tienen error de medición. En lugar de calcular una única trayectoria determinista, el proyecto muestrea miles de trayectorias plausibles consistentes con esa incertidumbre y cuenta cuántas de ellas intersectan la Tierra.

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

```
aᵢ = G · Σⱼ≠ᵢ  mⱼ · (rⱼ − rᵢ) / |rⱼ − rᵢ|³
```

Sol, Tierra y Luna se integran **entre sí** en cada paso de tiempo usando N-cuerpos completo (`step_verlet_bodies`). El asteroide se integra bajo la gravedad de los tres cuerpos actualizados (`step_verlet`).

### Integrador Velocity Verlet

Cada paso de tiempo `Δt = 1 día (86 400 s)` aplica el esquema simétrico de Verlet:

```
r(t+Δt) = r(t) + v(t)·Δt + ½·a(t)·Δt²     — nueva posición
a(t+Δt) = calcular con la nueva r            — nueva aceleración
v(t+Δt) = v(t) + ½·[a(t) + a(t+Δt)]·Δt    — nueva velocidad
```

El promedio de aceleraciones en el paso de velocidad es lo que distingue a Velocity Verlet de métodos más simples: **conserva la energía** a largo plazo sin disipación artificial.

### Condiciones iniciales de los cuerpos principales

Las velocidades orbitales se calculan analíticamente a partir de la condición de órbita circular `v = √(G·M/r)`:

| Cuerpo | Posición inicial | Velocidad inicial |
|--------|-----------------|-------------------|
| Sol    | (0, 0, 0) | (0, 0, 0) |
| Tierra | (1.496×10¹¹ m, 0, 0) | (0, ~29 784 m/s, 0) |
| Luna   | (Tierra + 3.844×10⁸ m, 0, 0) | (0, v_tierra + ~1 022 m/s, 0) |

Cada simulación Monte Carlo **resetea los cuerpos a estas condiciones iniciales**, garantizando que el sistema físico de fondo es idéntico para todas las muestras.

### Diagnóstico de conservación de energía

El programa evalúa el error relativo de energía total `|ΔE/E₀|` a lo largo de la simulación de referencia (muestra 0) y reporta:

- `< 0.01%` → **EXCELENTE** — el integrador funciona con alta precisión.
- `< 1.00%` → **ACEPTABLE** — deriva pequeña, esperada con Δt = 1 día.
- `> 1.00%` → **REVISAR** — el paso de tiempo es demasiado grande.

---

## Método Monte Carlo

La incertidumbre de medición del telescopio se modela como una distribución gaussiana sobre la posición y velocidad base del asteroide. Para cada simulación `i` de `N = 1000`:

```
posᵢ = base_pos + N(0, σ_pos)   [σ_pos = 1,000 km]
velᵢ = base_vel + N(0, σ_vel)   [σ_vel = 100 m/s]
```

Los números aleatorios se generan con el método **Box-Muller**:

```
z = √(−2 ln u₁) · cos(2π u₂)    donde u₁, u₂ ~ Uniforme(0,1)
```

Cada muestra se integra durante **365 pasos (1 año)**. Si la distancia entre el asteroide y la Tierra cae por debajo de `RADIUS_EARTH = 6.371×10⁶ m`, se registra un impacto.

La probabilidad final es:

```
P(impacto) = N_impactos / N_simulaciones
```


### Ejecutar

```bash
# Semilla aleatoria — resultados distintos en cada ejecución
./main

# Semilla fija — resultados 100% reproducibles (útil para comparar y presentar)
./main 42
```

### Limpiar binarios

```bash
make clean
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
║  Semilla usada       : 42             ║
║  Simulaciones        : 1000           ║
║  Duracion por sim.   : 365        dias ║
╠══════════════════════════════════════╣
║  Impactos detectados : 0              ║
║  Probabilidad        : 0.0000       % ║
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

El archivo [`visualizer.html`](visualizer.html) es un visualizador autocontenido que **no requiere el CSV** — implementa el mismo integrador Velocity Verlet N-cuerpos directamente en JavaScript y corre la simulación en el navegador.

```bash
open visualizer.html
```

### Características

| Funcionalidad | Descripción |
|---|---|
| N-cuerpos en JS | Sol, Tierra y Luna se mueven con física real |
| Enjambre Monte Carlo | Hasta 400 asteroides con incertidumbre gaussiana ajustable |
| Diagnósticos en vivo | Energía cinética K, potencial U, total E; Virial Q = 2K/\|U\| |
| Controles | Play/Pausa, Reset, velocidad (0.5×–20×), sliders de incertidumbre |
| Cámara | Scroll para zoom, drag para pan |
| Colores | Blanco = activo, Rojo = impacto, Verde = se alejó |



---

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
