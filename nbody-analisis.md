# Gravitación de N cuerpos: análisis de referencia

Fuente analizada: https://claude.ai/artifact/AUxpmcwxPTjZ6QGteTzFzM ("Gravitación de N cuerpos", notas de clase interactivas, artifact público de otra organización). Este documento resume y analiza esa fuente. El código de abajo es una implementación de referencia escrita aparte, no una copia.

Uso previsto: dar contexto a Claude Code para un proyecto de simulación. Léelo antes de proponer cambios al integrador, a las unidades o a las pruebas de conservación.

## 1. Qué contiene la fuente

Siete secciones, cada figura construye la siguiente:

1. Leyes de partida: segunda ley de Newton y gravitación universal.
2. Fuerza entre dos partículas (par acción-reacción).
3. Fuerza resultante sobre la partícula i con N partículas (superposición).
4. Cambio de posición: Δr_i = v_i Δt + ½ (F_i / m_i) Δt².
5. Trayectoria paso a paso (una partícula alrededor de una masa fija m_j, G = 1).
6. Condiciones iniciales en una caja M×M×M con distribución de Maxwell-Boltzmann.
7. Simulación de N cuerpos en 3D, más el escenario Sol, Tierra y Luna.

## 2. Modelo físico

Aceleración de la partícula i, con suavizado ε para evitar divergencias en encuentros cercanos:

```
a_i = -G * sum_{j != i} m_j * (r_i - r_j) / (|r_i - r_j|^2 + eps^2)^(3/2)
```

Condiciones iniciales de velocidad (Maxwell-Boltzmann): cada componente sigue una normal N(0, k_B T / m_i).
Rapidez cuadrática media: v_rms = sqrt(3 k_B T / m).

En la caja se usan unidades de simulación: G = 1, k_B = 1, masas de 1 a 10. La fuente ofrece fronteras reflejante, periódica y abierta.

## 3. Integradores

La fuente ofrece dos. Referencia en Python/NumPy (r, v: arreglos (N, 3); m: (N,)):

```python
import numpy as np

def accel(r, m, G, eps=0.0):
    d = r[:, None, :] - r[None, :, :]              # d[i, j] = r_i - r_j
    r2 = (d ** 2).sum(-1) + eps ** 2
    np.fill_diagonal(r2, np.inf)                   # sin autointeracción
    return -G * (m[None, :, None] * d / r2[..., None] ** 1.5).sum(axis=1)

def step_taylor(r, v, a, m, G, dt, eps=0.0):
    """Esquema 'nuestras ecuaciones' de la fuente: Δr = vΔt + ½aΔt², v += a_vieja·Δt."""
    r = r + v * dt + 0.5 * a * dt ** 2
    a_new = accel(r, m, G, eps)
    v = v + a * dt                                  # usa la aceleración VIEJA
    return r, v, a_new

def step_verlet(r, v, a, m, G, dt, eps=0.0):
    """Verlet de velocidades: promedia la aceleración antes y después del paso."""
    r = r + v * dt + 0.5 * a * dt ** 2
    a_new = accel(r, m, G, eps)
    v = v + 0.5 * (a + a_new) * dt
    return r, v, a_new
```

## 4. Diagnósticos que usa la fuente

- Energía total E = K + U. Con suavizado, U = -G * sum_{i<j} m_i m_j / sqrt(r_ij^2 + eps^2). Se vigila el error relativo ΔE / |E_0|.
- Razón virial Q = 2K / |U|. Q < 1: el sistema colapsa y oscila hacia el equilibrio virial (Q ≈ 1). Q muy grande: el gas llena la caja.
- Momento lineal total |P| (debería conservarse).
- Temperatura medida T = 2K / (f k_B), con f = 3N - 3 grados de libertad si se resta la velocidad del centro de masa.

## 5. Escenario Sol, Tierra y Luna

- Estado inicial: 21 de octubre de 2026, 17:53:27 UTC, efeméride DE421 de JPL.
- Posiciones en el plano de la eclíptica, medidas desde el centro de masa de los tres cuerpos.
- Unidades: UA, días y masas solares. G = k² = 2.959122e-4 UA³/(M☉·día²), con k = 0.01720209895 (constante gaussiana). Sin suavizado (ε = 0).
- Valores iniciales declarados: distancia Tierra-Luna 390 851 km, Tierra-Sol 0.99543 UA, rapidez orbital de la Tierra 29.91 km/s, de la Luna respecto a la Tierra 1.004 km/s.
- Comprobaciones esperadas: la Luna completa una vuelta en unos 27.3 días y la Tierra en unos 365.25. El Sol atrae a la Luna con más del doble de fuerza que la Tierra (mi cálculo da una razón cercana a 2.2).
- Limitación declarada: sin los demás planetas, la simulación se aparta poco a poco de las posiciones reales.

## 6. Hallazgos del análisis

Lo que verifiqué por lectura (no ejecuté la simulación):

- Las ecuaciones de las secciones 3, 4 y 7 son consistentes entre sí.
- La varianza de velocidades k_B T / m y la distribución de rapideces de Maxwell son correctas.
- El criterio virial Q ≈ 1 y la descripción de Verlet son correctos.
- Los datos orbitales del escenario real son plausibles. No comparé los vectores de estado exactos contra JPL Horizons.

Puntos a tener en cuenta:

1. **El esquema "nuestras ecuaciones" no conserva bien la energía.** Actualiza la posición con ½aΔt² pero la velocidad con la aceleración vieja (v += a·Δt). La energía se desvía con el tiempo. La fuente lo reconoce y recomienda bajar Δt o usar Verlet. Consecuencia: la afirmación de la sección 5 de que "con Δt pequeño la elipse se cierra" es, a mi juicio, aproximada. Probablemente la órbita precesa o se abre poco a poco (no lo probé).
2. **Notación de F₁₂.** La fuente la define como la fuerza que ejerce m₁ sobre m₂, con r̂₁₂ de 1 a 2. Es coherente internamente, pero muchos textos usan F₁₂ para la fuerza sobre 1 debida a 2. Conviene fijar una convención en el proyecto.
3. **El suavizado ε solo aplica a la caja.** En el escenario Sol-Tierra-Luna se fija ε = 0 porque los encuentros cercanos no ocurren.
4. **No verificado:** el comportamiento de las fronteras periódica y reflejante, y la interfaz.

## 7. Notas de mi proyecto (completar)

- Objetivo de la simulación:
- Unidades y sistema de referencia que voy a usar:
- Integrador elegido y por qué:
- Pruebas de conservación que quiero (energía, momento, virial):
- Datos externos (por ejemplo, condiciones iniciales de JPL Horizons):
