"""
visualize_results.py
Animación 2D del sistema Sol-Tierra-Luna usando datos generados por el programa C.

Uso:
    python3 visualize_results.py

Requiere:
    pip install matplotlib pandas numpy
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

CSV_PATH   = 'data/trajectories.csv'
AU         = 1.496e11          # metros por UA (para escalar ejes)
TRAIL_LEN  = 60                # días de estela visible por cuerpo
INTERVAL   = 20                # ms entre frames (≈50 fps)
STEP_SKIP  = 2                 # pasos del CSV por frame (velocidad de reproducción)

# ─────────────────────────────────────────────────────────────────────────────
#  CARGA DE DATOS
# ─────────────────────────────────────────────────────────────────────────────
df = pd.read_csv(CSV_PATH)

# Convertir posiciones de metros a UA para que la escala sea legible
sx = df['sx'].values / AU
sy = df['sy'].values / AU
ex = df['ex'].values / AU
ey = df['ey'].values / AU
mx = df['mx'].values / AU
my = df['my'].values / AU

n_frames = len(df)
steps    = df['step'].values

# ─────────────────────────────────────────────────────────────────────────────
#  FIGURA Y EJES
# ─────────────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 7),
                         facecolor='#04080f',
                         gridspec_kw={'width_ratios': [2, 1]})
fig.suptitle('Sistema Sol–Tierra–Luna  |  Integrador Velocity Verlet (N-cuerpos)',
             color='#c9d1d9', fontsize=13, y=0.97)

# Panel izquierdo: vista heliocéntrica (órbita terrestre completa)
ax_main = axes[0]
ax_main.set_facecolor('#04080f')
ax_main.set_aspect('equal')
ax_main.set_xlim(-1.6, 1.6)
ax_main.set_ylim(-1.6, 1.6)
ax_main.set_xlabel('X  [UA]', color='#8b949e')
ax_main.set_ylabel('Y  [UA]', color='#8b949e')
ax_main.set_title('Vista heliocéntrica', color='#8b949e', fontsize=10)
ax_main.tick_params(colors='#8b949e')
for spine in ax_main.spines.values():
    spine.set_edgecolor('#21262d')

# Panel derecho: sistema Tierra-Luna (posiciones relativas a la Tierra)
ax_zoom = axes[1]
ax_zoom.set_facecolor('#04080f')
ax_zoom.set_aspect('equal')
ax_zoom.set_title('Sistema Tierra–Luna  (referencial geocéntrico)', color='#8b949e', fontsize=10)
ax_zoom.set_xlabel('ΔX  [×10⁻³ UA]', color='#8b949e')
ax_zoom.set_ylabel('ΔY  [×10⁻³ UA]', color='#8b949e')
ax_zoom.tick_params(colors='#8b949e')
for spine in ax_zoom.spines.values():
    spine.set_edgecolor('#21262d')

# Escala del panel derecho en unidades de distancia Tierra-Luna
R_MOON_AU = 3.844e8 / AU   # ~0.00257 UA
PAD_Z     = R_MOON_AU * 2.0
ax_zoom.set_xlim(-PAD_Z * 1e3, PAD_Z * 1e3)   # en ×10⁻³ UA
ax_zoom.set_ylim(-PAD_Z * 1e3, PAD_Z * 1e3)

plt.tight_layout(rect=[0, 0, 1, 0.95])

# ─────────────────────────────────────────────────────────────────────────────
#  ELEMENTOS ANIMADOS
# ─────────────────────────────────────────────────────────────────────────────
# Estelas (líneas)
trail_sun_m,   = ax_main.plot([], [], color='#ffe066', lw=0.8, alpha=0.5)
trail_earth_m, = ax_main.plot([], [], color='#4d9de0', lw=0.8, alpha=0.6)
trail_moon_m,  = ax_main.plot([], [], color='#aaaaaa', lw=0.6, alpha=0.5)

# Panel derecho: Tierra fija en origen, Luna en coordenadas relativas (×10⁻³ UA)
trail_moon_z,  = ax_zoom.plot([], [], color='#aaaaaa', lw=0.6, alpha=0.7)

# Puntos actuales
pt_sun,     = ax_main.plot([], [], 'o', color='#ffe066', ms=14, zorder=5)
pt_earth_m, = ax_main.plot([], [], 'o', color='#4d9de0', ms=6,  zorder=5)
pt_moon_m,  = ax_main.plot([], [], 'o', color='#aaaaaa', ms=3,  zorder=5)

# Panel derecho: Tierra como punto fijo en (0,0), Luna orbitando
pt_earth_z, = ax_zoom.plot([0], [0], 'o', color='#4d9de0', ms=10, zorder=5)
pt_moon_z,  = ax_zoom.plot([], [],   'o', color='#aaaaaa', ms=5,  zorder=5)

# Contador de días
txt_day = ax_main.text(0.02, 0.97, '', transform=ax_main.transAxes,
                       color='#c9d1d9', fontsize=10, va='top')

# ─────────────────────────────────────────────────────────────────────────────
#  FUNCIONES DE ANIMACIÓN
# ─────────────────────────────────────────────────────────────────────────────
def init():
    for artist in (trail_sun_m, trail_earth_m, trail_moon_m, trail_moon_z,
                   pt_sun, pt_earth_m, pt_moon_m, pt_moon_z):
        artist.set_data([], [])
    pt_earth_z.set_data([0], [0])   # Tierra siempre fija en el origen del panel derecho
    txt_day.set_text('')
    return (trail_sun_m, trail_earth_m, trail_moon_m, trail_moon_z,
            pt_sun, pt_earth_m, pt_moon_m, pt_earth_z, pt_moon_z, txt_day)


def update(frame):
    idx = min(frame * STEP_SKIP, n_frames - 1)

    # Panel izquierdo:
    #   Tierra: estela acumulada desde el inicio
    #   Sol y Luna: ventana deslizante de TRAIL_LEN días
    t0_sl = max(0, idx - TRAIL_LEN)

    trail_sun_m.set_data(sx[t0_sl:idx+1],  sy[t0_sl:idx+1])
    trail_earth_m.set_data(ex[0:idx+1],    ey[0:idx+1])      # acumulada
    trail_moon_m.set_data(mx[t0_sl:idx+1], my[t0_sl:idx+1])

    pt_sun.set_data([sx[idx]],     [sy[idx]])
    pt_earth_m.set_data([ex[idx]], [ey[idx]])
    pt_moon_m.set_data([mx[idx]],  [my[idx]])

    # Panel derecho: posición de la Luna RELATIVA a la Tierra (×10⁻³ UA)
    # La Tierra queda fija en (0, 0) — referencial geocéntrico
    rel_mx = (mx[0:idx+1] - ex[0:idx+1]) * 1e3
    rel_my = (my[0:idx+1] - ey[0:idx+1]) * 1e3
    trail_moon_z.set_data(rel_mx, rel_my)
    pt_moon_z.set_data([(mx[idx] - ex[idx]) * 1e3],
                       [(my[idx] - ey[idx]) * 1e3])

    # Día actual
    year = steps[idx] // 365
    day  = steps[idx] %  365
    txt_day.set_text(f'Año {year}  Día {day}  (t = {steps[idx]} días)')

    return (trail_sun_m, trail_earth_m, trail_moon_m, trail_moon_z,
            pt_sun, pt_earth_m, pt_moon_m, pt_earth_z, pt_moon_z, txt_day)


total_frames = (n_frames + STEP_SKIP - 1) // STEP_SKIP

anim = animation.FuncAnimation(
    fig,
    update,
    frames=total_frames,
    init_func=init,
    interval=INTERVAL,
    blit=True
)

plt.show()
