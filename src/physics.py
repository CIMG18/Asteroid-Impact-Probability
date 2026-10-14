import numpy as np

# --- Constantes Astronómicas (S.I.: metros, kilogramos, segundos) ---
G = 6.67430e-11  # Constante de gravitación universal

# Masas
MASS_SUN = 1.989e30
MASS_EARTH = 5.972e24
MASS_MOON = 7.348e22

# Radios (para detección de colisiones)
RADIUS_EARTH = 6.371e6  # 6,371 km

def compute_acceleration(pos_meteorite, pos_bodies, mass_bodies):
    """
    Calcula la aceleración neta sobre el meteorito debido a múltiples cuerpos celestes.

    pos_meteorite: array [x, y, z]
    pos_bodies: array de arrays [[x1, y1, z1], [x2, y2, z2], ...]
    mass_bodies: array de masas [m1, m2, ...]
    """
    accel_total = np.array([0.0, 0.0, 0.0])

    for pos_body, mass_body in zip(pos_bodies, mass_bodies):
        # Vector desde el meteorito al cuerpo
        r_vector = pos_body - pos_meteorite
        distance = np.linalg.norm(r_vector)

        # Ley de Gravitación Universal: a = G * M / r^2
        # Direccion: r_vector / distance
        accel_mag = (G * mass_body) / (distance**2)
        accel_total += accel_mag * (r_vector / distance)

    return accel_total

def rk4_step(pos, vel, pos_bodies, mass_bodies, dt):
    """
    Implementa un paso del integrador Runge-Kutta de 4to orden (RK4)
    para actualizar la posición y velocidad del meteorito.
    """

    # k1
    v1 = vel
    a1 = compute_acceleration(pos, pos_bodies, mass_bodies)

    # k2
    v2 = vel + a1 * (dt / 2)
    a2 = compute_acceleration(pos + v1 * (dt / 2), pos_bodies, mass_bodies)

    # k3
    v3 = vel + a2 * (dt / 2)
    a3 = compute_acceleration(pos + v2 * (dt / 2), pos_bodies, mass_bodies)

    # k4
    v4 = vel + a3 * dt
    a4 = compute_acceleration(pos + v3 * dt, pos_bodies, mass_bodies)

    # Promedio ponderado
    new_pos = pos + (dt / 6.0) * (v1 + 2*v2 + 2*v3 + v4)
    new_vel = vel + (dt / 6.0) * (a1 + 2*a2 + 2*a3 + a4)

    return new_pos, new_vel
