import numpy as np
from .physics import rk4_step, RADIUS_EARTH, MASS_SUN, MASS_EARTH, MASS_MOON

def run_single_trajectory(initial_pos, initial_vel, bodies_config, dt, max_steps):
    """
    Simula la trayectoria de un solo meteorito.
    Retorna (impacted, final_pos, trajectory)
    """
    pos = np.array(initial_pos, dtype=float)
    vel = np.array(initial_vel, dtype=float)

    # Extraemos posiciones y masas de los cuerpos celestes
    # En una simulación simple, asumimos que los cuerpos se mueven poco
    # o tienen órbitas predefinidas. Para el MVP, los dejaremos estáticos
    # o con movimiento lineal simple.
    pos_bodies = np.array([b['pos'] for b in bodies_config])
    mass_bodies = np.array([b['mass'] for b in bodies_config])

    trajectory = []

    for step in range(max_steps):
        trajectory.append(pos.copy())

        # Actualizar posición y velocidad usando RK4
        pos, vel = rk4_step(pos, vel, pos_bodies, mass_bodies, dt)

        # Verificar colisión con la Tierra (Cuerpo índice 1)
        dist_to_earth = np.linalg.norm(pos - pos_bodies[1])
        if dist_to_earth <= RADIUS_EARTH:
            return True, pos, np.array(trajectory)

    return False, pos, np.array(trajectory)

def run_monte_carlo(n_simulations, base_pos, base_vel, uncertainty_pos, uncertainty_vel, bodies_config, dt, max_steps):
    """
    Ejecuta múltiples simulaciones variando las condiciones iniciales.
    """
    impacts = 0
    all_trajectories = []

    print(f"Iniciando simulación de Montecarlo: {n_simulations} trayectorias...")

    for i in range(n_simulations):
        # Muestreo aleatorio basado en incertidumbre (Distribución Normal)
        pos_sample = np.array(base_pos) + np.random.normal(0, uncertainty_pos, 3)
        vel_sample = np.array(base_vel) + np.random.normal(0, uncertainty_vel, 3)

        impacted, final_pos, traj = run_single_trajectory(pos_sample, vel_sample, bodies_config, dt, max_steps)

        if impacted:
            impacts += 1

        # Guardamos algunas trayectorias para visualizar (no todas si n es muy grande)
        if i % (n_simulations // 10 if n_simulations > 10 else 1) == 0:
            all_trajectories.append(traj)

        if (i + 1) % 100 == 0:
            print(f"Progreso: {i + 1}/{n_simulations} simulaciones completadas.")

    probability = impacts / n_simulations
    return probability, all_trajectories
