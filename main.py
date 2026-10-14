import numpy as np
import matplotlib.pyplot as plt
from src.simulation import run_monte_carlo
from src.physics import MASS_SUN, MASS_EARTH, MASS_MOON

def plot_results(trajectories, bodies_config):
    """
    Grafica las trayectorias en 3D.
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Dibujar cuerpos celestes
    colors = ['yellow', 'blue', 'gray'] # Sol, Tierra, Luna
    for i, body in enumerate(bodies_config):
        pos = body['pos']
        ax.scatter(pos[0], pos[1], pos[2], color=colors[i], s=100, label=body['name'])

    # Dibujar trayectorias
    for traj in trajectories:
        ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], alpha=0.5, linewidth=0.5)

    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    ax.set_title('Simulaciones de Trayectorias de Meteoritos')
    ax.legend()
    plt.show()

def main():
    # 1. Configuración de Cuerpos Celestes (Posiciones simplificadas en metros)
    # Colocamos al Sol en el origen
    bodies_config = [
        {'name': 'Sol',   'mass': MASS_SUN,   'pos': np.array([0, 0, 0])},
        {'name': 'Tierra', 'mass': MASS_EARTH, 'pos': np.array([1.496e11, 0, 0])}, # 1 AU
        {'name': 'Luna',   'mass': MASS_MOON,   'pos': np.array([1.496e11 + 3.844e8, 0, 0])},
    ]

    # 2. Parámetros del Meteorito (Condiciones Base)
    # Posición inicial: un poco antes de la órbita de la Tierra
    base_pos = np.array([1.4e11, 1e9, 0])
    # Velocidad inicial: aproximadamente velocidad orbital (~30km/s)
    base_vel = np.array([0, 30000, 0])

    # Incertidumbres (Desviación estándar)
    uncertainty_pos = 1e6 # 1,000 km
    uncertainty_vel = 100 # 100 m/s

    # 3. Parámetros de Simulación
    n_simulations = 1000
    dt = 3600 * 24 # Paso de tiempo: 1 día (86,400 s)
    max_steps = 100 # Simular por 100 días

    # Ejecutar Montecarlo
    prob, trajectories = run_monte_carlo(
        n_simulations, base_pos, base_vel,
        uncertainty_pos, uncertainty_vel,
        bodies_config, dt, max_steps
    )

    print(f"\n--- RESULTADOS ---")
    print(f"Probabilidad de impacto estimada: {prob * 100:.4f}%")

    # Visualización
    plot_results(trajectories, bodies_config)

if __name__ == "__main__":
    main()
