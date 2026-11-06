import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_trajectories(csv_path):
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error leyendo el archivo: {e}")
        return

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Cuerpos celestes (aproximados)
    bodies = {
        'Sol': [0, 0, 0],
        'Tierra': [1.496e11, 0, 0],
        'Luna': [1.496e11 + 3.844e8, 0, 0]
    }
    colors = {'Sol': 'yellow', 'Tierra': 'blue', 'Luna': 'gray'}

    for name, pos in bodies.items():
        ax.scatter(pos[0], pos[1], pos[2], color=colors[name], s=100, label=name)

    # Dibujar trayectorias
    for sim_id in df['sim_id'].unique():
        traj = df[df['sim_id'] == sim_id]
        ax.plot(traj['x'], traj['y'], traj['z'], alpha=0.5, linewidth=0.5)

    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    ax.set_title('Trayectorias de Meteoritos (Simuladas en C)')
    ax.legend()
    plt.show()

if __name__ == "__main__":
    plot_trajectories('data/trajectories.csv')
