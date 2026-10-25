#include "physics.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Configuración de la simulación ---
#define N_SIMULATIONS 1000
#define DT 86400.0        // 1 día en segundos
#define MAX_STEPS 100     // 100 días de simulación
#define UNCERTAINTY_POS 1e6 // 1,000 km
#define UNCERTAINTY_VEL 100 // 100 m/s

// Función para generar números aleatorios en un rango normal (aproximación simple de Box-Muller)
double random_normal(double mu, double sigma) {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return z0 * sigma + mu;
}

int main() {
    srand(time(NULL));

    // 1. Definir Cuerpos Celestes
    Body bodies[3] = {
        {"Sol",   MASS_SUN,   {0, 0, 0}},
        {"Tierra", MASS_EARTH, {1.496e11, 0, 0}},
        {"Luna",   MASS_MOON,   {1.496e11 + 3.844e8, 0, 0}}
    };

    // 2. Condiciones Base del Meteorito
    Vector3 base_pos = {1.4e11, 1e9, 0};
    Vector3 base_vel = {0, 30000, 0};

    int impacts = 0;

    // Abrir archivo para guardar trayectorias de algunos meteoritos (para visualización)
    FILE *fp = fopen("data/trajectories.csv", "w");
    if (fp == NULL) {
        perror("Error abriendo archivo de datos");
        return 1;
    }
    fprintf(fp, "sim_id,step,x,y,z,vx,vy,vz,ek,ep,etot\n");

    printf("Iniciando Montecarlo en C: %d simulaciones...\n", N_SIMULATIONS);

    for (int i = 0; i < N_SIMULATIONS; i++) {
        // Muestreo aleatorio de condiciones iniciales
        Vector3 pos = {
            base_pos.x + random_normal(0, UNCERTAINTY_POS),
            base_pos.y + random_normal(0, UNCERTAINTY_POS),
            base_pos.z + random_normal(0, UNCERTAINTY_POS)
        };
        Vector3 vel = {
            base_vel.x + random_normal(0, UNCERTAINTY_VEL),
            base_vel.y + random_normal(0, UNCERTAINTY_VEL),
            base_vel.z + random_normal(0, UNCERTAINTY_VEL)
        };

        bool impacted = false;
        double mass_meteorite = 1e10; // Masa arbitraria para el meteorito para cálculos de energía

        for (int step = 0; step < MAX_STEPS; step++) {
            // Guardar trayectoria de cada 10 simulaciones para no llenar el disco
            if (i % 10 == 0) {
                double ek = compute_kinetic_energy(vel, mass_meteorite);
                double ep = compute_potential_energy(pos, bodies, 3, mass_meteorite);
                fprintf(fp, "%d,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",
                        i, step, pos.x, pos.y, pos.z, vel.x, vel.y, vel.z, ek, ep, ek + ep);
            }

            step_verlet(&pos, &vel, bodies, 3, DT);

            // Verificar colisión con la Tierra (bodies[1])
            Vector3 diff = vec_sub(pos, bodies[1].pos);
            if (vec_norm(diff) <= RADIUS_EARTH) {
                impacted = true;
                break;
                // El 'break' ha sido eliminado para que el meteorito siga volando
            }
        }


        if (impacted) impacts++;

        if ((i + 1) % 100 == 0) {
            printf("Progreso: %d/%d simulaciones completadas.\n", i + 1, N_SIMULATIONS);
        }
    }

    fclose(fp);

    double probability = (double)impacts / N_SIMULATIONS;
    printf("\n--- RESULTADOS FINALES ---\n");
    printf("Impactos detectados: %d\n", impacts);
    printf("Probabilidad de impacto: %.4f%%\n", probability * 100.0);
    printf("Datos guardados en data/trajectories.csv\n");

    return 0;
}
