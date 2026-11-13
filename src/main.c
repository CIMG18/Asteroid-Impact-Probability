#include "physics.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Configuración de la simulación ---
#define N_SIMULATIONS 1000
#define DT 86400.0          // 1 día en segundos
#define MAX_STEPS 365       // 365 días de simulación (~1 año)
#define UNCERTAINTY_POS 1e6 // 1,000 km
#define UNCERTAINTY_VEL 100 // 100 m/s

// Función para generar números aleatorios con distribución normal (Box-Muller)
double random_normal(double mu, double sigma) {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return z0 * sigma + mu;
}

int main(int argc, char *argv[]) {
    // -------------------------------------------------------------------------
    // Semilla aleatoria:
    //   ./main          → semilla basada en el tiempo (resultados distintos cada vez)
    //   ./main 42       → semilla fija 42 (resultados 100% reproducibles)
    // -------------------------------------------------------------------------
    unsigned int seed;
    if (argc >= 2) {
        seed = (unsigned int)atoi(argv[1]);
        printf("Semilla fija: %u (resultados reproducibles)\n", seed);
    } else {
        seed = (unsigned int)time(NULL);
        printf("Semilla aleatoria: %u\n", seed);
    }
    srand(seed);

    // -------------------------------------------------------------------------
    // 1. Condiciones iniciales de los cuerpos principales (S.I.)
    //
    // Tierra: órbita circular alrededor del Sol a 1 UA.
    //   - velocidad orbital: v = sqrt(G*M_sol / r) ≈ 29 784 m/s en +Y
    //
    // Luna: posición relativa a la Tierra + 3.844e8 m en X,
    //   velocidad orbital relativa a la Tierra: v = sqrt(G*M_tierra / r) ≈ 1 022 m/s en +Y
    //   (se suma la velocidad de la Tierra para obtener velocidad heliocéntrica)
    // -------------------------------------------------------------------------
    double r_earth  = 1.496e11;               // distancia Sol-Tierra [m]
    double r_moon   = 3.844e8;                // distancia Tierra-Luna [m]
    double v_earth  = sqrt(G * MASS_SUN   / r_earth);   // ~29 784 m/s
    double v_moon   = sqrt(G * MASS_EARTH / r_moon);    // ~1 022 m/s

    Body ic[3] = {
        {"Sol",    MASS_SUN,   {0,                   0, 0}, {0,                   0, 0}},
        {"Tierra", MASS_EARTH, {r_earth,              0, 0}, {0,          v_earth,   0}},
        {"Luna",   MASS_MOON,  {r_earth + r_moon,     0, 0}, {0,  v_earth + v_moon,  0}}
    };

    // 2. Condiciones base del meteorito (aproximación a ~0.94 UA del Sol, acercándose)
    Vector3 base_pos = {1.4e11, 1e9, 0};
    Vector3 base_vel = {0, 30000, 0};

    int impacts = 0;

    // Variables para diagnóstico de conservación de energía (simulación 0 como referencia)
    double E0_ref   = 0.0; // energía total inicial de la sim. de referencia
    double E_max_err= 0.0; // máximo error relativo observado
    double E_fin    = 0.0; // energía total final de la sim. de referencia

    // Abrir archivo para guardar trayectorias de algunos meteoritos (para visualización)
    FILE *fp = fopen("data/trajectories.csv", "w");
    if (fp == NULL) {
        perror("Error abriendo archivo de datos");
        return 1;
    }
    // Formato CSV extendido: incluye posición de Tierra y Luna para el visualizador
    fprintf(fp, "sim_id,step,x,y,z,vx,vy,vz,ek,ep,etot,ex,ey,mx,my\n");

    printf("Iniciando Montecarlo en C: %d simulaciones, %d dias cada una...\n",
           N_SIMULATIONS, MAX_STEPS);

    for (int i = 0; i < N_SIMULATIONS; i++) {
        // Resetear cuerpos principales a condiciones iniciales en cada simulación
        Body bodies[3];
        for (int b = 0; b < 3; b++) bodies[b] = ic[b];

        // Muestreo aleatorio de condiciones iniciales del meteorito
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
        double mass_meteorite = 1e10; // Masa arbitraria para cálculos de energía

        // Energía inicial de este meteorito (para diagnóstico en sim. 0)
        double ek0 = compute_kinetic_energy(vel, mass_meteorite);
        double ep0 = compute_potential_energy(pos, bodies, 3, mass_meteorite);
        if (i == 0) E0_ref = ek0 + ep0;

        for (int step = 0; step < MAX_STEPS; step++) {
            // Guardar 1 de cada 10 simulaciones (asteroides) con posición de Tierra y Luna
            if (i % 10 == 0) {
                double ek = compute_kinetic_energy(vel, mass_meteorite);
                double ep = compute_potential_energy(pos, bodies, 3, mass_meteorite);
                fprintf(fp, "%d,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",
                        i, step,
                        pos.x, pos.y, pos.z,
                        vel.x, vel.y, vel.z,
                        ek, ep, ek + ep,
                        bodies[1].pos.x, bodies[1].pos.y,   // Tierra x,y
                        bodies[2].pos.x, bodies[2].pos.y);  // Luna  x,y

                // Seguimiento del error de energía en la simulación 0
                if (i == 0 && E0_ref != 0.0) {
                    double err = fabs((ek + ep - E0_ref) / E0_ref);
                    if (err > E_max_err) E_max_err = err;
                    if (step == MAX_STEPS - 1) E_fin = ek + ep;
                }
            }

            // Avanzar los cuerpos principales entre sí (N-cuerpos real)
            step_verlet_bodies(bodies, 3, DT);

            // Avanzar el meteorito bajo la gravedad de los cuerpos actualizados
            step_verlet(&pos, &vel, bodies, 3, DT);

            // Verificar colisión con la Tierra (bodies[1])
            Vector3 diff = vec_sub(pos, bodies[1].pos);
            if (vec_norm(diff) <= RADIUS_EARTH) {
                impacted = true;
                break;
            }
        }

        if (impacted) impacts++;

        if ((i + 1) % 100 == 0) {
            printf("Progreso: %d/%d simulaciones completadas.\n", i + 1, N_SIMULATIONS);
        }
    }

    fclose(fp);

    double probability = (double)impacts / N_SIMULATIONS;

    printf("\n╔══════════════════════════════════════╗\n");
    printf("║        RESULTADOS FINALES            ║\n");
    printf("╠══════════════════════════════════════╣\n");
    printf("║  Semilla usada       : %-14u ║\n", seed);
    printf("║  Simulaciones        : %-14d ║\n", N_SIMULATIONS);
    printf("║  Duracion por sim.   : %-10d dias ║\n", MAX_STEPS);
    printf("╠══════════════════════════════════════╣\n");
    printf("║  Impactos detectados : %-14d ║\n", impacts);
    printf("║  Probabilidad        : %-13.4f%% ║\n", probability * 100.0);
    printf("╠══════════════════════════════════════╣\n");
    printf("║  CONSERVACION DE ENERGIA (sim. 0)    ║\n");
    printf("║  E inicial : %+.6e J  ║\n", E0_ref);
    printf("║  E final   : %+.6e J  ║\n", E_fin);
    printf("║  Error max : %-13.6f%%  ║\n", E_max_err * 100.0);
    if (E_max_err < 1e-4)
        printf("║  Integrador: EXCELENTE (<0.01%%)     ║\n");
    else if (E_max_err < 1e-2)
        printf("║  Integrador: ACEPTABLE (<1.00%%)     ║\n");
    else
        printf("║  Integrador: REVISAR   (>1.00%%)     ║\n");
    printf("╚══════════════════════════════════════╝\n");
    printf("\nDatos guardados en data/trajectories.csv\n");

    return 0;
}
