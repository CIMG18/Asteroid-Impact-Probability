#include "physics.h"
#include <math.h>

// --- Configuración de la simulación ---
#define DT       86400.0  // paso de tiempo: 1 día en segundos
#define N_YEARS  10       // duración total de la simulación en años
#define MAX_STEPS (N_YEARS * 365)  // pasos totales

int main(void) {
    // -------------------------------------------------------------------------
    // Condiciones iniciales de los cuerpos (S.I.)
    //
    // Tierra: órbita ELÍPTICA real (e = 0.0167).
    //   El cuerpo se coloca en el PERIHELIO (punto más cercano al Sol).
    //   La velocidad en el perihelio se obtiene de la conservación del momento
    //   angular para una elipse kepleriana:
    //
    //     r_peri  = a * (1 - e)          a = semieje mayor = 1 UA
    //     v_peri  = sqrt(G*M * (1+e) / (a*(1-e)))
    //
    //   Con e = 0.0167:
    //     r_peri  ≈ 1.4710×10¹¹ m  (147.1 millones de km)
    //     r_aph   ≈ 1.5210×10¹¹ m  (152.1 millones de km)  — afelio opuesto
    //     v_peri  ≈ 30 287 m/s     (más rápido que en órbita circular)
    //
    // Luna: órbita circular alrededor de la Tierra (e_luna ≈ 0.0549 real,
    //   simplificada a circular para mantener el enfoque en la órbita terrestre).
    //   v = sqrt(G * M_tierra / r_luna) ≈ 1 022 m/s en +Y (relativa a Tierra)
    // -------------------------------------------------------------------------
    double e_earth  = 0.0167;               // excentricidad orbital de la Tierra
    double a_earth  = 1.496e11;             // semieje mayor [m] (~1 UA)
    double r_peri   = a_earth * (1.0 - e_earth);  // distancia en perihelio [m]
    double v_peri   = sqrt(G * MASS_SUN * (1.0 + e_earth) / (a_earth * (1.0 - e_earth)));

    double r_moon   = 3.844e8;              // distancia Tierra-Luna [m]
    double v_moon   = sqrt(G * MASS_EARTH / r_moon);  // ~1 022 m/s

    Body bodies[3] = {
        {"Sol",    MASS_SUN,   {0,                0, 0}, {0,                 0, 0}},
        {"Tierra", MASS_EARTH, {r_peri,            0, 0}, {0,        v_peri,   0}},
        {"Luna",   MASS_MOON,  {r_peri + r_moon,   0, 0}, {0, v_peri + v_moon,  0}}
    };

    // -------------------------------------------------------------------------
    // Diagnóstico de conservación de energía
    // -------------------------------------------------------------------------
    double E0_ref    = 0.0;
    double E_max_err = 0.0;

    // Energía inicial del sistema
    {
        double K = 0.0, U = 0.0;
        for (int i = 0; i < 3; i++) {
            double v2 = bodies[i].vel.x * bodies[i].vel.x
                      + bodies[i].vel.y * bodies[i].vel.y
                      + bodies[i].vel.z * bodies[i].vel.z;
            K += 0.5 * bodies[i].mass * v2;
        }
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                Vector3 d = vec_sub(bodies[j].pos, bodies[i].pos);
                U -= G * bodies[i].mass * bodies[j].mass / vec_norm(d);
            }
        }
        E0_ref = K + U;
    }

    // -------------------------------------------------------------------------
    // Abrir CSV
    // -------------------------------------------------------------------------
    FILE *fp = fopen("data/trajectories.csv", "w");
    if (fp == NULL) {
        perror("Error abriendo data/trajectories.csv");
        return 1;
    }
    fprintf(fp, "step,sx,sy,ex,ey,mx,my\n");
    // Columnas: paso | Sol x,y | Tierra x,y | Luna x,y

    printf("Simulando sistema Sol-Tierra-Luna: %d anios (%d pasos de 1 dia)...\n",
           N_YEARS, MAX_STEPS);

    // -------------------------------------------------------------------------
    // Loop de integración N-cuerpos
    // -------------------------------------------------------------------------
    for (int step = 0; step < MAX_STEPS; step++) {
        // Guardar posiciones actuales
        fprintf(fp, "%d,%e,%e,%e,%e,%e,%e\n",
                step,
                bodies[0].pos.x, bodies[0].pos.y,   // Sol
                bodies[1].pos.x, bodies[1].pos.y,   // Tierra
                bodies[2].pos.x, bodies[2].pos.y);  // Luna

        // Integrar un paso (todos los cuerpos entre sí)
        step_verlet_bodies(bodies, 3, DT);

        // Diagnóstico de energía
        double K = 0.0, U = 0.0;
        for (int i = 0; i < 3; i++) {
            double v2 = bodies[i].vel.x * bodies[i].vel.x
                      + bodies[i].vel.y * bodies[i].vel.y
                      + bodies[i].vel.z * bodies[i].vel.z;
            K += 0.5 * bodies[i].mass * v2;
        }
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                Vector3 d = vec_sub(bodies[j].pos, bodies[i].pos);
                U -= G * bodies[i].mass * bodies[j].mass / vec_norm(d);
            }
        }
        double E = K + U;
        double err = fabs((E - E0_ref) / E0_ref);
        if (err > E_max_err) E_max_err = err;
    }

    fclose(fp);

    printf("\n╔══════════════════════════════════════╗\n");
    printf("║   SISTEMA SOL-TIERRA-LUNA            ║\n");
    printf("╠══════════════════════════════════════╣\n");
    printf("║  Duracion        : %-10d anios  ║\n", N_YEARS);
    printf("║  Pasos           : %-14d ║\n", MAX_STEPS);
    printf("║  DT              : 1 dia (86400 s)  ║\n");
    printf("╠══════════════════════════════════════╣\n");
    printf("║  CONSERVACION DE ENERGIA             ║\n");
    printf("║  E inicial : %+.6e J  ║\n", E0_ref);
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
