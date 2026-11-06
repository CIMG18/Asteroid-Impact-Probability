#include "physics.h"

Vector3 vec_add(Vector3 a, Vector3 b) {
    return (Vector3){a.x + b.x, a.y + b.y, a.z + b.z};
}

Vector3 vec_sub(Vector3 a, Vector3 b) {
    return (Vector3){a.x - b.x, a.y - b.y, a.z - b.z};
}

Vector3 vec_mult(Vector3 a, double scalar) {
    return (Vector3){a.x * scalar, a.y * scalar, a.z * scalar};
}

double vec_norm(Vector3 a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

Vector3 compute_acceleration(Vector3 pos_meteorite, Body *bodies, int num_bodies) {
    Vector3 accel_total = {0, 0, 0};

    for (int i = 0; i < num_bodies; i++) {
        Vector3 r_vector = vec_sub(bodies[i].pos, pos_meteorite);
        double dist = vec_norm(r_vector);

        // Fórmula con suavizado: a = G*M / (dist^2 + eps^2)^(3/2) * r_vector
        double mag = (G * bodies[i].mass) / pow(dist * dist + EPSILON * EPSILON, 1.5);
        Vector3 accel_body = vec_mult(r_vector, mag);

        accel_total = vec_add(accel_total, accel_body);
    }

    return accel_total;
}

void step_verlet(Vector3 *pos, Vector3 *vel, Body *bodies, int num_bodies, double dt) {
    // 1. Calcular aceleración actual
    Vector3 a_current = compute_acceleration(*pos, bodies, num_bodies);

    // 2. Actualizar posición: r(t+dt) = r(t) + v(t)dt + 0.5*a(t)dt^2
    Vector3 term1 = vec_mult(*vel, dt);
    Vector3 term2 = vec_mult(a_current, 0.5 * dt * dt);
    *pos = vec_add(*pos, vec_add(term1, term2));

    // 3. Calcular nueva aceleración en la nueva posición
    Vector3 a_next = compute_acceleration(*pos, bodies, num_bodies);

    // 4. Actualizar velocidad: v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))dt
    Vector3 a_avg = vec_mult(vec_add(a_current, a_next), 0.5 * dt);
    *vel = vec_add(*vel, a_avg);
}

// Integra todos los cuerpos principales entre sí usando Velocity Verlet (N-cuerpos).
// Calcula aceleraciones actuales, avanza posiciones, recalcula aceleraciones y promedia
// para actualizar velocidades — igual que step_verlet pero aplicado a todos los pares.
void step_verlet_bodies(Body *bodies, int num_bodies, double dt) {
    // Arreglo temporal de aceleraciones actuales (tamaño máximo estático)
    // Usamos VLA o un máximo razonable; 16 cuerpos es más que suficiente para este proyecto
    Vector3 a_curr[16];
    Vector3 a_next[16];

    // 1. Calcular aceleración actual de cada cuerpo debida a los demás
    for (int i = 0; i < num_bodies; i++) {
        a_curr[i] = (Vector3){0, 0, 0};
        for (int j = 0; j < num_bodies; j++) {
            if (i == j) continue;
            Vector3 r = vec_sub(bodies[j].pos, bodies[i].pos);
            double dist = vec_norm(r);
            double mag = (G * bodies[j].mass) / (dist * dist * dist);
            a_curr[i] = vec_add(a_curr[i], vec_mult(r, mag));
        }
    }

    // 2. Avanzar posiciones: r(t+dt) = r(t) + v(t)dt + 0.5*a(t)dt^2
    for (int i = 0; i < num_bodies; i++) {
        Vector3 dr = vec_add(vec_mult(bodies[i].vel, dt),
                             vec_mult(a_curr[i], 0.5 * dt * dt));
        bodies[i].pos = vec_add(bodies[i].pos, dr);
    }

    // 3. Recalcular aceleraciones en las nuevas posiciones
    for (int i = 0; i < num_bodies; i++) {
        a_next[i] = (Vector3){0, 0, 0};
        for (int j = 0; j < num_bodies; j++) {
            if (i == j) continue;
            Vector3 r = vec_sub(bodies[j].pos, bodies[i].pos);
            double dist = vec_norm(r);
            double mag = (G * bodies[j].mass) / (dist * dist * dist);
            a_next[i] = vec_add(a_next[i], vec_mult(r, mag));
        }
    }

    // 4. Actualizar velocidades: v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))dt
    for (int i = 0; i < num_bodies; i++) {
        Vector3 a_avg = vec_mult(vec_add(a_curr[i], a_next[i]), 0.5 * dt);
        bodies[i].vel = vec_add(bodies[i].vel, a_avg);
    }
}

double compute_kinetic_energy(Vector3 vel, double mass) {
    double v2 = vel.x * vel.x + vel.y * vel.y + vel.z * vel.z;
    return 0.5 * mass * v2;
}

double compute_potential_energy(Vector3 pos_meteorite, Body *bodies, int num_bodies, double mass_meteorite) {
    double u_total = 0;
    for (int i = 0; i < num_bodies; i++) {
        Vector3 r_vector = vec_sub(bodies[i].pos, pos_meteorite);
        double dist = vec_norm(r_vector);
        u_total -= (G * bodies[i].mass * mass_meteorite) / dist;
    }
    return u_total;
}
