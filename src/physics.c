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
