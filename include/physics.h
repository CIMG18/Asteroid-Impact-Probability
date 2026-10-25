#ifndef PHYSICS_H
#define PHYSICS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// --- Constantes Astronómicas (S.I.) ---
#define G 6.67430e-11
#define MASS_SUN 1.989e30
#define MASS_EARTH 5.972e24
#define MASS_MOON 7.348e22
#define RADIUS_EARTH 6.371e6
#define EPSILON 1e-3 // Suavizado para evitar divergencias

typedef struct {
    double x, y, z;
} Vector3;

typedef struct {
    char name[32];
    double mass;
    Vector3 pos;
} Body;

// Operaciones básicas de vectores
Vector3 vec_add(Vector3 a, Vector3 b);
Vector3 vec_sub(Vector3 a, Vector3 b);
Vector3 vec_mult(Vector3 a, double scalar);
double vec_norm(Vector3 a);

// Motor Físico
Vector3 compute_acceleration(Vector3 pos_meteorite, Body *bodies, int num_bodies);
void step_verlet(Vector3 *pos, Vector3 *vel, Body *bodies, int num_bodies, double dt);
double compute_kinetic_energy(Vector3 vel, double mass);
double compute_potential_energy(Vector3 pos_meteorite, Body *bodies, int num_bodies, double mass_meteorite);


#endif
