# Asteroid-Impact-Probability
Real asteroid orbits are never known with perfect precision — telescope observations carry measurement error. Instead of predicting a single deterministic path, this project samples thousands of plausible trajectories consistent with that uncertainty and checks how many of them actually intersect Earth's effective gravitational cross-section.

## Data Sources 
1. JPL Small-Body Database (SBDB) API https://ssd-api.jpl.nasa.gov/sbdb.api Provides an asteroid's orbital elements and their 1-sigma uncertainties

2. JPL Sentry System API https://ssd-api.jpl.nasa.gov/sentry.api Lists real near-Earth objects currently monitored for potential impact, including relative velocity, absolute magnitude, estimated diameter, and JPL's own official cumulative impact probability. 

## Objetive Variable
- Probability of impact which is a fraction of simulated trajectories where the sampled impact parameter falls inside the effective impact radius. 