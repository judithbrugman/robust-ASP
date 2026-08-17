# Robust Appointment Scheduling for General Convex Uncertainty Sets

MATLAB code accompanying the paper  
**"Robust Appointment Scheduling for General Convex Uncertainty Sets"**.

This repository contains implementations of the proposed methods for polyhedral and ellipsoidal uncertainty sets, together with scripts for the numerical experiments.

## Requirements

The numerical experiments in the paper were run using:

- MATLAB R2023b
- YALMIP R20210331
- MOSEK 10.0.27
- Gurobi 9.5.1

YALMIP and the required solvers should be installed separately and available on the MATLAB path.

## Repository structure

- `src/`: implementation of the optimization methods
- `scripts/`: scripts for running the numerical experiments
- `results/`: stored numerical results

## Running the code

The main scripts are:

- `run_budget_base.m`: base instance with budget uncertainty
- `run_ball_base.m`: base instance with ball uncertainty
- `run_budget_scalability.m`: scalability experiment with budget uncertainty

The scripts automatically add the `src/` directory to the MATLAB path.

For example, the base budget instance can be run from MATLAB using:

```matlab
run('scripts/run_budget_base.m')
```

## Results

The `results/` directory contains MATLAB data files with stored outputs of the numerical experiments.

The stored scalability results cover instances up to `N = 80`. The paper reports additional scalability results up to `N = 100`.

Computation times may differ from those reported in the paper depending on hardware and software versions.

## Authors

See `AUTHORS`.

## License

See `LICENSE`.
