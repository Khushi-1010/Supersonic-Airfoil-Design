# Double Wedge Supersonic Airfoil Optimization

Find a double-wedge airfoil that maximizes L/D at Mach 2.6 and 5 deg AoA while meeting lift and pitch-moment constraints. Aerodynamics are computed with classical shock–expansion theory for 2D supersonic flow.

## Problem in one line
Search over vertex locations and thickness split to maximize `Cl/Cd` subject to `Cl ≥ 0.2`, `|Cm_LE| ≤ 0.05`, and `tu + tl = 0.1`.

## What this repo contains
- `optimize_double_wedge_airfoil.m` — exhaustive search with clear comments and a final plot of the best shape
- Helper functions you must provide on the MATLAB path:
  - `oswbeta`, `nswr_newmach`, `nswr_pressure`, `pm`, `stag2stat`

## How it works
- Geometry: two straight panels per side defined by `(xu, tu)` and `(xl, tl)`
- Physics: oblique shock at the leading edge if required, Prandtl–Meyer expansion at vertices, quasi-1D relations for pressure and Mach updates
- Objective: maximize `Cl/Cd` among feasible designs
- Selection: script reports the best design and plots the airfoil

## Quick start
1. Put the helper aerodynamic functions on your MATLAB path.
2. Open `optimize_double_wedge_airfoil.m`.
3. Adjust grid resolutions `tu_grid`, `xu_grid`, `xl_grid` for speed vs accuracy.
4. Run. The console prints the best parameters and L/D, and a figure shows the optimal shape.

## Tuning tips
- Start coarse, for example step of `0.01`, then refine to `0.002` or `0.001`.
- If no feasible designs are found, widen the grid or relax step size to avoid skipping narrow feasible regions.
- Keep `tu + tl = 0.1` enforced by construction as done in the script.

## Citation
This implementation follows the course assignment brief on supersonic double-wedge airfoil design using shock–expansion analysis. :contentReference[oaicite:0]{index=0}

