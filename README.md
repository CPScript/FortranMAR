# Projectile Motion Simulation with Air Resistance

> A physics project using projectile motion with quadratic air drag using Fortran for numerical simulation and R for data analysis and visualization.

## Overview

This simulates the trajectory of a projectile under the influence of gravity and air resistance. The simulation uses 4th-order Runge-Kutta integration to solve the coupled differential equation system and provides comprehensive analysis of the results.

<img width="720" height="631" alt="image" src="https://github.com/user-attachments/assets/d5d05b76-1727-48d9-9b10-6a91eb4dfd7e" />


## Files

- projectile_motion.f90 - Main Fortran simulation program
- trajectory_analysis.R - R script for data analysis and visualization
- Makefile - Build system configuration

## Requirements

- gfortran compiler (GNU Fortran)
- R statistical software
- R packages: ggplot2, dplyr, gridExtra, viridis, plotly, tidyr, htmlwidgets

## Building

To compile the Fortran program:

```
make
```

For debug build:

```
make debug
```

To install R dependencies:

```
make install-deps
```

## Running

Execute the simulation:

```
./projectile_motion
```

Run complete analysis pipeline:

```
make analyze
```

This will run the Fortran simulation followed by R analysis and generate all output files.

## Output Files

- trajectory_data.csv - Raw simulation data with time, position, and velocity
- processed_trajectory_data.csv - Enhanced dataset with derived calculations
- trajectory_analysis.png - Four-panel analysis plot
- trajectory_comparison.png - Comparison with and without air resistance
- interactive_trajectory.html - 3D interactive visualization

## Physics Model

The simulation implements:

- Quadratic air drag: F_drag = 0.5 * rho * Cd * A * v^2
- Gravitational acceleration: 9.81 m/s^2
- Baseball parameters: mass 0.145 kg, radius 0.037 m
- Initial conditions: velocity 45 m/s at 35 degree angle

## Analysis Features

- Trajectory plotting with velocity color coding
- Velocity component analysis over time
- Energy conservation analysis (kinetic, potential, total)
- Velocity angle tracking
- Comparative analysis with theoretical no-drag trajectory
- Statistical quantification of air resistance effects

## Build Targets

- make - Build release version
- make debug - Build with debugging symbols
- make run - Execute simulation
- make analyze - Run simulation and R analysis
- make clean - Remove build artifacts
- make distclean - Remove all generated files
- make help - Show all available targets

## Results

Typical simulation results show:

- Range reduction of approximately 15-20% due to air resistance
- Flight time reduction compared to vacuum trajectory
- Energy loss quantification due to drag forces
- Detailed velocity and position tracking throughout flight

## Compiler Support

Tested with:

- gfortran (GNU Fortran)
- Intel Fortran (ifort) with appropriate flag modifications

## Platform Compatibility

- Linux (primary development platform)
- macOS (with Homebrew gfortran)
- Windows (with MinGW-w64 or similar)
